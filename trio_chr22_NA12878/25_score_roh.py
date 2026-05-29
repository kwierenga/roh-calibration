"""
25_score_roh.py - reference implementation: score runs of homozygosity (ROH) for
recent autozygosity, calibrated per (locus, ancestry, platform, prior).

This is the deployable artifact behind the project's thesis that a fixed length
is not a fixed weight of evidence. For each ROH it reports the prior-free Bayes
factor (weight of evidence that the tract is recent IBD rather than the
individual's population background), the posterior at the chosen prior, the
locus-specific decisive length L*, and a FLAG / REVIEW / BACKGROUND decision.

CALIBRATION INPUTS (all produced earlier in this pipeline; nothing fitted here):
  trio_null_pchance.tsv         empirical per-population weight of evidence
                                log10 BF(L) = log10(c / p_background(L)), from the
                                leakage-free, cryptic-relatedness-screened trio
                                children null (script 21). This is the genome-wide
                                (median-rate) reference curve.
  cross_pop_hap_diversity.tsv   per 1-Mb-window deCODE recombination rate (cMperMb),
                                used for the per-locus adjustment.

MODEL (each step traceable to a prior script; the biostatistics partner should
vet the locus and platform adjustments):
  1. Weight of evidence scales with (recombination rate r) x (length L) (script 22):
     a length L at local rate r carries the evidence of an effective length
     L_eff = L * (r_locus / r_median) at the genome-wide median rate.
  2. Platform (script 23): sparser array markers give less evidence per Mb; an
     array ROH is scored at L_eff / array_penalty (WGS penalty = 1).
  3. BF read from the population curve at L_eff (interpolated); posterior at prior
     pi is  post = pi*BF / (pi*BF + (1 - pi)).
  4. Prior: --prior sets pi directly; --froh sets pi to the individual's genome-
     wide autozygosity fraction (F_ROH), operationalizing "condition on the
     individual's own ROH burden, not the superpopulation label."

Decision bands: FLAG post >= 0.95; REVIEW 0.50 <= post < 0.95; BACKGROUND < 0.50.

Usage:
  python 25_score_roh.py --demo
  python 25_score_roh.py --roh=chr2:1000000-4500000 --ancestry=EUR --prior=0.0625
  python 25_score_roh.py my_roh.bed --ancestry=SAS --platform=array --froh=0.03
    (BED/TSV: tab/space-separated chrom start end [...]; a header line is skipped)
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PCHANCE = HERE / "trio_null_pchance.tsv"
DIVTSV = HERE / "cross_pop_hap_diversity.tsv"

POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
GENO_ERR = 0.001
C_IBD = (1 - GENO_ERR) ** 1000          # emission constant, matches scripts 21-23
DEFAULT_PRIOR = 0.0625                   # first-cousin offspring
DEFAULT_ARRAY_PENALTY = 1.5             # array L* inflation vs dense WGS (script 23)
FLAG_POST, REVIEW_POST = 0.95, 0.50
PI_CLAMP = (1e-4, 0.5)

DEMO = [   # (chrom, start, end) - illustrative tracts spanning the decision bands
    ("chr1", 0, 600_000),               # 0.6 Mb, cold telomere      -> too short to call
    ("chr1", 3_000_000, 4_200_000),     # 1.2 Mb, recombination-rich -> credible despite length
    ("chr1", 0, 1_200_000),             # 1.2 Mb, recombination-poor -> SAME length, not decisive
    ("chr8", 5_000_000, 17_000_000),    # 12 Mb large tract          -> strongly credible
]


def load_curves():
    """Return (L grid array, {pop: log10BF array})."""
    with PCHANCE.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        bf_cols = {p: ix.get(f"{p}_log10BF_clean") for p in POPS}
        L, bf = [], {p: [] for p in POPS}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            L.append(float(f[0]))
            for p in POPS:
                c = bf_cols[p]
                bf[p].append(float(f[c]) if c is not None and c < len(f) else np.nan)
    L = np.asarray(L)
    bf = {p: np.asarray(v) for p, v in bf.items() if not np.all(np.isnan(bf[p]))}
    return L, bf


def load_rates():
    """Return (rate_by_window {chrom: {win_start: cMperMb}}, genome median rate)."""
    rate = {}
    allr = []
    with DIVTSV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                r = float(f[ix["cMperMb"]])
            except (ValueError, IndexError):
                continue
            ch = f[ix["chrom"]]
            ws = int(f[ix["window_start"]])
            rate.setdefault(ch, {})[ws] = r   # same rate across pops -> dedupe
            allr.append(r)
    pos = [r for r in allr if r > 0]
    return rate, float(np.median(pos)) if pos else 1.0


def locus_rate(rate, r_median, chrom, start, end):
    """Mean cMperMb over the 1-Mb windows the ROH overlaps; r_median if unknown."""
    wins = rate.get(chrom, {})
    if not wins:
        return r_median, False
    lo = (start // 1_000_000) * 1_000_000
    hi = ((end - 1) // 1_000_000) * 1_000_000
    vals = [wins[w] for w in range(lo, hi + 1, 1_000_000) if w in wins and wins[w] > 0]
    if not vals:
        return r_median, False
    return float(np.mean(vals)), True


def interp_log10bf(L_grid, bf_curve, L_eff):
    """Linear interpolation of log10 BF at L_eff, clamped to the grid ends."""
    if L_eff <= L_grid[0]:
        return float(bf_curve[0])
    if L_eff >= L_grid[-1]:
        return float(bf_curve[-1])
    return float(np.interp(L_eff, L_grid, bf_curve))


def posterior_from_bf(log10bf, pi):
    bf = 10.0 ** log10bf
    return pi * bf / (pi * bf + (1.0 - pi))


def lstar_eff(L_grid, bf_curve, pi):
    """Smallest effective length (Mb) where posterior >= FLAG_POST at prior pi."""
    post = posterior_from_bf(bf_curve, pi)
    hit = np.flatnonzero(post >= FLAG_POST)
    return float(L_grid[hit[0]]) if hit.size else float("inf")


def decide(post):
    return "FLAG" if post >= FLAG_POST else ("REVIEW" if post >= REVIEW_POST else "BACKGROUND")


def parse_roh_arg(s):
    chrom, rng = s.split(":")
    a, b = rng.split("-")
    return chrom, int(a), int(b)


def read_bed(path):
    out = []
    for raw in Path(path).read_text().splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        f = raw.replace(",", "").split()
        if len(f) < 3 or not f[1].isdigit() or not f[2].isdigit():
            continue   # header or malformed line
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def main():
    args = sys.argv[1:]

    def opt(name, default=None):
        hit = [a for a in args if a.startswith(f"--{name}=")]
        return hit[0].split("=", 1)[1] if hit else default

    ancestry = (opt("ancestry", "EUR") or "EUR").upper()
    platform = (opt("platform", "wgs") or "wgs").lower()
    array_penalty = float(opt("array-penalty", DEFAULT_ARRAY_PENALTY))
    froh = opt("froh")
    if froh is not None:
        pi = min(max(float(froh), PI_CLAMP[0]), PI_CLAMP[1])
        prior_src = f"F_ROH={froh} (individualized)"
    else:
        pi = float(opt("prior", DEFAULT_PRIOR))
        prior_src = f"prior={pi} (declared)"

    # collect ROH inputs
    rohs = []
    roharg = opt("roh")
    if "--demo" in args:
        rohs = list(DEMO)
    elif roharg:
        rohs = [parse_roh_arg(roharg)]
    else:
        files = [a for a in args if not a.startswith("--")]
        if files:
            rohs = read_bed(files[0])
    if not rohs:
        print(__doc__)
        sys.exit(0)

    L_grid, bf = load_curves()
    rate, r_median = load_rates()

    penalty = array_penalty if platform == "array" else 1.0
    if ancestry == "UNKNOWN":
        # most conservative population (largest decisive length) at this prior
        ancestry = max(bf, key=lambda p: lstar_eff(L_grid, bf[p], pi))
        anc_note = f"UNKNOWN -> conservative ({ancestry})"
    else:
        anc_note = ancestry
    if ancestry not in bf:
        sys.exit(f"ancestry {ancestry} not in calibration curves {sorted(bf)}")
    curve = bf[ancestry]
    lstar_e = lstar_eff(L_grid, curve, pi)

    print(f"# ROH autozygosity scoring (reference implementation v0)")
    print(f"# ancestry={anc_note}  platform={platform}(penalty x{penalty:g})  {prior_src}")
    print(f"# genome-wide median rate={r_median:.2f} cM/Mb; effective L* at median "
          f"rate={lstar_e:.2f} Mb; emission c={C_IBD:.3f}")
    print(f"# decision: FLAG post>={FLAG_POST}; REVIEW>={REVIEW_POST}; else BACKGROUND")
    print("chrom\tstart\tend\tL_Mb\tcMperMb\tL_eff_Mb\tlog10BF\tposterior\t"
          "Lstar_locus_Mb\tdecision")

    for chrom, start, end in rohs:
        L = (end - start) / 1e6
        r_loc, known = locus_rate(rate, r_median, chrom, start, end)
        L_eff = L * (r_loc / r_median) / penalty
        log10bf = interp_log10bf(L_grid, curve, L_eff)
        post = posterior_from_bf(log10bf, pi)
        # physical decisive length at this locus + platform
        lstar_locus = lstar_e * (r_median / r_loc) * penalty
        rate_str = f"{r_loc:.2f}" + ("" if known else "*")
        print(f"{chrom}\t{start}\t{end}\t{L:.2f}\t{rate_str}\t{L_eff:.2f}\t"
              f"{log10bf:.2f}\t{post:.3f}\t{lstar_locus:.2f}\t{decide(post)}")

    if any(not locus_rate(rate, r_median, c, s, e)[1] for c, s, e in rohs):
        print("# * recombination rate unknown for that locus; genome-wide median used")


if __name__ == "__main__":
    main()
