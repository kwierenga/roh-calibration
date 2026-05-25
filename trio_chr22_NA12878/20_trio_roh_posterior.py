"""
Trio ROH autozygosity tool (clinical-facing prototype).

Given a proband (optionally with both parents) it:
  1. Calls runs of homozygosity (ROH) in the proband over common SNPs, using the
     clinical rules established in this project (max-SNP-gap 1 Mb to avoid marker
     deserts; tolerate isolated genotyping errors; minimum reported length).
  2. For each ROH computes the per-locus calibrated readout: local recombination
     rate r and background homozygosity H-bar -> the closed-form minimum callable
     length, EMPIRICALLY CORRECTED (the analytic block-independent value is ~2x
     anti-conservative; see scripts 18-19), and a confident-autozygous CALL plus
     the (labelled) analytic posterior.
  3. Uses the PARENTS to confirm true autozygosity: counts Mendelian-inconsistent
     sites within each ROH (proband homozygous for an allele a parent lacks),
     flagging possible deletion/hemizygosity rather than autozygosity.
  4. Pedigree-F check: observed autozygosity (ROH burden) vs the value expected
     for the stated parental relationship.
  5. If candidate variant positions are supplied, reports whether each falls in a
     confidently-autozygous ROH (i.e., whether a homozygous variant there is
     plausibly autozygous/IBD).

This operationalizes the calibrated posterior for a single case. It is a
prototype: thresholds are the project's genome-wide values; platform-specific
(SNP-array) calibration and a per-locus empirical correction are future work.

Usage (defaults run the outbred NA12878 CEU trio on chr22 as a pipeline smoke test):
  python 20_trio_roh_posterior.py
  python 20_trio_roh_posterior.py --vcf=chr22_phased.vcf.gz --proband=NA12878 \
         --father=NA12891 --mother=NA12892 --pop=EUR --relationship=first_cousins \
         --variants=chr22:21000000,chr22:38000000
"""

import gzip
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DIVTSV = HERE / "cross_pop_hap_diversity.tsv"      # per-window Hbar + cMperMb (genome-wide H run)

# ---- clinical-calling knobs ----
MAF_MIN = 0.05
GAP_TOL = 1                       # tolerate this many isolated het (error) sites in a run
MAX_SNP_GAP_BP = 1_000_000        # break a run across a common-SNP desert (cf. PLINK --homozyg-gap)
MIN_ROH_MB = 0.5                  # report floor
WINDOW_BP = 1_000_000
BLOCK_CM = 0.5
H_FLOOR = 1e-4
GENO_ERR = 0.001
T_DECISION = 0.95
DEFAULT_PI = 0.0625
CALIB = 2.1                       # empirical correction to the analytic L* (gap1, LD heavy tail; scripts 18-19)
DELETION_FRAC = 0.03             # >this fraction of Mendel-inconsistent sites in a run -> flag possible deletion
POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
EXP_F = {"first_cousins": 0.0625, "second_cousins": 0.0156, "third_cousins": 0.0039,
         "avuncular": 0.125, "double_first_cousins": 0.125,
         "first_cousins_once_removed": 0.0313, "uncle_niece": 0.125,
         "incest_first_degree": 0.25, "unrelated": 0.0}


def args():
    a = {x.split("=", 1)[0].lstrip("-"): x.split("=", 1)[1]
         for x in sys.argv[1:] if "=" in x}
    a.setdefault("vcf", "chr22_phased.vcf.gz")
    a.setdefault("proband", "NA12878"); a.setdefault("father", "NA12891")
    a.setdefault("mother", "NA12892"); a.setdefault("pop", "EUR")
    a.setdefault("pi", str(DEFAULT_PI)); a.setdefault("relationship", "")
    a.setdefault("variants", "")
    return a


def load_div(pop):
    """{chrom:{window_start:(Hbar, cMperMb)}} for the chosen population."""
    out = {}
    if not DIVTSV.exists():
        return out
    with DIVTSV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["population"]] != pop:
                continue
            out.setdefault(f[ix["chrom"]], {})[int(f[ix["window_start"]])] = (
                float(f[ix["Hbar"]]), float(f[ix["cMperMb"]]))
    return out


def parse_trio(vcf, pop, ids):
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    col = {}
    for role, sid in ids.items():
        if sid and sid in samples:
            col[role] = samples.index(sid)
    if "proband" not in col:
        sys.exit(f"proband {ids['proband']} not in VCF")
    pre = f"AF_{pop}="
    chrom = None
    pos, P0, P1, F0, F1, M0, M1 = [], [], [], [], [], [], []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            af = None
            for kv in f[7].split(";"):
                if kv.startswith(pre):
                    try:
                        af = float(kv[len(pre):])
                    except ValueError:
                        af = None
                    break
            if af is None or min(af, 1 - af) < MAF_MIN:
                continue
            chrom = f[0]
            gts = f[9:]
            def al(role):
                c = col.get(role)
                if c is None:
                    return -1, -1
                g = gts[c]
                return (ord(g[0]) - 48, ord(g[2]) - 48)
            p = al("proband")
            if p[0] < 0:
                continue
            pos.append(int(f[1])); P0.append(p[0]); P1.append(p[1])
            fa = al("father"); mo = al("mother")
            F0.append(fa[0]); F1.append(fa[1]); M0.append(mo[0]); M1.append(mo[1])
    arr = lambda L: np.asarray(L, dtype=np.int8)
    return chrom, (np.asarray(pos, dtype=np.int64), arr(P0), arr(P1),
                   arr(F0), arr(F1), arr(M0), arr(M1)), col


def call_roh(pos, p0, p1):
    """Maximal homozygous runs in the proband (clinical rules). Returns list of
    (start_idx, end_idx) inclusive."""
    hom = (p0 == p1)
    m = hom.copy()
    if GAP_TOL > 0:                                  # bridge isolated het (error) sites
        pad = np.concatenate(([1], m.astype(np.int8), [1])); dd = np.diff(pad)
        hs = np.flatnonzero(dd == -1); he = np.flatnonzero(dd == 1)
        short = (he - hs) <= GAP_TOL
        if short.any():
            diff = np.zeros(m.size + 1, dtype=np.int32)
            np.add.at(diff, hs[short], 1); np.add.at(diff, he[short], -1)
            m = m | (np.cumsum(diff[:-1]) > 0)
    n = m.size
    intra = np.zeros(n, bool)
    intra[1:] = m[1:] & m[:-1] & ((pos[1:] - pos[:-1]) <= MAX_SNP_GAP_BP)
    starts = np.flatnonzero(m & ~intra)
    ends = m.copy(); ends[:-1] &= ~intra[1:]
    ends = np.flatnonzero(ends)
    return list(zip(starts, ends))


def min_callable(r, hbar, pi):
    b = max(hbar, H_FLOOR)
    if r <= 0 or b >= 1.0:
        return float("inf")
    num = pi * (1 - GENO_ERR) ** 1000
    rhs = num * (1 - T_DECISION) / (T_DECISION * (1 - pi))
    if rhs <= 0:
        return float("inf")
    return max(math.log(rhs) / math.log(b) * BLOCK_CM / r, BLOCK_CM / r)


def posterior(L, r, hbar, pi):
    if r <= 0 or not (hbar > 0):
        return pi
    pc = max(hbar, H_FLOOR) ** max(1.0, L * r / BLOCK_CM)
    num = pi * (1 - GENO_ERR) ** 1000
    return num / (num + (1 - pi) * pc)


def main():
    a = args()
    pop = a["pop"]; pi = float(a["pi"])
    vcf = HERE / a["vcf"]
    div = load_div(pop)
    chrom, D, col = parse_trio(vcf, pop, {"proband": a["proband"],
                                          "father": a["father"], "mother": a["mother"]})
    pos, p0, p1, f0, f1, m0, m1 = D
    has_parents = "father" in col and "mother" in col
    wins = div.get(chrom, {})
    def r_hbar(p):
        w = (int(p) // WINDOW_BP) * WINDOW_BP
        return wins.get(w, (float("nan"), float("nan")))

    runs = call_roh(pos, p0, p1)
    rows = []
    auto_mb = 0.0
    for s, e in runs:
        L = (pos[e] - pos[s]) / 1e6
        if L < MIN_ROH_MB:
            continue
        seg = slice(s, e + 1)
        rs = [r_hbar(p) for p in pos[seg]]
        rr = np.nanmean([x[1] for x in rs]); hh = np.nanmean([x[0] for x in rs])
        Lstar = CALIB * min_callable(rr, hh, pi)
        post = posterior(L, rr, hh, pi)
        call = "autozygous" if L >= Lstar else ("borderline" if L >= 0.7 * Lstar else "background")
        # trio Mendelian check within the run
        delflag = "-"
        share = "-"
        if has_parents:
            pa = p0[seg]                                   # proband allele (homozygous)
            fa_has = (f0[seg] == pa) | (f1[seg] == pa)
            mo_has = (m0[seg] == pa) | (m1[seg] == pa)
            valid = (f0[seg] >= 0) & (m0[seg] >= 0)
            inc = (~(fa_has & mo_has)) & valid
            frac = inc.mean() if valid.any() else float("nan")
            delflag = "POSSIBLE_DELETION" if frac > DELETION_FRAC else "ok"
            share = f"{1 - frac:.3f}" if valid.any() else "NA"
        if call == "autozygous" and delflag != "POSSIBLE_DELETION":
            auto_mb += L
        rows.append((chrom, int(pos[s]), int(pos[e]), round(L, 3), round(rr, 3),
                     round(hh, 4), round(Lstar, 2), round(post, 3), call, share, delflag))

    span_mb = (pos[-1] - pos[0]) / 1e6
    out = HERE / "trio_roh_report.tsv"
    with out.open("w", encoding="utf-8") as fh:
        fh.write("chrom\tstart\tend\tlength_Mb\tcMperMb\tHbar\tLstar_calib_Mb\t"
                 "posterior_analytic\tcall\tparental_sharing\tdeletion_flag\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"  proband={a['proband']} pop={pop} pi={pi}  chrom={chrom}  "
          f"common SNPs={pos.size}  scanned span={span_mb:.1f} Mb  "
          f"parents={'yes' if has_parents else 'no'}")
    print(f"  ROH (>= {MIN_ROH_MB} Mb): {len(rows)}  -> {out}")
    hdr = ["start", "end", "Mb", "cM/Mb", "Hbar", "L*calib", "post", "call",
           "share", "del?"]
    print("  " + "  ".join(f"{h:>9}" for h in hdr))
    for r in rows:
        print("  " + "  ".join(f"{str(x):>9}" for x in r[1:]))

    # pedigree-F (chromosome-scoped here; extrapolate over autosomes in production)
    F_obs = auto_mb / span_mb if span_mb else float("nan")
    print(f"\n  confidently-autozygous ROH burden (this chromosome): {auto_mb:.1f} Mb "
          f"of {span_mb:.1f} Mb  => F_obs(chr) = {F_obs:.4f}")
    if a["relationship"] in EXP_F:
        print(f"  expected F for '{a['relationship']}' = {EXP_F[a['relationship']]:.4f} "
              f"(genome-wide; compare against the all-autosome run)")
    if span_mb < 200:
        print("  NOTE: single-chromosome smoke scope — F is illustrative; run all "
              "autosomes for a real pedigree-F estimate.")

    # candidate variants
    if a["variants"]:
        print("\n  candidate variants:")
        for v in a["variants"].split(","):
            try:
                vc, vp = v.split(":"); vp = int(vp)
            except ValueError:
                continue
            hit = next((r for r in rows if r[0] == vc and r[1] <= vp <= r[2]), None)
            if hit:
                print(f"    {v}: inside ROH {hit[1]}-{hit[2]} ({hit[3]} Mb), "
                      f"call={hit[8]}, posterior~{hit[7]}, deletion_flag={hit[10]}")
            else:
                print(f"    {v}: not inside a called ROH (>= {MIN_ROH_MB} Mb)")


if __name__ == "__main__":
    main()
