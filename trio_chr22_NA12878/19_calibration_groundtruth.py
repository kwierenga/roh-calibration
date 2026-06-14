"""
Ground-truth calibration of the ROH autozygosity posterior (Tier-1 rigor item).

WHY: every result so far is internally consistent (the posterior's noise term H̄
was validated against an empirical chance-IBS measurement), but the posterior has
never been checked against KNOWN truth: does "posterior 0.95" actually mean 95%
truly autozygous? This script answers that with a simulation built from REAL
haplotypes (so real LD structure is preserved; no msprime needed).

DESIGN (per population, simulation-based calibration):
  - Split the panel's individuals into a REFERENCE set and a disjoint TEST set.
  - Estimate the noise term H̄ (per 0.5 cM block -> per 1 Mb window) from the
    REFERENCE set only  ->  the posterior is NOT circular w.r.t. the test data.
  - Generate N_TRIALS labeled trials, autozygous with probability = prior π:
      * autozygous (positive): both homologs = ONE test haplotype over an implanted
        tract (length ~ Exp(mean IMPLANT_MEAN), truncated); independent test
        haplotypes elsewhere; then i.i.d. genotyping error ε flips alleles.
      * non-IBD (negative): the two homologs are haplotypes from two DIFFERENT
        test individuals (real chance-IBS background) + error.
  - Call the homozygous run spanning the locus with the CLINICAL rules
    (common SNPs, max-SNP-gap 1 Mb, tolerate GAP_TOL isolated errors); keep trials
    where a run is actually observed (length >= MIN_CALL_MB), as in practice.
  - Compute the model posterior P(IBD | L_obs, locus, pop) and compare to the
    label: reliability diagram (binned posterior vs observed fraction autozygous),
    expected calibration error (ECE), and the realized false-discovery rate among
    calls with posterior >= 0.95 (should be ~5% if calibrated).

LIMITATIONS (documented; for the biostatistics partner): calibration is reported
under an assumed autozygous-tract-length distribution and prior π; the negative
class still derives from the same population panel (cross-individual pairs include
cryptic distant relatedness), so an independent coalescent null (msprime/SLiM) is
the planned cross-check. The emission term c is the model's flat (1-ε)^1000.

Output: calibration_reliability.tsv (per-bin), calibration_summary.txt.
Usage:  python 19_calibration_groundtruth.py [chrom ...]   (default chr22)
"""

import gzip
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"
PANEL = HERE / "samples_2504_pop.panel"
OUT_REL = HERE / "calibration_reliability.tsv"
OUT_SUMMARY = HERE / "calibration_summary.txt"

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]   # all 5 superpops (genome-wide); EUR/AFR were the v1 extremes
# ---- knobs ----
MAF_MIN = 0.05
BLOCK_CM = 0.5
MAX_SNPS_BLOCK = 25
N_REF = 150                            # individuals to estimate H̄
N_TEST = 150                          # disjoint individuals to build test diploids
N_TRIALS = 20000
PI = 0.0625
GENO_ERR = 0.001
GAP_TOL = 1                           # realistic operating point
MAX_SNP_GAP_BP = 1_000_000
MIN_CALL_MB = 0.10                    # a run must reach this to count as "an ROH"
IMPLANT_MEAN = 2.0                    # Mb, mean of (truncated) autozygous tract length
IMPLANT_LO, IMPLANT_HI = 0.2, 12.0    # Mb truncation
MAXSCAN_BP = 8_000_000              # local slice half-width for run finding
H_FLOOR = 1e-4
RNG_SEED = 17
WINDOW_BP = 1_000_000
AF_PREFIX = {p: f"AF_{p}=" for p in POPULATIONS}
C_IBD = (1.0 - GENO_ERR) ** 1000     # model's (flat) emission term, under test


def load_rate(path, chrom):
    out = {}
    with path.open() as fh:
        for line in fh:
            if line[0] == "#" or line.startswith("Chr"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == chrom:
                try:
                    out[int(f[1])] = float(f[3])
                except (IndexError, ValueError):
                    pass
    return out


def build_rate_cum(chrom):
    pat, mat = load_rate(PAT_MAP, chrom), load_rate(MAT_MAP, chrom)
    if not pat:
        return None
    nwin = max(pat) // WINDOW_BP + 2
    rate = np.zeros(nwin)
    for w in range(nwin):
        c = w * WINDOW_BP + 500_000
        for cc in (c, c - WINDOW_BP, c + WINDOW_BP):
            if cc in pat:
                rate[w] = 0.5 * (pat[cc] + mat.get(cc, pat[cc]))
                break
    cum = np.concatenate(([0.0], np.cumsum(rate)))   # cum cM at window starts
    return rate, cum


def cum_cm(pos, rate, cum):
    w = pos // WINDOW_BP
    if w + 1 >= len(cum):
        return cum[-1]
    return cum[w] + rate[w] * (pos - w * WINDOW_BP) / WINDOW_BP


def load_panel():
    pop_of = {}
    with PANEL.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] in POPULATIONS:
                pop_of[f[0]] = f[2]
    return pop_of


def parse(chrom, vcf, pop_of, rng):
    """Return per pop: (Href int8[nsnp,2*Nref], Htest int8[nsnp,2*Ntest], pos)."""
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    bypop = {p: [] for p in POPULATIONS}
    for j, s in enumerate(samples):
        if pop_of.get(s) in POPULATIONS:
            bypop[pop_of[s]].append(j)
    sel, ref_cols, test_cols = {}, {}, {}
    for p in POPULATIONS:
        idx = bypop[p][:]
        rng.shuffle(idx)
        ref = sorted(idx[:N_REF])
        test = sorted(idx[N_REF:N_REF + N_TEST])
        sel[p] = ref + test
        ref_cols[p] = list(range(0, 2 * len(ref)))
        test_cols[p] = list(range(2 * len(ref), 2 * len(ref) + 2 * len(test)))
    rows = {p: [] for p in POPULATIONS}
    pos = {p: [] for p in POPULATIONS}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            common = []
            for kv in f[7].split(";"):
                if kv[:3] != "AF_":
                    continue
                for p, pre in AF_PREFIX.items():
                    if kv.startswith(pre):
                        try:
                            af = float(kv[len(pre):])
                        except ValueError:
                            af = None
                        if af is not None and min(af, 1 - af) >= MAF_MIN:
                            common.append(p)
                        break
            if not common:
                continue
            gt = f[9].split("\t")
            for p in common:
                rows[p].append(bytes(ord(gt[c][k]) - 48 for c in sel[p] for k in (0, 2)))
                pos[p].append(int(f[1]))
    out = {}
    for p in POPULATIONS:
        if not rows[p]:
            continue
        n = len(rows[p]); width = 2 * len(sel[p])
        mat = np.frombuffer(b"".join(rows[p]), dtype=np.int8).reshape(n, width)
        out[p] = (mat[:, ref_cols[p]].copy(), mat[:, test_cols[p]].copy(),
                  np.asarray(pos[p], dtype=np.int64))
    return out


def hbar_per_window(Href, pos, rate, cum):
    """Nei unbiased homozygosity per 0.5 cM block (thinned), averaged per 1 Mb window."""
    blk = np.array([int(cum_cm(p, rate, cum) // BLOCK_CM) for p in pos])
    nwin = int(pos[-1] // WINDOW_BP) + 1
    sumH = np.zeros(nwin + 1); cnt = np.zeros(nwin + 1)
    nhap = Href.shape[1]
    order = np.argsort(blk, kind="stable")
    bsort = blk[order]
    bounds = np.flatnonzero(np.diff(bsort)) + 1
    starts = np.concatenate(([0], bounds)); ends = np.concatenate((bounds, [len(bsort)]))
    for s, e in zip(starts, ends):
        rowidx = order[s:e]
        if rowidx.size < 1:
            continue
        if rowidx.size > MAX_SNPS_BLOCK:
            step = rowidx.size / MAX_SNPS_BLOCK
            rowidx = rowidx[(np.arange(MAX_SNPS_BLOCK) * step).astype(int)]
        sub = Href[rowidx, :]                          # [nsnp_blk, nhap]
        # haplotype string per column -> Nei unbiased homozygosity
        keys = np.ascontiguousarray(sub.T).view([('', sub.dtype)] * sub.shape[0])
        _, counts = np.unique(keys, return_counts=True)
        H = (counts * (counts - 1)).sum() / (nhap * (nhap - 1))
        w = int(pos[rowidx[0]] // WINDOW_BP)
        sumH[w] += H; cnt[w] += 1
    hbar = np.full(nwin + 1, np.nan)
    nz = cnt > 0
    hbar[nz] = sumH[nz] / cnt[nz]
    return hbar


def observed_run_mb(a, b, pos, k0):
    """Length (Mb) of the homozygous run spanning local index k0 (clinical rules)."""
    m = a == b
    if GAP_TOL > 0:
        pad = np.concatenate(([1], m.astype(np.int8), [1]))
        dd = np.diff(pad)
        ms = np.flatnonzero(dd == -1); me = np.flatnonzero(dd == 1)
        short = (me - ms) <= GAP_TOL
        if short.any():
            diff = np.zeros(m.size + 1, dtype=np.int32)
            np.add.at(diff, ms[short], 1); np.add.at(diff, me[short], -1)
            m = m | (np.cumsum(diff[:-1]) > 0)
    n = m.size
    intra = np.zeros(n, bool)
    intra[1:] = m[1:] & m[:-1] & ((pos[1:] - pos[:-1]) <= MAX_SNP_GAP_BP)
    if not m[k0]:
        return 0.0
    s = k0
    while s > 0 and intra[s]:
        s -= 1
    e = k0
    while e + 1 < n and intra[e + 1]:
        e += 1
    return (pos[e] - pos[s]) / 1e6


def posterior(L_obs, r, hbar):
    if not (hbar > 0) or r <= 0:
        return PI
    n_eff = max(1.0, L_obs * r / BLOCK_CM)
    pc = max(hbar, H_FLOOR) ** n_eff
    num = PI * C_IBD
    return num / (num + (1 - PI) * pc)


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    global N_TRIALS
    args = sys.argv[1:]
    tf = [a for a in args if a.startswith("--trials=")]
    if tf:
        N_TRIALS = int(tf[0].split("=", 1)[1])
        args = [a for a in args if not a.startswith("--trials=")]
    chroms = args or ["chr22"]
    rng = np.random.default_rng(RNG_SEED)
    pop_of = load_panel()
    trials = {p: {"L": [], "lab": [], "post": []} for p in POPULATIONS}

    for chrom in chroms:
        vcf = HERE / "chr22_phased.vcf.gz" if chrom == "chr22" else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"
        rc = build_rate_cum(chrom)
        if rc is None or not vcf.exists():
            print(f"  [{chrom}] SKIP"); continue
        rate, cum = rc
        data = parse(chrom, vcf, pop_of, rng)
        for p in POPULATIONS:
            if p not in data:
                continue
            Href, Htest, pos = data[p]
            hbar = hbar_per_window(Href, pos, rate, cum)
            ntest_ind = Htest.shape[1] // 2
            nsnp = pos.size
            per_pop_trials = N_TRIALS
            for _ in range(per_pop_trials):
                k0 = int(rng.integers(0, nsnp))
                center = int(pos[k0])
                lo = np.searchsorted(pos, center - MAXSCAN_BP)
                hi = np.searchsorted(pos, center + MAXSCAN_BP)
                if hi - lo < 5:
                    continue
                kloc = k0 - lo
                ia = int(rng.integers(0, ntest_ind)); ha = int(rng.integers(0, 2))
                a = Htest[lo:hi, 2 * ia + ha].copy()
                ib = ia
                while ib == ia:
                    ib = int(rng.integers(0, ntest_ind))
                b = Htest[lo:hi, 2 * ib + int(rng.integers(0, 2))].copy()
                label = rng.random() < PI
                if label:                                  # implant autozygous tract
                    L = min(IMPLANT_HI, max(IMPLANT_LO, rng.exponential(IMPLANT_MEAN)))
                    half = int(L * 1e6 / 2)
                    seg = (pos[lo:hi] >= center - half) & (pos[lo:hi] <= center + half)
                    b[seg] = a[seg]
                a ^= (rng.random(a.size) < GENO_ERR)        # i.i.d. genotyping error
                b ^= (rng.random(b.size) < GENO_ERR)
                L_obs = observed_run_mb(a, b, pos[lo:hi], kloc)
                if L_obs < MIN_CALL_MB:
                    continue
                w = center // WINDOW_BP
                post = posterior(L_obs, rate[w] if w < len(rate) else 0.0,
                                 hbar[w] if w < len(hbar) else np.nan)
                T = trials[p]
                T["L"].append(L_obs); T["lab"].append(int(label)); T["post"].append(post)
            print(f"  [{chrom}/{p}] Hbar~{np.nanmean(hbar):.4f} trials scored "
                  f"({time.time()-t0:.0f}s)")
            sys.stdout.flush()

    def fdr_at(L, lab, thr):
        sel = L >= thr
        return float(1 - lab[sel].mean()) if sel.sum() else float("nan")

    def gt_lstar(L, lab, target=0.05):
        for thr in np.arange(0.1, 8.001, 0.05):
            sel = L >= thr
            if sel.sum() >= 20 and (1 - lab[sel].mean()) <= target:
                return float(thr)
        return float("inf")

    summary = {}
    with OUT_REL.open("w", encoding="utf-8") as fh, OUT_SUMMARY.open("w", encoding="utf-8") as fs:
        gw = len(chroms) >= 20
        fh.write("population\tbin\tmean_posterior\tobs_frac_autozygous\tn\n")
        fs.write("# Ground-truth calibration of the ROH autozygosity posterior\n")
        if gw:
            fs.write(f"# GENOME-WIDE: {len(chroms)} autosomes -- these supersede the "
                     "earlier chr22 prototype L* values.\n")
        else:
            fs.write("# PROTOTYPE: chr22 only -- a genome-wide rerun is pending and will "
                     "supersede the specific L* values below. Do NOT cite these as final.\n")
        fs.write(f"# chroms={','.join(chroms)} N_TRIALS={N_TRIALS}/pop PI={PI} "
                 f"GAP_TOL={GAP_TOL} eps={GENO_ERR} N_REF={N_REF} N_TEST={N_TEST}\n")
        fs.write(f"# emission c=(1-eps)^1000={C_IBD:.3f}; implant mean {IMPLANT_MEAN} Mb\n\n")
        fs.write("Per population: ground-truth minimum callable length L* (observed-"
                 "length threshold where realized FDR=5%), the realized FDR if one "
                 "instead used the analytic (~0.7 Mb) or empirical (~1.5 Mb) "
                 "thresholds, ECE of the analytic posterior, and FDR among "
                 "posterior>=0.95 calls.\n\n")
        for p in POPULATIONS:
            T = trials[p]
            L = np.asarray(T["L"]); lab = np.asarray(T["lab"]); post = np.asarray(T["post"])
            N = L.size
            if N == 0:
                continue
            ece = 0.0
            for i in range(10):
                m = (post >= i / 10) & (post < (i + 1) / 10) if i < 9 else (post >= 0.9)
                if m.sum():
                    mc = post[m].mean(); mo = lab[m].mean()
                    ece += m.sum() / N * abs(mo - mc)
                    fh.write(f"{p}\t{i/10:.1f}-{(i+1)/10:.1f}\t{mc:.3f}\t{mo:.3f}\t{int(m.sum())}\n")
            sel95 = post >= 0.95
            fdr95 = float(1 - lab[sel95].mean()) if sel95.sum() else float("nan")
            gtL = gt_lstar(L, lab)
            summary[p] = dict(N=N, gtL=gtL, fdr07=fdr_at(L, lab, 0.7),
                              fdr15=fdr_at(L, lab, 1.5), ece=ece, fdr95=fdr95)
            fs.write(f"{p}: scored={N}  ground-truth L*(FDR=5%)={gtL:.2f} Mb  | "
                     f"FDR@0.7Mb={summary[p]['fdr07']:.3f}  FDR@1.5Mb={summary[p]['fdr15']:.3f}  "
                     f"| analytic-posterior ECE={ece:.3f}  FDR@post>=0.95={fdr95:.3f}\n")
        if summary:
            gtLs = [s["gtL"] for s in summary.values() if not math.isnan(s["gtL"])]
            eces = [s["ece"] for s in summary.values()]
            f15s = [s["fdr15"] for s in summary.values() if not math.isnan(s["fdr15"])]
            f95s = [s["fdr95"] for s in summary.values() if not math.isnan(s["fdr95"])]
            scope = f"genome-wide ({len(chroms)} autosomes)" if gw else f"{','.join(chroms)}"
            lo_p = min(summary, key=lambda k: summary[k]["gtL"]) if gtLs else "?"
            hi_p = max(summary, key=lambda k: summary[k]["gtL"]) if gtLs else "?"
            fs.write(f"\nReading ({scope}): the ground-truth L*(FDR=5%) is population-"
                     f"dependent, ranging {min(gtLs):.2f}-{max(gtLs):.2f} Mb "
                     f"({lo_p} shortest, {hi_p} longest), bracketed by the analytic "
                     "(~0.7 Mb) and empirical (~1.5 Mb) thresholds. The empirical "
                     f"~1.5 Mb operating point gives realized FDR@1.5Mb = "
                     f"{min(f15s)*100:.1f}-{max(f15s)*100:.1f}% "
                     f"(target 5%), so it is {'CONSERVATIVE' if max(f15s) < 0.05 else 'near/above target'} "
                     "-- a defensible clinical choice trading sensitivity for specificity. "
                     "The analytic *posterior*, separately, is over-confident as a "
                     f"calibrated probability (ECE {min(eces):.2f}-{max(eces):.2f}; among "
                     f"post>=0.95 calls realized FDR {min(f95s)*100:.0f}-{max(f95s)*100:.0f}%), "
                     "so the closed-form law should be used for intuition, not for "
                     "clinical probability statements.\n")
        fs.write("\n(Calibration is conditional on the assumed tract-length "
                 "distribution and prior; negative class from same-panel pairs; an "
                 "independent coalescent null is the planned cross-check.)\n")
    print(f"\n  -> {OUT_REL}\n  -> {OUT_SUMMARY}\n  total {time.time()-t0:.0f}s")
    for p, s in summary.items():
        print(f"    {p}: ground-truth L*={s['gtL']:.2f} Mb | FDR@0.7={s['fdr07']:.3f} "
              f"FDR@1.5={s['fdr15']:.3f} | analytic ECE={s['ece']:.3f} "
              f"FDR@post>=.95={s['fdr95']:.3f}")


if __name__ == "__main__":
    main()
