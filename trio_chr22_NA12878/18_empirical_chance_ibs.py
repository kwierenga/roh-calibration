"""
Empirical chance-IBS validation of the H-bar noise term.

WHY: 16_haplotype_ibs_noise.py fixed the population ORDERING but is OVER-CONFIDENT
in magnitude -- it assumes the 0.5 cM blocks are independent (p_chance = H-bar **
n_eff), which ignores LD spilling across block boundaries, and its level is set by
the MAX_SNPS_BLOCK thinning knob. Result: an implausible ~0.3-0.4 Mb median
minimum-callable ROH length. This script measures the SAME quantity the posterior
needs -- P(a length-L run is identical-by-state by chance, in a non-IBD pair) --
directly and assumption-free, by walking random haplotype pairs from DIFFERENT
unrelated individuals within each superpopulation. Any IBS run between two such
haplotypes is chance IBS (the null the ROH posterior must beat).

For an outbred individual the two homologs are ~one random pair of population
haplotypes, so cross-individual IBS run lengths ARE the per-individual chance-ROH
background. No block-independence assumption, no thinning knob.

OUTPUTS
  empirical_chance_ibs_pchance.tsv   L grid x pop: empirical vs analytic p_chance
  empirical_chance_ibs_summary.txt   per-pop chance-IBS run-length percentiles +
                                     empirical vs analytic minimum-callable length

METHODS KNOBS (top of file, for the biostatistics partner):
  MAF_MIN        common-SNP threshold defining the haplotypes (matches script 16)
  N_SUBSAMPLE    unrelated individuals sampled per population
  N_PAIRS        random cross-individual haplotype pairs scanned per population
  GAP_SWEEP      list of genotyping-error tolerances (mismatches bridged within
                 an IBS run) evaluated in one parse pass; 0 = strict. Real ROH
                 callers tolerate a few, which only LENGTHENS chance runs, so 0 is
                 the conservative baseline.
  MAX_SNP_GAP_KB break a run across a common-SNP desert wider than this (kb), so
                 centromeres / assembly gaps cannot fabricate long runs.
  PI             prior used for the minimum-callable comparison

Usage:
  python 18_empirical_chance_ibs.py chr22   # smoke test
  python 18_empirical_chance_ibs.py         # all autosomes (heavier; background)
"""

import gzip
import random
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

OUT_PC = HERE / "empirical_chance_ibs_pchance.tsv"
OUT_SUMMARY = HERE / "empirical_chance_ibs_summary.txt"

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

# ---- methods knobs ----
MAF_MIN = 0.05
N_SUBSAMPLE = 100
N_PAIRS = 300            # override at CLI with --pairs=N
GAP_SWEEP = [0, 1, 2]    # genotyping-error tolerance(s) swept in ONE parse pass
                         # (0 = strict IBS). Override with --gap=N for a single value.
PI = 0.0625
MAX_SNP_GAP_KB = 1000    # break a run across a common-SNP desert wider than this
                         # (kb), so centromeres / assembly gaps cannot fake long
                         # runs (cf. PLINK --homozyg-gap). Raise high to disable.
MIN_KEEP_MB = 0.005      # drop runs below this from storage to bound memory; does
                         # NOT affect p_chance on L_GRID (min L=0.1 Mb), and true
                         # percentiles are reconstructed from the total run count
RNG_SEED = 17
# -----------------------

MAX_SNP_GAP_BP = MAX_SNP_GAP_KB * 1000
GENOTYPING_ERROR = 0.001
WINDOW_BP = 1_000_000
BLOCK_CM = 0.5
ACMG_THRESHOLD = 0.95
AF_PREFIX = {pop: f"AF_{pop}=" for pop in POPULATIONS}
L_GRID = np.round(np.arange(0.1, 20.001, 0.05), 3)          # Mb, fine grid
L_REPORT = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]             # Mb, tabulated

# H-bar from the verified chr22 smoke test (16_haplotype_ibs_noise.py), used only
# to draw the analytic over-confident curve side-by-side with the empirical one.
ANALYTIC_HBAR = {"EUR": 0.0295, "AFR": 0.0126, "EAS": 0.0313,
                 "SAS": 0.0247, "AMR": 0.0301}


def load_decode_rate(path, chrom):
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Chr"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            try:
                out[int(f[1])] = float(f[3])
            except (IndexError, ValueError):
                continue
    return out


def mean_cmpermb(chrom):
    pat = load_decode_rate(PAT_MAP, chrom)
    mat = load_decode_rate(MAT_MAP, chrom)
    if not pat:
        return None
    vals = [0.5 * (pat[p] + mat.get(p, pat[p])) for p in pat]
    return sum(vals) / len(vals)


def load_panel():
    pop_of = {}
    with PANEL.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] in POPULATIONS:
                pop_of[f[0]] = f[2]
    return pop_of


def select_columns(vcf_path, pop_of):
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    by_pop = {pop: [] for pop in POPULATIONS}
    for j, s in enumerate(samples):
        p = pop_of.get(s)
        if p:
            by_pop[p].append(j)
    rng = random.Random(RNG_SEED)
    sel = {}
    for pop in POPULATIONS:
        idx = by_pop[pop]
        if len(idx) > N_SUBSAMPLE:
            idx = sorted(rng.sample(idx, N_SUBSAMPLE))
        sel[pop] = idx
    return sel


def parse_chrom(chrom, vcf_path, pop_of):
    """Return {pop: (alleles int8 [n_snps, 2*Nsub], pos int64 [n_snps])}."""
    sel = select_columns(vcf_path, pop_of)
    rows = {pop: [] for pop in POPULATIONS}
    pos = {pop: [] for pop in POPULATIONS}
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            fields = line.split("\t", 9)
            ref, alt = fields[3], fields[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            info = fields[7]
            common = []
            for kv in info.split(";"):
                if kv[:3] != "AF_":
                    continue
                for pop, pre in AF_PREFIX.items():
                    if kv.startswith(pre):
                        try:
                            af = float(kv[len(pre):])
                        except ValueError:
                            af = None
                        if af is not None and min(af, 1 - af) >= MAF_MIN:
                            common.append(pop)
                        break
            if not common:
                continue
            p = int(fields[1])
            gtblob = fields[9].split("\t")
            for pop in common:
                cols = sel[pop]
                rows[pop].append(
                    "".join(gtblob[c][0] + gtblob[c][2] for c in cols).encode())
                pos[pop].append(p)
    out = {}
    for pop in POPULATIONS:
        if not rows[pop]:
            continue
        n = len(rows[pop])
        width = 2 * len(sel[pop])
        arr = np.frombuffer(b"".join(rows[pop]), dtype=np.int8).reshape(n, width)
        out[pop] = (arr, np.asarray(pos[pop], dtype=np.int64))
    return out


def ibs_seglens_mb(a, b, pos, gap_tol, max_snp_gap_bp):
    """Lengths (Mb) of maximal IBS runs between haplotype vectors a, b.
    gap_tol        = max consecutive mismatches bridged within a run (0 = strict).
    max_snp_gap_bp = a run is broken between two matching SNPs farther apart than
                     this, so common-SNP deserts (centromeres, assembly gaps)
                     cannot fabricate long runs (cf. PLINK --homozyg-gap)."""
    m = (a == b)
    if gap_tol > 0:
        # bridge maximal mismatch runs of length <= gap_tol (vectorized via a
        # difference array, so it is O(n) regardless of mismatch count)
        pad = np.concatenate(([1], m.astype(np.int8), [1]))
        dd = np.diff(pad)
        miss_starts = np.flatnonzero(dd == -1)        # mismatch-run starts (m idx)
        miss_ends = np.flatnonzero(dd == 1)           # exclusive ends
        short = (miss_ends - miss_starts) <= gap_tol
        if short.any():
            diff = np.zeros(m.size + 1, dtype=np.int32)
            np.add.at(diff, miss_starts[short], 1)
            np.add.at(diff, miss_ends[short], -1)
            m = m | (np.cumsum(diff[:-1]) > 0)
    n = m.size
    if n == 0:
        return np.empty(0)
    # a run continues from SNP k-1 to k only if both match AND are within the max
    # physical gap; otherwise k begins a new run (kills marker-desert artifacts).
    intra = np.zeros(n, dtype=bool)
    intra[1:] = m[1:] & m[:-1] & ((pos[1:] - pos[:-1]) <= max_snp_gap_bp)
    starts_mask = m & ~intra
    ends_mask = m.copy()
    ends_mask[:-1] &= ~intra[1:]
    start_idx = np.flatnonzero(starts_mask)
    end_idx = np.flatnonzero(ends_mask)
    if start_idx.size == 0:
        return np.empty(0)
    return (pos[end_idx] - pos[start_idx]) / 1e6


def posterior(p_chance, pi):
    p_ibd = (1.0 - GENOTYPING_ERROR) ** 1000
    num = pi * p_ibd
    return num / (num + (1.0 - pi) * p_chance)


def min_callable_from_pc(L_fine, pc_fine, pi):
    """Smallest L (Mb) where posterior(p_chance(L)) >= ACMG_THRESHOLD."""
    post = posterior(np.maximum(pc_fine, 1e-300), pi)
    hit = np.flatnonzero(post >= ACMG_THRESHOLD)
    return float(L_fine[hit[0]]) if hit.size else float("inf")


def main():
    t0 = time.time()
    global N_PAIRS
    args = sys.argv[1:]
    pflag = [a for a in args if a.startswith("--pairs=")]
    if pflag:
        N_PAIRS = int(pflag[0].split("=", 1)[1])
        args = [a for a in args if not a.startswith("--pairs=")]
    gflag = [a for a in args if a.startswith("--gap=")]
    gaps = [int(gflag[0].split("=", 1)[1])] if gflag else list(GAP_SWEEP)
    args = [a for a in args if not a.startswith("--gap=")]
    pop_of = load_panel()
    chroms = args or [f"chr{n}" for n in range(1, 23)]
    rng = np.random.default_rng(RNG_SEED)

    seglens = {g: {pop: [] for pop in POPULATIONS} for g in gaps}
    total_runs = {g: {pop: 0 for pop in POPULATIONS} for g in gaps}
    exposure = {pop: 0.0 for pop in POPULATIONS}  # gap-independent (span per pair)
    r_means = []

    for chrom in chroms:
        vp = HERE / "chr22_phased.vcf.gz" if chrom == "chr22" else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"
        if not vp.exists():
            print(f"  [{chrom}] SKIP (no VCF)")
            continue
        r = mean_cmpermb(chrom)
        if r is None:
            print(f"  [{chrom}] SKIP (no map)")
            continue
        r_means.append(r)
        ts = time.time()
        data = parse_chrom(chrom, vp, pop_of)
        for pop in POPULATIONS:
            if pop not in data:
                continue
            arr, pos = data[pop]
            n_indiv = arr.shape[1] // 2
            span_mb = (pos[-1] - pos[0]) / 1e6
            for _ in range(N_PAIRS):
                i, j = rng.integers(0, n_indiv, size=2)
                while j == i:
                    j = rng.integers(0, n_indiv)
                a = arr[:, 2 * i + int(rng.integers(0, 2))]
                b = arr[:, 2 * j + int(rng.integers(0, 2))]
                for g in gaps:
                    sl = ibs_seglens_mb(a, b, pos, g, MAX_SNP_GAP_BP)
                    total_runs[g][pop] += int(sl.size)
                    sl = sl[sl > MIN_KEEP_MB]
                    if sl.size:
                        seglens[g][pop].append(sl.astype(np.float32))
                exposure[pop] += span_mb
        print(f"  [{chrom}] r={r:.3f} cM/Mb  "
              + " ".join(f"{p}:{data[p][0].shape[0] if p in data else 0}snp"
                         for p in POPULATIONS)
              + f"  ({time.time()-ts:.1f}s)")
        sys.stdout.flush()

    r_mean = sum(r_means) / len(r_means) if r_means else 1.0
    thr_pc = PI * (1.0 - GENOTYPING_ERROR) ** 1000 * (1.0 - ACMG_THRESHOLD) / (
        ACMG_THRESHOLD * (1.0 - PI))

    # empirical p_chance(L) = P(random L-window fully inside an IBS segment)
    #                       = sum_segments max(0, seg - L) / exposure
    # computed in O(n log n) via sorted suffix sums, separately for each gap.
    emp_pc = {g: {} for g in gaps}
    s_sorted = {g: {} for g in gaps}
    ana_pc = {}
    for pop in POPULATIONS:
        if exposure[pop] > 0:
            ana_pc[pop] = ANALYTIC_HBAR[pop] ** (L_GRID * r_mean / BLOCK_CM)
    for g in gaps:
        for pop in POPULATIONS:
            if not seglens[g][pop] or exposure[pop] == 0:
                continue
            s = np.sort(np.concatenate(seglens[g][pop]).astype(np.float64))
            s_sorted[g][pop] = s
            prefix = np.concatenate(([0.0], np.cumsum(s)))
            idx = np.searchsorted(s, L_GRID, side="right")
            sum_gt = prefix[-1] - prefix[idx]
            cnt_gt = s.size - idx
            emp_pc[g][pop] = (sum_gt - L_GRID * cnt_gt) / exposure[pop]

    base = gaps[0]
    bpc = emp_pc[base]

    with OUT_PC.open("w") as fh:
        fh.write(f"# emp = empirical p_chance at GAP_TOL={base}\n")
        fh.write("L_Mb\t" + "\t".join(
            f"{p}_emp\t{p}_analytic" for p in POPULATIONS if p in bpc) + "\n")
        for k, L in enumerate(L_GRID):
            cells = []
            for p in POPULATIONS:
                if p in bpc:
                    cells.append(f"{bpc[p][k]:.3e}\t{ana_pc[p][k]:.3e}")
            fh.write(f"{L:.2f}\t" + "\t".join(cells) + "\n")

    with OUT_SUMMARY.open("w") as fh:
        fh.write("# Empirical chance-IBS validation of the H-bar noise term\n")
        fh.write(f"# Chromosomes: {','.join(chroms)}\n")
        fh.write(f"# knobs: MAF_MIN={MAF_MIN} N_SUBSAMPLE={N_SUBSAMPLE} "
                 f"N_PAIRS={N_PAIRS} GAP_SWEEP={gaps} "
                 f"MAX_SNP_GAP_KB={MAX_SNP_GAP_KB} PI={PI}\n")
        fh.write(f"# mean cM/Mb={r_mean:.3f}  wall clock={time.time()-t0:.1f}s\n\n")

        fh.write("Minimum callable ROH length (Mb), posterior>="
                 f"{ACMG_THRESHOLD}, pi={PI}, by genotyping-error tolerance "
                 "(GAP_TOL) vs the analytic H-bar prediction:\n")
        fh.write("population\t" + "\t".join(f"gap={g}" for g in gaps)
                 + "\tanalytic_Hbar\n")
        for pop in POPULATIONS:
            if pop not in ana_pc:
                continue
            cells = [(f"{min_callable_from_pc(L_GRID, emp_pc[g][pop], PI):.2f}"
                      if pop in emp_pc[g] else "NA") for g in gaps]
            ana_L = min_callable_from_pc(L_GRID, ana_pc[pop], PI)
            fh.write(f"{pop}\t" + "\t".join(cells) + f"\t{ana_L:.2f}\n")

        fh.write(f"\nChance-IBS run-length distribution (Mb), GAP_TOL={base} "
                 "baseline\n")
        fh.write("(percentiles over ALL runs via total count; runs <= "
                 f"{MIN_KEEP_MB} Mb dropped from storage only):\n")
        fh.write("population\tn_runs_all\tp99\tp99.9\tp99.99\tmax\n")
        for pop in POPULATIONS:
            if pop not in s_sorted[base]:
                continue
            s = s_sorted[base][pop]
            tot = total_runs[base][pop]

            def pct(q):  # true q-th percentile over all `tot` runs
                rank = q / 100.0 * tot
                pos_in_kept = rank - (tot - s.size)  # dropped runs are the smallest
                if pos_in_kept <= 0:
                    return 0.0
                return float(s[min(s.size - 1, int(np.ceil(pos_in_kept)) - 1)])

            fh.write(f"{pop}\t{tot}\t{pct(99):.3f}\t{pct(99.9):.3f}\t"
                     f"{pct(99.99):.4f}\t{s.max():.3f}\n")

        fh.write(f"\np_chance(L) at tabulated lengths -- emp(gap={base}) / "
                 "analytic:\n")
        fh.write("L_Mb\t" + "\t".join(bpc) + "\n")
        for L in L_REPORT:
            k = int(np.argmin(np.abs(L_GRID - L)))
            cells = [f"{bpc[p][k]:.1e}/{ana_pc[p][k]:.1e}" for p in bpc]
            fh.write(f"{L}\t" + "\t".join(cells) + "\n")
        fh.write(f"\n(threshold p_chance for posterior>={ACMG_THRESHOLD}: {thr_pc:.2e})\n")

    print(f"\n  total wall clock: {time.time()-t0:.1f}s")
    print(f"  -> {OUT_PC}")
    print(f"  -> {OUT_SUMMARY}")
    print("  empirical min-callable ROH (Mb) by gap | analytic:")
    for pop in POPULATIONS:
        if pop not in ana_pc:
            continue
        cells = " ".join(
            f"g{g}={min_callable_from_pc(L_GRID, emp_pc[g][pop], PI):.2f}"
            for g in gaps if pop in emp_pc[g])
        ana_L = min_callable_from_pc(L_GRID, ana_pc[pop], PI)
        print(f"    {pop}: {cells} | analytic {ana_L:.2f}")


if __name__ == "__main__":
    main()
