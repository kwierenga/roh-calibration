"""
LD-aware haplotype-IBS noise term for the per-locus ROH posterior.

WHY: 15_cross_population.py's noise term is mean per-site 2pq, which tracks the
allele-frequency SPECTRUM rather than haplotype diversity. It inverted the
population axis (AFR appeared LEAST calibrated despite its shorter LD / higher
haplotype diversity). This module replaces that term with an empirical,
LD-aware quantity:

    Hbar(window, pop) = mean over 0.5 cM sub-blocks of the
                        sample-size-UNBIASED haplotype homozygosity
                        = P(two random haplotypes are identical-by-state
                          across a 0.5 cM block) in unrelated individuals.

This is exactly the per-block "chance match" probability the posterior needs,
and it plugs into the existing geometry unchanged: n_eff = L*cMperMb/0.5cM is
already a count of 0.5 cM blocks, so p_chance = Hbar ** n_eff. Lower haplotype
diversity (founder/bottleneck) -> higher Hbar -> more background homozygosity;
higher diversity (AFR) -> lower Hbar. This should restore the expected
population ordering.

METHODS KNOBS (for the biostatistics partner to scrutinize) -- all at top:
  - MAF_MIN          common-SNP threshold that defines the haplotypes
  - BLOCK_CM         genetic sub-block size (must match n_eff's block in 15)
  - MAX_SNPS_BLOCK   SNPs per block are thinned to this, so Hbar measures
                     diversity at a controlled resolution rather than exact
                     match over hundreds of SNPs (which collapses to ~0)
  - N_SUBSAMPLE      individuals per population (unbiased estimator is valid
                     at any N; subsampling bounds the genotype-parsing cost)
  - H_FLOOR          floor on Hbar so p_chance stays finite
  - unbiased H = sum_h c_h(c_h-1) / (N(N-1))   [Nei's homozygosity, unbiased]

Inputs (all on disk): per-chr full 3,202-sample phased VCFs; deCODE maps;
samples_2504_pop.panel (unrelated subset -> superpopulation).

Outputs:
  cross_pop_hap_diversity.tsv   per-window per-pop Hbar + diagnostics
  cross_pop_hap_summary.txt     Hbar-based calibration matrix, side-by-side
                                with the 2pq-based one (does AFR flip?)

Usage:
  python 16_haplotype_ibs_noise.py chr22   # smoke test
  python 16_haplotype_ibs_noise.py         # all autosomes (background)
"""

import gzip
import math
import random
import statistics
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"
PANEL = HERE / "samples_2504_pop.panel"
MASTER_2PQ = HERE / "cross_pop_master_lookup.tsv"

OUT_DIV = HERE / "cross_pop_hap_diversity.tsv"
OUT_SUMMARY = HERE / "cross_pop_hap_summary.txt"

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
PRIORS = [(0.0156, "2nd_cousin"), (0.0625, "1st_cousin"),
          (0.125, "avuncular_or_double_1c"), (0.25, "incest_or_sibling")]
DEFAULT_PI = 0.0625

# ---- methods knobs ----
MAF_MIN = 0.05
BLOCK_CM = 0.5
MAX_SNPS_BLOCK = 25
N_SUBSAMPLE = 200
H_FLOOR = 1e-4
RNG_SEED = 17
# -----------------------

GENOTYPING_ERROR = 0.001
WINDOW_BP = 1_000_000
CONV_LENGTH = 10        # conventional clinical-lab comparison length (Mb). NB the
                        # ACMG-2021 standard is >3-5 Mb, not 10 Mb; 10 Mb is the
                        # common lab-practice operating point, used here only as a
                        # fixed-length reference for the calibration-fraction test.
POST_THRESHOLD = 0.95   # posterior decision threshold
# Length sweep: at 10 Mb the H-bar posterior saturates to ~1.0 in every pop
# (a 10 Mb run is essentially never chance-IBS), so the population signal that
# now lives in H-bar is invisible there. Sweep shorter lengths to expose it.
SWEEP_LENGTHS_MB = [1.5, 2.0, 3.0, 5.0, 10.0]
AF_PREFIX = {pop: f"AF_{pop}=" for pop in POPULATIONS}


def load_decode_map(path, chrom):
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Chr"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            out[int(f[1])] = float(f[3])
    return out


def build_cm(pat, mat, max_pos):
    """Sex-averaged cM/Mb per 1 Mb window + cumulative cM at window starts."""
    nwin = max_pos // WINDOW_BP + 2
    rate = [0.0] * nwin

    def r_at(center):
        for c in (center, center - WINDOW_BP, center + WINDOW_BP):
            if c in pat:
                return 0.5 * (pat[c] + mat.get(c, pat[c]))
        return 0.0

    for w in range(nwin):
        rate[w] = r_at(w * WINDOW_BP + 500_000)
    cum = [0.0] * (nwin + 1)
    for w in range(nwin):
        cum[w + 1] = cum[w] + rate[w]  # each window is 1 Mb wide
    return rate, cum


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
    """Return {pop: [blob-column indices]} for a subsample of unrelated samples."""
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


def unbiased_homozygosity(hap_strings):
    n = len(hap_strings)
    if n < 2:
        return None
    counts = Counter(hap_strings)
    num = sum(c * (c - 1) for c in counts.values())
    return num / (n * (n - 1)), len(counts)


def finalize_block(cur, sel):
    """cur[pop] = list of per-SNP allele strings (len 2*Nsub each).
    Returns {pop: (H, n_distinct, n_hap)} after thinning to MAX_SNPS_BLOCK."""
    res = {}
    for pop in POPULATIONS:
        snps = cur[pop]
        if not snps:
            continue
        if len(snps) > MAX_SNPS_BLOCK:
            stride = len(snps) / MAX_SNPS_BLOCK
            snps = [snps[int(i * stride)] for i in range(MAX_SNPS_BLOCK)]
        nhap = len(snps[0])
        haps = ["".join(snp[j] for snp in snps) for j in range(nhap)]
        hh = unbiased_homozygosity(haps)
        if hh is not None:
            res[pop] = (hh[0], hh[1], nhap)
    return res


def process_chromosome(chrom, vcf_path, pop_of):
    pat = load_decode_map(PAT_MAP, chrom)
    mat = load_decode_map(MAT_MAP, chrom)
    if not pat:
        return None
    max_pos = max(pat) + WINDOW_BP
    rate, cum = build_cm(pat, mat, max_pos)
    sel = select_columns(vcf_path, pop_of)

    def cum_cm(pos):
        w = pos // WINDOW_BP
        if w + 1 >= len(cum):
            return cum[-1]
        return cum[w] + rate[w] * (pos - w * WINDOW_BP) / WINDOW_BP

    # aggregation: agg[pop][w_start] = [sumH, cnt, sum_distinct, sum_nhap]
    agg = {pop: {} for pop in POPULATIONS}
    cur = {pop: [] for pop in POPULATIONS}
    cur_block = None
    cur_block_pos = None

    def flush():
        if cur_block is None:
            return
        w = (cur_block_pos // WINDOW_BP) * WINDOW_BP
        for pop, (H, ndist, nhap) in finalize_block(cur, sel).items():
            a = agg[pop].setdefault(w, [0.0, 0, 0, 0])
            a[0] += H
            a[1] += 1
            a[2] += ndist
            a[3] += nhap
        for pop in POPULATIONS:
            cur[pop].clear()

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
            pos = int(fields[1])
            blk = int(cum_cm(pos) // BLOCK_CM)
            if blk != cur_block:
                flush()
                cur_block = blk
                cur_block_pos = pos
            gtblob = fields[9].split("\t")
            for pop in common:
                cols = sel[pop]
                cur[pop].append("".join(gtblob[c][0] + gtblob[c][2] for c in cols))
    flush()

    rows = []
    for pop in POPULATIONS:
        for w, (sH, cnt, sdist, snhap) in agg[pop].items():
            rows.append((chrom, w, pop, sH / cnt, cnt, sdist / cnt, snhap / cnt))
    return rows, rate


def posterior_hap(L_mb, r, p_block, pi):
    n_eff = max(1.0, (L_mb * r) / BLOCK_CM)
    p_chance = max(p_block, H_FLOOR) ** n_eff
    p_ibd = (1.0 - GENOTYPING_ERROR) ** 1000
    num = pi * p_ibd
    return num / (num + (1.0 - pi) * p_chance)


def min_callable_length(r, p_block, pi, threshold=POST_THRESHOLD):
    """Smallest ROH length (Mb) whose H-bar posterior reaches `threshold`.
    Analytic inverse of posterior_hap: solve b**(L*r/BLOCK_CM) <= RHS for L,
    with b = max(p_block, H_FLOOR), RHS = pi*p_ibd*(1-T)/(T*(1-pi))."""
    b = max(p_block, H_FLOOR)
    if r <= 0 or b >= 1.0:
        return float("inf")
    p_ibd = (1.0 - GENOTYPING_ERROR) ** 1000
    num = pi * p_ibd
    rhs = num * (1.0 - threshold) / (threshold * (1.0 - pi))
    if rhs <= 0:
        return float("inf")
    n_eff_req = math.log(rhs) / math.log(b)        # both logs < 0 -> positive
    return max(n_eff_req * BLOCK_CM / r, BLOCK_CM / r)  # respect n_eff >= 1


def main():
    t0 = time.time()
    pop_of = load_panel()
    argv = sys.argv[1:]
    chroms = argv if argv else [f"chr{n}" for n in range(1, 23)]

    # Hbar[(chrom,w,pop)] = (Hbar, nblk, mean_distinct, mean_nhap)
    Hbar = {}
    rate_by_chrom = {}
    for chrom in chroms:
        vp = HERE / "chr22_phased.vcf.gz" if chrom == "chr22" else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"
        if not vp.exists():
            print(f"  [{chrom}] SKIP (no VCF)")
            continue
        ts = time.time()
        out = process_chromosome(chrom, vp, pop_of)
        if out is None:
            print(f"  [{chrom}] SKIP (no map)")
            continue
        rows, rate = out
        rate_by_chrom[chrom] = rate
        for (c, w, pop, H, nblk, mdist, mnhap) in rows:
            Hbar[(c, w, pop)] = (H, nblk, mdist, mnhap)

        def qcal(pop):
            nc = nr = 0
            for (c, w, p), (H, *_ ) in Hbar.items():
                if c != chrom or p != pop:
                    continue
                r = rate[w // WINDOW_BP]
                if r <= 0:
                    continue
                nr += 1
                if posterior_hap(CONV_LENGTH, r, H, DEFAULT_PI) >= POST_THRESHOLD:
                    nc += 1
            return nc, nr

        ec, er = qcal("EUR")
        ac, ar = qcal("AFR")
        mh_e = sum(v[0] for k, v in Hbar.items() if k[0] == chrom and k[2] == "EUR")
        mh_a = sum(v[0] for k, v in Hbar.items() if k[0] == chrom and k[2] == "AFR")
        ne = sum(1 for k in Hbar if k[0] == chrom and k[2] == "EUR")
        na = sum(1 for k in Hbar if k[0] == chrom and k[2] == "AFR")
        print(f"  [{chrom}] EUR Hbar~{mh_e/ne if ne else 0:.4f} cal {ec}/{er} "
              f"| AFR Hbar~{mh_a/na if na else 0:.4f} cal {ac}/{ar} "
              f"({time.time()-ts:.1f}s)")
        sys.stdout.flush()

    # diversity table
    with OUT_DIV.open("w") as fh:
        fh.write("chrom\twindow_start\tpopulation\tHbar\tn_blocks\tmean_distinct_hap\tmean_n_hap\tcMperMb\n")
        for (c, w, pop), (H, nblk, mdist, mnhap) in sorted(Hbar.items()):
            r = rate_by_chrom[c][w // WINDOW_BP]
            fh.write(f"{c}\t{w}\t{pop}\t{H:.6f}\t{nblk}\t{mdist:.1f}\t{mnhap:.0f}\t{r:.4f}\n")

    # 2pq-based calibration for the SAME windows, for side-by-side comparison
    old_m2pq = {}
    if MASTER_2PQ.exists():
        with MASTER_2PQ.open() as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {n: i for i, n in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                key = (f[ix["chrom"]], int(f[ix["window_start"]]))
                for pop in POPULATIONS:
                    old_m2pq[(key[0], key[1], pop)] = float(f[ix[f"mean_2pq_{pop}"]])

    def cal_fraction(pi, use_hap, length=CONV_LENGTH):
        out = {}
        for pop in POPULATIONS:
            nc = nr = 0
            for (c, w, p), (H, *_ ) in Hbar.items():
                if p != pop:
                    continue
                r = rate_by_chrom[c][w // WINDOW_BP]
                if r <= 0:
                    continue
                nr += 1
                pblock = H if use_hap else (1.0 - old_m2pq.get((c, w, pop), 0.0))
                if posterior_hap(length, r, pblock, pi) >= POST_THRESHOLD:
                    nc += 1
            out[pop] = (nc, nr)
        return out

    with OUT_SUMMARY.open("w") as fh:
        fh.write("# Haplotype-IBS (LD-aware) noise term vs per-site 2pq\n")
        fh.write(f"# Chromosomes: {','.join(chroms)}\n")
        fh.write(f"# knobs: MAF_MIN={MAF_MIN} BLOCK_CM={BLOCK_CM} MAX_SNPS_BLOCK={MAX_SNPS_BLOCK} "
                 f"N_SUBSAMPLE={N_SUBSAMPLE} H_FLOOR={H_FLOOR}\n")
        fh.write(f"# wall clock: {time.time()-t0:.1f}s\n\n")

        gw = {}
        for pop in POPULATIONS:
            vals = [v[0] for k, v in Hbar.items() if k[2] == pop]
            gw[pop] = sum(vals) / len(vals) if vals else 0.0
        fh.write("Genome-wide mean Hbar (P[2 haplotypes IBS over 0.5 cM]) by population:\n")
        for pop in sorted(POPULATIONS, key=lambda p: -gw[p]):
            fh.write(f"  {pop}\t{gw[pop]:.4f}\n")
        fh.write("  (expect founder/bottleneck HIGH, AFR LOW)\n\n")

        fh.write("Conventional 10 Mb >= 0.95 calibration fraction -- HAPLOTYPE-IBS noise:\n"
                 "(10 Mb = common lab-practice point; ACMG-2021 standard is >3-5 Mb)\n")
        fh.write("population\t" + "\t".join(f"pi={pi}" for pi, _ in PRIORS) + "\n")
        for pop in POPULATIONS:
            cells = []
            for pi, _ in PRIORS:
                nc, nr = cal_fraction(pi, True)[pop]
                cells.append(f"{nc/nr:.3f}" if nr else "NA")
            fh.write(pop + "\t" + "\t".join(cells) + "\n")

        fh.write("\nCalibration fraction vs ROH length (HAPLOTYPE-IBS, "
                 f"pi={DEFAULT_PI}) -- the population signal that 10 Mb hides:\n")
        sweep = {L: cal_fraction(DEFAULT_PI, True, length=L) for L in SWEEP_LENGTHS_MB}
        fh.write("population\t" + "\t".join(f"{L}Mb" for L in SWEEP_LENGTHS_MB) + "\n")
        for pop in POPULATIONS:
            cells = [(f"{sweep[L][pop][0]/sweep[L][pop][1]:.3f}"
                      if sweep[L][pop][1] else "NA") for L in SWEEP_LENGTHS_MB]
            fh.write(pop + "\t" + "\t".join(cells) + "\n")

        fh.write(f"\nMedian minimum callable ROH length (Mb), posterior>="
                 f"{POST_THRESHOLD}, pi={DEFAULT_PI} -- CLINICAL HEADLINE\n"
                 "(higher H-bar / lower diversity -> needs a longer run):\n")
        med_min = {}
        for pop in POPULATIONS:
            Ls = []
            for (c, w, p), (H, *_) in Hbar.items():
                if p != pop:
                    continue
                r = rate_by_chrom[c][w // WINDOW_BP]
                if r <= 0:
                    continue
                Ls.append(min_callable_length(r, H, DEFAULT_PI))
            med_min[pop] = statistics.median(Ls) if Ls else float("nan")
        fh.write("population\tmedian_min_callable_Mb\n")
        for pop in sorted(POPULATIONS, key=lambda p: med_min[p]):
            fh.write(f"{pop}\t{med_min[pop]:.2f}\n")

        fh.write("\n[RETIRED] per-site 2pq noise (script 15) -- kept only as a "
                 "foil; note the inverted AFR ordering that motivated retirement:\n")
        fh.write("population\t" + "\t".join(f"pi={pi}" for pi, _ in PRIORS) + "\n")
        for pop in POPULATIONS:
            cells = []
            for pi, _ in PRIORS:
                nc, nr = cal_fraction(pi, False)[pop]
                cells.append(f"{nc/nr:.3f}" if nr else "NA")
            fh.write(pop + "\t" + "\t".join(cells) + "\n")

    print()
    print(f"  total wall clock: {time.time()-t0:.1f}s")
    print(f"  -> {OUT_DIV}")
    print(f"  -> {OUT_SUMMARY}")
    print("  genome-wide mean Hbar (high=more background homozygosity):")
    for pop in sorted(POPULATIONS, key=lambda p: -gw[p]):
        print(f"    {pop}: {gw[pop]:.4f}")


if __name__ == "__main__":
    main()
