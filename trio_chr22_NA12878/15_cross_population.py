"""
================================ RETIRED (2026-05-23) ========================
SUPERSEDED by 16_haplotype_ibs_noise.py. The per-site 2pq noise term used here
is dimensionally wrong -- it raises (1 - mean_2pq), a PER-SITE homozygosity, to
n_eff, a count of 0.5 cM BLOCKS -- and it inverts the population axis (AFR comes
out LEAST calibrated despite its high haplotype diversity / short LD). The
adopted noise term is the LD-aware haplotype-IBS quantity H-bar (script 16).
This script is kept only for provenance and as the documented foil reproduced in
cross_pop_hap_summary.txt. Do NOT use its outputs as the cross-population result.
=============================================================================

Cross-population per-locus IBD posterior framework, all 22 autosomes.

Extends 11_run_all_autosomes.py from a single population (EUR) to all five
1000G superpopulations (EUR, AFR, EAS, SAS, AMR) in a SINGLE parse pass:
all five AF_{POP} keys live on the same VCF line, so the per-population
diversity term (mean_2pq) costs no extra I/O over the original EUR run.

The recombination input (deCODE / Palsson 2024 sex-averaged cM/Mb), the
1 Mb windowing, and the Bayesian posterior are identical to script 11/14 ---
only the allele-frequency source is iterated over populations and priors.

Inputs (all already on disk):
  - chr22_phased.vcf.gz + all_autosomes/chr1..21_phased.vcf.gz  (3,202 samples,
    per-population AF in INFO)
  - Palsson 2024 deCODE maps in external/

Outputs:
  - cross_pop_master_lookup.tsv          per-window cM/Mb + per-pop n & mean_2pq
  - cross_pop_{POP}_pi_{pi}_{name}.tsv    20 lookup tables (5 pop x 4 prior)
  - cross_pop_summary.txt                 calibration-fraction matrix (the headline)

Usage:
  python 15_cross_population.py            # all 22 autosomes (~1 h, run in bg)
  python 15_cross_population.py chr22      # single-chromosome smoke test
"""

import gzip
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"

OUT_MASTER = HERE / "cross_pop_master_lookup.tsv"
OUT_SUMMARY = HERE / "cross_pop_summary.txt"

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
PRIORS = [
    (0.0156, "2nd_cousin"),
    (0.0625, "1st_cousin"),
    (0.125,  "avuncular_or_double_1c"),
    (0.25,   "incest_or_sibling"),
]
DEFAULT_PI = 0.0625  # 1st-cousin offspring; the headline prior

GENOTYPING_ERROR = 0.001
BLOCK_CM = 0.5
TRACT_LENGTHS_MB = [1, 2, 3, 5, 7, 10, 15]
THRESHOLDS = [0.50, 0.90, 0.95, 0.99]
WINDOW_BP = 1_000_000
ACMG_LENGTH = 10
ACMG_THRESHOLD = 0.95

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


def stream_chr_af_multi(vcf_path):
    """Yield (pos, {pop: af}) for biallelic SNPs, all populations in one pass."""
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            afs = {}
            for kv in fields[7].split(";"):
                if kv[:3] != "AF_":
                    continue
                for pop, pre in AF_PREFIX.items():
                    if kv.startswith(pre):
                        try:
                            afs[pop] = float(kv[len(pre):])
                        except ValueError:
                            pass
                        break
            yield pos, afs


def posterior(L_mb, r_cmpermb, mean_2pq, pi, eps=GENOTYPING_ERROR,
              block_cm=BLOCK_CM, n_snp_in_tract=1000):
    n_eff = max(1.0, (L_mb * r_cmpermb) / block_cm)
    p_hom = max(1e-12, 1.0 - mean_2pq)
    p_chance = p_hom ** n_eff
    p_ibd_data = (1.0 - eps) ** n_snp_in_tract
    num = pi * p_ibd_data
    denom = num + (1.0 - pi) * p_chance
    return num / denom if denom > 0 else 1.0


def length_for_posterior(target, r, pq2, pi):
    lo, hi = 0.1, 50.0
    if posterior(hi, r, pq2, pi) < target:
        return None
    if posterior(lo, r, pq2, pi) >= target:
        return lo
    for _ in range(60):
        mid = (lo + hi) / 2
        if posterior(mid, r, pq2, pi) >= target:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def vcf_path_for(chrom):
    if chrom == "chr22":
        return HERE / "chr22_phased.vcf.gz"
    return DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"


def process_chromosome(chrom):
    vp = vcf_path_for(chrom)
    if not vp.exists():
        return chrom, None, f"VCF not found: {vp}"
    pat = load_decode_map(PAT_MAP, chrom)
    mat = load_decode_map(MAT_MAP, chrom)
    if not pat:
        return chrom, None, f"no deCODE entries for {chrom}"

    def r_at(pos):
        idx = pos // WINDOW_BP
        center = idx * WINDOW_BP + 500_000
        for c in [center, center - WINDOW_BP, center + WINDOW_BP]:
            if c in pat:
                return 0.5 * (pat[c] + mat.get(c, pat[c]))
        return 0.0

    # windows[w_start][pop] = [sum_2pq, n_var]
    windows = {}
    for pos, afs in stream_chr_af_multi(vp):
        w = (pos // WINDOW_BP) * WINDOW_BP
        slot = windows.get(w)
        if slot is None:
            slot = {pop: [0.0, 0] for pop in POPULATIONS}
            windows[w] = slot
        for pop, af in afs.items():
            if 0.0 < af < 1.0:
                slot[pop][0] += 2.0 * af * (1.0 - af)
                slot[pop][1] += 1

    rows = []
    for w_start in sorted(windows):
        r = r_at(w_start + WINDOW_BP // 2)
        per_pop = {}
        any_var = False
        for pop in POPULATIONS:
            s, n = windows[w_start][pop]
            if n == 0:
                per_pop[pop] = (0, 0.0)
                continue
            any_var = True
            per_pop[pop] = (n, s / n)
        if not any_var:
            continue
        rows.append({
            "chrom": chrom, "w_start": w_start, "w_end": w_start + WINDOW_BP,
            "cMperMb": r, "per_pop": per_pop,
        })
    return chrom, rows, None


def write_per_pop_prior_tables(all_rows):
    """20 lookup tables (pop x prior), schema matching all_autosomes_pi_*.tsv."""
    # stats[(pop, pi_name)] = [n_calibrated, n_with_recomb]
    stats = {}
    for pop in POPULATIONS:
        for pi, name in PRIORS:
            out_path = HERE / f"cross_pop_{pop}_pi_{pi:.4f}_{name}.tsv"
            n_cal = 0
            n_rec = 0
            with out_path.open("w") as fh:
                cols = ["chrom", "window_start", "window_end", "n_variants",
                        "mean_2pq", "cMperMb"] + \
                       [f"post_L{L}Mb" for L in TRACT_LENGTHS_MB] + \
                       [f"L_for_post_{thr}_Mb" for thr in THRESHOLDS]
                fh.write("\t".join(cols) + "\n")
                for r in all_rows:
                    n_var, m2pq = r["per_pop"][pop]
                    if n_var == 0:
                        continue
                    rr = r["cMperMb"]
                    if rr > 0:
                        n_rec += 1
                    vals = [r["chrom"], str(r["w_start"]), str(r["w_end"]),
                            str(n_var), f"{m2pq:.5f}", f"{rr:.4f}"]
                    for L in TRACT_LENGTHS_MB:
                        p = posterior(L, rr, m2pq, pi)
                        vals.append(f"{p:.5f}")
                        if L == ACMG_LENGTH and rr > 0 and p >= ACMG_THRESHOLD:
                            n_cal += 1
                    for thr in THRESHOLDS:
                        L = length_for_posterior(thr, rr, m2pq, pi)
                        vals.append("NA" if L is None else f"{L:.2f}")
                    fh.write("\t".join(vals) + "\n")
            stats[(pop, name)] = [n_cal, n_rec, pi]
    return stats


def main():
    t0 = time.time()
    argv_chroms = sys.argv[1:]
    chroms = argv_chroms if argv_chroms else [f"chr{n}" for n in range(1, 23)]

    all_rows = []
    for chrom in chroms:
        ts = time.time()
        try:
            c, rows, err = process_chromosome(chrom)
        except Exception as e:
            print(f"  [{chrom}] CRASHED ({type(e).__name__}: {e})")
            sys.stdout.flush()
            continue
        if err:
            print(f"  [{chrom}] SKIPPED ({err})")
            sys.stdout.flush()
            continue
        # quick per-chrom readout: AFR vs EUR calibration at default prior
        def cal_frac(pop):
            nc = nr = 0
            for r in rows:
                n_var, m2pq = r["per_pop"][pop]
                if n_var == 0 or r["cMperMb"] <= 0:
                    continue
                nr += 1
                if posterior(ACMG_LENGTH, r["cMperMb"], m2pq, DEFAULT_PI) >= ACMG_THRESHOLD:
                    nc += 1
            return nc, nr
        ec, er = cal_frac("EUR")
        ac, ar = cal_frac("AFR")
        print(f"  [{chrom}] {len(rows)} windows "
              f"| EUR cal {ec}/{er} ({100*ec/er if er else 0:.0f}%) "
              f"| AFR cal {ac}/{ar} ({100*ac/ar if ar else 0:.0f}%) "
              f"({time.time()-ts:.1f}s)")
        sys.stdout.flush()
        all_rows.extend(rows)

    if not all_rows:
        sys.exit("no rows produced")

    # master per-window table
    with OUT_MASTER.open("w") as fh:
        cols = ["chrom", "window_start", "window_end", "cMperMb"]
        for pop in POPULATIONS:
            cols += [f"n_{pop}", f"mean_2pq_{pop}"]
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            vals = [r["chrom"], str(r["w_start"]), str(r["w_end"]),
                    f"{r['cMperMb']:.4f}"]
            for pop in POPULATIONS:
                n_var, m2pq = r["per_pop"][pop]
                vals += [str(n_var), f"{m2pq:.5f}"]
            fh.write("\t".join(vals) + "\n")

    stats = write_per_pop_prior_tables(all_rows)

    # genome-wide mean diversity per population (variant-weighted)
    div = {}
    for pop in POPULATIONS:
        s = n = 0.0
        for r in all_rows:
            nv, m2pq = r["per_pop"][pop]
            s += m2pq * nv
            n += nv
        div[pop] = s / n if n else 0.0

    with OUT_SUMMARY.open("w") as fh:
        fh.write("# Cross-population per-locus IBD posterior framework\n")
        fh.write(f"# Chromosomes: {','.join(chroms)}\n")
        fh.write(f"# deCODE map: Palsson 2024 (sex-averaged); window {WINDOW_BP//10**6} Mb; "
                 f"block {BLOCK_CM} cM; eps {GENOTYPING_ERROR}\n")
        fh.write(f"# ACMG test: posterior(IBD | {ACMG_LENGTH} Mb) >= {ACMG_THRESHOLD}\n")
        fh.write(f"# Total wall clock: {time.time()-t0:.1f}s\n\n")

        fh.write("Genome-wide mean diversity (variant-weighted mean 2pq) by population:\n")
        for pop in POPULATIONS:
            fh.write(f"  {pop}\t{div[pop]:.4f}\n")
        fh.write("\n")

        fh.write("=" * 78 + "\n")
        fh.write("HEADLINE: fraction of autosomal windows where ACMG 10 Mb >= 0.95 posterior\n")
        fh.write("(rows = population, cols = prior pi / relationship)\n")
        fh.write("=" * 78 + "\n\n")
        header = ["population"] + [f"pi={pi}({name})" for pi, name in PRIORS]
        fh.write("\t".join(header) + "\n")
        for pop in POPULATIONS:
            line = [pop]
            for pi, name in PRIORS:
                nc, nr, _ = stats[(pop, name)]
                line.append(f"{nc}/{nr}={nc/nr:.3f}" if nr else "NA")
            fh.write("\t".join(line) + "\n")

    print()
    print(f"  total wall clock: {time.time()-t0:.1f}s")
    print(f"  -> {OUT_MASTER}")
    print(f"  -> {OUT_SUMMARY}")
    print(f"  -> 20 per-(pop x prior) lookup tables: cross_pop_<POP>_pi_<pi>_<name>.tsv")
    print()
    print(f"  diversity (mean 2pq): " +
          ", ".join(f"{pop}={div[pop]:.3f}" for pop in POPULATIONS))
    print(f"  ACMG-calibrated fraction @ pi={DEFAULT_PI} (1st-cousin):")
    for pop in POPULATIONS:
        nc, nr, _ = stats[(pop, "1st_cousin")]
        print(f"    {pop}: {nc}/{nr} = {100*nc/nr if nr else 0:.1f}%")


if __name__ == "__main__":
    main()
