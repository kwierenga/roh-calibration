"""
Run the per-locus IBD posterior framework on all 22 autosomes.

Inputs:
  - chr22_phased.vcf.gz (already on disk in this directory)
  - chr1..chr21 phased VCFs (in all_autosomes/ subdirectory after 10_download)
  - Palsson 2024 deCODE maps in external/

Output:
  - all_autosomes_posterior_lookup.tsv  (per-window posterior across all 22)
  - all_autosomes_summary.txt           (per-chromosome summary + ACMG-threshold scrutiny)

Re-uses the per-chromosome logic from 09_posterior_prototype.py with chrom as
a parameter. Fully scripted; no prompts.
"""

import gzip
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"

OUT_LOOKUP = HERE / "all_autosomes_posterior_lookup.tsv"
OUT_SUMMARY = HERE / "all_autosomes_summary.txt"

POPULATION = "EUR"
PRIOR_PI = 0.0625
GENOTYPING_ERROR = 0.001
BLOCK_CM = 0.5
TRACT_LENGTHS_MB = [1, 2, 3, 5, 7, 10, 15]
THRESHOLDS = [0.50, 0.90, 0.95, 0.99]
WINDOW_BP = 1_000_000


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


def stream_chr_maf(vcf_path, population=POPULATION):
    af_key = f"AF_{population}="
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
            info = fields[7]
            af = None
            for kv in info.split(";"):
                if kv.startswith(af_key):
                    try:
                        af = float(kv[len(af_key):])
                    except ValueError:
                        af = None
                    break
            if af is None or af <= 0.0 or af >= 1.0:
                continue
            yield pos, 2.0 * af * (1.0 - af)


def posterior(L_mb, r_cmpermb, mean_2pq, pi=PRIOR_PI, eps=GENOTYPING_ERROR,
              block_cm=BLOCK_CM, n_snp_in_tract=1000):
    n_eff = max(1.0, (L_mb * r_cmpermb) / block_cm)
    p_hom = max(1e-12, 1.0 - mean_2pq)
    p_chance = p_hom ** n_eff
    p_ibd_data = (1.0 - eps) ** n_snp_in_tract
    num = pi * p_ibd_data
    denom = num + (1.0 - pi) * p_chance
    return num / denom if denom > 0 else 1.0


def length_for_posterior(target, r, pq2):
    lo, hi = 0.1, 50.0
    if posterior(hi, r, pq2) < target:
        return None
    if posterior(lo, r, pq2) >= target:
        return lo
    for _ in range(60):
        mid = (lo + hi) / 2
        if posterior(mid, r, pq2) >= target:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def vcf_path_for(chrom):
    if chrom == "chr22":
        return HERE / "chr22_phased.vcf.gz"
    else:
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

    windows = {}
    for pos, pq2 in stream_chr_maf(vp, POPULATION):
        w = (pos // WINDOW_BP) * WINDOW_BP
        if w not in windows:
            windows[w] = [0.0, 0]
        windows[w][0] += pq2
        windows[w][1] += 1

    rows = []
    for w_start in sorted(windows):
        sum_pq, n_var = windows[w_start]
        if n_var == 0:
            continue
        mean_2pq = sum_pq / n_var
        r = r_at(w_start + WINDOW_BP // 2)
        post_by_L = {L: posterior(L, r, mean_2pq) for L in TRACT_LENGTHS_MB}
        length_at_thr = {thr: length_for_posterior(thr, r, mean_2pq) for thr in THRESHOLDS}
        rows.append({
            "chrom": chrom, "w_start": w_start, "w_end": w_start + WINDOW_BP,
            "n_variants": n_var, "mean_2pq": mean_2pq, "cMperMb": r,
            "posteriors_by_L": post_by_L, "length_at_thr": length_at_thr,
        })
    return chrom, rows, None


def main():
    t0 = time.time()
    all_rows = []
    per_chrom_stats = []
    chroms = [f"chr{n}" for n in range(1, 23)]

    for chrom in chroms:
        ts = time.time()
        try:
            c, rows, err = process_chromosome(chrom)
        except Exception as e:
            print(f"  [{chrom}] CRASHED ({type(e).__name__}: {e})")
            per_chrom_stats.append((chrom, 0, 0, 0, f"crash: {e}"))
            sys.stdout.flush()
            continue
        if err:
            print(f"  [{chrom}] SKIPPED ({err})")
            per_chrom_stats.append((chrom, 0, 0, 0, err))
            sys.stdout.flush()
            continue
        n_w = len(rows)
        n_v = sum(r["n_variants"] for r in rows)
        n_acmg_calibrated = sum(1 for r in rows if r["cMperMb"] > 0 and r["posteriors_by_L"][10] >= 0.95)
        n_w_with_recomb = sum(1 for r in rows if r["cMperMb"] > 0)
        print(f"  [{chrom}] {n_w} windows, {n_v:,} variants, "
              f"{n_acmg_calibrated}/{n_w_with_recomb} windows have P(IBD|10Mb)>=0.95 "
              f"({time.time() - ts:.1f}s)")
        sys.stdout.flush()
        per_chrom_stats.append((chrom, n_w, n_v, n_acmg_calibrated, n_w_with_recomb))
        all_rows.extend(rows)

    # Write big lookup
    with OUT_LOOKUP.open("w") as fh:
        cols = ["chrom", "window_start", "window_end", "n_variants", "mean_2pq", "cMperMb"] + \
               [f"post_L{L}Mb" for L in TRACT_LENGTHS_MB] + \
               [f"L_for_post_{thr}_Mb" for thr in THRESHOLDS]
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            vals = [r["chrom"], str(r["w_start"]), str(r["w_end"]),
                    str(r["n_variants"]), f"{r['mean_2pq']:.5f}", f"{r['cMperMb']:.4f}"]
            for L in TRACT_LENGTHS_MB:
                vals.append(f"{r['posteriors_by_L'][L]:.5f}")
            for thr in THRESHOLDS:
                L = r["length_at_thr"][thr]
                vals.append("NA" if L is None else f"{L:.2f}")
            fh.write("\t".join(vals) + "\n")

    # Summary
    with OUT_SUMMARY.open("w") as fh:
        fh.write(f"# All-autosomes per-locus IBD posterior framework\n")
        fh.write(f"# Population AF: {POPULATION}; Prior pi: {PRIOR_PI}; Block cM: {BLOCK_CM}\n")
        fh.write(f"# Total wall clock: {time.time() - t0:.1f}s\n\n")
        fh.write("chrom\tn_windows\tn_variants\tn_windows_acmg_calibrated\tn_windows_with_recomb\tfraction_calibrated\n")
        total_calibrated = 0
        total_with_recomb = 0
        for c, nw, nv, nc, nr_or_err in per_chrom_stats:
            if isinstance(nr_or_err, str):
                fh.write(f"{c}\tSKIPPED\t-\t-\t-\t-\n")
                continue
            nr = nr_or_err
            frac = nc / nr if nr else 0.0
            fh.write(f"{c}\t{nw}\t{nv}\t{nc}\t{nr}\t{frac:.3f}\n")
            total_calibrated += nc
            total_with_recomb += nr
        if total_with_recomb:
            fh.write(f"TOTAL\t-\t-\t{total_calibrated}\t{total_with_recomb}\t{total_calibrated/total_with_recomb:.3f}\n")

    print()
    print(f"  total wall clock: {time.time() - t0:.1f}s")
    print(f"  -> {OUT_LOOKUP}")
    print(f"  -> {OUT_SUMMARY}")
    if total_with_recomb:
        print(f"  HEADLINE: ACMG 10 Mb threshold delivers >=0.95 posterior at "
              f"{total_calibrated}/{total_with_recomb} = "
              f"{100*total_calibrated/total_with_recomb:.1f}% of autosomal windows")


if __name__ == "__main__":
    main()
