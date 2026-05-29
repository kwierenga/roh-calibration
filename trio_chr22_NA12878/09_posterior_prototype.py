"""
chr22 prototype of the per-locus IBD-vs-chance posterior framework.

This is the project's first recognizable deliverable artifact: a function
that takes (chromosome, position, tract length L, population, prior pi)
and returns P(IBD | observation). Implemented end-to-end on chr22 using
data already on disk:

  Inputs (already downloaded):
    - 1000G high-coverage 20220422 phased VCF for chr22 (with per-variant
      AF_EUR, AF_AFR, AF_EAS, AF_SAS, AF_AMR in INFO).
    - Palsson 2024 deCODE recombination maps (maps.pat.tsv + maps.mat.tsv).

  Math (per planning doc section 2 + methodology log entry):
    P(L | IBD) = (1 - epsilon)^N_snp  ~  1.0 at WGS density
    P(L | chance, locus) = (1 - <2pq>(locus))^N_eff(L, locus)
      where N_eff(L, locus) = L * r(locus) / d_block_cM
        r(locus)  = local recombination rate from deCODE (cM/Mb)
        d_block_cM = effective genetic distance per independent block
                     (default 0.5 cM, empirical convention).
    P(IBD | L) = pi * P(L|IBD) / [pi * P(L|IBD) + (1-pi) * P(L|chance,locus)]

  Output:
    - Per-window lookup table (TSV): for each 1 Mb window on chr22, the
      posterior at multiple tract lengths.
    - "Length-to-threshold" table: for each window, the tract length at
      which the posterior crosses {0.50, 0.90, 0.95, 0.99}.
    - Summary table at 3 representative loci (pericentromeric, mid-arm,
      subtelomeric) showing the differential.

This is the figure-1-of-the-paper artifact: same physical tract length
gives different clinical confidence at different loci.
"""

import gzip
import math
import sys
from pathlib import Path

HERE = Path(__file__).parent
VCF_PATH = HERE / "chr22_phased.vcf.gz"
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"

OUT_LOOKUP = HERE / "chr22_posterior_lookup.tsv"
OUT_LENGTHS = HERE / "chr22_length_to_threshold.tsv"
OUT_SUMMARY = HERE / "chr22_posterior_summary.txt"

# parameters
CHROM = "chr22"
POPULATION = "EUR"                  # which AF to use (EUR, AFR, EAS, SAS, AMR, or "" for all)
PRIOR_PI = 0.0625                   # first-cousin offspring (Wright's coefficient)
GENOTYPING_ERROR = 0.001
BLOCK_CM = 0.5                      # genetic-distance per effective independent block
TRACT_LENGTHS_MB = [1, 2, 3, 5, 7, 10, 15]
THRESHOLDS = [0.50, 0.90, 0.95, 0.99]
WINDOW_BP = 1_000_000               # 1 Mb to match deCODE map's native resolution


def load_decode_map(path, chrom):
    """Return list of (pos, cMperMb) for the given chromosome."""
    out = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Chr"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            pos = int(f[1])
            cm_per_mb = float(f[3])
            out.append((pos, cm_per_mb))
    return out


def interp_recomb_rate(pat_map, mat_map):
    """
    Build a function pos -> sex-averaged cM/Mb.
    deCODE supplies rates at 1 Mb-window centers (500000, 1500000, ...);
    we do nearest-window lookup (simpler than linear interp; window-scale
    is what we operate at anyway).
    """
    pat_dict = dict(pat_map)
    mat_dict = dict(mat_map)
    centers = sorted(pat_dict.keys())

    def lookup(pos):
        # find nearest window center
        # window centers are at 500000, 1500000, ... => window of (n-0.5)Mb to (n+0.5)Mb
        idx_kb = pos // WINDOW_BP
        center = idx_kb * WINDOW_BP + 500_000
        if center not in pat_dict:
            # try neighbors
            for c in [center, center - WINDOW_BP, center + WINDOW_BP]:
                if c in pat_dict:
                    center = c
                    break
            else:
                return 0.0
        return 0.5 * (pat_dict[center] + mat_dict.get(center, pat_dict[center]))

    return lookup, centers


def stream_chr22_maf(vcf_path, population=POPULATION):
    """
    Stream the chr22 VCF and yield (pos, 2pq, snp_count_implicit_one).
    Skips multi-allelic and indels. Uses per-population AF from INFO.
    """
    af_key = f"AF_{population}=" if population else "AF="
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
            pq2 = 2.0 * af * (1.0 - af)
            yield pos, pq2


def posterior(L_mb, r_cmpermb, mean_2pq, pi=PRIOR_PI, eps=GENOTYPING_ERROR,
              block_cm=BLOCK_CM, n_snp_in_tract=1000):
    """
    Compute P(IBD | tract of length L_mb at a locus with recombination
    rate r_cmpermb and average heterozygosity mean_2pq).
    """
    # Effective number of independent observations within the tract
    genetic_distance_cm = L_mb * r_cmpermb
    n_eff = max(1.0, genetic_distance_cm / block_cm)

    # Probability of homozygosity at one independent block under HWE+unrelated
    p_hom_per_block = max(1e-12, 1.0 - mean_2pq)

    # Joint probability of observing a homozygous tract by chance
    p_chance = p_hom_per_block ** n_eff

    # P(L | IBD) ~ (1-eps)^N_snp, ~1 at WGS density
    p_ibd_data = (1.0 - eps) ** n_snp_in_tract

    num = pi * p_ibd_data
    denom = num + (1.0 - pi) * p_chance
    if denom <= 0:
        return 1.0
    return num / denom


def length_for_posterior(target_post, r_cmpermb, mean_2pq, **kwargs):
    """
    Numerically find the tract length L (Mb) at which the posterior crosses
    target_post. Uses bisection on a log-uniform grid 0.1 .. 50 Mb.
    """
    lo, hi = 0.1, 50.0
    # if even L=50 doesn't reach target, return None
    if posterior(hi, r_cmpermb, mean_2pq, **kwargs) < target_post:
        return None
    # if L=0.1 already exceeds, return 0.1
    if posterior(lo, r_cmpermb, mean_2pq, **kwargs) >= target_post:
        return lo
    for _ in range(60):
        mid = (lo + hi) / 2
        p = posterior(mid, r_cmpermb, mean_2pq, **kwargs)
        if p >= target_post:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def main():
    if not VCF_PATH.exists():
        sys.exit(f"missing: {VCF_PATH}")
    if not PAT_MAP.exists() or not MAT_MAP.exists():
        sys.exit(f"missing deCODE maps in {DECODE_DIR}")

    print(f"  chromosome: {CHROM}")
    print(f"  population (AF source): {POPULATION}")
    print(f"  prior pi (declared IBD probability): {PRIOR_PI}")
    print(f"  block_cm (effective independent block size): {BLOCK_CM} cM")
    print()

    # Load deCODE maps
    pat_map = load_decode_map(PAT_MAP, CHROM)
    mat_map = load_decode_map(MAT_MAP, CHROM)
    print(f"  deCODE chr22 windows: paternal {len(pat_map)}, maternal {len(mat_map)}")
    r_lookup, decode_centers = interp_recomb_rate(pat_map, mat_map)

    # Stream chr22 MAFs and aggregate by 1 Mb window
    print(f"  scanning {VCF_PATH.name} for per-variant AF_{POPULATION} ...")
    windows = {}   # window_start -> [sum_2pq, n_variants]
    n_total = 0
    n_used = 0
    for pos, pq2 in stream_chr22_maf(VCF_PATH, POPULATION):
        n_total += 1
        w = (pos // WINDOW_BP) * WINDOW_BP
        if w not in windows:
            windows[w] = [0.0, 0]
        windows[w][0] += pq2
        windows[w][1] += 1
        n_used += 1
    print(f"    variants seen: {n_total:,}; used (AF_{POPULATION} valid biallelic SNV): {n_used:,}")
    print(f"    1 Mb windows with data: {len(windows)}")

    # Compute posterior table
    rows = []
    for w_start in sorted(windows):
        sum_pq, n_var = windows[w_start]
        if n_var == 0:
            continue
        mean_2pq = sum_pq / n_var
        r = r_lookup(w_start + WINDOW_BP // 2)
        posteriors_by_L = {}
        for L in TRACT_LENGTHS_MB:
            posteriors_by_L[L] = posterior(L, r, mean_2pq)
        length_at_thr = {}
        for thr in THRESHOLDS:
            length_at_thr[thr] = length_for_posterior(thr, r, mean_2pq)
        rows.append({
            "w_start": w_start,
            "w_end": w_start + WINDOW_BP,
            "n_variants": n_var,
            "mean_2pq": mean_2pq,
            "cMperMb": r,
            "posteriors_by_L": posteriors_by_L,
            "length_at_thr": length_at_thr,
        })

    # Write lookup TSV
    with OUT_LOOKUP.open("w") as fh:
        cols = ["window_start", "window_end", "n_variants", "mean_2pq", "cMperMb"] + \
               [f"post_L{L}Mb" for L in TRACT_LENGTHS_MB]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = [str(r["w_start"]), str(r["w_end"]), str(r["n_variants"]),
                    f"{r['mean_2pq']:.5f}", f"{r['cMperMb']:.4f}"]
            for L in TRACT_LENGTHS_MB:
                vals.append(f"{r['posteriors_by_L'][L]:.5f}")
            fh.write("\t".join(vals) + "\n")

    # Write length-to-threshold TSV
    with OUT_LENGTHS.open("w") as fh:
        cols = ["window_start", "window_end", "cMperMb", "mean_2pq"] + \
               [f"L_for_post_{thr}_Mb" for thr in THRESHOLDS]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = [str(r["w_start"]), str(r["w_end"]),
                    f"{r['cMperMb']:.4f}", f"{r['mean_2pq']:.5f}"]
            for thr in THRESHOLDS:
                L = r["length_at_thr"][thr]
                vals.append("NA" if L is None else f"{L:.2f}")
            fh.write("\t".join(vals) + "\n")

    # Summary: 3 representative loci
    summary_lines = []
    summary_lines.append("=" * 75)
    summary_lines.append(f"chr22 per-locus IBD posterior, prior pi = {PRIOR_PI} (1st cousin), pop = {POPULATION}")
    summary_lines.append("=" * 75)
    summary_lines.append("")

    # Pick representative loci: lowest r, median r, highest r among data-bearing windows
    sorted_by_r = sorted([r for r in rows if r["cMperMb"] > 0], key=lambda x: x["cMperMb"])
    if len(sorted_by_r) >= 3:
        lo_locus = sorted_by_r[0]
        mid_locus = sorted_by_r[len(sorted_by_r) // 2]
        hi_locus = sorted_by_r[-1]
        representatives = [("LOWEST r  ", lo_locus), ("MEDIAN r  ", mid_locus), ("HIGHEST r ", hi_locus)]
        for label, r in representatives:
            summary_lines.append(f"{label} chr22:{r['w_start']:,}-{r['w_end']:,}  "
                                 f"r={r['cMperMb']:.2f} cM/Mb  <2pq>={r['mean_2pq']:.3f}  "
                                 f"n_var={r['n_variants']:,}")
            line = "    P(IBD | L) ="
            for L in TRACT_LENGTHS_MB:
                line += f"  L={L}Mb: {r['posteriors_by_L'][L]:.3f}"
            summary_lines.append(line)
            line = "    Length to reach posterior: "
            for thr in THRESHOLDS:
                L = r["length_at_thr"][thr]
                line += f"  {thr}: {('NA' if L is None else f'{L:.1f}Mb')}"
            summary_lines.append(line)
            summary_lines.append("")

    # The wedge: same length, different posterior across loci
    summary_lines.append("=" * 75)
    summary_lines.append("THE CLINICAL WEDGE: same 5 Mb tract, different loci")
    summary_lines.append("=" * 75)
    for L_focus in [2, 5, 10]:
        summary_lines.append(f"  tract length {L_focus} Mb:")
        for label, r in representatives:
            p = r["posteriors_by_L"][L_focus]
            summary_lines.append(f"    {label} (r={r['cMperMb']:.2f}): P(IBD) = {p:.3f}")

    # conventional 10 Mb comparison point
    summary_lines.append("")
    summary_lines.append("=" * 75)
    summary_lines.append("CONVENTIONAL 10 Mb CLINICAL-LAB POINT: locus-by-locus posterior")
    summary_lines.append("=" * 75)
    summary_lines.append("  (The ACMG-2021 standard counts segments >3-5 Mb as likely IBD; many")
    summary_lines.append("   clinical labs operationally use 5/7/10 Mb. This table shows whether")
    summary_lines.append("   the field's '10 Mb = confident IBD' actually holds across chr22 loci.)")
    summary_lines.append("")
    summary_lines.append("    locus_range          r(cM/Mb)   <2pq>    P(IBD|10Mb)")
    for r in rows:
        if r["cMperMb"] == 0:
            continue
        p10 = r["posteriors_by_L"][10]
        flag = "**" if p10 < 0.95 else "  "
        summary_lines.append(f"    {flag}chr22:{r['w_start']:>11,}-{r['w_end']:>11,}  "
                             f"{r['cMperMb']:6.3f}    {r['mean_2pq']:.3f}      {p10:.3f}")

    summary_lines.append("")
    summary_lines.append("  ** = posterior at 10 Mb < 0.95 (the field's nominal 'confident' threshold)")

    out_text = "\n".join(summary_lines)
    print(out_text)
    with OUT_SUMMARY.open("w") as fh:
        fh.write(out_text + "\n")

    print()
    print(f"  -> {OUT_LOOKUP}")
    print(f"  -> {OUT_LENGTHS}")
    print(f"  -> {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
