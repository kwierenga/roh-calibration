"""
Recompute the per-locus IBD posterior at multiple priors, using the
per-window stats (cM/Mb, mean_2pq, variant count) already on disk from
11_run_all_autosomes.py. The heavy lifting (VCF parsing) is already done;
this is pure recomputation of the Bayesian posterior at different priors.

Priors:
  pi = 0.0625  (1st-cousin offspring)
  pi = 0.125   (double-1st-cousin / avuncular / half-sibling offspring)
  pi = 0.25    (sibling-incest offspring)

Output:
  per-prior TSV (full lookup table)
  combined summary showing calibration fraction at each prior

This makes the prior-dependence of the conventional 10 Mb clinical-lab
comparison point explicit and clinically actionable. (The ACMG-2021 standard
counts segments >3-5 Mb as likely IBD; 10 Mb is common lab practice, not the
standard.)
"""

import math
import sys
from pathlib import Path

HERE = Path(__file__).parent
IN_LOOKUP = HERE / "all_autosomes_posterior_lookup.tsv"
OUT_SUMMARY = HERE / "all_autosomes_multi_prior_summary.txt"

PRIORS = [
    (0.0156, "2nd_cousin"),
    (0.0625, "1st_cousin"),
    (0.125,  "avuncular_or_double_1c"),
    (0.25,   "incest_or_sibling"),
]

GENOTYPING_ERROR = 0.001
BLOCK_CM = 0.5
TRACT_LENGTHS_MB = [1, 2, 3, 5, 7, 10, 15]
THRESHOLDS = [0.50, 0.90, 0.95, 0.99]
CONV_LENGTH = 10        # conventional clinical-lab comparison length (Mb); the
                        # ACMG-2021 standard is >3-5 Mb, not 10 Mb (10 Mb = lab practice).
POST_THRESHOLD = 0.95


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


def main():
    if not IN_LOOKUP.exists():
        sys.exit(f"missing: {IN_LOOKUP}. Run 11_run_all_autosomes.py first.")

    print(f"  loading per-window stats from {IN_LOOKUP.name} ...")
    rows = []
    with IN_LOOKUP.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append({
                "chrom":     f[idx["chrom"]],
                "w_start":   int(f[idx["window_start"]]),
                "w_end":     int(f[idx["window_end"]]),
                "n_variants": int(f[idx["n_variants"]]),
                "mean_2pq":  float(f[idx["mean_2pq"]]),
                "cMperMb":   float(f[idx["cMperMb"]]),
            })
    print(f"    {len(rows):,} windows loaded across all 22 autosomes")
    print()

    # per-prior tables + summary stats
    per_prior_stats = {}
    for pi, name in PRIORS:
        out_path = HERE / f"all_autosomes_pi_{pi:.4f}_{name}.tsv"
        n_conv_calibrated = 0
        n_with_recomb = 0
        per_chrom = {}
        with out_path.open("w") as fh:
            cols = ["chrom", "window_start", "window_end", "n_variants",
                    "mean_2pq", "cMperMb"] + \
                   [f"post_L{L}Mb" for L in TRACT_LENGTHS_MB] + \
                   [f"L_for_post_{thr}_Mb" for thr in THRESHOLDS]
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                if r["cMperMb"] > 0:
                    n_with_recomb += 1
                    per_chrom.setdefault(r["chrom"], [0, 0])
                    per_chrom[r["chrom"]][1] += 1
                vals = [r["chrom"], str(r["w_start"]), str(r["w_end"]),
                        str(r["n_variants"]),
                        f"{r['mean_2pq']:.5f}", f"{r['cMperMb']:.4f}"]
                for L in TRACT_LENGTHS_MB:
                    p = posterior(L, r["cMperMb"], r["mean_2pq"], pi)
                    vals.append(f"{p:.5f}")
                    if L == CONV_LENGTH and p >= POST_THRESHOLD and r["cMperMb"] > 0:
                        n_conv_calibrated += 1
                        per_chrom[r["chrom"]][0] += 1
                for thr in THRESHOLDS:
                    L = length_for_posterior(thr, r["cMperMb"], r["mean_2pq"], pi)
                    vals.append("NA" if L is None else f"{L:.2f}")
                fh.write("\t".join(vals) + "\n")
        per_prior_stats[(pi, name)] = {
            "out_path": out_path,
            "n_conv_calibrated": n_conv_calibrated,
            "n_with_recomb": n_with_recomb,
            "per_chrom": per_chrom,
        }
        frac = n_conv_calibrated / n_with_recomb if n_with_recomb else 0
        print(f"  pi = {pi:.4f}  ({name:25s})  "
              f"conventional-10Mb-calibrated windows: {n_conv_calibrated:4d}/{n_with_recomb:4d} = {100*frac:5.1f}%")

    # summary file
    with OUT_SUMMARY.open("w") as fh:
        fh.write("# Per-locus IBD posterior framework — calibration at the conventional 10 Mb clinical-lab point\n")
        fh.write("# (ACMG-2021 standard is >3-5 Mb; 10 Mb is common lab practice, not the standard)\n")
        fh.write(f"# Threshold: posterior(IBD | 10 Mb) >= {POST_THRESHOLD}\n")
        fh.write(f"# Genotyping error eps: {GENOTYPING_ERROR}\n")
        fh.write(f"# Effective independent block size: {BLOCK_CM} cM\n")
        fh.write(f"# AF source: 1000G high-coverage 20220422, EUR superpopulation\n\n")

        fh.write("=" * 80 + "\n")
        fh.write("HEADLINE: fraction of autosomal windows where conventional 10 Mb >= 0.95 posterior\n")
        fh.write("=" * 80 + "\n\n")
        fh.write("pi (Wright's F)\trelationship\tcalibrated/total\tfraction\n")
        for (pi, name) in PRIORS:
            s = per_prior_stats[(pi, name)]
            frac = s["n_conv_calibrated"] / s["n_with_recomb"]
            fh.write(f"{pi:.4f}\t{name}\t{s['n_conv_calibrated']}/{s['n_with_recomb']}\t{frac:.3f}\n")

        fh.write("\n" + "=" * 80 + "\n")
        fh.write("Per-chromosome breakdown at each prior (fraction calibrated)\n")
        fh.write("=" * 80 + "\n\n")

        all_chroms = sorted(set().union(*[set(per_prior_stats[k]["per_chrom"].keys()) for k in per_prior_stats]),
                            key=lambda c: int(c.replace("chr", "")))
        prior_labels = [f"pi={pi}" for pi, _ in PRIORS]
        fh.write("chrom\tn_windows\t" + "\t".join(prior_labels) + "\n")
        for chrom in all_chroms:
            line = [chrom]
            n_w = per_prior_stats[PRIORS[0]]["per_chrom"].get(chrom, [0, 0])[1]
            line.append(str(n_w))
            for pi, name in PRIORS:
                nc, nr = per_prior_stats[(pi, name)]["per_chrom"].get(chrom, [0, 0])
                frac = nc / nr if nr else 0
                line.append(f"{frac:.3f}")
            fh.write("\t".join(line) + "\n")

    print()
    print(f"  -> {OUT_SUMMARY}")
    for (pi, name), s in per_prior_stats.items():
        print(f"  -> {s['out_path']}")


if __name__ == "__main__":
    main()
