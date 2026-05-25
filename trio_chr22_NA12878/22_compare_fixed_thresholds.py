"""
Comparator: a fixed ROH length is NOT a fixed weight of evidence (vs ACMG/PLINK
fixed-Mb thresholds). Reviewer-requested comparison against existing fixed cutoffs.

From the length-evidence law, the weight of evidence for autozygosity scales with
the product (recombination rate r) x (ROH length L): log10 BF is ~ proportional to
r*L*log10(1/Hbar). Within a population H-bar is ~constant across loci, so at a
FIXED length the evidence varies directly with the local recombination rate, which
ranges several-fold across the genome. Hence a single Mb cutoff (e.g., ACMG 10 Mb)
conflates strong and weak evidence, and a short ROH at a recombination-rich locus
can carry the same evidence as a long ROH at a recombination-poor locus.

This script quantifies that using the per-window deCODE rate (read from
cross_pop_hap_diversity.tsv, cMperMb column) and reports, for a reference 10 Mb
ROH at the median-rate locus, the EQUIVALENT-evidence ROH length at each
recombination percentile.

Output: fixed_threshold_comparison.txt
Usage:  python 22_compare_fixed_thresholds.py
"""

from pathlib import Path
import numpy as np

HERE = Path(__file__).parent
DIV = HERE / "cross_pop_hap_diversity.tsv"
OUT = HERE / "fixed_threshold_comparison.txt"
REF_L = 10.0   # reference ROH length (Mb) at the median-rate locus (ACMG-style cutoff)


def main():
    seen = {}
    with DIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[ix["chrom"]], f[ix["window_start"]])
            try:
                r = float(f[ix["cMperMb"]])
            except (ValueError, IndexError):
                continue
            if r > 0:
                seen[key] = r            # same cMperMb across populations -> dedupe
    r = np.array(list(seen.values()))
    pct = {p: float(np.percentile(r, p)) for p in (10, 25, 50, 75, 90)}
    rmed = pct[50]
    # equal evidence <=> equal r*L  =>  L_equiv(p) = rmed*REF_L / r_p
    lines = []
    lines.append("# Fixed length is not fixed evidence: ROH lengths of EQUAL weight of")
    lines.append(f"# evidence to a {REF_L:.0f} Mb ROH at the median-rate locus.")
    lines.append(f"# deCODE sex-averaged rate over {r.size} autosomal 1-Mb windows; "
                 f"median {rmed:.2f} cM/Mb, range {r.min():.2f}-{r.max():.2f}.\n")
    lines.append("recombination percentile\tcM/Mb\tequal-evidence ROH length (Mb)")
    for p in (10, 25, 50, 75, 90):
        lines.append(f"{p}th\t{pct[p]:.2f}\t{rmed*REF_L/pct[p]:.1f}")
    fold = pct[90] / pct[10]
    lines.append(f"\nAcross the 10th-90th percentile of recombination rate, the same "
                 f"evidence is reached at ROH lengths spanning ~{fold:.0f}-fold "
                 f"({rmed*REF_L/pct[90]:.1f} Mb at hot loci vs {rmed*REF_L/pct[10]:.1f} "
                 f"Mb at cold loci).")
    lines.append("Implication: a single Mb threshold (ACMG 2/5/10 Mb; PLINK fixed "
                 "--homozyg-kb) treats these as equivalent; the per-locus weight of "
                 "evidence does not. A short ROH at a recombination-rich locus can be "
                 "more diagnostic than a long ROH at a recombination-poor locus.")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
