"""
Robust crossover detection on NA12878 chr22, with two methodological refinements
beyond 03_detect_crossovers.py:

1. Mendelian-deterministic sites only. For paternal crossover detection we use
   only sites where father is heterozygous AND mother is homozygous (so the
   transmitted paternal allele is unambiguously determinable from the child's
   genotype). Symmetric for maternal. This removes ambiguity from the both-
   heterozygous case where statistical phasing is the only source of phase.

2. Biological-plausibility filter. Crossover interference forbids real meiotic
   crossovers from being closer than ~50 cM. chr22's q arm has ~37 cM (male) /
   ~58 cM (female), so the per-meiosis maximum is realistically 1 paternal and
   1-2 maternal. Among our raw candidates, we keep at most the top-N by
   support score, where N=1 for paternal and N=2 for maternal. Calls beyond
   that biological maximum are flagged as candidate phasing artifacts.

Output also includes a "verdict" file that contrasts the raw candidate count
with the biologically plausible top-N, making clear which calls survive both
the statistical-support filter and the biological-plausibility filter.
"""

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).parent
IN_TSV = HERE / "trio_chr22.tsv.gz"
OUT_ROBUST = HERE / "crossovers_chr22_robust.tsv"
OUT_VERDICT = HERE / "verdict_chr22.txt"

MIN_RUN = 500          # informative-site run length for a block to be "stable"
MAX_PAT_PER_CHR = 1    # biological max for chr22 paternal
MAX_MAT_PER_CHR = 2    # biological max for chr22 maternal (slightly more lenient)


def parse_gt(gt):
    if "|" not in gt:
        return None
    a, b = gt.split("|", 1)
    try:
        return int(a), int(b)
    except ValueError:
        return None


def load_trio(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            father = parse_gt(f[idx["FATHER"]])
            mother = parse_gt(f[idx["MOTHER"]])
            child = parse_gt(f[idx["CHILD"]])
            if father is None or mother is None or child is None:
                continue
            rows.append((int(f[idx["POS"]]), father, mother, child))
    return rows


def trace_deterministic(rows, parent_side):
    """
    Build a haplotype-state trace using Mendelian-deterministic sites only.

    parent_side = "paternal": use sites where father is het AND mother is hom.
                              The transmitted paternal allele is then
                              unambiguous and we tag which paternal haplotype
                              label (0 or 1) was transmitted.
    parent_side = "maternal": symmetric.

    Returns (trace, n_used, n_total_informative)
      trace: list of (pos, hap_label_0_or_1) at deterministic sites only.
      n_used: number of deterministic sites where the call succeeded.
      n_total_informative: number of het-parent sites total (the larger
                           pool that the previous algorithm used).
    """
    trace = []
    n_total = 0
    n_used = 0
    for pos, father, mother, child in rows:
        if parent_side == "paternal":
            parent = father
            other = mother
            child_idx = 0  # paternal-inherited allele under convention A
        else:
            parent = mother
            other = father
            child_idx = 1
        if parent[0] == parent[1]:
            continue  # parent hom: uninformative for this parent's crossovers
        n_total += 1
        if other[0] != other[1]:
            continue  # other parent is also het: ambiguous transmission, skip
        # Now: parent is het, other parent is hom. Mendelian-deterministic.
        child_allele = child[child_idx]
        if child_allele == parent[0]:
            trace.append((pos, 0))
            n_used += 1
        elif child_allele == parent[1]:
            trace.append((pos, 1))
            n_used += 1
        # else: Mendelian inconsistency, skip silently.
    return trace, n_used, n_total


def long_runs(trace, min_run):
    """Return list of (start_idx, end_idx_exclusive, hap_value) for runs >= min_run."""
    if not trace:
        return []
    runs = []
    i = 0
    while i < len(trace):
        v = trace[i][1]
        j = i + 1
        while j < len(trace) and trace[j][1] == v:
            j += 1
        if (j - i) >= min_run:
            runs.append((i, j, v))
        i = j
    return runs


def call_candidate_crossovers(trace, min_run):
    runs = long_runs(trace, min_run)
    if len(runs) < 2:
        return []
    out = []
    for k in range(1, len(runs)):
        prev = runs[k - 1]
        new = runs[k]
        if prev[2] == new[2]:
            continue
        prev_pos = trace[prev[1] - 1][0]
        new_pos = trace[new[0]][0]
        # support score: harmonic mean of the two run lengths
        sb = prev[1] - prev[0]
        sa = new[1] - new[0]
        support = 2 * sb * sa / (sb + sa)
        out.append({
            "prev_pos": prev_pos,
            "new_pos": new_pos,
            "prev_hap": prev[2],
            "new_hap": new[2],
            "support_before": sb,
            "support_after": sa,
            "support_score": support,
        })
    return out


def apply_interference_filter(candidates, max_n):
    """
    Keep up to max_n candidates with highest support_score, sorted by position.
    Annotate each with whether it survived the filter.
    """
    sorted_by_support = sorted(candidates, key=lambda c: c["support_score"], reverse=True)
    keepers = set(id(c) for c in sorted_by_support[:max_n])
    for c in candidates:
        c["survives_interference"] = id(c) in keepers
    return candidates


def main():
    if not IN_TSV.exists():
        sys.exit(f"input not found: {IN_TSV}. Run 02_extract_trio.py first.")

    print(f"  loading {IN_TSV} ...")
    rows = load_trio(IN_TSV)
    print(f"  loaded {len(rows):,} variants")

    pat_trace, pat_used, pat_total_het = trace_deterministic(rows, "paternal")
    mat_trace, mat_used, mat_total_het = trace_deterministic(rows, "maternal")

    print()
    print(f"  PATERNAL (father het + mother hom): {pat_used:,} deterministic sites")
    print(f"    (vs {pat_total_het:,} total father-het sites used in 03_detect_crossovers.py)")
    print(f"    fraction usable: {pat_used/pat_total_het:.1%}")
    print(f"  MATERNAL (mother het + father hom): {mat_used:,} deterministic sites")
    print(f"    (vs {mat_total_het:,} total mother-het sites used in 03_detect_crossovers.py)")
    print(f"    fraction usable: {mat_used/mat_total_het:.1%}")

    pat_cands = call_candidate_crossovers(pat_trace, MIN_RUN)
    mat_cands = call_candidate_crossovers(mat_trace, MIN_RUN)

    pat_cands = apply_interference_filter(pat_cands, MAX_PAT_PER_CHR)
    mat_cands = apply_interference_filter(mat_cands, MAX_MAT_PER_CHR)

    def fmt(cand, label):
        check = "[KEEP]" if cand["survives_interference"] else "[drop]"
        return (f"    {check} {label}  chr22:{cand['prev_pos']:,} <-> chr22:{cand['new_pos']:,}  "
                f"({cand['new_pos']-cand['prev_pos']:,} bp gap; "
                f"hap{cand['prev_hap']}->hap{cand['new_hap']}; "
                f"support {cand['support_before']} / {cand['support_after']} sites, "
                f"score {cand['support_score']:.0f})")

    print()
    print(f"  PATERNAL candidates (raw count {len(pat_cands)}, biological max {MAX_PAT_PER_CHR}):")
    if not pat_cands:
        print("    (none)")
    for c in pat_cands:
        print(fmt(c, "PAT"))

    print()
    print(f"  MATERNAL candidates (raw count {len(mat_cands)}, biological max {MAX_MAT_PER_CHR}):")
    if not mat_cands:
        print("    (none)")
    for c in mat_cands:
        print(fmt(c, "MAT"))

    # write outputs
    with OUT_ROBUST.open("w") as fh:
        fh.write("parent\tsurvives_interference\tprev_pos\tnew_pos\tgap_bp\tprev_hap\tnew_hap\tsupport_before\tsupport_after\tsupport_score\n")
        for c in pat_cands:
            fh.write(f"paternal\t{c['survives_interference']}\t{c['prev_pos']}\t{c['new_pos']}\t{c['new_pos']-c['prev_pos']}\t{c['prev_hap']}\t{c['new_hap']}\t{c['support_before']}\t{c['support_after']}\t{c['support_score']:.1f}\n")
        for c in mat_cands:
            fh.write(f"maternal\t{c['survives_interference']}\t{c['prev_pos']}\t{c['new_pos']}\t{c['new_pos']-c['prev_pos']}\t{c['prev_hap']}\t{c['new_hap']}\t{c['support_before']}\t{c['support_after']}\t{c['support_score']:.1f}\n")

    with OUT_VERDICT.open("w") as fh:
        fh.write("NA12878 chr22 robust crossover detection verdict\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Method: Mendelian-deterministic sites only (one parent het, other hom)\n")
        fh.write(f"        Long-run threshold: {MIN_RUN} consecutive sites on both sides\n")
        fh.write(f"        Interference filter: top {MAX_PAT_PER_CHR} paternal, top {MAX_MAT_PER_CHR} maternal by support score\n\n")
        fh.write(f"Deterministic sites used: {pat_used:,} paternal, {mat_used:,} maternal\n\n")
        fh.write(f"Raw candidates: {len(pat_cands)} paternal, {len(mat_cands)} maternal\n")
        fh.write(f"Surviving interference filter: "
                 f"{sum(1 for c in pat_cands if c['survives_interference'])} paternal, "
                 f"{sum(1 for c in mat_cands if c['survives_interference'])} maternal\n\n")
        fh.write("CAVEAT\n------\n")
        fh.write("This analysis is performed on the 1000G high-coverage release 20220422\n")
        fh.write("phased VCF, which uses SHAPEIT-family phasing across the 3,202-sample\n")
        fh.write("panel with trio constraints. The raw candidate count is high (5+5 with\n")
        fh.write("MIN_RUN=500) vs. the deCODE-predicted ~1 expected for chr22 in a single\n")
        fh.write("meiosis. The interference filter selects the top candidates by support\n")
        fh.write("score, but cannot definitively distinguish:\n")
        fh.write("  (a) real meiotic crossovers, and\n")
        fh.write("  (b) parental-haplotype block-level label swaps in the released phasing.\n")
        fh.write("Both produce identical trace signatures and identical Mendelian consistency.\n\n")
        fh.write("Definitive validation requires one of:\n")
        fh.write("  - Grandparent genotypes (CEPH 1463 has NA12877/NA12879/NA12889/NA12890;\n")
        fh.write("    not present in 1000G 3,202-sample panel)\n")
        fh.write("  - Independent re-phasing with strict trio mode (SHAPEIT5)\n")
        fh.write("  - Platinum Genomes phased VCFs (Eberle 2017, 17-member pedigree phasing)\n")
        fh.write("  - deCODE recombination map (Halldorsson 2019) as a Bayesian prior on position\n\n")
        fh.write("Of these, Platinum Genomes is the lowest-friction next step.\n")

    print()
    print(f"  -> {OUT_ROBUST}")
    print(f"  -> {OUT_VERDICT}")


if __name__ == "__main__":
    main()
