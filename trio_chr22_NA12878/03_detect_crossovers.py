"""
Detect crossovers on NA12878's chr22 haplotypes from the trio TSV.

Method
------
At each variant we have phased genotypes for father (F0|F1), mother (M0|M1),
and child (C0|C1). Under the SHAPEIT/Eagle trio-aware phasing convention,
C0 is the paternally-inherited allele, C1 is the maternally-inherited one.
We empirically verify this convention by counting Mendelian-consistent vs
-inconsistent assignments under each ordering and picking the consistent one.

For each "paternal-informative" site (father heterozygous, F0 != F1) we
identify which paternal haplotype was transmitted by checking whether
C0 == F0 or C0 == F1. We get a sequence of 0/1 calls along the chromosome.
A meiotic crossover is a switch in this sequence. To suppress phasing-switch
errors (which look like singleton flips), we require MIN_RUN consecutive
informative sites in the new state before calling a switch a crossover.

The crossover position is reported as the interval between the last variant
of the previous stable block and the first variant of the new stable block.
"""

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).parent
IN_TSV = HERE / "trio_chr22.tsv.gz"
OUT_CROSSOVERS = HERE / "crossovers_chr22.tsv"
OUT_SUMMARY = HERE / "summary_chr22.txt"

# Real meiotic crossovers produce haplotype-state changes that persist for
# megabases (= thousands of informative sites). Phasing-switch errors produce
# short flips of 1-100 sites. We require a stable run of at least MIN_RUN
# consecutive informative sites on BOTH sides of a switch before calling it
# a crossover. 500 sites at chr22 density (~32K informative sites in ~50 Mb)
# corresponds to roughly 750 kb of persistence, which is conservative for
# real meiotic crossovers (typically separated by ~10s of Mb) and stringent
# against phasing errors.
MIN_RUN = 2000


def parse_gt(gt: str):
    """Parse 'a|b' -> (int(a), int(b)). Returns None if unparseable."""
    if "|" not in gt:
        return None
    a, b = gt.split("|", 1)
    try:
        return int(a), int(b)
    except ValueError:
        return None


def load_trio(path: Path):
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos = int(f[idx["POS"]])
            father = parse_gt(f[idx["FATHER"]])
            mother = parse_gt(f[idx["MOTHER"]])
            child = parse_gt(f[idx["CHILD"]])
            if father is None or mother is None or child is None:
                continue
            rows.append((pos, father, mother, child))
    return rows


def check_phasing_convention(rows):
    """
    Test which child haplotype is paternally inherited.
    Under convention A (child[0]=paternal, child[1]=maternal):
      child[0] must be in {father[0], father[1]}
      child[1] must be in {mother[0], mother[1]}
    Under convention B (reversed), the opposite.
    Returns ('A', frac_consistent_A) or ('B', frac_consistent_B), whichever wins.
    """
    a_ok = b_ok = 0
    total = 0
    for _pos, f, m, c in rows:
        # Mendelian: under convention A
        ok_a = (c[0] in f) and (c[1] in m)
        ok_b = (c[1] in f) and (c[0] in m)
        if ok_a:
            a_ok += 1
        if ok_b:
            b_ok += 1
        total += 1
    return ("A", a_ok / total, b_ok / total) if a_ok >= b_ok else ("B", a_ok / total, b_ok / total)


def trace_parent(rows, parent_idx_in_row, child_allele_idx):
    """
    Walk along the chromosome. At each site where the relevant parent is
    heterozygous, record which parental haplotype was transmitted.

    parent_idx_in_row: 1 for father, 2 for mother in our row tuple.
    child_allele_idx: 0 if we believe child[0]=this-parent's-contribution, else 1.

    Returns a list of (pos, parental_hap_index_0_or_1) for informative sites only.
    """
    trace = []
    skip = 0
    for row in rows:
        pos = row[0]
        parent = row[parent_idx_in_row]
        child = row[3]
        if parent[0] == parent[1]:
            continue  # parent homozygous: site is not informative
        child_allele = child[child_allele_idx]
        if child_allele == parent[0]:
            trace.append((pos, 0))
        elif child_allele == parent[1]:
            trace.append((pos, 1))
        else:
            # Mendelian-inconsistent (genotyping error). Skip.
            skip += 1
    return trace, skip


def call_crossovers(trace, min_run=MIN_RUN):
    """
    Identify stable haplotype-state blocks of at least min_run consecutive
    informative sites, then report crossovers as transitions between
    successive stable blocks of DIFFERENT haplotype values.

    Short flips between long stable blocks are interpreted as phasing
    errors and ignored. If two successive stable blocks happen to have
    the same haplotype value (i.e., the chromosome was on hap A for a
    long stretch, briefly flipped to hap B for >= min_run sites due to
    a clustered phasing-error region, then resumed hap A), that is
    suspicious and reported separately as a "candidate flip-region"
    rather than two crossovers.

    Returns
    -------
    crossovers : list of (prev_pos, new_pos, prev_hap, new_hap, support_before, support_after)
    """
    if len(trace) < 2 * min_run:
        return []

    # Build all maximal runs of consecutive same-hap calls.
    runs = []  # (start_idx, end_idx_exclusive, hap)
    i = 0
    while i < len(trace):
        v = trace[i][1]
        j = i + 1
        while j < len(trace) and trace[j][1] == v:
            j += 1
        runs.append((i, j, v))
        i = j

    # Keep only long runs (>= min_run consecutive sites).
    long_runs = [r for r in runs if (r[1] - r[0]) >= min_run]

    if len(long_runs) < 2:
        return []

    # Crossovers = transitions between successive long runs of DIFFERENT haplotype.
    crossovers = []
    for k in range(1, len(long_runs)):
        prev = long_runs[k - 1]
        new = long_runs[k]
        if prev[2] == new[2]:
            # Two long runs of same hap separated by a flipped region.
            # Could be a clustered phasing-error region; not a crossover.
            continue
        prev_pos = trace[prev[1] - 1][0]
        new_pos = trace[new[0]][0]
        crossovers.append((
            prev_pos,
            new_pos,
            prev[2],
            new[2],
            prev[1] - prev[0],   # support before (long-run length)
            new[1] - new[0],     # support after (long-run length)
        ))
    return crossovers


def main():
    if not IN_TSV.exists():
        sys.exit(f"input not found: {IN_TSV}. Run 02_extract_trio.py first.")

    print(f"  loading {IN_TSV} ...")
    rows = load_trio(IN_TSV)
    print(f"  loaded {len(rows):,} variants with complete phased trio genotypes")

    conv, frac_a, frac_b = check_phasing_convention(rows)
    print(f"  Mendelian-consistency under convention A (child[0]=paternal): {frac_a:.4f}")
    print(f"  Mendelian-consistency under convention B (child[0]=maternal): {frac_b:.4f}")
    print(f"  -> using convention {conv}")

    if conv == "A":
        pat_child_idx, mat_child_idx = 0, 1
    else:
        pat_child_idx, mat_child_idx = 1, 0

    print()
    print("  tracing paternal haplotype...")
    pat_trace, pat_skip = trace_parent(rows, parent_idx_in_row=1, child_allele_idx=pat_child_idx)
    print(f"    paternal-informative sites (father het): {len(pat_trace):,}")
    print(f"    Mendelian-inconsistent (skipped): {pat_skip:,}")

    print("  tracing maternal haplotype...")
    mat_trace, mat_skip = trace_parent(rows, parent_idx_in_row=2, child_allele_idx=mat_child_idx)
    print(f"    maternal-informative sites (mother het): {len(mat_trace):,}")
    print(f"    Mendelian-inconsistent (skipped): {mat_skip:,}")

    pat_cos = call_crossovers(pat_trace, min_run=MIN_RUN)
    mat_cos = call_crossovers(mat_trace, min_run=MIN_RUN)

    print()
    print(f"  PATERNAL crossovers (NA12891 spermatogenesis -> NA12878):")
    for prev, new, ph, nh, sb, sa in pat_cos:
        gap = new - prev
        print(f"    between chr22:{prev:,} and chr22:{new:,}  ({gap:,} bp gap;  hap{ph}->hap{nh};  support {sb} -> {sa} sites)")
    if not pat_cos:
        print("    (none)")

    print()
    print(f"  MATERNAL crossovers (NA12892 oogenesis -> NA12878):")
    for prev, new, ph, nh, sb, sa in mat_cos:
        gap = new - prev
        print(f"    between chr22:{prev:,} and chr22:{new:,}  ({gap:,} bp gap;  hap{ph}->hap{nh};  support {sb} -> {sa} sites)")
    if not mat_cos:
        print("    (none)")

    # write outputs
    with OUT_CROSSOVERS.open("w") as fh:
        fh.write("parent\tprev_pos\tnew_pos\tgap_bp\tprev_hap\tnew_hap\tsupport_before\tsupport_after\n")
        for prev, new, ph, nh, sb, sa in pat_cos:
            fh.write(f"paternal\t{prev}\t{new}\t{new-prev}\t{ph}\t{nh}\t{sb}\t{sa}\n")
        for prev, new, ph, nh, sb, sa in mat_cos:
            fh.write(f"maternal\t{prev}\t{new}\t{new-prev}\t{ph}\t{nh}\t{sb}\t{sa}\n")

    with OUT_SUMMARY.open("w") as fh:
        fh.write("# NA12878 chr22 crossover detection summary\n")
        fh.write(f"variants_total\t{len(rows)}\n")
        fh.write(f"phasing_convention\t{conv}\n")
        fh.write(f"mendel_consistency_A\t{frac_a:.4f}\n")
        fh.write(f"mendel_consistency_B\t{frac_b:.4f}\n")
        fh.write(f"paternal_informative_sites\t{len(pat_trace)}\n")
        fh.write(f"maternal_informative_sites\t{len(mat_trace)}\n")
        fh.write(f"paternal_mendel_inconsistent\t{pat_skip}\n")
        fh.write(f"maternal_mendel_inconsistent\t{mat_skip}\n")
        fh.write(f"paternal_crossovers\t{len(pat_cos)}\n")
        fh.write(f"maternal_crossovers\t{len(mat_cos)}\n")
        fh.write(f"min_run_threshold\t{MIN_RUN}\n")

    print()
    print(f"  -> {OUT_CROSSOVERS}")
    print(f"  -> {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
