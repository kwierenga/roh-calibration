"""
Compare NA12878 chr22 phasing between two sources:
  A) 1000G high-coverage 20220422 release (statistical phasing of 3,202 samples
     with trio constraints), extracted into trio_chr22.tsv.gz.
  B) Illumina Platinum Genomes 2017-1.0 release (pedigree-phasing against full
     CEPH 1463 pedigree, with NA12891 and NA12892 explicitly declared as
     parents in the VCF PEDIGREE header).

Both use convention "haplotype 0 = paternal-inherited". A disagreement at a
heterozygous site means one of the two phasings has assigned the alleles to
the wrong parental side.

If the Platinum pedigree-phasing is correct (the strongest pedigree
constraints in this comparison), then any disagreement is a 1000G phasing
artifact. Long blocks of consistent disagreement are *block-level label
swaps* in the 1000G statistical phasing — exactly the methodology issue we
hypothesized was inflating our crossover-detection candidate count.

Output
------
- per-position disagreement TSV (chr22:pos, REF, ALT, 1000g_call, platinum_call,
  agree True/False)
- a windowed summary: in each 250 kb window, what fraction of comparable
  heterozygous sites agree between the two phasings?
- written into ASCII chromosome trace showing agreement pattern, so
  long-stretch disagreements (= block swaps) are visible at a glance
"""

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).parent
TRIO_TSV = HERE / "trio_chr22.tsv.gz"
PLATINUM_VCF = HERE / "external" / "platinum_NA12878_hg38.vcf.gz"
OUT_DISAGREE = HERE / "phasing_comparison_chr22.tsv.gz"
OUT_WINDOW = HERE / "phasing_comparison_windowed.tsv"
OUT_ASCII = HERE / "phasing_comparison_ascii.txt"

WINDOW_BP = 250_000


def load_thousand_g_na12878():
    """
    Load NA12878's chr22 phased genotypes from the trio TSV.
    Returns dict: pos -> (ref, alt, gt_str)  where gt_str is 'a|b' with ints.
    """
    out = {}
    with gzip.open(TRIO_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos = int(f[idx["POS"]])
            ref = f[idx["REF"]]
            alt = f[idx["ALT"]]
            child_gt = f[idx["CHILD"]]
            out[pos] = (ref, alt, child_gt)
    return out


def load_platinum_na12878():
    """
    Stream the Platinum NA12878 VCF (chr22 only) and return dict
    pos -> (ref, alt, gt_str).
    """
    out = {}
    with gzip.open(PLATINUM_VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    if cols[-1] != "NA12878":
                        print(f"  warning: last sample column is {cols[-1]}, expected NA12878")
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != "chr22":
                continue
            chrom, pos_s, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            # skip multi-allelic / indels
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            gt = f[-1].split(":")[0]
            if "|" not in gt:
                continue
            if "." in gt:
                continue
            out[int(pos_s)] = (ref, alt, gt)
    return out


def main():
    if not TRIO_TSV.exists() or not PLATINUM_VCF.exists():
        sys.exit("required inputs not found.")

    print("  loading NA12878 chr22 from 1000G trio TSV ...")
    g1000 = load_thousand_g_na12878()
    print(f"    {len(g1000):,} phased biallelic SNV positions")

    print("  loading NA12878 chr22 from Platinum Genomes VCF ...")
    plat = load_platinum_na12878()
    print(f"    {len(plat):,} phased biallelic SNV positions")

    # intersection of positions
    common = sorted(set(g1000.keys()) & set(plat.keys()))
    print(f"  intersection: {len(common):,} shared positions")

    # at each shared site, check if it's heterozygous in both and if so,
    # whether the phasing agrees on which allele is on hap 0.
    n_total = 0
    n_both_het = 0
    n_agree = 0
    n_disagree = 0
    n_ref_alt_mismatch = 0
    disagreements = []
    rows_for_output = []

    for pos in common:
        ref1, alt1, gt1 = g1000[pos]
        ref2, alt2, gt2 = plat[pos]
        n_total += 1
        if ref1 != ref2 or alt1 != alt2:
            n_ref_alt_mismatch += 1
            continue

        a1, b1 = (int(x) for x in gt1.split("|"))
        a2, b2 = (int(x) for x in gt2.split("|"))

        is_het_1 = a1 != b1
        is_het_2 = a2 != b2
        if not (is_het_1 and is_het_2):
            continue
        n_both_het += 1

        # compare which allele is on hap 0
        if a1 == a2:
            n_agree += 1
            rows_for_output.append((pos, ref1, alt1, gt1, gt2, True))
        else:
            n_disagree += 1
            disagreements.append(pos)
            rows_for_output.append((pos, ref1, alt1, gt1, gt2, False))

    print()
    print(f"  comparison summary:")
    print(f"    shared sites: {n_total:,}")
    print(f"    REF/ALT mismatched (skipped): {n_ref_alt_mismatch:,}")
    print(f"    both-heterozygous sites (the comparable set): {n_both_het:,}")
    print(f"    phasings AGREE on hap0 allele: {n_agree:,} ({n_agree/n_both_het*100:.2f}%)")
    print(f"    phasings DISAGREE on hap0 allele: {n_disagree:,} ({n_disagree/n_both_het*100:.2f}%)")

    # write per-site TSV
    with gzip.open(OUT_DISAGREE, "wt") as fh:
        fh.write("pos\tref\talt\t1000g_gt\tplatinum_gt\tagree\n")
        for pos, ref, alt, gt1, gt2, agree in rows_for_output:
            fh.write(f"{pos}\t{ref}\t{alt}\t{gt1}\t{gt2}\t{agree}\n")

    # windowed agreement
    windows = {}
    for pos, ref, alt, gt1, gt2, agree in rows_for_output:
        w = (pos // WINDOW_BP) * WINDOW_BP
        if w not in windows:
            windows[w] = [0, 0]  # [agree, disagree]
        if agree:
            windows[w][0] += 1
        else:
            windows[w][1] += 1

    with OUT_WINDOW.open("w") as fh:
        fh.write("window_start\twindow_end\tn_agree\tn_disagree\tagree_frac\n")
        for w in sorted(windows):
            n_a, n_d = windows[w]
            tot = n_a + n_d
            frac = n_a / tot if tot else 0
            fh.write(f"{w}\t{w+WINDOW_BP}\t{n_a}\t{n_d}\t{frac:.3f}\n")

    # ASCII chromosome trace: one char per 250 kb window
    # 'A' = >=95% agree;  'a' = 80-95%;  'x' = mixed;  'D' = >=95% disagree
    def code(n_a, n_d):
        tot = n_a + n_d
        if tot < 5:
            return "."
        f = n_a / tot
        if f >= 0.95:
            return "A"
        elif f >= 0.80:
            return "a"
        elif f <= 0.05:
            return "D"
        elif f <= 0.20:
            return "d"
        else:
            return "x"

    sorted_w = sorted(windows)
    trace = "".join(code(*windows[w]) for w in sorted_w)

    with OUT_ASCII.open("w") as fh:
        fh.write("# NA12878 chr22 phasing-agreement trace (1000G vs Platinum)\n")
        fh.write(f"# each char = {WINDOW_BP/1000:.0f} kb window\n")
        fh.write("# 'A' = >=95% of het sites agree (consistent phasing)\n")
        fh.write("# 'a' = 80-95% agree\n")
        fh.write("# 'x' = 20-80% (mixed)\n")
        fh.write("# 'd' = 5-20% agree\n")
        fh.write("# 'D' = <=5% agree (== >=95% disagree, i.e. systematic label SWAP in 1000G phasing)\n")
        fh.write("# '.' = fewer than 5 comparable sites\n\n")
        WRAP = 100
        for chunk_start in range(0, len(trace), WRAP):
            chunk_end = min(chunk_start + WRAP, len(trace))
            pos_start = sorted_w[chunk_start]
            pos_end = sorted_w[chunk_end - 1] + WINDOW_BP
            fh.write(f"\nchr22:{pos_start:,}-{pos_end:,}\n")
            fh.write(f"  {trace[chunk_start:chunk_end]}\n")

    print()
    print("CHROMOSOME-WIDE AGREEMENT TRACE:")
    print("  'A'=agree, 'D'=systematic SWAP, 'x'=mixed, '.'=sparse")
    print()
    print(f"  {trace}")
    print()
    print(f"  -> {OUT_DISAGREE}")
    print(f"  -> {OUT_WINDOW}")
    print(f"  -> {OUT_ASCII}")


if __name__ == "__main__":
    main()
