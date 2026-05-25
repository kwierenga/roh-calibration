"""
Stream the chr22 phased VCF, extract just the NA12878/NA12891/NA12892 trio,
and write a compact TSV with phased genotypes. Discards all other 3,199 samples.

Output columns (per phased-VCF convention from trio-aware phasing):
  CHROM, POS, REF, ALT, FATHER_h0|h1, MOTHER_h0|h1, CHILD_h0|h1
  where for the CHILD: h0 = paternal-inherited allele, h1 = maternal-inherited
  (by SHAPEIT trio convention; we verify Mendelian consistency in step 3).

Skips multi-allelic sites and non-SNV variants for simplicity (chr22 phased
file is biallelic-decomposed already, but we double-check).
"""

import gzip
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
IN_VCF = HERE / "chr22_phased.vcf.gz"
OUT_TSV = HERE / "trio_chr22.tsv.gz"

CHILD, FATHER, MOTHER = "NA12878", "NA12891", "NA12892"


def main():
    if not IN_VCF.exists():
        sys.exit(f"input not found: {IN_VCF}. Run 01_download.py first.")

    t0 = time.time()
    n_lines = 0
    n_kept = 0
    n_skipped_multi = 0
    n_skipped_missing = 0

    with gzip.open(IN_VCF, "rt") as fin, gzip.open(OUT_TSV, "wt", encoding="utf-8") as fout:
        col_child = col_father = col_mother = None

        fout.write("CHROM\tPOS\tID\tREF\tALT\tFATHER\tMOTHER\tCHILD\n")

        for line in fin:
            n_lines += 1
            if n_lines % 100_000 == 0:
                elapsed = time.time() - t0
                rate = n_lines / max(elapsed, 1e-9)
                sys.stdout.write(f"\r  read {n_lines:,} lines, kept {n_kept:,}  ({rate:,.0f} lines/s, {elapsed:.0f}s)")
                sys.stdout.flush()

            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                # samples start at column 9 (after FORMAT)
                try:
                    col_child = header.index(CHILD)
                    col_father = header.index(FATHER)
                    col_mother = header.index(MOTHER)
                except ValueError as e:
                    sys.exit(f"sample not in header: {e}")
                print(f"  sample columns: father={col_father} mother={col_mother} child={col_child}")
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            # skip multi-allelic and non-SNV (we want clean biallelic SNPs for crossover detection)
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                n_skipped_multi += 1
                continue

            gt_father = fields[col_father].split(":")[0]
            gt_mother = fields[col_mother].split(":")[0]
            gt_child = fields[col_child].split(":")[0]

            # require phased and complete for all three
            if "|" not in gt_father or "|" not in gt_mother or "|" not in gt_child:
                n_skipped_missing += 1
                continue
            if "." in gt_father or "." in gt_mother or "." in gt_child:
                n_skipped_missing += 1
                continue

            fout.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{gt_father}\t{gt_mother}\t{gt_child}\n")
            n_kept += 1

    print()
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s")
    print(f"  total VCF lines: {n_lines:,}")
    print(f"  kept (biallelic SNV, all phased): {n_kept:,}")
    print(f"  skipped multi-allelic / indel: {n_skipped_multi:,}")
    print(f"  skipped unphased / missing: {n_skipped_missing:,}")
    print(f"  -> {OUT_TSV}")


if __name__ == "__main__":
    main()
