"""
Subset a gnomAD v3.1.2 HGDP+1KG per-chrom VCF (~55-275 GB) down to an analysis-
ready slice we actually need:

  - Samples : the 828 HGDP samples only (column name starts with "HGDP")
  - Variants: biallelic SNVs, single-base REF and ALT, no commas in ALT
  - Filter  : minor allele frequency (computed on the HGDP subset) >= MAF_MIN
  - Columns : strip INFO/FORMAT to just `AF_HGDP=<float>` and `GT`

Result: ~50-200 MB per chrom (vs 50-275 GB raw).

Usage:
  python 29_subset_hgdp_vcf.py <chromN>
       e.g.   python 29_subset_hgdp_vcf.py 22

Reads:   hgdp_tgp/chr{N}.hgdp_tgp.vcf.bgz
Writes:  hgdp_tgp/chr{N}.hgdp_subset.vcf.gz
"""

import gzip
import sys
import time
from pathlib import Path

MAF_MIN = 0.05
HERE = Path(__file__).parent
DDIR = HERE / "hgdp_tgp"


def subset(chrom: str) -> dict:
    src = DDIR / f"chr{chrom}.hgdp_tgp.vcf.bgz"
    dst = DDIR / f"chr{chrom}.hgdp_subset.vcf.gz"
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    if not src.exists():
        raise FileNotFoundError(src)

    t0 = time.time()
    n_in = n_pass = n_kept = 0
    hgdp_cols: list[int] = []
    n_hgdp = 0

    with gzip.open(src, "rt") as fin, gzip.open(tmp, "wt", compresslevel=6) as fout:
        # 1. header pass-through (only ##fileformat + #CHROM with HGDP cols)
        for line in fin:
            if line.startswith("##"):
                # only forward fileformat + contig + reference lines (keep tiny)
                if (line.startswith("##fileformat")
                        or line.startswith("##contig=<ID=chr" + chrom + ",")
                        or line.startswith("##reference")):
                    fout.write(line)
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                hgdp_cols = [9 + i for i, s in enumerate(samples) if s.startswith("HGDP")]
                n_hgdp = len(hgdp_cols)
                hgdp_names = [samples[i - 9] for i in hgdp_cols]
                fout.write("##INFO=<ID=AF_HGDP,Number=1,Type=Float,"
                           "Description=\"Alt allele frequency in HGDP samples in this VCF\">\n")
                fout.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                           + "\t".join(hgdp_names) + "\n")
                break

        if not hgdp_cols:
            raise RuntimeError("no HGDP columns found in VCF")

        # 2. variant pass
        for line in fin:
            n_in += 1
            # split lazily; cheap pre-filter on tab boundaries to skip multiallelic
            # and indels before splitting full line
            # find first 5 tab positions
            t1 = line.find("\t")
            t2 = line.find("\t", t1 + 1)
            t3 = line.find("\t", t2 + 1)
            t4 = line.find("\t", t3 + 1)
            t5 = line.find("\t", t4 + 1)
            if t5 == -1:
                continue
            ref = line[t3 + 1:t4]
            alt = line[t4 + 1:t5]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            n_pass += 1

            f = line.rstrip("\n").split("\t")
            fmt = f[8]
            # find GT index in FORMAT
            fmt_parts = fmt.split(":")
            try:
                gt_idx = fmt_parts.index("GT")
            except ValueError:
                continue

            # collect GTs from HGDP columns, compute AC/AN
            ac = an = 0
            gts_out = []
            for c in hgdp_cols:
                cell = f[c]
                # extract GT subfield
                if gt_idx == 0:
                    colon = cell.find(":")
                    gt = cell if colon == -1 else cell[:colon]
                else:
                    parts = cell.split(":")
                    gt = parts[gt_idx] if len(parts) > gt_idx else "."
                gts_out.append(gt)
                # parse
                if len(gt) == 3 and (gt[1] == "/" or gt[1] == "|"):
                    a, b = gt[0], gt[2]
                    if a in "01" and b in "01":
                        an += 2
                        ac += (a == "1") + (b == "1")
                elif len(gt) == 1 and gt in "01":
                    # haploid (shouldn't happen on autosomes, but defensive)
                    an += 1
                    ac += (gt == "1")
                # else: missing, ignore
            if an == 0:
                continue
            af = ac / an
            maf = af if af <= 0.5 else 1.0 - af
            if maf < MAF_MIN:
                continue

            n_kept += 1
            out_line = "\t".join((
                f[0], f[1], f[2], ref, alt, f[5], f[6],
                f"AF_HGDP={af:.6g}", "GT", *gts_out,
            )) + "\n"
            fout.write(out_line)

            if n_kept % 50_000 == 0:
                dt = time.time() - t0
                print(f"  [chr{chrom}] kept={n_kept:,} bial_snv={n_pass:,} read={n_in:,} "
                      f"  {dt:.0f}s  ({n_in/dt:.0f} var/s)")
                sys.stdout.flush()

    tmp.replace(dst)
    dt = time.time() - t0
    out_size = dst.stat().st_size
    print(f"  [chr{chrom}] DONE: read={n_in:,}  biallelic_snv={n_pass:,}  "
          f"kept(MAF>={MAF_MIN})={n_kept:,}")
    print(f"  [chr{chrom}] samples={n_hgdp}  out={dst.name} ({out_size/1e6:.1f} MB)  "
          f"  wall={dt:.0f}s")
    return {"chrom": chrom, "n_in": n_in, "n_bial": n_pass, "n_kept": n_kept,
            "n_hgdp": n_hgdp, "out_bytes": out_size, "seconds": dt}


def main():
    if len(sys.argv) < 2:
        print("usage: python 29_subset_hgdp_vcf.py <chromN> [chromN ...]")
        sys.exit(2)
    for c in sys.argv[1:]:
        subset(c)


if __name__ == "__main__":
    main()
