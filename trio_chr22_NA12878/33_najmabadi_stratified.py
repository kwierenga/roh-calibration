"""
33_najmabadi_stratified.py - locus-rate enrichment of homozygosity-mapped causative-
candidate variants from Najmabadi et al. Nature 2011 (PMID 21937992, doi 10.1038/
nature10423), stratified by the length of the homozygous linkage interval that
contains each variant.

CENTRAL FINDING: when stratified by linkage-interval length, the locus-rate
distribution of causal-candidate variants from a consanguineous AR-intellectual-
disability cohort behaves exactly as the locus-aware ROH calibration predicts:

  - LONG  intervals (> 10 Mb) -> locus distribution indistinguishable from a
    random genome-wide sample (~19% at r >= 1.5x median, vs 25% baseline; NS).
  - SHORT intervals (<= 5 Mb) -> ~54% of variants sit at r >= 1.5x median
    (5.5x10^-6) and ~37% at r >= 2x (2.3x10^-6). Strong enrichment at hot loci.

Mechanism (biology): generations of recombination trim the ancestral ROH down
faster at high-recombination loci than at cold loci, so the short causative-ROH
cases that still localize are disproportionately at hot loci. The locus-aware
scorer (script 25) exploits exactly this in its calibrated weight of evidence.

Sources (in `docs/`):
  - Table S1 (linkage intervals): 41586_2011_BFnature10423_MOESM219_ESM.xls
    columns: Family ID | ... | Chr | Start marker | End marker |
             Length [Mbp] | LOD score | Degree of Consanguinity | ...
  - Table S2 (mutations):         41586_2011_BFnature10423_MOESM220_ESM.xls
    columns: Family ID | Base Change (chr:pos) | Protein Change | ...

Method: join Table S1 with Table S2 on (Family ID, Chromosome). For each
mutation, take the matching family+chrom interval's Length [Mbp] as the
linkage-interval length. Where multiple intervals exist for a (family, chrom)
pair (54 of 412 ~ 13%), take the smallest (conservative). Drop chrX/chrY.

Outputs:
  najmabadi_stratified.tsv  per-variant fam, gene, chrom, pos, ROH_Mb, cMperMb
  najmabadi_stratified.txt  stratified enrichment summary
"""

from math import comb
from pathlib import Path
import re
import statistics
import sys

import xlrd

HERE = Path(__file__).parent
DOCS = HERE.parent / "docs"
DIV = HERE / "cross_pop_hap_diversity.tsv"
S1 = DOCS / "41586_2011_BFnature10423_MOESM219_ESM.xls"
S2 = DOCS / "41586_2011_BFnature10423_MOESM220_ESM.xls"
OUT_TSV = HERE / "najmabadi_stratified.tsv"
OUT_TXT = HERE / "najmabadi_stratified.txt"

STRATA = [
    ("ROH <= 3 Mb",     lambda L: L <= 3),
    ("ROH <= 5 Mb",     lambda L: L <= 5),
    ("ROH 5-10 Mb",     lambda L: 5 < L <= 10),
    ("ROH > 10 Mb",     lambda L: L > 10),
]


def load_rates():
    rate = {}
    with DIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        seen = set()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            k = (f[ix["chrom"]], int(f[ix["window_start"]]))
            if k in seen:
                continue
            seen.add(k)
            try:
                r = float(f[ix["cMperMb"]])
            except (ValueError, IndexError):
                continue
            if r > 0:
                rate.setdefault(f[ix["chrom"]], {})[k[1]] = r
    return rate


def load_intervals():
    """Table S1 -> {(family, chrom_int): [Length_Mbp, ...]}."""
    wb = xlrd.open_workbook(str(S1))
    s = wb.sheet_by_index(0)
    intervals = {}
    last_fam = ""
    for r in range(2, s.nrows):
        fam = str(s.cell_value(r, 0)).strip()
        if fam:
            last_fam = fam
        try:
            ch = int(float(s.cell_value(r, 7)))
        except (ValueError, TypeError):
            continue
        try:
            L = float(s.cell_value(r, 10))
        except (ValueError, TypeError):
            continue
        intervals.setdefault((last_fam, ch), []).append(L)
    return intervals


def load_variants():
    """Table S2 -> list of (family, chrom_str, pos_int, gene)."""
    wb = xlrd.open_workbook(str(S2))
    s = wb.sheet_by_index(0)
    out = []
    last_fam = ""
    for r in range(2, s.nrows):
        fam = str(s.cell_value(r, 0)).strip()
        if fam:
            last_fam = fam
        bc = str(s.cell_value(r, 1)).strip()
        pc = str(s.cell_value(r, 2)).strip()
        m = re.match(r"chr([0-9XY]+):([0-9,]+)", bc)
        if not m:
            continue
        ch = m.group(1)
        if ch in ("X", "Y"):
            continue
        gene = pc.split(":")[0].strip()
        out.append((last_fam, ch, int(m.group(2).replace(",", "")), gene))
    return out


def binom_ge(k, n, p):
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    rate = load_rates()
    all_r = sorted(r for ws in rate.values() for r in ws.values())
    median = statistics.median(all_r)
    n_w = len(all_r)
    frac15 = sum(1 for r in all_r if r >= 1.5 * median) / n_w
    frac20 = sum(1 for r in all_r if r >= 2.0 * median) / n_w

    intervals = load_intervals()
    variants = load_variants()

    # Join: per variant, look up smallest interval for (family, chrom_int)
    joined = []
    multi = missed = 0
    for fam, ch, pos, gene in variants:
        ch_int = int(ch)
        ivs = intervals.get((fam, ch_int))
        if not ivs:
            missed += 1
            continue
        if len(ivs) > 1:
            multi += 1
        L = min(ivs)
        chrom_str = f"chr{ch}"
        w = (pos // 1_000_000) * 1_000_000
        r = rate.get(chrom_str, {}).get(w)
        if r is None:
            continue
        joined.append((fam, gene, chrom_str, pos, L, r))

    rows_tsv = ["family\tgene\tchrom\tpos\tROH_length_Mb\tcMperMb\tratio_to_median"]
    for fam, gene, ch, pos, L, r in joined:
        rows_tsv.append(f"{fam}\t{gene}\t{ch}\t{pos}\t{L:.3f}\t{r:.3f}\t{r/median:.3f}")
    OUT_TSV.write_text("\n".join(rows_tsv) + "\n", encoding="utf-8")

    lines = [
        "# Najmabadi 2011 Nature stratified locus-rate analysis",
        f"# n_variants_joined={len(joined)}  (multi-interval={multi}, no-S1-match={missed})",
        f"# genome-wide median rate={median:.2f} cM/Mb  n_windows={n_w}",
        f"# baseline: {100*frac15:.1f}% windows at r>=1.5x; {100*frac20:.1f}% at >=2x",
        "",
        f"{'stratum':16s} {'n':>4s}  "
        f"{'hot>=1.5x':>10s} {'%':>5s} {'p_one_sided':>13s}  "
        f"{'hot>=2.0x':>10s} {'%':>5s} {'p_one_sided':>13s}",
    ]
    for name, pred in STRATA:
        sub = [x for x in joined if pred(x[4])]
        n = len(sub)
        if n == 0:
            lines.append(f"{name:16s} {n:>4d}  (no cases)")
            continue
        h15 = sum(1 for x in sub if x[5] >= 1.5 * median)
        h20 = sum(1 for x in sub if x[5] >= 2.0 * median)
        p15 = binom_ge(h15, n, frac15) if h15 else 1.0
        p20 = binom_ge(h20, n, frac20) if h20 else 1.0
        lines.append(f"{name:16s} {n:>4d}  "
                     f"{h15:>10d} {100*h15/n:>4.1f}%  {p15:>11.4g}    "
                     f"{h20:>10d} {100*h20/n:>4.1f}%  {p20:>11.4g}")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
