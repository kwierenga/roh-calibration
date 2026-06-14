"""
32_locus_rate_enrichment.py - test whether the loci of published short-ROH
recessive cases cluster at high-recombination regions of the genome, as
predicted by the locus-aware calibration (evidence ~ r x L).

For each published causative or candidate AR locus where the ROH is short
(<=10 Mb) we look up the sex-averaged deCODE/Palsson 2024 cM/Mb in the
overlapping 1 Mb window(s) and compare the rate distribution to the
genome-wide baseline. Test: are these case loci over-represented at
r >= k x median?

Sources:
  - Hildebrandt 2009 Table 1, 5 outbred solved cases (NPHP4, NPHS2 x2,
    NPHP5/IQCB1 x2). PMID 19165332.
  - Prasad 2018 Table 1, 6 homozygous pathogenic variants (GJB2, TPP1,
    SLC25A15, TYR, PCCB, NDUFV2). PMID 29554876.
  - Schuurs-Hoeijmakers 2011 Supp Table 3, 12 per-family SROHs from 10
    outbred ID families (these are SHARED-ROH candidate loci, not
    confirmed causative). PMID 21248743.
  - Schuurs-Hoeijmakers 2011 Supp Table 2, 5 SROHs that overlap previously
    reported MRT loci (MRT7-MRT11).

Build note: Schuurs coordinates are GRCh37/hg19 (2011 era). The deCODE map
is GRCh38. Per-1-Mb-window shifts between builds are typically <1 window;
we use the published midpoints as approximate GRCh38 positions and flag
this in the output. For Hildebrandt and Prasad we use GRCh38 gene midpoints.

Outputs:
  locus_rate_enrichment.tsv  per-case locus, rate, percentile
  locus_rate_enrichment.txt  summary table + binomial tests
"""

from math import comb
from pathlib import Path
import statistics

HERE = Path(__file__).parent
DIV = HERE / "cross_pop_hap_diversity.tsv"
OUT_TSV = HERE / "locus_rate_enrichment.tsv"
OUT_TXT = HERE / "locus_rate_enrichment.txt"

# (source, case_id, gene_or_locus, chrom, midpoint_bp, ROH_Mb, build_note)
CASES = [
    # ----- Hildebrandt 2009, outbred solved (Table 1) -----
    ("Hildebrandt2009", "F30-1",   "NPHP4",        "chr1",   5_988_000, 2.10, "GRCh38 gene midpoint"),
    ("Hildebrandt2009", "A1730",   "NPHS2",        "chr1", 179_534_000, 2.70, "GRCh38 gene midpoint"),
    ("Hildebrandt2009", "F399-1",  "NPHP5/IQCB1",  "chr3", 121_750_000, 8.25, "GRCh38 gene midpoint"),
    ("Hildebrandt2009", "F408",    "NPHP5/IQCB1",  "chr3", 121_750_000, 5.18, "GRCh38 gene midpoint"),
    ("Hildebrandt2009", "A237/A825", "NPHS2 R138Q founder", "chr1", 179_534_000, 2.30, "GRCh38 gene midpoint"),
    # ----- Prasad 2018, homozygous pathogenic (6 genes named) -----
    ("Prasad2018", "ROH31",  "NDUFV2",    "chr18",   9_100_000, None, "GRCh38 gene midpoint"),
    ("Prasad2018", "—",      "GJB2",      "chr13",  20_700_000, None, "GRCh38 gene midpoint"),
    ("Prasad2018", "—",      "TPP1",      "chr11",   6_600_000, None, "GRCh38 gene midpoint"),
    ("Prasad2018", "—",      "SLC25A15",  "chr13",  40_800_000, None, "GRCh38 gene midpoint"),
    ("Prasad2018", "—",      "TYR",       "chr11",  89_000_000, None, "GRCh38 gene midpoint"),
    ("Prasad2018", "—",      "PCCB",      "chr3",  136_000_000, None, "GRCh38 gene midpoint"),
    # ----- Schuurs 2011 Supp Table 3, per-family SROHs (10 outbred ID) -----
    # using SROH midpoint, GRCh37 -> approx GRCh38 (~no liftover; <1 Mb shift)
    ("Schuurs2011-T3", "ARMR1",  "chr16 78-80Mb",       "chr16",  79_500_000,  2.6, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR1",  "chr19 39-50Mb (MRT11)", "chr19", 44_000_000, 11.0, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR4",  "chr4 32-35Mb",        "chr4",   33_500_000,  2.5, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR4",  "chr6 27-29Mb (HLA)",  "chr6",   28_000_000,  2.7, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR4",  "chr11 48-50Mb",       "chr11",  49_200_000,  2.2, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR7",  "chr6 26-29Mb (HLA)",  "chr6",   27_500_000,  2.4, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR7",  "chr7 117-120Mb",      "chr7",  118_600_000,  2.4, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR8",  "chr6 130-139Mb",      "chr6",  134_700_000,  8.4, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR8",  "chr9 131-134Mb",      "chr9",  132_600_000,  2.2, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR8",  "chr11 48-50Mb",       "chr11",  49_000_000,  2.4, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR9",  "chr6 62-65Mb",        "chr6",   63_300_000,  2.6, "GRCh37 SROH midpoint"),
    ("Schuurs2011-T3", "ARMR9",  "chr8 47-50Mb",        "chr8",   48_500_000,  2.9, "GRCh37 SROH midpoint"),
]


def load_rates():
    rate = {}
    with DIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        seen = set()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[ix["chrom"]], int(f[ix["window_start"]]))
            if key in seen:
                continue
            seen.add(key)
            try:
                r = float(f[ix["cMperMb"]])
            except (ValueError, IndexError):
                continue
            rate.setdefault(f[ix["chrom"]], {})[key[1]] = r
    return rate


def binom_ge(k, n, p):
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def main():
    rate = load_rates()
    # genome-wide distribution
    all_r = sorted(r for ws in rate.values() for r in ws.values() if r > 0)
    median = statistics.median(all_r)
    n_w = len(all_r)
    thr15 = 1.5 * median
    thr20 = 2.0 * median
    frac15 = sum(1 for r in all_r if r >= thr15) / n_w
    frac20 = sum(1 for r in all_r if r >= thr20) / n_w

    # score each case
    rows = []
    for src, cid, locus, ch, pos, roh, note in CASES:
        w = (pos // 1_000_000) * 1_000_000
        r = rate.get(ch, {}).get(w, float("nan"))
        ratio = r / median if r > 0 else float("nan")
        pct = (sum(1 for rr in all_r if rr <= r) / n_w * 100) if r > 0 else float("nan")
        rows.append((src, cid, locus, ch, pos, roh, r, ratio, pct, note))

    # tests, by data source
    sources = ["Hildebrandt2009", "Prasad2018", "Schuurs2011-T3"]
    tally = {s: [0, 0, 0, 0] for s in sources}   # n_total, n_hot_1.5x, n_hot_2x, n_short<=5Mb
    for src, cid, locus, ch, pos, roh, r, ratio, pct, note in rows:
        if src not in tally:
            continue
        tally[src][0] += 1
        if r >= thr15: tally[src][1] += 1
        if r >= thr20: tally[src][2] += 1
        if roh is not None and roh <= 5.0: tally[src][3] += 1

    # combined
    n_all = sum(t[0] for t in tally.values())
    n15_all = sum(t[1] for t in tally.values())
    n20_all = sum(t[2] for t in tally.values())
    p15 = binom_ge(n15_all, n_all, frac15)
    p20 = binom_ge(n20_all, n_all, frac20)

    # output
    lines = [
        "# Locus-rate enrichment test: do published short-ROH AR cases cluster at",
        "# high-recombination loci, as predicted by the locus-aware calibration?",
        f"# deCODE/Palsson 2024 sex-averaged rate map; n_windows={n_w}; "
        f"genome-wide median = {median:.2f} cM/Mb.",
        f"# baseline: {100*frac15:.1f}% of 1 Mb windows have r >= 1.5x median ({thr15:.2f}); "
        f"{100*frac20:.1f}% at >= 2.0x ({thr20:.2f}).",
        "\nPer-source tally (case loci at r >= k x genome-wide median):\n",
        f"{'source':18s} {'n':>3s} {'n_hot >=1.5x':>13s} {'n_hot >=2.0x':>13s}  {'n_ROH<=5Mb':>10s}",
    ]
    for s in sources:
        t = tally[s]
        lines.append(f"{s:18s} {t[0]:>3d} {t[1]:>13d} {t[2]:>13d}  {t[3]:>10d}")
    lines.append(f"{'COMBINED':18s} {n_all:>3d} {n15_all:>13d} {n20_all:>13d}")

    lines.append("\nBinomial tests (one-sided, vs genome-wide baseline):")
    lines.append(f"  P(X >= {n15_all} | n={n_all}, p={frac15:.3f}) at r >= 1.5x median = {p15:.4f}")
    lines.append(f"  P(X >= {n20_all} | n={n_all}, p={frac20:.3f}) at r >= 2.0x median = {p20:.4f}")

    lines.append("\nPer-locus detail:\n")
    lines.append(f"{'source':18s} {'case':12s} {'locus':28s} {'chrom':6s} "
                 f"{'midpoint':>11s} {'ROH(Mb)':>8s} {'r(cM/Mb)':>9s} "
                 f"{'ratio':>6s} {'pct':>6s}  note")
    for src, cid, locus, ch, pos, roh, r, ratio, pct, note in rows:
        roh_s = f"{roh:.2f}" if roh is not None else "  NA"
        lines.append(f"{src:18s} {cid:12s} {locus:28s} {ch:6s} "
                     f"{pos:>11d} {roh_s:>8s} {r:>9.2f} "
                     f"{ratio:>6.2f} {pct:>6.1f}  {note}")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    # also TSV
    tsv = ["source\tcase\tlocus\tchrom\tmidpoint_bp\tROH_Mb\tcMperMb\tratio_to_median\tpercentile\tbuild_note"]
    for src, cid, locus, ch, pos, roh, r, ratio, pct, note in rows:
        tsv.append(f"{src}\t{cid}\t{locus}\t{ch}\t{pos}\t"
                   f"{roh if roh is not None else 'NA'}\t{r:.4f}\t"
                   f"{ratio:.3f}\t{pct:.2f}\t{note}")
    OUT_TSV.write_text("\n".join(tsv) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
