"""
35_literature_cases.py - consolidated master list of published patients with
HOMOZYGOUS causative mutations in autosomal-recessive disease genes where the
surrounding region/run of homozygosity (ROH) size is reported, drawn from a
PubMed deep dive (2026-05-31).

For each case we record:
  source paper, case ID, gene, chromosome, GRCh38 (or marked GRCh37) midpoint,
  ROH length (Mb), background population class (outbred / consanguineous /
  population isolate / UPD / founder), and the local deCODE recombination rate
  at the locus.

Headline use cases:
  - Section "ROH descriptive statistics" of the manuscript: mean / SD / median /
    fraction sub-cutoff over the outbred + isolate subset.
  - Section "Locus-aware exemplars" of the manuscript: per-case verdicts under
    the calibrated scorer (script 25 / 26).
  - Anchor table for the Discussion's literature-meta paragraph.

The Najmabadi 2011 *Nature* cohort is processed separately by 33_najmabadi_
stratified.py (n=356 candidate variants joined to Table-S1 linkage intervals)
because it is large enough for an independent stratified statistical test.
This script holds the case-by-case literature dataset (Hildebrandt 2009 outbred
Table 1; Schuurs-Hoeijmakers 2011 outbred MR SROHs; Gomez-Diaz 2017 CAPN3
Mexican isolate; Bayes 2017 DDD PIGT autozygous family).

Outputs:
  literature_cases.tsv  one row per case with locus rate + classification
  literature_cases.txt  human-readable summary tables
"""

from math import comb
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).parent
DIV = HERE / "cross_pop_hap_diversity.tsv"
OUT_TSV = HERE / "literature_cases.tsv"
OUT_TXT = HERE / "literature_cases.txt"

# ---------------------------------------------------------------- master case list
# (source, case_id, gene, chrom, midpoint_bp, ROH_Mb, class, build, note)
# class: outbred / consanguineous_1st-cousin / isolate / UPD / founder / candidate
CASES = [
    # ---- Hildebrandt et al. 2009 PLoS Genet (PMID 19165332) Table 1, outbred
    ("Hildebrandt2009", "F30-1",  "NPHP4",        "chr1",   5_988_000,  2.10, "outbred",            "GRCh38", "chr1p36.31 telomere; parents related 10 generations back"),
    ("Hildebrandt2009", "A1730-1","NPHS2",        "chr1", 179_534_000,  2.70, "outbred",            "GRCh38", "chr1q25.2; sib pair"),
    ("Hildebrandt2009", "A1730-2","NPHS2",        "chr1", 179_534_000,  2.70, "outbred",            "GRCh38", "chr1q25.2; sib pair"),
    ("Hildebrandt2009", "F408",   "NPHP5/IQCB1",  "chr3", 121_750_000,  5.18, "outbred (Swiss founder)", "GRCh38", "chr3q13.33"),
    ("Hildebrandt2009", "F399-1", "NPHP5/IQCB1",  "chr3", 121_750_000,  8.25, "outbred",            "GRCh38", "chr3q13.33"),
    ("Hildebrandt2009", "A237",   "NPHS2 R138Q",  "chr1", 179_534_000,  2.30, "founder (European)", "GRCh38", "chr1q25.2 founder mutation"),
    ("Hildebrandt2009", "A825",   "NPHS2 R138Q",  "chr1", 179_534_000,  2.30, "founder (European)", "GRCh38", "chr1q25.2 founder mutation"),
    # ---- Schuurs-Hoeijmakers et al. 2011 EJHG (PMID 21248743) Supp Table 3, outbred MR (candidate loci)
    ("Schuurs2011",    "ARMR1",  "chr16 78-80Mb",        "chr16",  79_500_000,  2.6,  "outbred (candidate)", "GRCh37", "shared sib SROH"),
    ("Schuurs2011",    "ARMR1",  "chr19 39-50Mb MRT11",  "chr19",  44_000_000, 11.0,  "outbred (candidate)", "GRCh37", "MRT11 locus overlap"),
    ("Schuurs2011",    "ARMR4",  "chr4 32-35Mb",         "chr4",   33_500_000,  2.5,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR4",  "chr6 27-29Mb HLA",     "chr6",   28_000_000,  2.7,  "outbred (candidate)", "GRCh37", "HLA region"),
    ("Schuurs2011",    "ARMR4",  "chr11 48-50Mb",        "chr11",  49_200_000,  2.2,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR7",  "chr6 26-29Mb HLA",     "chr6",   27_500_000,  2.4,  "outbred (candidate)", "GRCh37", "HLA region"),
    ("Schuurs2011",    "ARMR7",  "chr7 117-120Mb",       "chr7",  118_600_000,  2.4,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR8",  "chr6 130-139Mb",       "chr6",  134_700_000,  8.4,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR8",  "chr9 131-134Mb",       "chr9",  132_600_000,  2.2,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR8",  "chr11 48-50Mb",        "chr11",  49_000_000,  2.4,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR9",  "chr6 62-65Mb",         "chr6",   63_300_000,  2.6,  "outbred (candidate)", "GRCh37", ""),
    ("Schuurs2011",    "ARMR9",  "chr8 47-50Mb",         "chr8",   48_500_000,  2.9,  "outbred (candidate)", "GRCh37", ""),
    # ---- Gomez-Diaz et al. 2017 PLoS One (PMID 28103310) Mexican LGMD2A isolate
    ("GomezDiaz2017",  "Tlaxcala isolate", "CAPN3",      "chr15",  42_400_000,  6.6,  "isolate (older founder)", "GRCh37", "Mexican Nahuatl + mestizo isolate; pooled DNA of 3 affecteds"),
    # ---- Bayes/Pagnamenta 2017 EJHG (PMID 28327575) DDD PIGT consanguineous family
    ("Bayes2017",      "PIGT fam2 270250", "PIGT",       "chr20",  44_100_000, 11.0,  "consanguineous_1st-cousin (Afghani)", "GRCh38", "chr20q13; F=1/15-1/19 estimated"),
    # ---- Prasad 2018 BMC Med Genet (PMID 29554876) gene loci only (no per-ROH lengths)
    ("Prasad2018",     "ROH04", "TYR",      "chr11",  89_000_000, None, "consanguineous (declared)", "GRCh38", "chr11q14.3"),
    ("Prasad2018",     "ROH22", "PCCB",     "chr3",  136_000_000, None, "outbred-declared",         "GRCh38", "chr3q22.3"),
    ("Prasad2018",     "ROH26", "SLC25A15", "chr13",  40_800_000, None, "outbred",                  "GRCh38", "chr13q14.11"),
    ("Prasad2018",     "ROH31", "NDUFV2",   "chr18",   9_100_000, None, "2nd-cousin",               "GRCh38", "chr18p11.22 subtelomere"),
    ("Prasad2018",     "ROH44", "TPP1",     "chr11",   6_600_000, None, "1st-cousin",               "GRCh38", "chr11p15.4 telomere"),
    ("Prasad2018",     "ROH52", "GJB2",     "chr13",  20_700_000, None, "outbred",                  "GRCh38", "chr13q12.11 DFNB1"),
    # ---- Pengelly 2018 Hum Mutat (PMID 29573052) DDD PIGH 1st-cousin Pakistani
    ("Pengelly2018",   "DDD 265247", "PIGH",       "chr14",  68_060_000, 25.2,  "consanguineous_1st-cousin (Pakistani)", "GRCh38", "chr14q24"),
]
PRIORS_TXT = ("(class abbreviations: outbred = no consanguinity declared; "
              "isolate = old founder event in a non-recently-consanguineous "
              "population; UPD = uniparental disomy; founder = known specific "
              "founder mutation; consanguineous = 1st-cousin or closer)")


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

    rows = []
    for src, cid, gene, ch, pos, L, klass, build, note in CASES:
        w = (pos // 1_000_000) * 1_000_000
        r = rate.get(ch, {}).get(w)
        rratio = (r / median) if r else None
        rows.append((src, cid, gene, ch, pos, L, klass, build, r, rratio, note))

    # write TSV
    tsv = ["source\tcase\tgene\tchrom\tmidpoint_bp\tROH_Mb\tclass\tbuild\t"
           "cMperMb\tratio_to_median\tnote"]
    for src, cid, gene, ch, pos, L, klass, build, r, rratio, note in rows:
        L_s = f"{L:.2f}" if L is not None else "NA"
        r_s = f"{r:.3f}" if r else "NA"
        rr_s = f"{rratio:.3f}" if rratio else "NA"
        tsv.append(f"{src}\t{cid}\t{gene}\t{ch}\t{pos}\t{L_s}\t{klass}\t{build}\t"
                   f"{r_s}\t{rr_s}\t{note}")
    OUT_TSV.write_text("\n".join(tsv) + "\n", encoding="utf-8")

    lines = [
        "# Consolidated literature cases for the ROH calibration manuscript",
        "# (excludes Najmabadi 2011 n=356; that cohort handled separately by 33_najmabadi_stratified.py)",
        f"# baseline: {100*frac15:.1f}% of 1 Mb windows at r>=1.5x median; "
        f"{100*frac20:.1f}% at >=2.0x. genome-wide median={median:.2f} cM/Mb.",
        f"# {PRIORS_TXT}\n",
    ]
    # full per-case table
    lines.append(f"{'source':16s} {'case':16s} {'gene':18s} {'chrom':6s} "
                 f"{'ROH_Mb':>7s} {'r(cM/Mb)':>9s} {'ratio':>6s}  class")
    for src, cid, gene, ch, pos, L, klass, build, r, rratio, note in rows:
        L_s = f"{L:6.2f}" if L is not None else "    NA"
        r_s = f"{r:>9.2f}" if r else "       NA"
        rr_s = f"{rratio:>6.2f}" if rratio else "    NA"
        lines.append(f"{src:16s} {cid[:16]:16s} {gene[:18]:18s} {ch:6s} "
                     f"{L_s} {r_s} {rr_s}  {klass}")

    # outbred / isolate subset descriptive stats (ROH lengths where known)
    outbred_classes = {"outbred", "outbred (Swiss founder)", "outbred (candidate)",
                       "founder (European)", "isolate (older founder)"}
    outbred_rows = [r_ for r_ in rows if r_[6] in outbred_classes and r_[5] is not None]
    if outbred_rows:
        Ls = [r_[5] for r_ in outbred_rows]
        rs_known = [r_[8] for r_ in outbred_rows if r_[8] is not None]
        n = len(Ls)
        lines.append(f"\n=== Outbred/isolate subset, n={n} cases with ROH length ===")
        lines.append(f"  mean ROH = {statistics.mean(Ls):.2f} Mb")
        lines.append(f"  SD       = {statistics.stdev(Ls):.2f}" if n > 1 else "")
        lines.append(f"  median   = {statistics.median(Ls):.2f} Mb")
        lines.append(f"  range    = [{min(Ls):.2f}, {max(Ls):.2f}] Mb")
        for thr in (3.0, 5.0, 10.0):
            k = sum(1 for L in Ls if L <= thr)
            lines.append(f"  ROH <= {thr:.0f} Mb: {k}/{n} ({100*k/n:.0f}%)")
        # hot-locus fraction
        if rs_known:
            h15 = sum(1 for r_ in rs_known if r_ >= 1.5 * median)
            h20 = sum(1 for r_ in rs_known if r_ >= 2.0 * median)
            nr = len(rs_known)
            p15 = binom_ge(h15, nr, frac15) if h15 else 1.0
            p20 = binom_ge(h20, nr, frac20) if h20 else 1.0
            lines.append(f"\n  locus enrichment (outbred/isolate subset only):")
            lines.append(f"  {h15}/{nr} = {100*h15/nr:.1f}% at r>=1.5x median "
                         f"(baseline {100*frac15:.1f}%, p={p15:.3g})")
            lines.append(f"  {h20}/{nr} = {100*h20/nr:.1f}% at r>=2.0x median "
                         f"(baseline {100*frac20:.1f}%, p={p20:.3g})")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
