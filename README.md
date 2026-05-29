# ROH calibration: a recombination- and diversity-aware weight of evidence for autozygosity

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20387797.svg)](https://doi.org/10.5281/zenodo.20387797)

Calibrating runs of homozygosity (ROH) to a **prior-free Bayes factor** (weight of
evidence) that an observed ROH reflects *recent autozygosity* rather than an
individual's population background — resolved per locus (recombination rate) and
per population (background haplotype structure). Proof-of-principle in public
reference data; the manuscript targets *Genetics in Medicine*.

## Headline results
- ROH length is a weight of evidence that rises ~linearly: `log10 BF(L) ≈ (λ/ln10)·L`,
  with `λ = 2r·ln(1/H̄)` per Mb (r = cM/Mb, H̄ = LD-aware background homozygosity).
- Three independent nulls (cross-individual pairs, simulated known-truth tracts,
  and a leakage-free trio-children background) converge on a **decisive length
  ~1.5–1.6 Mb** at a first-cousin prior; the closed-form block-independent score is
  **~2–3× over-confident** and is used only as a mechanistic ranking.
- The per-individual evidence law is **ancestry-robust** (~1.6 Mb across five
  superpopulations); the ancestry-dependent quantity is the **prevalence of
  elevated background autozygosity** (7% AFR → 68% AMR).
- The decisive length is **stable under resampling**: bootstrapping the trio
  children (`24_bootstrap_lstar.py`, B=1000) gives 95% CIs EUR [1.50, 1.60],
  AFR [1.60, 1.70], EAS [1.50, 1.65], SAS [1.55, 1.70], AMR [1.60, 1.80] Mb —
  overlapping across all five superpopulations.

## Threshold framing (ACMG-2021 vs lab practice)
The ACMG-2021 technical standard (Gonzales et al., *Genet Med* 2022;24(2):255–261)
recommends counting homozygous segments **>3 to 5 Mb** as likely IBD (Pemberton
2012; Hildebrandt 2009; Kearney 2011 — empirical length-distribution evidence).
Clinical laboratories operationally apply **5, 7, or 10 Mb** without consensus.
The "10 Mb" used as a comparison point in scripts 09/11/14/15/16/22 is that
**conventional lab-practice operating point**, not the ACMG recommendation. This
project supplies the per-locus quantitative resolution the standard's IBS-vs-IBD
rationale calls for.

## Pipeline (in `trio_chr22_NA12878/`)
| Script | Purpose |
|---|---|
| `15_cross_population.py` | RETIRED per-site 2pq noise term (documented foil) |
| `16_haplotype_ibs_noise.py` | LD-aware H̄ noise term (mechanistic score), genome-wide |
| `17_recomb_per_arm.py` | Per-arm recombination summary (deCODE maps) |
| `18_empirical_chance_ibs.py` | Empirical chance-IBS null (cross-individual pairs) + max-SNP-gap guard |
| `19_calibration_groundtruth.py` | Calibration vs simulated known-truth tracts (reliability, FDR) |
| `20_trio_roh_posterior.py` | Clinical trio tool: ROH calling + per-locus evidence + pedigree-F + deletion check |
| `21_trio_background_null.py` | Leakage-free trio-children background + cryptic-relatedness screen |
| `22_compare_fixed_thresholds.py` | A fixed length is not a fixed weight of evidence (per-locus span vs fixed cutoffs) |
| `24_bootstrap_lstar.py` | Bootstrap 95% CI on the decisive length L\* (resampling unit = trio child) |
| `25_score_roh.py` | **Reference implementation**: score any ROH per (locus, ancestry, platform, prior/F_ROH) → Bayes factor, posterior, FLAG/REVIEW/BACKGROUND |

Figure and manuscript generators (pure-stdlib HTML + inline SVG → headless-Edge
PDF/PNG) are released with the published manuscript, not in this code+data deposit.

## Data (NOT in this repository — re-downloadable)
- 1000 Genomes Project high-coverage release (3,202 samples, GRCh38), per-population AF.
- deCODE recombination maps: Palsson G, et al. *Complete human recombination maps.*
  Nature. 2025. doi:10.1038/s41586-024-08450-5.

Raw VCFs, the deCODE map archive, and `all_autosomes/` are excluded via
`.gitignore` (storage discipline: ~28 GB, trivially re-downloadable). Derived
look-up tables and summaries are versioned here; figures and the in-preparation
manuscript are released with publication.

## Reproduction (outline)
1. Download the 1000G high-coverage phased VCFs and the deCODE maps into
   `trio_chr22_NA12878/` (and `all_autosomes/`, `external/`).
2. Run scripts 15 → 22 (each prints usage; `chrom` args optional, default genome-wide).
3. `24_bootstrap_lstar.py` adds bootstrap CIs on the decisive length.
4. `25_score_roh.py` is the user-facing scorer: e.g.
   `python 25_score_roh.py --demo --ancestry=EUR`, or
   `python 25_score_roh.py my_roh.bed --ancestry=SAS --platform=array --froh=0.03`.
Pure Python + NumPy; Python 3.14.

## Provenance and AI assistance
Portions of the analysis code, figures, and drafting were produced with an AI
coding assistant under the author's direction and verification; the full, runnable
pipeline is provided here for audit and independent re-execution. See
`METHODOLOGY_LOG.md` for the dated decision record.

## License
MIT (see `LICENSE`).
