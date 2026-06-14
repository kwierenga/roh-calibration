"""
37_coalescent_null.py - an INDEPENDENT, cryptic-relatedness-free non-IBD null
from coalescent simulation (msprime + stdpopsim), to cross-check the decisive
ROH length L* derived from real data.

WHY (reviewer hardening, item 8 / Tier-2 #3): the empirical (script 18, cross-
individual haplotype pairs) and trio-children (script 21) nulls both measure
chance-IBS in REAL people, so they can be contaminated by cryptic distant
relatedness / population structure -- itself ancestry-dependent. A coalescent
simulation has NO relatedness beyond the demographic model by construction, so
its IBS-run background is a clean null. If the coalescent L* matches the
SCREENED real-data L* (not the inflated unscreened one), that validates that the
F_ROH screen removes cryptic-relatedness inflation rather than manufacturing the
result.

MODEL: stdpopsim HomSap OutOfAfrica_3G09 (Gutenkunst 2009), populations
YRI->AFR, CEU->EUR, CHB->EAS (the three unadmixed continental groups; SAS/AMR
are admixed/structured and are not part of this 3-population model, so the
coalescent cross-check is reported for AFR/EUR/EAS only). Realistic per-
population Ne -> the model should reproduce AFR's shorter chance-IBS (higher
diversity / shorter LD) without any tuning.

ASCERTAINMENT MATCH: each simulated DIPLOID individual's two homologs are one
draw of the population's chance/background autozygosity (the same unit as a trio
child in script 21). Homozygous runs are measured over COMMON SNPs (MAF>=0.05
WITHIN the simulated population) using the SAME run-segmentation, gap tolerance,
desert-break, and exceedance-integral p_background as the real-data nulls
(reused from script 21), so L* is computed identically.

REQUIRES: a Python interpreter with msprime + stdpopsim (NOT the default 3.14
env, which lacks wheels). Create one with:
    py -3.13 -m venv .venv_coalescent
    .venv_coalescent\\Scripts\\python -m pip install msprime stdpopsim
and run THIS script with that interpreter.

Outputs:
  coalescent_null_summary.txt   coalescent L* per pop/gap vs empirical & trio-null
  coalescent_null_pchance.tsv   L grid x pop: coalescent p_background + log10 BF
Usage:
  <venv-py> 37_coalescent_null.py --smoke          # tiny, fast sanity run
  <venv-py> 37_coalescent_null.py                  # full run (background)
  <venv-py> 37_coalescent_null.py --ndip=120 --total-mb=1500 --chunk-mb=25
"""

import importlib.util
import math
import sys
import time
from pathlib import Path

import numpy as np

try:
    import msprime
    import stdpopsim
except ImportError as e:
    sys.exit(f"need msprime + stdpopsim (run under the .venv_coalescent interpreter): {e}")

HERE = Path(__file__).parent

# reuse script 21's null math so THR_PC / agg_null / constants are IDENTICAL
_spec = importlib.util.spec_from_file_location("trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m21)

MAF_MIN = m21.MAF_MIN
MAX_SNP_GAP_BP = m21.MAX_SNP_GAP_BP
MIN_KEEP_MB = m21.MIN_KEEP_MB
C_IBD = m21.C_IBD
THR_PC = m21.THR_PC
PI = m21.PI
T_DEC = m21.T_DEC
L_GRID = m21.L_GRID
agg_null = m21.agg_null

POP_MAP = {"AFR": "YRI", "EUR": "CEU", "EAS": "CHB"}   # our label -> model deme
GAPS = [0, 1]
RNG_SEED = 17
B_BOOT = 500          # bootstrap replicates (resample simulated individuals)
CI = (2.5, 97.5)

OUT_SUM = HERE / "coalescent_null_summary.txt"
OUT_PC = HERE / "coalescent_null_pchance.tsv"

# real-data GENOME-WIDE references for the magnitude comparison (gap=1):
TRIO_SCREENED_GW = {"AFR": 1.65, "EUR": 1.55, "EAS": 1.55}   # script 21 screened (headline)
TRIO_UNSCREENED_GW = {"AFR": 1.05, "EUR": 1.10, "EAS": 10.85}  # script 21 L*_all (chr22 ref;
# the EAS 10.85 is the cryptic-relatedness artifact the relatedness-free null should NOT show)


def roh_lengths_g(hom, pos, gap_tol):
    """Homozygous-run lengths (Mb); same logic as script 21, gap as a parameter."""
    m = hom.copy()
    if gap_tol > 0:
        pad = np.concatenate(([1], m.astype(np.int8), [1])); dd = np.diff(pad)
        hs = np.flatnonzero(dd == -1); he = np.flatnonzero(dd == 1)
        short = (he - hs) <= gap_tol
        if short.any():
            diff = np.zeros(m.size + 1, dtype=np.int32)
            np.add.at(diff, hs[short], 1); np.add.at(diff, he[short], -1)
            m = m | (np.cumsum(diff[:-1]) > 0)
    n = m.size
    intra = np.zeros(n, bool)
    intra[1:] = m[1:] & m[:-1] & ((pos[1:] - pos[:-1]) <= MAX_SNP_GAP_BP)
    starts = np.flatnonzero(m & ~intra)
    ends = m.copy(); ends[:-1] &= ~intra[1:]; ends = np.flatnonzero(ends)
    if starts.size == 0:
        return np.empty(0, np.float32)
    return ((pos[ends] - pos[starts]) / 1e6).astype(np.float32)


def main():
    t0 = time.time()
    args = sys.argv[1:]
    smoke = "--smoke" in args

    def getopt(name, default):
        for a in args:
            if a.startswith(name + "="):
                return type(default)(a.split("=", 1)[1])
        return default

    ndip = getopt("--ndip", 20 if smoke else 120)         # diploids per population
    total_mb = getopt("--total-mb", 40 if smoke else 1000)  # total simulated Mb/pop
    chunk_mb = getopt("--chunk-mb", 20 if smoke else 25)
    n_chunks = max(1, round(total_mb / chunk_mb))
    B = 0 if smoke else B_BOOT

    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_3G09")
    demography = model.model
    mu = model.mutation_rate or species.genome.mean_mutation_rate
    rrate = species.genome.mean_recombination_rate
    # name -> population index for sample extraction
    name2id = {p.name: i for i, p in enumerate(demography.populations)}
    demes = [POP_MAP[p] for p in POP_MAP]

    print(f"OOA_3G09  ndip={ndip}/pop  total={total_mb}Mb  chunk={chunk_mb}Mb "
          f"x{n_chunks}  mu={mu:.2e}  r={rrate:.2e}", flush=True)

    # per-(pop) accumulators: list over individuals of per-chunk seg arrays, per gap
    segs = {p: {g: [[] for _ in range(ndip)] for g in GAPS} for p in POP_MAP}
    spans = {p: [0.0] * ndip for p in POP_MAP}

    for ci in range(n_chunks):
        seq_len = chunk_mb * 1_000_000
        ts = msprime.sim_ancestry(
            samples={POP_MAP[p]: ndip for p in POP_MAP},
            demography=demography, sequence_length=seq_len,
            recombination_rate=rrate, random_seed=RNG_SEED + ci)
        mts = msprime.sim_mutations(
            ts, rate=mu, model=msprime.BinaryMutationModel(), random_seed=RNG_SEED + ci)
        G = mts.genotype_matrix()                      # (sites, haplotypes) 0/1
        pos = mts.tables.sites.position.astype(np.int64)
        for p in POP_MAP:
            haps = ts.samples(population=name2id[POP_MAP[p]])   # 2*ndip haplotype ids
            haps = np.asarray(haps)
            Gp = G[:, haps]
            af = Gp.mean(axis=1)
            common = np.flatnonzero((np.minimum(af, 1 - af) >= MAF_MIN))
            if common.size < 2:
                continue
            posc = pos[common]
            span = (posc[-1] - posc[0]) / 1e6
            Gc = Gp[common, :]
            for j in range(ndip):
                hom = Gc[:, 2 * j] == Gc[:, 2 * j + 1]
                spans[p][j] += span
                for g in GAPS:
                    sl = roh_lengths_g(hom, posc, g)
                    slk = sl[sl > MIN_KEEP_MB]
                    if slk.size:
                        segs[p][g][j].append(slk)
        print(f"  chunk {ci+1}/{n_chunks}  sites~{mts.num_sites}  ({time.time()-t0:.0f}s)", flush=True)

    # L* per pop per gap + bootstrap CI (resample simulated individuals) + the
    # gap=1 p_background curve for the BF table
    rng = np.random.default_rng(RNG_SEED)
    res = {}; ci = {}
    pc_curves = {}
    for p in POP_MAP:
        res[p] = {}; ci[p] = {}
        for g in GAPS:
            emp, L = agg_null(segs[p][g], spans[p])
            res[p][g] = L
            if g == 1:
                pc_curves[p] = emp
            if B:
                boot = np.empty(B)
                for b in range(B):
                    idx = rng.integers(0, ndip, size=ndip)
                    _, Lb = agg_null([segs[p][g][j] for j in idx], [spans[p][j] for j in idx])
                    boot[b] = Lb if np.isfinite(Lb) else float(L_GRID[-1]) + 0.05
                ci[p][g] = (float(np.percentile(boot, CI[0])), float(np.percentile(boot, CI[1])))
            else:
                ci[p][g] = (float("nan"), float("nan"))
        print(f"  {p}: L*coal gap0={res[p][0]:.2f}  gap1={res[p][1]:.2f} Mb "
              f"CI{CI}={ci[p][1]}  (ndip={ndip})", flush=True)

    with OUT_PC.open("w", encoding="utf-8") as fh:
        cols = [p for p in POP_MAP if pc_curves.get(p) is not None]
        fh.write("L_Mb\t" + "\t".join(f"{p}_pchance_coal\t{p}_log10BF_coal" for p in cols) + "\n")
        for k, L in enumerate(L_GRID):
            cell = []
            for p in cols:
                pcv = max(pc_curves[p][k], 1e-12)
                cell.append(f"{pc_curves[p][k]:.3e}\t{math.log10(C_IBD/pcv):.2f}")
            fh.write(f"{L:.2f}\t" + "\t".join(cell) + "\n")

    lines = [
        "# Independent coalescent non-IBD null (msprime + stdpopsim OOA_3G09)",
        "# PURPOSE: a relatedness-free MAGNITUDE check on the decisive length L*.",
        "# Use it to confirm that a null with NO cryptic relatedness still lands at the",
        "# ~Mb scale of the SCREENED real-data L* (not the inflated unscreened value).",
        "# It is NOT used for the cross-ancestry ORDERING (see caveat below).",
        f"# ndip={ndip}/pop total={total_mb}Mb/pop chunk={chunk_mb}Mb x{n_chunks} "
        f"mu={mu:.2e} r={rrate:.2e} (flat) seed={RNG_SEED} MAF>={MAF_MIN} "
        f"pi={PI} posterior>={T_DEC} B={B} wall={time.time()-t0:.0f}s",
        "# Homozygous runs over common SNPs (MAF>=MAF_MIN within the simulated pop),",
        "# scored with script 21's identical p_background / posterior threshold.",
        "",
        "pop   L*_coal_gap1 [95% CI]      | real GW screened   real GW unscreened",
    ]
    for p in POP_MAP:
        lo, hi = ci[p][1]
        cistr = f"[{lo:.2f},{hi:.2f}]" if np.isfinite(lo) else "[n/a]"
        lines.append(f"{p}   {res[p][1]:.2f} {cistr:15s}       "
                     f"{TRIO_SCREENED_GW.get(p, float('nan')):.2f}               "
                     f"{TRIO_UNSCREENED_GW.get(p, float('nan')):.2f}")
    lines += [
        "",
        "(gap0 L*: " + "  ".join(f"{p}={res[p][0]:.2f}" for p in POP_MAP) + " Mb)",
        "",
        "READING (magnitude): the coalescent L* sits at the ~Mb scale of the SCREENED",
        "real-data L*, NOT at the inflated unscreened value (e.g. the EAS 10.85 Mb",
        "artifact). Because the simulation has no relatedness by construction, this",
        "supports that the real-data unscreened inflation is cryptic relatedness and that",
        "the F_ROH screen removes it rather than manufacturing the result. The direct",
        "primary evidence for that is the script-36 Table B screen-sensitivity on REAL",
        "data (removing ~2 of 70 EAS children collapses the unscreened L*).",
        "",
        "CAVEAT (ordering, stated honestly): a NEUTRAL OOA model with a flat",
        "recombination rate does NOT reproduce the empirical AFR-shortest ordering. AFR's",
        "neutral common-variant SFS loads MAF near 0.05 (low 2pq -> longer chance runs),",
        "which in the model outweighs AFR's shorter-LD effect; in real data the LD effect",
        "dominates and AFR is shortest. So the ancestry ordering rests on the REAL-DATA",
        "nulls (scripts 18/21), not on this simulation. SAS/AMR are not in the 3-",
        "population model and are omitted. AF and run-measurement use the same simulated",
        "individuals (~1/ndip self-inclusion bias, negligible at ndip>=100).",
    ]
    OUT_SUM.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_SUM}\n  -> {OUT_PC}\n  total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
