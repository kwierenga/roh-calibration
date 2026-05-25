"""
HMM-based crossover detection. The principled replacement for the heuristic
long-run detector in 03_detect_crossovers.py and 05_robust_detection.py.

Model
-----
Hidden state at each Mendelian-deterministic site (father het + mother hom, or
the symmetric maternal case): which parental haplotype was transmitted to the
child at that site. Two states: 0 or 1.

Observations: at each site, the child's parental-inherited allele (0 or 1
where 0=REF, 1=ALT) and the parent's two alleles. Under convention A,
child[0] = paternal-inherited allele. The observation tells us, conditional
on the parent's two alleles at this site, which parental haplotype was
transmitted (up to genotyping error epsilon).

Emission probabilities
- E(observe | state s) = 1 - eps  if the observed child allele matches the
                                  parent's allele at haplotype-s.
                       = eps      otherwise (genotyping/phasing error).
  Default eps = 0.001.

Transition probabilities
- Between consecutive deterministic sites at positions p_i, p_{i+1}, the
  genetic distance is integrated from the local recombination rate map.
  If a deCODE-style map is provided: integrate the sex-specific rate.
  If not: use a chromosome-average sex-specific rate (chr22: 0.7 cM/Mb male,
                                                       1.2 cM/Mb female).
- Haldane's formula: P(switch) = (1 - exp(-2 * D_morgans)) / 2.

Inference
- Forward-backward in log space for per-site posterior P(state | data).
- Viterbi for the most likely state path.
- Crossovers = transitions in the Viterbi path. The posterior probability
  of each transition is reported as a confidence score.

Output
- per-site posterior TSV (for plotting + downstream uses)
- crossover candidate TSV with posterior support per call
"""

import gzip
import math
import sys
from pathlib import Path

HERE = Path(__file__).parent
IN_TSV = HERE / "trio_chr22.tsv.gz"
OUT_POSTERIOR = HERE / "hmm_posterior_chr22.tsv.gz"
OUT_CROSSOVERS = HERE / "hmm_crossovers_chr22.tsv"
OUT_SUMMARY = HERE / "hmm_summary_chr22.txt"

# parameters
GENOTYPING_ERROR = 0.001
PATERNAL_AVG_CM_PER_MB_CHR22 = 0.7
MATERNAL_AVG_CM_PER_MB_CHR22 = 1.2
LOG_HALF = math.log(0.5)


def parse_gt(gt):
    if "|" not in gt:
        return None
    a, b = gt.split("|", 1)
    try:
        return int(a), int(b)
    except ValueError:
        return None


def load_deterministic_sites(path, parent_side):
    """
    Load Mendelian-deterministic sites only.
    Returns list of (pos, parent_alleles_tuple, observed_child_allele).
    """
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            father = parse_gt(f[idx["FATHER"]])
            mother = parse_gt(f[idx["MOTHER"]])
            child = parse_gt(f[idx["CHILD"]])
            if None in (father, mother, child):
                continue
            if parent_side == "paternal":
                parent = father
                other = mother
                child_idx = 0
            else:
                parent = mother
                other = father
                child_idx = 1
            if parent[0] == parent[1] or other[0] != other[1]:
                continue
            rows.append((int(f[idx["POS"]]), parent, child[child_idx]))
    return rows


def emission_logprob(observed, parent_alleles, eps=GENOTYPING_ERROR):
    """log P(observed allele | state s) for s in {0, 1}."""
    log_correct = math.log(1.0 - eps)
    log_error = math.log(eps)
    e0 = log_correct if observed == parent_alleles[0] else log_error
    e1 = log_correct if observed == parent_alleles[1] else log_error
    return e0, e1


def transition_logprob(distance_bp, cm_per_mb):
    """
    Probability of haplotype switch between two consecutive sites,
    Haldane's mapping function.
    """
    morgans = distance_bp * (cm_per_mb / 1e8)
    # Haldane: P(switch) = (1 - exp(-2*M)) / 2
    p_switch = 0.5 * (1.0 - math.exp(-2.0 * morgans))
    # clamp away from 0 and 1
    p_switch = max(min(p_switch, 0.5 - 1e-12), 1e-15)
    p_stay = 1.0 - p_switch
    return math.log(p_stay), math.log(p_switch)


def hmm_forward_backward(sites, cm_per_mb):
    """
    Run forward-backward on the HMM. Returns:
      log_fwd, log_bwd, log_post, log_evidence
    All in log space. log_fwd[i] = log P(obs_1..i, state_i = s); shape (N, 2).
    log_post[i][s] = log P(state_i = s | obs_1..N).
    """
    n = len(sites)
    if n == 0:
        return [], [], [], -math.inf

    # initial state distribution: uniform
    log_fwd = [[0.0, 0.0] for _ in range(n)]
    e0, e1 = emission_logprob(sites[0][2], sites[0][1])
    log_fwd[0][0] = LOG_HALF + e0
    log_fwd[0][1] = LOG_HALF + e1

    for i in range(1, n):
        prev_pos = sites[i - 1][0]
        cur_pos = sites[i][0]
        d = cur_pos - prev_pos
        log_stay, log_switch = transition_logprob(d, cm_per_mb)
        e0, e1 = emission_logprob(sites[i][2], sites[i][1])
        for s in (0, 1):
            # log P(state_i = s | obs_1..i) up to constant
            from_0 = log_fwd[i - 1][0] + (log_stay if s == 0 else log_switch)
            from_1 = log_fwd[i - 1][1] + (log_switch if s == 0 else log_stay)
            log_fwd[i][s] = log_sum_exp(from_0, from_1) + (e0 if s == 0 else e1)

    log_evidence = log_sum_exp(log_fwd[-1][0], log_fwd[-1][1])

    # backward
    log_bwd = [[0.0, 0.0] for _ in range(n)]
    log_bwd[n - 1][0] = 0.0
    log_bwd[n - 1][1] = 0.0
    for i in range(n - 2, -1, -1):
        next_pos = sites[i + 1][0]
        cur_pos = sites[i][0]
        d = next_pos - cur_pos
        log_stay, log_switch = transition_logprob(d, cm_per_mb)
        e0_next, e1_next = emission_logprob(sites[i + 1][2], sites[i + 1][1])
        for s in (0, 1):
            to_0 = (log_stay if s == 0 else log_switch) + e0_next + log_bwd[i + 1][0]
            to_1 = (log_switch if s == 0 else log_stay) + e1_next + log_bwd[i + 1][1]
            log_bwd[i][s] = log_sum_exp(to_0, to_1)

    # posterior
    log_post = [[0.0, 0.0] for _ in range(n)]
    for i in range(n):
        unnorm = [log_fwd[i][s] + log_bwd[i][s] for s in (0, 1)]
        norm = log_sum_exp(unnorm[0], unnorm[1])
        log_post[i][0] = unnorm[0] - norm
        log_post[i][1] = unnorm[1] - norm

    return log_fwd, log_bwd, log_post, log_evidence


def hmm_viterbi(sites, cm_per_mb):
    """Most likely state sequence. Returns list of states (0/1) of length n."""
    n = len(sites)
    if n == 0:
        return []
    log_delta = [[0.0, 0.0] for _ in range(n)]
    psi = [[0, 0] for _ in range(n)]
    e0, e1 = emission_logprob(sites[0][2], sites[0][1])
    log_delta[0][0] = LOG_HALF + e0
    log_delta[0][1] = LOG_HALF + e1
    for i in range(1, n):
        d = sites[i][0] - sites[i - 1][0]
        log_stay, log_switch = transition_logprob(d, cm_per_mb)
        e0, e1 = emission_logprob(sites[i][2], sites[i][1])
        for s in (0, 1):
            from_0 = log_delta[i - 1][0] + (log_stay if s == 0 else log_switch)
            from_1 = log_delta[i - 1][1] + (log_switch if s == 0 else log_stay)
            if from_0 >= from_1:
                log_delta[i][s] = from_0 + (e0 if s == 0 else e1)
                psi[i][s] = 0
            else:
                log_delta[i][s] = from_1 + (e0 if s == 0 else e1)
                psi[i][s] = 1
    # backtrack
    states = [0] * n
    states[n - 1] = 0 if log_delta[n - 1][0] >= log_delta[n - 1][1] else 1
    for i in range(n - 2, -1, -1):
        states[i] = psi[i + 1][states[i + 1]]
    return states


def log_sum_exp(a, b):
    if a > b:
        return a + math.log1p(math.exp(b - a))
    else:
        return b + math.log1p(math.exp(a - b))


def call_crossovers(sites, viterbi_states, log_post, min_posterior=0.95):
    """
    Crossovers = positions where the Viterbi path switches.
    Each crossover gets a confidence = posterior probability of the new state
    at the first site of the new run.
    """
    crossovers = []
    for i in range(1, len(viterbi_states)):
        if viterbi_states[i] != viterbi_states[i - 1]:
            prev_pos = sites[i - 1][0]
            new_pos = sites[i][0]
            new_state = viterbi_states[i]
            post = math.exp(log_post[i][new_state])
            if post >= min_posterior:
                crossovers.append({
                    "prev_pos": prev_pos,
                    "new_pos": new_pos,
                    "gap_bp": new_pos - prev_pos,
                    "prev_state": viterbi_states[i - 1],
                    "new_state": new_state,
                    "posterior": post,
                })
    return crossovers


def run_one_side(side, cm_per_mb):
    print(f"\n  === {side.upper()} ===")
    sites = load_deterministic_sites(IN_TSV, side)
    print(f"    Mendelian-deterministic sites: {len(sites):,}")
    if len(sites) < 100:
        print("    (too few sites; aborting)")
        return [], [], []
    log_fwd, log_bwd, log_post, log_evidence = hmm_forward_backward(sites, cm_per_mb)
    viterbi = hmm_viterbi(sites, cm_per_mb)
    crossovers = call_crossovers(sites, viterbi, log_post)
    print(f"    log evidence: {log_evidence:.2f}")
    print(f"    Viterbi-path switches: {sum(1 for i in range(1, len(viterbi)) if viterbi[i] != viterbi[i-1])}")
    print(f"    high-confidence crossovers (post >= 0.95): {len(crossovers)}")
    for c in crossovers:
        print(f"      chr22:{c['prev_pos']:,} <-> chr22:{c['new_pos']:,}  "
              f"({c['gap_bp']:,} bp gap; {c['prev_state']}->{c['new_state']}; posterior {c['posterior']:.4f})")
    return sites, log_post, crossovers


def main():
    if not IN_TSV.exists():
        sys.exit(f"input not found: {IN_TSV}. Run 02_extract_trio.py first.")

    print(f"  HMM crossover detection on chr22, using chromosome-average")
    print(f"  recombination rates (deCODE map integration TODO).")
    print(f"  paternal rate: {PATERNAL_AVG_CM_PER_MB_CHR22} cM/Mb")
    print(f"  maternal rate: {MATERNAL_AVG_CM_PER_MB_CHR22} cM/Mb")
    print(f"  genotyping error eps: {GENOTYPING_ERROR}")

    pat_sites, pat_post, pat_cos = run_one_side("paternal", PATERNAL_AVG_CM_PER_MB_CHR22)
    mat_sites, mat_post, mat_cos = run_one_side("maternal", MATERNAL_AVG_CM_PER_MB_CHR22)

    # write posterior
    with gzip.open(OUT_POSTERIOR, "wt") as fh:
        fh.write("parent\tpos\tP_state_0\tP_state_1\n")
        for (sites, post, label) in [(pat_sites, pat_post, "paternal"),
                                       (mat_sites, mat_post, "maternal")]:
            for site, lp in zip(sites, post):
                fh.write(f"{label}\t{site[0]}\t{math.exp(lp[0]):.6f}\t{math.exp(lp[1]):.6f}\n")

    # write crossovers
    with OUT_CROSSOVERS.open("w") as fh:
        fh.write("parent\tprev_pos\tnew_pos\tgap_bp\tprev_state\tnew_state\tposterior\n")
        for label, cos in [("paternal", pat_cos), ("maternal", mat_cos)]:
            for c in cos:
                fh.write(f"{label}\t{c['prev_pos']}\t{c['new_pos']}\t{c['gap_bp']}\t{c['prev_state']}\t{c['new_state']}\t{c['posterior']:.6f}\n")

    # summary
    with OUT_SUMMARY.open("w") as fh:
        fh.write("HMM crossover detection summary, chr22 NA12878\n")
        fh.write("=" * 50 + "\n\n")
        fh.write(f"Paternal: {len(pat_sites):,} deterministic sites, {len(pat_cos)} high-confidence crossovers\n")
        for c in pat_cos:
            fh.write(f"  chr22:{c['prev_pos']:,}-{c['new_pos']:,}  (post={c['posterior']:.4f})\n")
        fh.write(f"\nMaternal: {len(mat_sites):,} deterministic sites, {len(mat_cos)} high-confidence crossovers\n")
        for c in mat_cos:
            fh.write(f"  chr22:{c['prev_pos']:,}-{c['new_pos']:,}  (post={c['posterior']:.4f})\n")

    print(f"\n  -> {OUT_POSTERIOR}")
    print(f"  -> {OUT_CROSSOVERS}")
    print(f"  -> {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
