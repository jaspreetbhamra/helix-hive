# utils/pssm.py
from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Dict, List

import numpy as np

# Canonical amino-acid order (20aa)
AA = "ACDEFGHIKLMNPQRSTVWY"

# Industry-standard BLOSUM62 background frequencies (Henikoff & Henikoff, 1992)
BLOSUM62_BG = {
    "A": 0.074,
    "R": 0.052,
    "N": 0.045,
    "D": 0.054,
    "C": 0.025,
    "Q": 0.034,
    "E": 0.054,
    "G": 0.074,
    "H": 0.026,
    "I": 0.068,
    "L": 0.099,
    "K": 0.058,
    "M": 0.025,
    "F": 0.047,
    "P": 0.039,
    "S": 0.057,
    "T": 0.051,
    "W": 0.013,
    "Y": 0.032,
    "V": 0.073,
}


def _column_counts(msa: List[str], col: int) -> Counter:
    """Count non-gap residues in a given column of the MSA."""
    return Counter([seq[col] for seq in msa if seq[col] != "-"])


def _dirichlet_smooth(counts: Counter, alpha: float = 1.0) -> Dict[str, float]:
    """
    Dirichlet smoothing with background priors:
      p_i = (c_i + alpha * BG_i) / (N + alpha)
    where N = sum c_i ; alpha ~ 1 is a safe default.
    """
    N = sum(counts.values())
    denom = N + alpha
    probs: Dict[str, float] = {}
    for aa in AA:
        c = counts.get(aa, 0)
        probs[aa] = (
            (c + alpha * BLOSUM62_BG[aa]) / denom if denom > 0 else BLOSUM62_BG[aa]
        )
    return probs


# def sequence_weights(msa: List[str], threshold: float = 0.8) -> np.ndarray:
#     """
#     Compute sequence weights to reduce redundancy in MSA.
#     Simple scheme: each sequence gets weight = 1 / (number of sequences >= threshold similarity).
#     """
#     n = len(msa)
#     L = len(msa[0])
#     msa_arr = np.array([list(seq) for seq in msa])
#     weights = np.ones(n)

#     for i in range(n):
#         sim_count = 0
#         for j in range(n):
#             if i == j:
#                 continue
#             matches = np.sum(msa_arr[i] == msa_arr[j])
#             ident = matches / L
#             if ident >= threshold:
#                 sim_count += 1
#         weights[i] = 1.0 / (1 + sim_count)

#     # Normalize weights so they sum to n
#     weights *= n / np.sum(weights)
#     return weights


def sequence_identity(seq1: str, seq2: str) -> float:
    """Fraction identity ignoring gaps."""
    matches, length = 0, 0
    for a, b in zip(seq1, seq2, strict=False):
        if a == "-" or b == "-":
            continue
        length += 1
        if a == b:
            matches += 1
    return matches / length if length > 0 else 0.0


def sequence_weights(msa: List[str], threshold: float = 0.8) -> np.ndarray:
    """
    Compute redundancy-aware sequence weights.
    Each sequence weight = 1 / (# sequences with identity >= threshold).
    Returns: array of weights (length N).
    """
    N = len(msa)
    weights = np.zeros(N, dtype=float)

    for i in range(N):
        sim_count = 0
        for j in range(N):
            if sequence_identity(msa[i], msa[j]) >= threshold:
                sim_count += 1
        weights[i] = 1.0 / sim_count if sim_count > 0 else 0.0
        print(f"[DEBUG] Seq {i}: sim_count={sim_count}, weight={weights[i]:.3f}")

    Meff = float(np.sum(weights))
    print(f"[DEBUG] Effective depth (Meff)={Meff:.3f}")
    return weights


def adaptive_alpha(alpha: float, Meff: float) -> float:
    """
    Scale Dirichlet pseudocounts based on effective depth.
    - Large Meff → trust MSA more → keep alpha small.
    - Small Meff → increase alpha (stronger prior).
    """
    if Meff > 50:
        return alpha
    else:
        return alpha * (10.0 / (Meff + 1e-6))


def _weighted_column_counts(
    msa: List[str], col: int, weights: np.ndarray
) -> Dict[str, float]:
    counts = defaultdict(float)
    for seq, w in zip(msa, weights, strict=False):
        res = seq[col]
        if res != "-":
            counts[res] += w
    return counts


def effective_msa_depth(weights: np.ndarray) -> float:
    return float(np.sum(weights))


def compute_pssm_from_msa(
    msa: List[str], alpha: float = 1.0, ident_thresh: float = 0.8
):
    if not msa:
        raise ValueError("Empty MSA")
    L = len(msa[0])
    for s in msa:
        if len(s) != L:
            raise ValueError("All MSA sequences must have equal length")

    # Compute sequence weights
    weights = sequence_weights(msa, threshold=ident_thresh)
    Meff = effective_msa_depth(weights)

    positions, probs_cols, scores_cols = [], [], []

    for i in range(L):
        counts = _weighted_column_counts(msa, i, weights)

        # scale alpha automatically
        adj_alpha = adaptive_alpha(alpha, Meff)

        probs = _dirichlet_smooth(counts, alpha=adj_alpha)
        scores = {aa: math.log2(probs[aa] / BLOSUM62_BG[aa]) for aa in AA}
        positions.append(i + 1)
        probs_cols.append(probs)
        scores_cols.append(scores)

    print(f"[DEBUG] Final Meff={Meff:.2f}, using alpha={adj_alpha:.3f}")
    return positions, probs_cols, scores_cols


def debug_msa_columns(msa: list[str], positions: list[int]):
    """Print raw counts for selected alignment columns (1-based positions)."""
    for pos in positions:
        col = pos - 1  # convert to 0-based
        counts = _column_counts(msa, col)
        print(f"DEBUG Column {pos}: {counts}, total={sum(counts.values())}")
