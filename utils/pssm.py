# utils/pssm.py
from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Tuple

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


def compute_pssm_from_msa(
    msa: List[str], alpha: float = 1.0
) -> Tuple[List[int], List[Dict[str, float]], List[Dict[str, float]]]:
    """
    Compute PSSM (log2 odds) from an aligned MSA (list of equal-length strings).
    Returns:
      positions: 1-based alignment positions
      probs: per-column smoothed probabilities dict {aa: p}
      scores: per-column log-odds dict {aa: log2(p/bg)}
    """
    if not msa:
        raise ValueError("Empty MSA")
    L = len(msa[0])
    for s in msa:
        if len(s) != L:
            raise ValueError("All MSA sequences must have equal length")

    positions = []
    probs_cols: List[Dict[str, float]] = []
    scores_cols: List[Dict[str, float]] = []

    print("Fuck you!!")

    for i in range(L):
        counts = _column_counts(msa, i)
        probs = _dirichlet_smooth(counts, alpha=alpha)
        scores = {aa: math.log2(probs[aa] / BLOSUM62_BG[aa]) for aa in AA}
        positions.append(i + 1)
        probs_cols.append(probs)
        scores_cols.append(scores)

    print("GTH!!!")

    # DEBUG: start
    if i in [33, 36, 44]:  # check columns around 34, 37, 45
        print("DEBUG checking column", i + 1)
        debug_msa_columns(msa, [i + 1])
    # DEBUG: end

    return positions, probs_cols, scores_cols


def debug_msa_columns(msa: list[str], positions: list[int]):
    """Print raw counts for selected alignment columns (1-based positions)."""
    for pos in positions:
        col = pos - 1  # convert to 0-based
        counts = _column_counts(msa, col)
        print(f"DEBUG Column {pos}: {counts}, total={sum(counts.values())}")
