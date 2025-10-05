# utils/pssm_utils.py
import json
from typing import Dict


def load_pssm(path: str) -> Dict:
    with open(path) as f:
        return json.load(f)


def mutation_score(pssm: Dict, mutation: str) -> float:
    """
    Compute ΔPSSM score for a mutation like "E35Q".
    Uses 1-based positions from the PSSM.
    Returns: new_score - old_score (log-odds difference, bits).
    """
    import re

    m = re.match(r"^([A-Z])(\d+)([A-Z])$", mutation)
    if not m:
        raise ValueError(f"Invalid mutation format: {mutation} (expected e.g. E35Q)")
    wt, pos_str, mut = m.groups()
    pos = int(pos_str)

    try:
        idx = pssm["positions"].index(pos)
    except ValueError:
        raise KeyError(f"Position {pos} not found in PSSM")

    scores = pssm["pssm"][idx]
    old_score = scores.get(wt)
    new_score = scores.get(mut)
    if old_score is None or new_score is None:
        raise KeyError(f"Residue {wt}->{mut} not in alphabet")

    return new_score - old_score


# Usage
# from utils.pssm_utils import load_pssm, mutation_score

# pssm = load_pssm("data/cache/lysozyme_pssm.json")
# delta = mutation_score(pssm, "E35Q")
# print(f"E35Q ΔPSSM = {delta:.2f} bits")
