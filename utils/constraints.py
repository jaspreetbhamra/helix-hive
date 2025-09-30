# utils/constraints.py
import re
from typing import Dict, List, Literal

import numpy as np

HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}


def _choose_conserved_sites(
    conservation: List[Dict],
    mode: Literal["threshold", "percentile", "topk"] = "percentile",
    entropy_thresh: float | None = None,
    percentile: float = 5.0,  # lowest-entropy X%
    topk: int | None = None,
    min_majority: float = 0.9,  # >= 90% same residue among non-gaps
    min_occupancy: float = 0.8,  # >= 80% non-gap occupancy in column
) -> List[int]:
    cols = [
        c
        for c in conservation
        if c.get("occupancy", 0.0) >= min_occupancy
        and c.get("majority", 0.0) >= min_majority
    ]
    if not cols:
        return []

    if mode == "threshold":
        if entropy_thresh is None:
            entropy_thresh = 0.1
        kept = [c for c in cols if c["entropy"] <= entropy_thresh]
        kept.sort(key=lambda x: (x["entropy"], -x["majority"]))
        return [c["position"] for c in kept]

    ent = np.array([c["entropy"] for c in cols], dtype=float)
    if mode == "percentile":
        cut = np.percentile(ent, percentile)
        kept = [c for c in cols if c["entropy"] <= cut]
        kept.sort(key=lambda x: (x["entropy"], -x["majority"]))
        return [c["position"] for c in kept]
    elif mode == "topk":
        if topk is None:
            topk = 10
        cols.sort(key=lambda x: (x["entropy"], -x["majority"]))
        return [c["position"] for c in cols[:topk]]
    else:
        return []


def extract_constraints(
    sequence: str,
    conservation: List[Dict],
    *,
    conserved_mode: str = "percentile",
    entropy_thresh: float | None = None,
    percentile: float = 5.0,
    topk: int | None = None,
    min_majority: float = 0.9,
    min_occupancy: float = 0.8,
) -> Dict:
    """
    Extract constraints; conserved_sites now configurable to avoid huge lists.
    """
    constraints = {
        "conserved_sites": _choose_conserved_sites(
            conservation,
            mode=conserved_mode,
            entropy_thresh=entropy_thresh,
            percentile=percentile,
            topk=topk,
            min_majority=min_majority,
            min_occupancy=min_occupancy,
        ),
        "cysteines": [],
        "glycosylation_sites": [],
        "signal_peptide": False,
        "transmembrane_spans": [],
    }

    # cysteines
    constraints["cysteines"] = [
        i for i, aa in enumerate(sequence, start=1) if aa == "C"
    ]

    # N-glycosylation motifs (N-X-[ST], X != P)
    for m in re.finditer(r"N[^P][ST]", sequence):
        start = m.start() + 1
        constraints["glycosylation_sites"].append((start, start + 2))

    # signal peptide heuristic (first 25 aa)
    window = sequence[:25]
    if window:
        avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in window) / len(window)
        if avg_hydro > 1.5:
            constraints["signal_peptide"] = True

    # simple TM scan (window 19, avg hydropathy > 1.6)
    w = 19
    for i in range(len(sequence) - w + 1):
        win = sequence[i : i + w]
        avg_hydro = sum(HYDROPATHY.get(aa, 0.0) for aa in win) / w
        if avg_hydro > 1.6:
            constraints["transmembrane_spans"].append((i + 1, i + w))

    return constraints
