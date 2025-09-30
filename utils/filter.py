# utils/filter.py
import re
from typing import List, Tuple

from Bio.Align import PairwiseAligner


def deduplicate_sequences(
    seqs: List[Tuple[str, str]], identity_thresh: float = 0.95
) -> List[Tuple[str, str]]:
    """
    Remove near-duplicate sequences (> identity_thresh) using pairwise identity.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    filtered = []
    for sid, seq in seqs:
        is_dup = False
        for _, kept in filtered:
            matches = aligner.score(seq, kept)
            ident = matches / max(len(seq), len(kept))
            if ident >= identity_thresh:
                is_dup = True
                break
        if not is_dup:
            filtered.append((sid, seq))
    return filtered


def balance_taxa(
    seqs: List[Tuple[str, str]],
    max_per_taxon: int = 2,
) -> List[Tuple[str, str]]:
    """
    Balance sequences by taxon inferred from ID (prefix before '_').
    Example IDs: LYSC_HUMAN, LYSC_MOUSE, ...
    """
    groups = {}
    for sid, seq in seqs:
        m = re.match(r"^[^_]+_([A-Za-z0-9]+)", sid)
        taxon = m.group(1) if m else "UNK"
        groups.setdefault(taxon, []).append((sid, seq))

    balanced = []
    for taxon, items in groups.items():
        balanced.extend(items[:max_per_taxon])  # truncate
    return balanced
