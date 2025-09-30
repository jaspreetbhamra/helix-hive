# utils/conservation.py
import math
from collections import Counter
from typing import Dict, List


def compute_conservation(msa: List[str]) -> List[Dict]:
    """
    Compute per-position conservation from an MSA.
    Adds: entropy, majority freq, occupancy (non-gap fraction), non-gap depth.
    """
    n_seqs = len(msa)
    aln_len = len(msa[0])
    results: List[Dict] = []

    for i in range(aln_len):
        column = [seq[i] for seq in msa]
        residues = [aa for aa in column if aa != "-"]
        counts = Counter(residues)
        total = sum(counts.values())

        if total > 0:
            freqs = {aa: c / total for aa, c in counts.items()}
            entropy = -sum(p * math.log2(p) for p in freqs.values())
            majority = max(freqs.values())
        else:
            freqs, entropy, majority = {}, 0.0, 0.0

        occupancy = total / max(1, n_seqs)

        results.append(
            {
                "position": i + 1,  # 1-based in alignment space
                "frequencies": freqs,
                "entropy": float(entropy),
                "majority": float(majority),
                "occupancy": float(occupancy),
                "depth_non_gap": int(total),
            }
        )
    return results
