# agents/design.py
from typing import List, Tuple
import random
from Bio.Align import substitution_matrices


class DesignAgent:
    def __init__(
        self, seed: int = 42, score_cutoff: int = 1, use_full_matrix: bool = False
    ):
        """
        seed: random seed for reproducibility
        score_cutoff: min BLOSUM62 score to consider substitution conservative
        use_full_matrix: if True, load full Biopython matrix, else fallback to PoC subset
        """
        random.seed(seed)
        self.score_cutoff = score_cutoff
        self.use_full_matrix = use_full_matrix

        if use_full_matrix:
            # Load official Biopython substitution matrix
            self.matrix = substitution_matrices.load("BLOSUM62")
        else:
            # PoC lightweight subset (faster to demo, simpler logic)
            self.matrix = {
                "A": {"V": 0, "L": -1, "I": -1, "G": 0},
                "V": {"A": 0, "L": 1, "I": 3, "M": 1},
                "L": {"I": 2, "V": 1, "M": 2},
                "I": {"L": 2, "V": 3, "M": 1},
                "K": {"R": 2, "E": 1},
                "R": {"K": 2, "H": 0},
                "E": {"D": 2, "Q": 2},
                "D": {"E": 2, "N": 1},
            }

    def _suggest_substitutions(self, residue: str) -> List[str]:
        """Suggest substitutions based on cutoff score."""
        subs = []
        if self.use_full_matrix:
            # Biopython matrix lookup
            for aa in self.matrix.alphabet:
                if aa == residue:
                    continue
                score = self.matrix[residue, aa]
                if score >= self.score_cutoff:
                    subs.append(aa)
        else:
            # PoC subset lookup
            if residue not in self.matrix:
                return []
            for aa, score in self.matrix[residue].items():
                if score >= self.score_cutoff:
                    subs.append(aa)
        return subs

    def run(
        self, sequence: str, context: List[str] | None = None, n_mutants: int = 5
    ) -> List[Tuple[str, str]]:
        """
        Propose candidate mutations.
        Returns list of (mutation_notation, mutated_sequence).
        Example: [("A10V", "MKT...V..."), ...]
        """
        seq = list(sequence)
        length = len(seq)
        variants: List[Tuple[str, str]] = []

        alphabet = "ACDEFGHIKLMNPQRSTVWY"

        for _ in range(n_mutants):
            pos = random.randint(0, length - 1)
            orig = seq[pos]

            subs = self._suggest_substitutions(orig)
            if subs:
                new_res = random.choice(subs)
            else:
                # fallback random pick
                choices = [aa for aa in alphabet if aa != orig]
                new_res = random.choice(choices)

            mutant = seq.copy()
            mutant[pos] = new_res
            mut_str = "".join(mutant)
            notation = f"{orig}{pos + 1}{new_res}"
            variants.append((notation, mut_str))

        return variants


# CLI test
if __name__ == "__main__":
    seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWQ"

    print("=== PoC mode ===")
    agent_poc = DesignAgent(use_full_matrix=False)
    for n, m in agent_poc.run(seq, n_mutants=5):
        print(n, m)

    print("\n=== Full BLOSUM62 mode ===")
    agent_full = DesignAgent(use_full_matrix=True)
    for n, m in agent_full.run(seq, n_mutants=5):
        print(n, m)
