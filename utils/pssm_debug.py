import matplotlib.pyplot as plt
import numpy as np

from utils.pssm import effective_msa_depth, sequence_weights


def analyze_msa(msa_file: str, ident_thresh: float = 0.8):
    # Load FASTA sequences
    sequences = []
    with open(msa_file, "r", encoding="utf-8") as f:
        seq = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            sequences.append(seq)

    N = len(sequences)
    weights = sequence_weights(sequences, threshold=ident_thresh)
    Meff = effective_msa_depth(weights)

    print(f"[INFO] MSA size (N)={N}, Effective depth (Meff)={Meff:.2f}")
    print(f"[INFO] Diversity ratio Meff/N={Meff / N:.2f}")

    return N, Meff


def plot_meff_vs_threshold(msa_file: str):
    thresholds = np.linspace(0.5, 0.95, 10)  # test a range of cutoffs
    Ns, Meffs, ratios = [], [], []

    for t in thresholds:
        N, Meff = analyze_msa(msa_file, ident_thresh=t)
        Ns.append(N)
        Meffs.append(Meff)
        ratios.append(Meff / N)

    plt.figure(figsize=(7, 5))
    plt.plot(thresholds, ratios, marker="o")
    plt.title("Meff/N vs Sequence Identity Threshold")
    plt.xlabel("Identity threshold")
    plt.ylabel("Effective diversity (Meff/N)")
    plt.ylim(0, 1.05)
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    msa_file = "data/cache/lysozyme_alignment.fasta"
    plot_meff_vs_threshold(msa_file)
