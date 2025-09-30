# utils/visualize.py
import json
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np


def plot_conservation_logo(
    json_path: str,
    out_png: str = "data/cache/conservation_logo.png",
    top_n: int = 5,
    constraints_json: str | None = None,
):
    """
    Plot a sequence logo (residue frequencies), Shannon entropy, and optional constraints.
    Args:
        json_path: path to conservation.json (from analyze_msa)
        out_png: output image path
        top_n: number of most frequent residues to plot per position
        constraints_json: optional path to constraints.json
    """
    with open(json_path) as f:
        data: List[Dict] = json.load(f)

    positions = [d["position"] for d in data]
    n_pos = len(positions)

    # Collect residues and frequencies
    residue_set = set()
    for d in data:
        residue_set.update(d["frequencies"].keys())
    residues = sorted(residue_set)

    freq_matrix = {aa: [d["frequencies"].get(aa, 0.0) for d in data] for aa in residues}
    entropy = [d["entropy"] for d in data]

    # Sort residues by mean frequency, take top_n for readability
    mean_freqs = {aa: np.mean(vals) for aa, vals in freq_matrix.items()}
    sorted_residues = sorted(mean_freqs, key=mean_freqs.get, reverse=True)[:top_n]

    # Colors for top residues
    colors = plt.cm.tab20(np.linspace(0, 1, len(sorted_residues)))

    # --- Plot ---
    fig, (ax_logo, ax_entropy) = plt.subplots(
        2,
        1,
        figsize=(max(10, n_pos // 3), 6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    # Stacked bar chart for logo
    bottom = np.zeros(n_pos)
    for aa, color in zip(sorted_residues, colors, strict=False):
        vals = freq_matrix[aa]
        ax_logo.bar(positions, vals, bottom=bottom, label=aa, color=color, width=1.0)
        bottom += vals

    ax_logo.set_ylabel("Frequency")
    ax_logo.set_title("Sequence Logo (top residues)")
    ax_logo.legend(ncol=min(10, len(sorted_residues)), fontsize="small")

    # Entropy line plot
    ax_entropy.plot(positions, entropy, color="black", linewidth=1.2)
    ax_entropy.set_ylabel("Entropy")
    ax_entropy.set_xlabel("Position")
    ax_entropy.set_title("Shannon Entropy + Constraints")

    # Overlay constraints if available
    if constraints_json:
        with open(constraints_json) as f:
            constraints = json.load(f)

        # conserved sites
        if constraints.get("conserved_sites"):
            ax_entropy.scatter(
                constraints["conserved_sites"],
                [0] * len(constraints["conserved_sites"]),
                color="red",
                marker="o",
                s=40,
                label="Conserved",
            )
        # cysteines
        if constraints.get("cysteines"):
            ax_entropy.scatter(
                constraints["cysteines"],
                [0.1] * len(constraints["cysteines"]),
                color="orange",
                marker="^",
                s=50,
                label="Cysteines",
            )
        # glycosylation sites
        for start, end in constraints.get("glycosylation_sites", []):
            ax_entropy.axvspan(
                start, end, color="green", alpha=0.3, label="N-glycosylation"
            )
        # signal peptide
        if constraints.get("signal_peptide"):
            ax_entropy.axvspan(1, 25, color="blue", alpha=0.2, label="Signal peptide")
        # TM spans
        for start, end in constraints.get("transmembrane_spans", []):
            ax_entropy.axvspan(start, end, color="purple", alpha=0.2, label="TM span")

        ax_entropy.legend(fontsize="x-small", ncol=3)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"[Visualizer] Conservation + constraints plot saved to {out_png}")
