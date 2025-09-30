# tools/compare_conservation.py
import json
import argparse
import matplotlib.pyplot as plt


def load_conservation(path):
    with open(path) as f:
        return json.load(f)


def load_constraints(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def main():
    parser = argparse.ArgumentParser("Compare dedup vs raw conservation")
    parser.add_argument(
        "--dedup", required=True, help="Conservation JSON with deduplication"
    )
    parser.add_argument(
        "--raw", required=True, help="Conservation JSON without deduplication"
    )
    parser.add_argument("--dedup-constraints", help="Constraints JSON with dedup")
    parser.add_argument("--raw-constraints", help="Constraints JSON without dedup")
    parser.add_argument("--out", default="data/cache/conservation_compare.png")
    args = parser.parse_args()

    cons_dedup = load_conservation(args.dedup)
    cons_raw = load_conservation(args.raw)
    c_dedup = load_constraints(args.dedup_constraints) if args.dedup_constraints else {}
    c_raw = load_constraints(args.raw_constraints) if args.raw_constraints else {}

    # extract entropy
    pos = [c["position"] for c in cons_dedup]
    ent_dedup = [c["entropy"] for c in cons_dedup]
    ent_raw = [c["entropy"] for c in cons_raw]

    fig, ax = plt.subplots(figsize=(max(10, len(pos) // 3), 4))

    ax.plot(pos, ent_dedup, label="Dedup + balanced", color="blue", linewidth=1.5)
    ax.plot(
        pos, ent_raw, label="Raw (no filtering)", color="red", linewidth=1.5, alpha=0.7
    )

    # Overlay conserved sites (optional)
    if "conserved_sites" in c_dedup:
        ax.scatter(
            c_dedup["conserved_sites"],
            [0] * len(c_dedup["conserved_sites"]),
            color="blue",
            marker="o",
            s=20,
            label="Conserved (dedup)",
        )
    if "conserved_sites" in c_raw:
        ax.scatter(
            c_raw["conserved_sites"],
            [0.1] * len(c_raw["conserved_sites"]),
            color="red",
            marker="x",
            s=20,
            label="Conserved (raw)",
        )

    ax.set_xlabel("Position")
    ax.set_ylabel("Entropy")
    ax.set_title("Conservation comparison: dedup vs raw")
    ax.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"[Compare] Figure saved to {args.out}")


if __name__ == "__main__":
    main()
