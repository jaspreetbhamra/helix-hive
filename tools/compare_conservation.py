# tools/compare_conservation.py
import argparse
import json

import matplotlib.pyplot as plt


def load_json(path):
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

    cons_dedup = load_json(args.dedup)
    cons_raw = load_json(args.raw)
    c_dedup = load_json(args.dedup_constraints) if args.dedup_constraints else {}
    c_raw = load_json(args.raw_constraints) if args.raw_constraints else {}

    # Extract entropy curves
    pos = [c["position"] for c in cons_dedup]
    ent_dedup = [c["entropy"] for c in cons_dedup]
    ent_raw = [c["entropy"] for c in cons_raw]

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(max(10, len(pos) // 3), 4))
    ax.plot(pos, ent_dedup, label="Dedup + balanced", color="blue", linewidth=1.5)
    ax.plot(
        pos, ent_raw, label="Raw (no filtering)", color="red", linewidth=1.5, alpha=0.7
    )

    # Overlay conserved sites
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

    # --- Numerical Summary ---
    n_cons_dedup = len(c_dedup.get("conserved_sites", []))
    n_cons_raw = len(c_raw.get("conserved_sites", []))
    cyst_dedup = len(c_dedup.get("cysteines", []))
    cyst_raw = len(c_raw.get("cysteines", []))

    print("\n=== Numerical Summary ===")
    print(f"Conserved sites (dedup+balance): {n_cons_dedup}")
    print(f"Conserved sites (raw):           {n_cons_raw}")
    print(f"Cysteines (dedup):               {cyst_dedup}")
    print(f"Cysteines (raw):                 {cyst_raw}")

    if n_cons_raw > n_cons_dedup:
        print(
            f"→ Filtering reduced apparent conservation by {n_cons_raw - n_cons_dedup} positions"
        )
    elif n_cons_raw < n_cons_dedup:
        print(
            f"→ Deduplication unexpectedly increased conserved calls by {n_cons_dedup - n_cons_raw}"
        )
    else:
        print("→ Number of conserved sites unchanged (rare, but possible)")


if __name__ == "__main__":
    main()
