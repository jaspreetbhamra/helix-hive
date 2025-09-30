from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
from Bio import AlignIO
from Bio.Align import PairwiseAligner  # biopython

from models.esm import embed_sequence_mean
from utils.conservation import compute_conservation
from utils.constraints import extract_constraints
from utils.filter import balance_taxa, deduplicate_sequences
from utils.io import read_fasta, validate_protein_sequence
from utils.visualize import plot_conservation_logo


@dataclass
class RetrievalHit:
    id: str
    sequence: str
    similarity: float
    identity: float | None = None


class RetrieverAgent:
    def __init__(self, index_path: str | Path = "data/cache/retriever_index.npz"):
        self.index_path = Path(index_path)
        self._ids: np.ndarray | None = None
        self._seqs: np.ndarray | None = None
        self._embs: np.ndarray | None = None

    # ---------- Index ----------
    def build_index(
        self, fasta_path: str | Path, metadata_out: str | Path | None = None
    ) -> None:
        fasta_path = Path(fasta_path)
        items = read_fasta(fasta_path)
        if not items:
            raise ValueError(f"No sequences found in {fasta_path}")

        ids, seqs, embs = [], [], []
        for sid, seq in items:
            seq = validate_protein_sequence(seq)
            vec = embed_sequence_mean(seq)
            ids.append(sid)
            seqs.append(seq)
            embs.append(vec.astype(np.float32))

        ids_np = np.array(ids, dtype=object)
        seqs_np = np.array(seqs, dtype=object)
        embs_np = np.stack(embs, axis=0).astype(np.float32)

        self.index_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(self.index_path, ids=ids_np, seqs=seqs_np, embs=embs_np)

        meta = {
            "source_fasta": str(fasta_path),
            "n_sequences": len(ids),
            "embedding_dim": int(embs_np.shape[1]),
        }
        if metadata_out is None:
            metadata_out = self.index_path.with_suffix(".metadata.json")
        with Path(metadata_out).open("w") as f:
            json.dump(meta, f, indent=2)
        print(f"[Retriever] Indexed {len(ids)} seqs → {self.index_path}")

    def _ensure_loaded(self) -> None:
        if self._ids is None:
            if not self.index_path.exists():
                raise FileNotFoundError(
                    f"Index not found: {self.index_path}. Run build_index first."
                )
            data = np.load(self.index_path, allow_pickle=True)
            self._ids = data["ids"]
            self._seqs = data["seqs"]
            self._embs = data["embs"].astype(np.float32)
            # L2-normalize for cosine
            norms = np.linalg.norm(self._embs, axis=1, keepdims=True) + 1e-8
            self._embs = self._embs / norms

    # ---------- Query ----------
    def query(
        self, sequence: str, top_k: int = 10, with_identity: bool = True
    ) -> List[RetrievalHit]:
        self._ensure_loaded()
        assert (
            self._ids is not None and self._seqs is not None and self._embs is not None
        )

        seq = validate_protein_sequence(sequence)
        q = embed_sequence_mean(seq).astype(np.float32)
        q = q / (np.linalg.norm(q) + 1e-8)

        sims = self._embs @ q  # cosine similarity
        idxs = np.argsort(-sims)[:top_k]

        hits: List[RetrievalHit] = []

        if with_identity:
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 1.0
            aligner.mismatch_score = 0.0
            aligner.open_gap_score = 0.0
            aligner.extend_gap_score = 0.0

        for j, i in enumerate(idxs):
            sid = str(self._ids[i])
            sseq = str(self._seqs[i])
            sim = float(sims[i])
            ident = None
            if with_identity:
                # identity = matches / max(lenA, lenB)
                matches = aligner.score(seq, sseq)
                ident = float(matches) / float(max(len(seq), len(sseq)))
            hits.append(
                RetrievalHit(id=sid, sequence=sseq, similarity=sim, identity=ident)
            )

        return hits

    # ---------- MSA ----------
    def build_msa(
        self,
        query: str,
        hits: List[RetrievalHit],
        out_fasta: str = "data/cache/alignment.fasta",
        tool: str = "clustalo",
        dedup: bool = True,
        balance: bool = True,
    ) -> str:
        """
        Build an MSA using Clustal Omega or MAFFT, with optional deduplication & taxonomic balancing.
        """
        seqs = [("query", query)] + [(h.id, h.sequence) for h in hits]

        if dedup:
            seqs = deduplicate_sequences(seqs, identity_thresh=0.95)
            print(f"[Retriever] Deduplicated → {len(seqs)} sequences")

        if balance:
            seqs = balance_taxa(seqs, max_per_taxon=2)
            print(f"[Retriever] Taxonomic balancing → {len(seqs)} sequences")

        tmp_in = Path("data/cache/tmp_input.fasta")
        out_fasta = Path(out_fasta)

        with tmp_in.open("w") as f:
            for sid, s in seqs:
                f.write(f">{sid}\n{s}\n")

        if tool == "clustalo":
            cmd = ["clustalo", "-i", str(tmp_in), "-o", str(out_fasta), "--force"]
            subprocess.run(cmd, check=True)
        elif tool == "mafft":
            cmd = ["mafft", str(tmp_in)]
            with out_fasta.open("w") as fout:
                subprocess.run(cmd, check=True, stdout=fout)
        else:
            raise ValueError(f"Unsupported MSA tool: {tool}")

        print(f"[Retriever] MSA written to {out_fasta} using {tool}")
        return str(out_fasta)

    def analyze_msa(
        self, msa_path: str, out_json: str = "data/cache/conservation.json"
    ) -> str:
        """
        Compute conservation metrics (entropy + frequencies) from an MSA FASTA file.
        """
        alignment = AlignIO.read(msa_path, "fasta")
        msa = [str(rec.seq) for rec in alignment]

        results = compute_conservation(msa)

        out_json = Path(out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with out_json.open("w") as f:
            json.dump(results, f, indent=2)

        print(f"[Retriever] Conservation metrics saved to {out_json}")
        return str(out_json)

    def extract_constraints_from_msa(
        self,
        msa_path: str,
        sequence: str,
        conservation_json: str = "data/cache/conservation.json",
        out_json: str = "data/cache/constraints.json",
    ) -> str:
        """
        Given an MSA + conservation analysis, extract biological constraints.
        """
        import json

        with open(conservation_json) as f:
            conservation = json.load(f)

        constraints = extract_constraints(sequence, conservation)

        out_json = Path(out_json)
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with out_json.open("w") as f:
            json.dump(constraints, f, indent=2)

        print(f"[Retriever] Constraints saved to {out_json}")
        return str(out_json)


# ---------- CLI ----------
def _cli():
    parser = argparse.ArgumentParser("RetrieverAgent CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("analyze", help="Analyze MSA for conservation metrics")
    a.add_argument("--msa", required=True, help="Input aligned FASTA")
    a.add_argument(
        "--out", default="data/cache/conservation.json", help="Output JSON file"
    )

    b = sub.add_parser("build", help="Build embedding index from FASTA")
    b.add_argument("--fasta", required=True, help="Path to reference FASTA")
    b.add_argument(
        "--out", default="data/cache/retriever_index.npz", help="Output index path"
    )

    c = sub.add_parser(
        "constraints", help="Extract biological constraints from MSA + conservation"
    )
    c.add_argument("--seq", required=True, help="Original protein sequence")
    c.add_argument("--msa", required=True, help="Aligned FASTA file")
    c.add_argument("--conservation", default="data/cache/conservation.json")
    c.add_argument("--out", default="data/cache/constraints.json")
    # new knobs
    c.add_argument(
        "--conserved-mode",
        choices=["threshold", "percentile", "topk"],
        default="percentile",
    )
    c.add_argument("--entropy-thresh", type=float, default=None)
    c.add_argument("--percentile", type=float, default=5.0)
    c.add_argument("--topk", type=int, default=None)
    c.add_argument("--min-majority", type=float, default=0.9)
    c.add_argument("--min-occupancy", type=float, default=0.8)

    q = sub.add_parser("query", help="Query top-k homologs")
    q.add_argument("--seq", required=True, help="Protein sequence (letters only)")
    q.add_argument("--k", type=int, default=5, help="Top-k")
    q.add_argument(
        "--index", default="data/cache/retriever_index.npz", help="Index path"
    )

    m = sub.add_parser("msa", help="Build MSA for query + top-k hits")
    m.add_argument("--seq", required=True, help="Protein sequence (letters only)")
    m.add_argument("--k", type=int, default=5, help="Top-k")
    m.add_argument(
        "--index", default="data/cache/retriever_index.npz", help="Index path"
    )
    m.add_argument(
        "--out", default="data/cache/alignment.fasta", help="Output aligned FASTA"
    )
    m.add_argument(
        "--tool",
        choices=["clustalo", "mafft"],
        default="clustalo",
        help="MSA tool to use",
    )
    m.add_argument("--no-dedup", action="store_true", help="Disable deduplication")
    m.add_argument(
        "--no-balance", action="store_true", help="Disable taxonomic balancing"
    )

    v = sub.add_parser("visualize", help="Visualize conservation logo from JSON")
    v.add_argument("--json", required=True, help="Input conservation.json")
    v.add_argument(
        "--out", default="data/cache/conservation_logo.png", help="Output PNG"
    )
    v.add_argument("--top", type=int, default=5, help="Top N residues to display")

    args = parser.parse_args()
    if args.cmd == "build":
        agent = RetrieverAgent(index_path=args.out)
        agent.build_index(args.fasta)
    elif args.cmd == "query":
        agent = RetrieverAgent(index_path=args.index)
        hits = agent.query(args.seq, top_k=args.k, with_identity=True)
        for h in hits:
            ident_pct = f"{h.identity * 100:.1f}%" if h.identity is not None else "n/a"
            print(f"{h.id}\tsim={h.similarity:.3f}\tident={ident_pct}")
    elif args.cmd == "msa":
        agent = RetrieverAgent(index_path=args.index)
        hits = agent.query(args.seq, top_k=args.k, with_identity=False)
        agent.build_msa(
            args.seq,
            hits,
            out_fasta=args.out,
            tool=args.tool,
            dedup=not args.no_dedup,
            balance=not args.no_balance,
        )

    elif args.cmd == "analyze":
        agent = RetrieverAgent()
        agent.analyze_msa(args.msa, out_json=args.out)
    elif args.cmd == "visualize":
        plot_conservation_logo(args.json, out_png=args.out, top_n=args.top)
    elif args.cmd == "constraints":
        agent = RetrieverAgent()
        # load conservation json
        import json

        with open(args.conservation) as f:
            conservation = json.load(f)
        cons = extract_constraints(
            args.seq,
            conservation,
            conserved_mode=args.conserved_mode,
            entropy_thresh=args.entropy_thresh,
            percentile=args.percentile,
            topk=args.topk,
            min_majority=args.min_majority,
            min_occupancy=args.min_occupancy,
        )
        outp = Path(args.out)
        outp.parent.mkdir(parents=True, exist_ok=True)
        with outp.open("w") as f:
            json.dump(cons, f, indent=2)
        print(f"[Retriever] Constraints saved to {outp}")
    else:
        raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    _cli()
