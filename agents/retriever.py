# agents/retriever.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List
import json
import argparse
import numpy as np

from utils.io import read_fasta, validate_protein_sequence
from models.esm import embed_sequence_mean
from Bio.Align import PairwiseAligner  # biopython


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
        # identities: List[float | None] = [None] * len(idxs)

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


# ---------- CLI ----------
def _cli():
    parser = argparse.ArgumentParser("RetrieverAgent CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="Build embedding index from FASTA")
    b.add_argument("--fasta", required=True, help="Path to reference FASTA")
    b.add_argument(
        "--out", default="data/cache/retriever_index.npz", help="Output index path"
    )

    q = sub.add_parser("query", help="Query top-k homologs")
    q.add_argument("--seq", required=True, help="Protein sequence (letters only)")
    q.add_argument("--k", type=int, default=5, help="Top-k")
    q.add_argument(
        "--index", default="data/cache/retriever_index.npz", help="Index path"
    )

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


if __name__ == "__main__":
    _cli()
