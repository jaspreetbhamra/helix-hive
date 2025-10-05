from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import AlignIO  # for mapping alignment columns to query positions


@dataclass
class DesignProposal:
    mutation: str  # e.g., "A123V"
    position: int  # 1-based in the *ungapped query* sequence
    old_res: str
    new_res: str
    delta_pssm: float  # log2 odds (new - old)
    entropy: float  # per-column Shannon entropy
    majority: float  # consensus frequency among non-gaps (0..1)
    reason: str  # "consensus"
    notes: str  # any flags/reasons


class DesignAgent:
    def __init__(
        self,
        respect_constraints: bool = True,
        skip_conserved_sites: bool = True,
        skip_cysteines: bool = True,
        skip_glyco: bool = True,
        skip_signal_peptide: bool = True,
        skip_tm_spans: bool = True,
        min_occupancy: float = 0.7,  # require enough non-gap depth
        min_majority: float = 0.6,  # require clear consensus
    ):
        self.respect_constraints = respect_constraints
        self.skip_conserved_sites = skip_conserved_sites
        self.skip_cysteines = skip_cysteines
        self.skip_glyco = skip_glyco
        self.skip_signal_peptide = skip_signal_peptide
        self.skip_tm_spans = skip_tm_spans
        self.min_occupancy = min_occupancy
        self.min_majority = min_majority

    # ----------------- helpers -----------------
    @staticmethod
    def _load_json(path: str | Path):
        with Path(path).open() as f:
            return json.load(f)

    @staticmethod
    def _alignment_map(msa_fasta: str | Path) -> Tuple[List[str], Dict[int, int]]:
        """
        Returns (query_aln_seq, col_to_query_pos):
          - query_aln_seq: the aligned sequence string for 'query'
          - col_to_query_pos: mapping from alignment column (0-based) to ungapped query position (1-based)
        """
        aln = AlignIO.read(str(msa_fasta), "fasta")
        # find the row labeled "query"
        try:
            qrec = next(rec for rec in aln if rec.id == "query")
        except StopIteration:
            raise ValueError(
                "MSA must contain a sequence with id 'query' as first row in our pipeline."
            )
        qseq = str(qrec.seq)

        col_to_pos: Dict[int, int] = {}
        running = 0
        for i, ch in enumerate(qseq):
            if ch != "-":
                running += 1
                col_to_pos[i] = running
        return [qseq], col_to_pos

    @staticmethod
    def _in_any_span(pos: int, spans: List[List[int] | Tuple[int, int]]) -> bool:
        for s, e in spans:
            if s <= pos <= e:
                return True
        return False

    # ----------------- core -----------------
    def run(
        self,
        query_sequence: str,
        msa_fasta: str | Path,
        conservation_json: str | Path,
        pssm_json: str | Path,
        constraints_json: Optional[str | Path] = None,
        top_n: int = 20,
    ) -> List[DesignProposal]:
        """
        Generate ranked consensus-based mutation proposals respecting constraints and PSSM.
        """
        # Load data
        conservation = self._load_json(conservation_json)  # list per alignment column
        pssm = self._load_json(
            pssm_json
        )  # {"positions":[...], "pssm":[{aa:score}, ...], "alphabet":[...]}
        constraints = (
            self._load_json(constraints_json)
            if (constraints_json and Path(constraints_json).exists())
            else {}
        )

        # Build mapping from alignment columns -> ungapped query positions
        _, col_to_qpos = self._alignment_map(msa_fasta)

        # Quick lookups for constraints
        conserved_sites = set(constraints.get("conserved_sites", []))
        cysteines = set(constraints.get("cysteines", []))
        glyco_spans = constraints.get("glycosylation_sites", []) or []
        tm_spans = constraints.get("transmembrane_spans", []) or []
        has_signal_pep = bool(constraints.get("signal_peptide", False))

        # Iterate alignment columns; propose consensus→query diffs
        proposals: List[DesignProposal] = []

        for col_idx, col_info in enumerate(
            conservation
        ):  # col_info has: position, frequencies, entropy, (maybe majority/occupancy in your latest version)
            # map column → query position (ungapped). If gap at query, skip.
            qpos = col_to_qpos.get(col_idx)
            if qpos is None:
                continue  # query has gap here

            # Depth / filters
            freqs: Dict[str, float] = col_info.get("frequencies", {})
            entropy: float = float(col_info.get("entropy", 0.0))
            majority: float = float(
                col_info.get("majority", max(freqs.values()) if freqs else 0.0)
            )
            occupancy: float = float(
                col_info.get("occupancy", 1.0)
            )  # if not present, assume high occupancy

            if occupancy < self.min_occupancy:
                continue
            if majority < self.min_majority:
                continue

            # consensus residue
            if not freqs:
                continue
            consensus_res = max(freqs.items(), key=lambda kv: kv[1])[0]

            old_res = query_sequence[qpos - 1]
            if old_res == consensus_res:
                continue  # already consensus

            # Constraint guardrails
            notes = []
            if self.respect_constraints:
                if self.skip_conserved_sites and (qpos in conserved_sites):
                    notes.append("skip:conserved_site")
                if self.skip_cysteines and (
                    qpos in cysteines or old_res == "C" or consensus_res == "C"
                ):
                    notes.append("skip:cysteine")
                if self.skip_glyco and self._in_any_span(qpos, glyco_spans):
                    notes.append("skip:glyco_motif")
                if self.skip_tm_spans and self._in_any_span(qpos, tm_spans):
                    notes.append("skip:tm_span")
                if self.skip_signal_peptide and has_signal_pep and qpos <= 25:
                    notes.append("skip:signal_peptide_region")

            if any(n.startswith("skip:") for n in notes):
                continue

            # PSSM delta score (alignment positions are 1-based in pssm["positions"])
            try:
                aln_pos_1b = col_info["position"]  # alignment space 1-based
                pssm_idx = pssm["positions"].index(aln_pos_1b)
            except Exception:
                # if positions mismatch, skip safely
                continue

            scores_at_col: Dict[str, float] = pssm["pssm"][pssm_idx]
            new_score = float(scores_at_col.get(consensus_res, 0.0))
            old_score = float(scores_at_col.get(old_res, 0.0))
            delta = new_score - old_score

            mutation = f"{old_res}{qpos}{consensus_res}"
            proposals.append(
                DesignProposal(
                    mutation=mutation,
                    position=qpos,
                    old_res=old_res,
                    new_res=consensus_res,
                    delta_pssm=delta,
                    entropy=entropy,
                    majority=majority,
                    reason="consensus",
                    notes=";".join(notes) if notes else "",
                )
            )

        # Rank: higher ΔPSSM first, then lower entropy (more conserved carries risk; prefer variable), then higher majority
        proposals.sort(
            key=lambda p: (-p.delta_pssm, p.entropy, -p.majority, p.position)
        )
        return proposals[:top_n]


# --------------- CLI ---------------
def _cli():
    import argparse

    parser = argparse.ArgumentParser("DesignAgent (consensus + constraints + PSSM)")
    parser.add_argument(
        "--seq", required=True, help="Ungapped query protein sequence (letters only)"
    )
    parser.add_argument(
        "--msa", required=True, help="Aligned FASTA (must contain entry id='query')"
    )
    parser.add_argument("--conservation", required=True, help="conservation.json")
    parser.add_argument("--pssm", required=True, help="pssm.json")
    parser.add_argument("--constraints", default=None, help="constraints.json")
    parser.add_argument("--top", type=int, default=20, help="Top-N proposals")
    parser.add_argument(
        "--out", default="data/cache/design_proposals.json", help="Output JSON"
    )

    # knobs
    parser.add_argument("--min-occupancy", type=float, default=0.7)
    parser.add_argument("--min-majority", type=float, default=0.6)
    parser.add_argument("--no-constraints", action="store_true")
    parser.add_argument("--keep-conserved", action="store_true")
    parser.add_argument("--allow-cys", action="store_true")
    parser.add_argument("--allow-glyco", action="store_true")
    parser.add_argument("--allow-sp", action="store_true")
    parser.add_argument("--allow-tm", action="store_true")

    args = parser.parse_args()

    agent = DesignAgent(
        respect_constraints=not args.no_constraints,
        skip_conserved_sites=not args.keep_conserved,
        skip_cysteines=not args.allow_cys,
        skip_glyco=not args.allow_glyco,
        skip_signal_peptide=not args.allow_sp,
        skip_tm_spans=not args.allow_tm,
        min_occupancy=args.min_occupancy,
        min_majority=args.min_majority,
    )
    proposals = agent.run(
        query_sequence=args.seq,
        msa_fasta=args.msa,
        conservation_json=args.conservation,
        pssm_json=args.pssm,
        constraints_json=args.constraints,
        top_n=args.top,
    )

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("w") as f:
        json.dump([asdict(p) for p in proposals], f, indent=2)

    # nice console preview
    print(f"[Design] Wrote {len(proposals)} proposals → {outp}")
    for p in proposals[: min(10, len(proposals))]:
        print(
            f"{p.mutation}\tΔPSSM={p.delta_pssm:+.2f}\tentropy={p.entropy:.2f}\tmaj={p.majority:.2f}\t{p.reason}"
        )


if __name__ == "__main__":
    _cli()
