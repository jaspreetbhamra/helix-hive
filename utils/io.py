# utils/io.py
from pathlib import Path
from typing import List, Tuple


def read_fasta(path: str | Path) -> List[Tuple[str, str]]:
    """Return list of (id, seq) from a FASTA file."""
    path = Path(path)
    records: List[Tuple[str, str]] = []
    if not path.exists():
        raise FileNotFoundError(f"FASTA not found: {path}")

    with path.open("r") as f:
        seq_id = None
        seq_parts: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    records.append((seq_id, "".join(seq_parts)))
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append("".join(ch for ch in line if ch.isalpha()).upper())
        if seq_id is not None:
            records.append((seq_id, "".join(seq_parts)))
    return records


def validate_protein_sequence(seq: str) -> str:
    """Uppercase, remove non-letters, basic amino-acid filter."""
    seq = "".join(ch for ch in seq if ch.isalpha()).upper()
    if not seq:
        raise ValueError("Empty protein sequence.")
    return seq
