import subprocess
from pathlib import Path
from typing import List, Literal


def run_alignment(
    sequences: List[tuple[str, str]],
    out_fasta: str | Path,
    tool: Literal["clustalo", "mafft"] = "clustalo",
) -> Path:
    """
    Run MSA using chosen tool (clustalo or mafft).
    """
    tmp_in = Path("data/cache/tmp_input.fasta")
    out_fasta = Path(out_fasta)

    with tmp_in.open("w") as f:
        for sid, seq in sequences:
            f.write(f">{sid}\n{seq}\n")

    if tool == "clustalo":
        cmd = ["clustalo", "-i", str(tmp_in), "-o", str(out_fasta), "--force"]
    elif tool == "mafft":
        cmd = ["mafft", str(tmp_in)]
        with out_fasta.open("w") as fout:
            subprocess.run(cmd, check=True, stdout=fout)
        return out_fasta
    else:
        raise ValueError(f"Unsupported tool: {tool}")

    subprocess.run(cmd, check=True)
    return out_fasta
