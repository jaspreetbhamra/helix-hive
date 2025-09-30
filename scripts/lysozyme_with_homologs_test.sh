# Build index
uv run python -m agents.retriever build --fasta data/examples/lysozyme_homologs.fasta --out data/cache/lysozyme_index.npz


# Query with chicken lysozyme (same as one sequence above, but simulates lookup)
uv run python -m agents.retriever query --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV --index data/cache/lysozyme_index.npz --k 5


# Build MSA
uv run python -m agents.retriever msa --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV --k 5 --index data/cache/lysozyme_index.npz --out data/cache/lysozyme_alignment.fasta --tool clustalo


# Analyze conservation
uv run python -m agents.retriever analyze --msa data/cache/lysozyme_alignment.fasta --out data/cache/lysozyme_conservation.json


# Extract constraints
uv run python -m agents.retriever constraints --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV --msa data/cache/lysozyme_alignment.fasta --conservation data/cache/lysozyme_conservation.json --out data/cache/lysozyme_constraints.json


# Visualize
uv run python -m agents.retriever visualize --json data/cache/lysozyme_conservation.json 


# What You Should See:
# Entropy plot: flat regions (conserved catalytic Glu35 & Asp52).
# Cysteine annotations: multiple conserved cysteines → disulfide bonds.
# Conserved sites: most of the catalytic cleft residues.
# Variable loops: entropy peaks → mutation candidates.