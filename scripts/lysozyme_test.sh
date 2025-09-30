#Build index (if you have a homolog dataset FASTA, even just 10–20 sequences).
uv run python -m agents.retriever build --fasta data/examples/lysozyme.fasta


#Query with the lysozyme sequence.
uv run python -m agents.retriever query --seq KVFGRCELAAAMKRHGLDNYR... --k 5


#MSA with Clustal Omega.
uv run python -m agents.retriever msa --seq KVFGRCELAAAMKRHGLDNYR... --k 5 --tool clustalo


#Conservation analysis.
uv run python -m agents.retriever analyze --msa data/cache/alignment.fasta --out data/cache/conservation.json


#Constraint extraction.
uv run python -m agents.retriever constraints --seq KVFGRCELAAAMKRHGLDNYR... --msa data/cache/alignment.fasta --conservation data/cache/conservation.json --out data/cache/constraints.json


#Visualization (logo + entropy + constraints).
uv run python -m agents.retriever visualize --json data/cache/conservation.json --out data/cache/lysozyme_plot.png --top 6 --constraints data/cache/constraints.json



# What You’ll See for Lysozyme

# Conserved catalytic Glu35 and Asp52 → very low entropy, marked as conserved.
# Eight conserved cysteines → disulfide bonds.
# No signal peptide (since mature protein).
# Entropy peaks in loops → candidate mutation sites.