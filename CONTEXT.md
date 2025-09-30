## Big Picture

You have a target protein (a string over the 20 amino acids). You want to change its sequence so some property improves—stability, activity, binding, solubility, expression, etc.—while keeping it folded and functional. Your system is building an AI co-pilot for protein engineering. Steps:
- Retrieve knowledge from nature (homologs, conserved positions, known motifs/constraints, structure hints).
- Design sequence edits (mutations) that are plausible given biology and aligned to your goal.
- Evaluate & rank candidates with cheap in-silico filters; then hand the best to the wet lab.

## Background Info
- Sequence → Structure → Function. A protein is a polymer of amino acids. The sequence (primary structure) folds into a 3D conformation (secondary/tertiary). Function emerges from that fold.
- Amino acids differ in size, charge, hydrophobicity, and special roles (e.g., Cys makes disulfide bonds; Pro bends backbones; Gly is flexible; His often catalyzes).
- Conservation encodes importance. Positions that barely change across evolution are likely critical (active site, structural core). Positions that vary are more tolerant to mutation.
- Fitness landscape & epistasis. Each mutation moves you on a rugged landscape; effects can interact (A is good alone, bad with B). Start with single mutants and keep combinations small unless you have strong prior.
- Trade-offs are real. Stabilizing mutations can hurt activity; solubility can fight binding; expression can oppose secretion efficiency. Know your primary KPI.

## Data

- Homologous sequences & multiple sequence alignments (MSA). Nature’s A/B tests. Show which residues are conserved vs variable; reveal consensus residues and co-varying positions.
- Substitution matrices (e.g., BLOSUM/PAM). Empirical log-odds for swapping residue X→Y based on what evolution tolerates in aligned positions. High score = conservative (safe-ish).
- Structure / structure confidence. Experimental (PDB) or predicted (e.g., AlphaFold). Lets you reason about core vs surface, loops vs helices, interfaces, disulfides, active sites.
- Annotations/motifs. Domains, catalytic residues, signal peptides, transmembrane spans, glycosylation motifs (N-X-S/T), low-complexity or aggregation-prone regions.


## Workflow - Agents in the Multi-Agent System

### Retrieval Agent

**Goal:** Pull in “priors” from biology so design isn’t blind.

1. **Homolog search** (BLAST/HHblits-like) for your target sequence.
2. **Build an MSA** (deduplicate, filter close/remote hits, balance clades).
3. **Per-position statistics**:

   * Frequency (f_{i,a}) of amino acid (a) at position (i)
   * **Conservation** (Shannon entropy) (H_i=-\sum_a f_{i,a}\log f_{i,a})
   * **Consensus** residue (c_i=\arg\max_a f_{i,a})
4. **Constraint extraction**: mark high-conservation sites, catalytic motifs, cystines/disulfides, glycosylation sites, signal peptides, TM spans, interface residues (from structure, if available).
5. **Optionally** compute **PSSM** (position-specific scoring) and/or co-variation hints to avoid breaking pairwise contacts.

**Outputs:** an “evidence pack” your designer can use:

* MSA + per-position conservation/consensus
* A “do-not-touch” mask (active site, disulfides, glyco, catalytic H/D/E, metal-binding)
* Region annotations (core vs surface, loop vs helix) if structure is present
* Optional thermostable homolog shortlist (if your KPI is stability)



### Build MSA in Retrieval Agent 
Options to build MSA using either clustal or mafft
```bash
# Build index
uv run python -m agents.retriever build --fasta data/examples/small.fasta

# Run query
uv run python -m agents.retriever query --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5

# Build MSA with Clustal Omega (default)
uv run python -m agents.retriever msa --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5

# Build MSA with MAFFT
uv run python -m agents.retriever msa --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5 --tool mafft

```

### Conservation analysis (entropy + per-position frequencies)

Ccompute conservation metrics:
- Per-position amino acid frequencies
- Shannon entropy (lower = more conserved)

This gives us the raw stats that later agents (Design, Evaluator) can use for constraint-aware mutations.


```bash
# Build index
uv run python -m agents.retriever build --fasta data/examples/small.fasta

# Build MSA (Clustal Omega)
uv run python -m agents.retriever msa --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5 --out data/cache/alignment.fasta

# Analyze conservation
uv run python -m agents.retriever analyze --msa data/cache/alignment.fasta --out data/cache/conservation.json
```

Interpretation
- Entropy = 0 → strictly conserved (likely critical residues).
- High entropy → variable region (safe to mutate).
- Frequencies → guide consensus-based mutations later.


### Visualize conservation patterns after MSA analysis
```bash
# Build MSA
uv run python -m agents.retriever msa --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5

# Analyze conservation
uv run python -m agents.retriever analyze --msa data/cache/alignment.fasta --out data/cache/conservation.json

# Visualize conservation
uv run python -m agents.retriever visualize --json data/cache/conservation.json --out data/cache/conservation_logo.png --top 6

```

Interpretation:
- Can see conserved positions at a glance (one letter dominates).
- Variable positions will show multiple residues stacked.
- Helps decide which residues are safe to mutate in the Design Agent.


### Constraint Extraction 
So you don’t just get statistics, but also biological insights that guide design

We’ll implement these constraints locally (all lightweight):
1. Highly conserved sites
- From entropy: mark positions with entropy ≤ threshold (e.g. 0.1).
- Likely functionally critical — avoid mutating. 
2. Cysteines / Disulfides
- Detect Cys (C) positions.
- If multiple conserved cysteines → candidate disulfide bonds.
3. N-linked Glycosylation motifs (N-X-[ST])
- Regex in sequence: N[^P][ST]
- Mark motif positions.
4. Signal peptide heuristic
- N-terminal hydrophobic stretch (~first 20–30 residues, ≥70% hydrophobic).
- Use Kyte–Doolittle hydropathy index.
5. Transmembrane spans (simple heuristic)
- Scan windows of 19 residues, mark if hydrophobic index ≥ threshold.
- (True TM predictors = heavy, but this gives a PoC).

```bash
# 1. Build MSA
uv run python -m agents.retriever msa --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 5

# 2. Analyze conservation
uv run python -m agents.retriever analyze --msa data/cache/alignment.fasta --out data/cache/conservation.json

# 3. Extract constraints
uv run python -m agents.retriever constraints --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --msa data/cache/alignment.fasta --conservation data/cache/conservation.json --out data/cache/constraints.json

```


* Run lysozyme_with_homologs_test.sh - to test the pipeline so far
* Compare against examples/expected_lysozyme_constraints.json
- What’s in there (and why)
   - conserved_sites: the catalytic pair E35 and D52 are known to be highly conserved in lysozymes; even with a tiny homolog set, your entropy should flag these positions as low-entropy (≈0).
   - cysteines: lysozymes carry 8 cysteines (disulfide bonds); positions above are for the exact chicken lysozyme sequence you’re using (1-based indexing).
   - glycosylation_sites: no N-X-[ST] motif in the query sequence provided → [].
   - signal_peptide: the sequence you’re using is the mature enzyme, not the prepro-peptide → false under our hydropathy heuristic.
   - transmembrane_spans: lysozyme is secreted/soluble → none expected.

```bash
# Generate actual constraints
uv run python -m agents.retriever constraints \
  --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV \
  --msa data/cache/lysozyme_alignment.fasta \
  --conservation data/cache/lysozyme_conservation.json \
  --out data/cache/lysozyme_constraints.json

```
* Evaluate Constraints
- `uv run python tools/validate_constraints.py data/examples/expected_lysozyme_constraints.json data/cache/lysozyme_constraints.json`


### MSA Deduplication and Taxonomic Balancing

🧹 1. Sequence De-duplication (before MSA)
- Why
   - Removes exact duplicates or sequences that are > 95% identical.
   - Keeps your alignment diverse, avoids overweighting one species.
- How
   - Use pairwise identity with Biopython.
   - Keep first occurrence, drop near-duplicates.

🌍 2. Taxonomic Balancing
- Why
   - Prevents bias toward one clade (e.g. all mammals).
   - Ensures you have a representative mix across taxonomic groups.
- How
   - If your FASTA headers include species info (>LYSC_HUMAN), group by species or genus prefix.
   - Subsample so each group contributes at most N sequences.

```bash
# With deduplication + balancing (default)
uv run python -m agents.retriever msa \
  --seq KVFGRCELAAAMKRHGL... \
  --k 10 \
  --tool clustalo

# Disable deduplication
uv run python -m agents.retriever msa \
  --seq KVFGRCELAAAMKRHGL... \
  --k 10 \
  --tool clustalo --no-dedup

```

* Why This Matters

- Deduplication prevents overcounting of near-identical sequences.
- Taxonomic balancing ensures you don’t bias conservation toward one species.
- Together → make entropy + consensus analysis biologically realistic.
- Without these, your Design Agent might wrongly think “everything is conserved” just because you loaded 20 nearly identical chicken sequences.

* Experiment
TODO: Duplicate one of the homologs n times to make the effect clear eg chicken 5x
```bash
# 1. Build index
uv run python -m agents.retriever build \
  --fasta data/examples/lysozyme_homologs_with_dupes.fasta \
  --out data/cache/lysozyme_index.npz

# 2. Run MSA with default (dedup + balance on)
uv run python -m agents.retriever msa \
  --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV \
  --index data/cache/lysozyme_index.npz \
  --k 10 \
  --out data/cache/lysozyme_msa_dedup.fasta \
  --tool clustalo

# Analyze conservation
uv run python -m agents.retriever analyze \
  --msa data/cache/lysozyme_msa_dedup.fasta \
  --out data/cache/lysozyme_cons_dedup.json

# Extract constraints
uv run python -m agents.retriever constraints \
  --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV \
  --msa data/cache/lysozyme_msa_dedup.fasta \
  --conservation data/cache/lysozyme_cons_dedup.json \
  --out data/cache/lysozyme_constraints_dedup.json \
  --conserved-mode percentile --percentile 5

# Visualize
uv run python -m agents.retriever visualize \
  --json data/cache/lysozyme_cons_dedup.json \
  --constraints data/cache/lysozyme_constraints_dedup.json \
  --out data/cache/lysozyme_plot_dedup.png \
  --top 6

# 3. Run MSA without deduplication / balancing
uv run python -m agents.retriever msa \
  --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQNRDVRQYVQGCGV \
  --index data/cache/lysozyme_index.npz \
  --k 10 \
  --out data/cache/lysozyme_msa_nodedup.fasta \
  --tool clustalo \
  --no-dedup --no-balance

# Analyze conservation
uv run python -m agents.retriever analyze \
  --msa data/cache/lysozyme_msa_nodedup.fasta \
  --out data/cache/lysozyme_cons_nodedup.json

# Extract constraints
uv run python -m agents.retriever constraints \
  --seq KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCQYVQGCGV \
  --msa data/cache/lysozyme_msa_nodedup.fasta \
  --conservation data/cache/lysozyme_cons_nodedup.json \
  --out data/cache/lysozyme_constraints_nodedup.json \
  --conserved-mode percentile --percentile 5

# Visualize
uv run python -m agents.retriever visualize \
  --json data/cache/lysozyme_cons_nodedup.json \
  --constraints data/cache/lysozyme_constraints_nodedup.json \
  --out data/cache/lysozyme_plot_nodedup.png \
  --top 6

```

🔍 What to Compare
- With dedup + balance:
   - Entropy plot: only truly conserved positions (like catalytic residues) stay flat.
   - Constraints: fewer “false” conserved sites, clearer biological signals.

- Without dedup/balance (especially if you added chicken duplicates):
   - Entropy plot: almost the whole sequence looks conserved.
   - Constraints: huge list of “conserved sites” (like what you saw earlier).
   - Misleading for Design Agent.

🌐 Why This Test is Important
This experiment demonstrates why bioinformatics pipelines spend so much time on redundancy filtering and clade balancing. Without it, the statistical signal is biased — and your downstream Design/Evaluator Agents make bad calls.

After running the scripts above, run the following to visualize
```bash
uv run python tools/compare_conservation.py \
  --dedup data/cache/lysozyme_cons_dedup.json \
  --raw data/cache/lysozyme_cons_nodedup.json \
  --dedup-constraints data/cache/lysozyme_constraints_dedup.json \
  --raw-constraints data/cache/lysozyme_constraints_nodedup.json \
  --out data/cache/lysozyme_compare.png
```


# Archived Explanation to get back to later - Generated by GPT


### Design Agent
**Goal:** Propose plausible mutations with guardrails.

Two design modes:

1. **Consensus / conservative edits (fast, safe baseline).**

   * For low-conservation positions ((H_i) high), suggest moving the target residue toward **consensus** (c_i).
   * Use **BLOSUM/PAM** to restrict to **conservative** substitutions (positive log-odds).
   * Skip protected positions (active site, disulfides, glyco motifs, TM core).
   * Great for **stability/expressibility** improvements with low risk.

2. **Model-guided proposals (smarter, goal-aware).**

   * **Protein language models** (e.g., ProtGPT-style) score the likelihood of residue choices in context; high-likelihood edits tend to be fold-compatible.
   * **Structure-aware designers** (e.g., ProteinMPNN-style) propose sidechain identities that pack well into a backbone—useful for stabilizing cores or interfaces.

**Ranking/filters your designer should apply now (cheap and effective):**

* Penalize changes at **highly conserved** or **annotated** residues.
* Favor **BLOSUM-positive** swaps; downrank radical physicochemical jumps (e.g., charged→hydrophobe in core).
* Prefer **surface** over **core** edits for first passes (less risk).
* Respect **sequence motifs** (N-X-S/T glyco; signal peptides; catalytic triads).
* Keep initial library to **single mutants** or **very small combos**.

# how the agents collaborate (the loop you already have)

1. **Retrieval**: “Here’s the MSA, consensus per position, conservation, structural masks, and sensitive sites.”
2. **Design**: “Given those priors + your KPI, here are N mutation suggestions, each with a rationale and safety checks passed.”
3. **(Soon) Evaluation**: “Score ΔΔG/aggregation/solubility/LLM likelihood; triage to a short list.”
4. **(Optional) Active learning** after assays: feed measured fitness back into the design policy.

# what you’ve effectively accomplished so far

* **You’re not designing in a vacuum.** Proposals are anchored in evolution (consensus + substitution matrices) and guided by “what nature already tried.”
* **You’ve added guardrails.** Catalytic/disulfide/glyco/TM constraints prevent catastrophic edits.
* **You can generate sane, diverse first-pass libraries** (dozens to a few hundred single-point suggestions) that are cheap to synthesize/test and likely to keep the fold intact.
* **You’re set up to plug in smarter scoring** (language-model or structure-aware) without changing the outer loop.

# a bit more math (so the heuristics feel principled)

* **Conservation (entropy):** low (H_i) ⇒ critical site; high (H_i) ⇒ tolerant.
* **PSSM** score for proposing (a) at (i): (S_{i,a}=\log \frac{f_{i,a}}{p_a}) (vs background).
* **BLOSUM log-odds** (L(x\to y)): positive = commonly tolerated swap; negative = historically disfavored.
* A simple **composite score** you can implement now per mutation ((i:x\to y)):
  [
  \text{Score} = w_1 \cdot L(x\to y) + w_2 \cdot S_{i,y} - w_3 \cdot \mathbb{1}[\text{conserved/annotated}] - w_4 \cdot \Delta\text{physchem}
  ]
  (tune (w_k); set the indicator penalty very large for “do-not-touch”.)

# practical constraints your designer should already respect

* **Disulfides:** never mutate Cys in a Cys–Cys pair unless you redesign both.
* **Catalytic/ligand sites:** keep the triad/metal chelators intact (commonly Ser/His/Asp; His/Glu/Asp for metals).
* **Glycosylation motif (N-X-S/T):** don’t break it if glycosylation matters—or *do* remove if your KPI is “no glyco.”
* **Signal peptides / TM segments:** hydrophobicity patterns are essential; avoid polar swaps in TM cores.
* **Stop codons & restriction sites** (downstream DNA step): avoid introducing them inadvertently.

# example: if your KPI is “more stable & better expressed”

* Retrieval: enrich MSA with **thermophilic homologs**; compute consensus; map **core vs surface** from structure.
* Design: prefer **consensus or conservative** swaps on **surface loops**, avoid fully conserved positions, propose hydrophobic→polar on solvent-exposed sticky patches to reduce aggregation; propose small-to-medium core packing fixes only if structure-aware scoring is available.
* Output: a ranked list of single-point edits with rationales like “A123S: matches consensus; BLOSUM +1; surface loop; reduces local hydrophobic patch.”

# common pitfalls (so you don’t chase ghosts)

* **Dirty MSAs** (paralogs, fragments) can fake “conservation” or “variability.” Filter length/identity and ensure domain boundaries align.
* **Domain boundaries**: editing across a linker vs inside a catalytic pocket is night-and-day risk.
* **Over-combining single hits**: additive is not guaranteed—watch epistasis.
* **Objective creep**: mixing multiple KPIs without a scoring rubric leads to noisy libraries.

# sensible next steps (when you’re ready)

1. **Evaluation agent (cheap first):**

   * Sequence-only scores: LM likelihood / perplexity deltas, heuristic aggregation/solubility predictors, isoelectric point shifts, net charge at pH.
   * Rule-based vetoes: motif breaks, Cys handling, TM polarity, glyco disruptions.
2. **Structure-aware scoring (if you have a model/backbone):**

   * Rapid ΔΔG or pack-quality proxies; interface clash checks if it’s a complex.
3. **Design “policies” by KPI:** Stability policy, solubility policy, expression policy—each with different weights/filters.
4. **Reporting & export:** FASTA of variants, CSV rationale table, optional DNA codon optimization for your expression system.

---

If you want, I can turn this into **concrete scoring/ranking code** that plugs right into your existing retrieval/design agents: read MSA stats + BLOSUM, apply masks, generate candidates, compute the composite score, and emit a tidy DataFrame ready for triage.
