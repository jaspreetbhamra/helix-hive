# helix-hive
Agentic Collaborative Protein Engineer

## Steps to Setup the Repo

1. Install uv - `brew install uv` (or the equivalent for Linux/Windows)
2. `uv venv --python $(which python)`
3. Ensure that the `.python-version` file has the appropriate python version - might cause failures if the global python version is different from the python version for the current venv (mentioned in the pyproject.toml) file
    - Additionally, torchtext isn't currently supported with py3.13. Hence python<=3.12 is required
4. `uv sync`
5. `brew install muscle`

## Testing components one by one

### Retrieval Agent
1. `make index`
    - build a sample index from `data/examples/small.fasta`
2. `make query`

### Design Agent
1. `make design`