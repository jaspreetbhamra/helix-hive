# Project Makefile for helix_hive

# === Environment ===
.PHONY: setup sync clean

setup:
	uv sync

sync:
	uv sync --upgrade

clean:
	rm -rf __pycache__ .pytest_cache .ruff_cache .mypy_cache .uv .coverage

# === Code Quality ===
.PHONY: lint format check

lint:
	uv run ruff check .

format:
	uv run ruff check . --fix

check:
	uv run ruff check . --exit-zero

# === Testing ===
.PHONY: test

test:
	uv run pytest -q

# === Build a retriever index ===
.PHONY: index
index:
	uv run python -m agents.retriever build --fasta data/examples/small.fasta --out data/cache/retriever_index.npz

# === Quick query (example; edit the seq) ===
.PHONY: query
query:
	uv run python -m agents.retriever query --seq MKTAYIAKQRQISFVKSHFSRQDILDLWQ --k 3 --index data/cache/retriever_index.npz

# === Test design agent ===
.PHONY: design
design:
	uv run python -m agents.design


# === Run App ===
.PHONY: dashboard

dashboard:
	uv run streamlit run ui/dashboard.py
