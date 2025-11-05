# GEMINI — You are here: `wave-crispr-signal`

**Purpose**: A framework for signal-theoretic analysis of DNA mutations, encoding DNA as a complex waveform and interrogating it with FFT-based metrics plus geodesic curvature weighting.

## TL;DR Build/Run
```bash
# Install dependencies
pip install -r requirements.txt

# Run quick validation demo
python proof_pack/quick_validation_demo.py
```

## Common Tasks (copy/paste)

*   Run full test suite: `python -m pytest -q`
*   Design CRISPR guides: `python applications/crispr_cli.py design "ATGCTGCGGA..." -n 5 -o guides.json`
*   Score an existing guide: `python applications/crispr_cli.py score "GACGATCGATCGATCGATCG"`

## Key Files & Dirs (minimal map)

*   `applications/` — CRISPR guide designer and related tools.
*   `proof_pack/` — Validation scripts and reproducible proofs.
*   `docs/` — Detailed documentation, including topological analysis.
*   `tests/` — Unit tests.
*   `notebooks/` — Jupyter notebooks for analysis and comparison.

## Pointers

*   **Parent map**: [`../GEMINI.md`](../GEMINI.md)
*   **Siblings**: [`../unified-framework/GEMINI.md`](../unified-framework/GEMINI.md), [`../z-sandbox/GEMINI.md`](../z-sandbox/GEMINI.md)

## Ask Gemini (do this, not that)

*   “Explain the geodesic-topological bridge in this project.”
*   “Show me how the framework encodes DNA as a complex waveform.”
*   “Find the scripts for running the full validation suite.”
