# Experiments

This directory now serves two purposes:

- Active CRISPR experiment work (for example `signal_theoretic_crispr/`).
- Compatibility wrappers for historical non-CRISPR modules that were relocated to `legacy/experiments/`.

## Active Scope

Primary active experiment area:

- `signal_theoretic_crispr/` — CRISPR-focused spectral pipeline and validation.

## Legacy Compatibility

Several `experiments/*.py` modules are thin wrappers that import from `legacy/experiments/*`.
This preserves historical imports/tests while keeping non-CRISPR exploratory work out of mainline scope.

## Usage

```bash
python -m pytest -q
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```
