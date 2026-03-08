# Repository Policy

## Purpose

`wave-crispr-signal` is a **CRISPR signal analytics toolkit**. Mainline code must directly support CRISPR scoring/design/validation outcomes.

## Scope Rules

- **Mainline (`core`/`support`)**: CRISPR package code, apps, tests, proof/repro tooling, datasets, and docs for current scope.
- **Legacy (`research-legacy`)**: non-CRISPR exploratory work must live under `legacy/`.
- **Generated artifacts (`delete/regenerate`)**: must not be committed to repo root.

## Top-Level Layout (Allowed)

- `.github/`, `applications/`, `bin/`, `configs/`, `data/`, `docs/`, `experiments/`, `legacy/`, `modules/`, `notebooks/`, `proof_pack/`, `scripts/`, `static/`, `templates/`, `test_data/`, `tests/`, `tools/`, `wave_crispr_signal/`, `web_apps/`

Allowed root files:
- `.gitignore`, `LICENSE`, `Makefile`, `README.md`, `pyproject.toml`, `repo_metadata.yaml`, `requirements.txt`, `targets.yaml`

## Guardrails

- No generated artifact files at root (for example `acf_scale_*.png`, `null_dist_output.txt`, root falsification outputs).
- No non-CRISPR domain modules in core package paths (`wave_crispr_signal/`, `applications/`).
- `docs/INDEX.md` must exist and reference both `docs/INVENTORY_CLASSIFICATION.md` and `legacy/README.md`.

## Compatibility

- Public CRISPR-facing package/API/CLI contracts must remain stable.
- Historical non-CRISPR import paths may be preserved through explicit wrappers in `experiments/` that forward to `legacy/`.

## Enforcement

- CI executes `python tools/check_repo_policy.py`.
- Structure policy tests run in `tests/test_repo_structure_policy.py`.
