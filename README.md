# wave-crispr-signal

Signal-theoretic CRISPR analytics toolkit for guide scoring, disruption analysis, and reproducible validation.

## Project Charter (Scope Contract)

### Mission
Build and maintain a focused CRISPR signal-analysis toolkit that turns sequence-level signals into reproducible guide-quality and disruption metrics.

### Intended Users
- Computational biologists evaluating gRNA candidates.
- Research engineers building CRISPR scoring/validation pipelines.
- Reproducibility reviewers validating published CRISPR signal claims.

### Primary Outputs
- CRISPR scoring APIs and CLIs (`applications/`, `wave_crispr_signal/`).
- Reproducible validation/proof scripts (`proof_pack/`, `tests/`).
- Core methodology documentation (`docs/`).

### Non-Goals
- This repo is **not** a general cross-domain Z-framework lab.
- Non-CRISPR exploratory work (MRI/FUS/pain and related side projects) is retained under `legacy/` for archival reproducibility.

## Keep/Archive Rules

- `core`: CRISPR scoring/design/validation code and tests stay in mainline.
- `support`: build/test/tooling/data helpers stay in mainline.
- `research-legacy`: non-CRISPR exploratory work is kept under `legacy/` (read-only by default).
- `delete/regenerate`: generated artifacts and machine-specific files are removed from root and ignored.

## Entrypoint Contract

- Primary package and CRISPR-facing CLI behaviors remain stable.
- Existing imports from moved non-CRISPR modules continue to work through compatibility wrappers in `experiments/`.
- New contributions must map to CRISPR outcomes or be added under `legacy/`.

## Repository Layout

- `wave_crispr_signal/`: core CRISPR feature and scoring package.
- `applications/`: user-facing CRISPR applications and CLIs.
- `experiments/`: active CRISPR experiments plus compatibility wrappers.
- `proof_pack/`: reproducible validation and proof scripts.
- `data/` and `test_data/`: datasets and fixtures.
- `docs/`: canonical docs for current scope.
- `legacy/`: archived non-CRISPR research and side projects.
- `tools/`: repo maintenance and policy tooling.
- `tests/`: regression and policy tests.

## Contributor Decision Tree

1. Does this change advance CRISPR scoring/design/validation?
   - Yes: implement in mainline (`wave_crispr_signal/`, `applications/`, `experiments/`, `proof_pack/`, `tests/`).
   - No: place under `legacy/`.
2. Is the output generated or environment-specific?
   - Yes: keep out of versioned root; add ignore rule if needed.
3. Does the change alter public CRISPR interfaces?
   - Yes: include migration note and compatibility test updates.

## Quick Start

```bash
python -m pip install -r requirements.txt
python -m pytest -q

# Example: phase-weighted scorecard
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

## Canonical Docs

- `docs/INDEX.md`: documentation index and navigation.
- `docs/INVENTORY_CLASSIFICATION.md`: current structure classification (`core/support/research-legacy/delete-regenerate`).
- `legacy/README.md`: archived scope and retention policy.

## License

MIT
