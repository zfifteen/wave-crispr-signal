# wave-crispr-signal

Signal-theoretic CRISPR analytics toolkit for guide scoring, disruption analysis, and reproducible validation.

## Project Charter

### Mission
Build and maintain a focused CRISPR signal-analysis toolkit that turns sequence-level signals into reproducible guide-quality and disruption metrics.

### Intended Users
- Computational biologists evaluating gRNA candidates.
- Research engineers building CRISPR scoring and validation pipelines.
- Reproducibility reviewers validating CRISPR signal claims.

### Primary Outputs
- CRISPR scoring APIs and CLIs (`applications/`, `wave_crispr_signal/`).
- Reproducible analysis scripts (`proof_pack/`).
- Methodology and usage documentation (`docs/`).

### Non-Goals
- This repository is not a general cross-domain experiment lab.

## Working Model

- Keep active CRISPR work in mainline directories.
- Keep generated artifacts and machine-specific files out of versioned root.
- Introduce tests only when a specific add/remove/alter change requires them.
- Historical tests are not a blocking source of truth.
- CI is intentionally absent for now and will be reintroduced later as a minimal, change-scoped gate.

## Repository Layout

- `wave_crispr_signal/`: core CRISPR feature and scoring package.
- `applications/`: user-facing CRISPR applications and CLIs.
- `experiments/`: active CRISPR experiments.
- `proof_pack/`: reproducible validation and proof scripts.
- `data/` and `test_data/`: datasets and fixtures.
- `docs/`: canonical docs for current scope.
- `tools/`: repository utility scripts.

## Contributor Decision Tree

1. Does this change advance CRISPR scoring/design/validation?
   - Yes: implement in mainline.
   - No: do not add it here without explicit direction.
2. Is the output generated or environment-specific?
   - Yes: keep it out of versioned root and add ignore rules if needed.
3. Does the change alter behavior?
   - Yes: add change-coupled validation (test or explicit manual validation steps).

## Quick Start

```bash
python -m pip install -r requirements.txt

# Example: phase-weighted scorecard
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
```

## Canonical Docs

- `docs/INDEX.md`: documentation index and navigation.
- `docs/CONTRIBUTING.md`: contribution workflow.

## License

MIT
