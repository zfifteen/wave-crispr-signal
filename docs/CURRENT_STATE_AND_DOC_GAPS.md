# Current State and Documentation Alignment

Date: 2026-03-08

## Active Runtime Surfaces

- `wave_crispr_signal/`
- `applications/phase_weighted_scorecard_cli.py`
- `applications/crispr_cli.py`
- `applications/genomic_disruption_api.py`
- `web_apps/demo_mvp/app.py`

## Current Functional Truth

- `phase_weighted_scorecard_cli.py`: implemented `score`, `batch`, `analyze` commands.
- `crispr_cli.py`: implemented `design`, `score`, `batch-score` commands.
- `genomic_disruption_api.py` CLI:
  - implemented: `score`, `design`
  - not implemented in CLI execution path: `batch`, `offtarget`
- `web_apps/demo_mvp/app.py`: importable FastAPI app with `/`, `/health`, `/info`, `/analyze`.

## Documentation Status

The active documentation set has been hard-pruned and rewritten to match implemented behavior.

Known implemented limitation that is intentionally documented:
- `genomic_disruption_api.py` exposes `batch` and `offtarget` parser options, but those two commands are not wired in the CLI dispatcher.

No additional documentation/implementation mismatches are currently tracked in the retained docs.
