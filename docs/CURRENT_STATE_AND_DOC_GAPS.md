# Current State and Documentation Alignment

Date: 2026-03-09

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

## Validation Gate Status (On-Target Comparator)

- Active gate implementation: `validation/ontarget_gate/scripts/run_gate_v3.py`
- Active protocol: `validation/ontarget_gate/PROTOCOL_V3.md`
- Active funding-governance addendum: `validation/ontarget_gate/BOARD_FUNDING_ADDENDUM.md`
- Active closeout record: `validation/ontarget_gate/CLOSEOUT_NO_GO.md`
- Comparator lane design:
  - decision-grade clean manifest (`decision_split_manifest_v3_clean.csv`) drives `GO/NO-GO`,
  - exploratory raw manifest (`exploratory_split_manifest_v3_raw.csv`) is report-only context.
- Current authoritative result:
  - `NO-GO` on clean holdouts vs `baseline_c_crisprpred` fixed-delta criteria.
- Current thesis state:
  - closed (`NO-GO`) for incremental continuation under this model approach.
- Current raw manifest remains overlap-contaminated and is explicitly non-authoritative.
