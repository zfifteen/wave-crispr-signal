# wave-crispr-signal

Focused CRISPR signal-analysis toolkit for guide scoring and disruption analysis.

## Active Runtime Surfaces

- `wave_crispr_signal/` (core sequence, spectral, and stats utilities)
- `applications/phase_weighted_scorecard_cli.py`
- `applications/crispr_cli.py`
- `applications/genomic_disruption_api.py`
- `web_apps/demo_mvp/app.py`

## Quick Start

```bash
python3 -m pip install -r requirements.txt

python3 applications/phase_weighted_scorecard_cli.py --help
python3 applications/crispr_cli.py --help
python3 applications/genomic_disruption_api.py --help
```

Current CLI truth for `applications/genomic_disruption_api.py`:
- `score` is implemented.
- `design` is implemented.
- `batch` and `offtarget` are parser options, but CLI execution returns `Command <name> not yet implemented in CLI`.

## Validation Model

- Validation is change-scoped.
- For behavior changes, default expectation is:
  - one targeted automated test, and
  - one manual verification command with observed outcome.

## Validation Milestone (Gate v3)

- External comparator gate is active under `validation/ontarget_gate/`.
- Gate v3 uses two manifest lanes:
  - decision-grade clean manifest (authoritative),
  - raw exploratory manifest (non-authoritative).
- Current authoritative Gate v3 outcome on clean holdouts: `NO-GO` against `baseline_c_crisprpred`.
- Current model thesis status: **closed (NO-GO)** for incremental continuation.
- See:
  - `validation/ontarget_gate/README.md`
  - `validation/ontarget_gate/PROTOCOL_V3.md`
  - `validation/ontarget_gate/BOARD_FUNDING_ADDENDUM.md`
  - `validation/ontarget_gate/CLOSEOUT_NO_GO.md`
  - `validation/ontarget_gate/outputs/gate_report_v3.md`

## Canonical Docs

- `docs/INDEX.md`
- `docs/CURRENT_STATE_AND_DOC_GAPS.md`
- `docs/CONTRIBUTING.md`
- `docs/PRIMARY_CLI_VALIDATION.md`
- `docs/PHASE_WEIGHTED_SCORECARD.md`
- `docs/MODEL_NOVELTY_STATEMENT.md`

## License

MIT
