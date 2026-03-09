# Documentation Index

This index contains the active documentation set for the current application surface.

## Core Documents

- `CURRENT_STATE_AND_DOC_GAPS.md` - current runtime state and documented limitations.
- `CONTRIBUTING.md` - contribution workflow and validation contract.
- `PRIMARY_CLI_VALIDATION.md` - manual validation commands for active CLIs and demo app import surface.
- `PHASE_WEIGHTED_SCORECARD.md` - authoritative technical spec for phase-weighted scoring implementation.
- `MODEL_NOVELTY_STATEMENT.md` - detailed novelty boundary and model-level claim framing.
- `../validation/ontarget_gate/README.md` - on-target comparator gate usage and outputs.
- `../validation/ontarget_gate/PROTOCOL_V3.md` - decision-grade comparator gate protocol.
- `../validation/ontarget_gate/BOARD_FUNDING_ADDENDUM.md` - board conditional-approval controls for the 4-week recovery cycle.

## Runtime Entry Points

```bash
python3 applications/phase_weighted_scorecard_cli.py --help
python3 applications/crispr_cli.py --help
python3 applications/genomic_disruption_api.py --help
python3 -c "from web_apps.demo_mvp.app import app; print(sorted([r.path for r in app.router.routes]))"
```
