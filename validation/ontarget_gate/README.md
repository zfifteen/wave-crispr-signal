# Fast Go/No-Go On-Target Validation Gate

The active external-comparator protocol is `PROTOCOL_V3.md`.
Board-facing conditional approval controls are in `BOARD_FUNDING_ADDENDUM.md`.

- v3 adds required `baseline_c_crisprpred` comparator wiring.
- v1/v2 outputs are historical and non-authoritative for comparator claims.

## Run (v3)

```bash
cd /Users/velocityworks/IdeaProjects/wave-crispr-signal
python3 validation/ontarget_gate/scripts/run_gate_v3.py
```

## Required Comparator Setup

```bash
cd /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/comparators/crisprpred
./setup_venv.sh
```

Gate v3 is fail-closed if comparator preconditions are not met.

## Outputs (v3)

Written to `validation/ontarget_gate/outputs/`:

- `locked_dataset_manifest_v3.json`
- `locked_split_manifest_v3.csv`
- `decision_split_manifest_v3_clean.csv`
- `exploratory_split_manifest_v3_raw.csv`
- `locked_schema_manifest_v3.csv`
- `gate_results_v3.json`
- `gate_report_v3.md`
- `week1_checkpoint_v3.json`

## Comparator Policy (v3)

- Required comparator: `baseline_c_crisprpred`
- Hard decision thresholds:
  - primary holdout: `delta(model - baseline_c) >= +0.01`
  - external holdout: `delta(model - baseline_c) >= +0.01`
- CIs are diagnostic-only in this phase.
- Authoritative decisions use the clean manifest only.
- Raw manifest results are exploratory and non-authoritative.
- `gate_results_v3.json` now includes decision/exploratory diagnostics blocks for ranking health, GC coupling, and subgroup deltas.

## Funding Recovery Governance (Current)

- Recovery cycle is conditionally approved under pre-registered controls.
- Hard caps apply (time, team, ablation count, compute envelope).
- Week-1 checkpoint requires explicit `CONTINUE` or `STOP`.
- Any short-cycle exception requires re-approval and is capped at 10 business days.
