# Model Thesis Closeout (NO-GO)

Date: 2026-03-09  
Status: Closed

## Decision

The current wave-crispr-signal model thesis is closed as `NO-GO` for continued incremental investment.

This closeout applies to the current model approach and recovery variants evaluated under locked Gate v3 conditions.

## Basis

- Authoritative Gate v3 clean-lane decision: `NO-GO` against `baseline_c_crisprpred`.
- Follow-on recovery experiments did not produce dual-holdout positive lift.
- Short-cycle exception criteria were not met.

Reference artifacts:

- `outputs/gate_results_v3.json`
- `outputs/gate_report_v3.md`
- `outputs/week1_checkpoint_v3.json`
- `outputs/recovery_ablation_v1.json`
- `outputs/recovery_param_sweep_v1.json`
- `outputs/recovery_model_family_v1.json`

## What Is Closed vs Retained

Closed:

- Further incremental tuning of the current model thesis as a funded continuation track.

Retained:

- Gate v3 protocol and comparator framework.
- Overlap-safe split logic and decision/reporting artifacts.
- CLI/runtime codebase for non-thesis work.

## Re-entry Rule

Any future attempt must be opened as a new hypothesis track with:

1. a fresh written hypothesis,
2. explicit success criteria locked before runs,
3. a separate experiment ID/plan, and
4. re-approval for funding.
