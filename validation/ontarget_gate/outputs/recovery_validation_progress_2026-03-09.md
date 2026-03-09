# Recovery Validation Progress (2026-03-09)

## Objective
Run immediate, execution-first validation experiments on locked Gate v3 clean-lane data to test whether practical in-scope model adjustments can recover holdout ranking lift vs `baseline_c`.

## Experiments Executed

1. Recovery Ablation v1 (`run_recovery_ablation_v1.py`)
- `A0`: current model
- `A1`: z-blend (model + baseline_a)
- `A2`: GC-debias + z-blend
- `A3`: stacked ridge meta-model

2. Core hyperparameter sweep (inline run)
- `k` in `[0.1, 0.6]`
- `phi_period` in `[19, 23]`
- evaluated on dev + holdouts against `baseline_c`

3. Alternative model-family stackers (inline run)
- ridge / gradient boosting / random forest stackers over `[model_score, baseline_a, gc, homopolymer]`

## Results Snapshot

### Ablation v1
- Best dev candidate: `A2_gc_debias_plus_zblend`
- Deltas vs baseline_c (Spearman):
  - dev: `+0.0484`
  - primary: `-0.0119`
  - external: `-0.0172`
- Conclusion: dev lift did not translate to holdouts.

### Hyperparameter sweep
- Best strict holdout config found: `k=0.1`, `phi_period=19..23`
- Deltas vs baseline_c:
  - dev: `+0.0057`
  - primary: `-0.0144`
  - external: `-0.0149`
- Conclusion: parameter tuning alone cannot clear holdout gap.

### Model-family stackers
- Best primary-holdout performance: tree-based stackers (primary ~ `+0.029`)
- External holdout remained strongly negative for all variants (about `-0.054` to `-0.059`).
- Conclusion: variants improve one holdout while degrading the other; no robust generalization.

## Current Decision Signal
No tested in-scope variant currently satisfies dual-holdout lift criteria. Failure mode appears to be generalization mismatch across holdout domains rather than simple scalar tuning.

## Artifacts
- `validation/ontarget_gate/outputs/recovery_ablation_v1.json`
- `validation/ontarget_gate/outputs/recovery_ablation_v1.md`
- `validation/ontarget_gate/outputs/recovery_param_sweep_v1.json`
- `validation/ontarget_gate/outputs/recovery_model_family_v1.json`

