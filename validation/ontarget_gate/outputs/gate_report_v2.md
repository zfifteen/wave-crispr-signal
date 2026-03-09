# Fast Go/No-Go Validation Report v2 (On-Target)

Generated: 2026-03-09T06:08:16.936497+00:00

## Decision: **GO**

## Split Integrity
- status: PASS
- dev_train: 1173
- dev_val: 414
- primary_holdout: 2217
- external_holdout: 4256
- aux_train_weak: 62834

## Per-Split Metrics (Spearman / MSE)
### dev_val (doench2016_hg19, n=414)
- baseline_a: spearman=-0.0423, mse=2.7625
- baseline_b: spearman=0.0146, mse=4.6956
- model: spearman=-0.0402, mse=5.4846
- delta(model - baseline_a) spearman: 0.0020 (95% CI [-0.1074, 0.1091])
- delta(model - baseline_b) spearman: -0.0548

### primary_holdout (doench2016_hg19, n=2217)
- baseline_a: spearman=-0.0010, mse=4.3547
- baseline_b: spearman=0.0647, mse=5.7667
- model: spearman=0.0711, mse=6.4274
- delta(model - baseline_a) spearman: 0.0721 (95% CI [0.0184, 0.1295])
- delta(model - baseline_b) spearman: 0.0064

### external_holdout (hart2016_hela_lib1_avg, n=4256)
- baseline_a: spearman=-0.0435, mse=9.2447
- baseline_b: spearman=-0.0042, mse=4.1856
- model: spearman=-0.0019, mse=3.6672
- delta(model - baseline_a) spearman: 0.0416 (95% CI [0.0031, 0.0820])
- delta(model - baseline_b) spearman: 0.0023

## Aggregated Holdout Delta (No Mixed-Scale Pooled Spearman)
- weighted delta(model - baseline_a): 0.0521 (95% CI [0.0197, 0.0864])

## Rule Evaluation
- hard_precondition_dev_not_materially_negative: PASS
- primary_delta_min: PASS
- primary_ci_excludes_zero: PASS
- external_delta_min: PASS
- external_ci_excludes_zero: PASS
- directional_consistency: PASS
- no_underperform_vs_baseline_b_primary: PASS
- no_underperform_vs_baseline_b_external: PASS
- aggregate_delta_positive: PASS
- aggregate_ci_excludes_zero: PASS

## One-Page Decision Summary
- outcome: GO
- rationale: derived exclusively from locked v2 criteria in `PROTOCOL_V2.md`.
