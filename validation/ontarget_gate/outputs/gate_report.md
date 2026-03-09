# Fast Go/No-Go Validation Report (On-Target Efficiency)

Generated: 2026-03-09T05:48:18.560336+00:00

## Decision: **INCONCLUSIVE**

## Dataset Counts
- `doench2016_gecko1_lenticrispr_hg19` (dev_pool): 62834 rows
- `doench2016_hg19` (primary_holdout): 3804 rows
- `hart2016_hela_lib1_avg` (external_holdout): 4256 rows

## Split Counts
- dev_train: 50030
- dev_val: 12804
- primary_holdout: 3804
- external_holdout: 4256

## Metrics (Spearman / MSE)
### dev_val (n=12804)
- baseline_a: spearman=0.2293, mse=0.3391
- baseline_b: spearman=0.0201, mse=1.8562
- model: spearman=0.0099, mse=2.5249
- delta(model - baseline_a) spearman: -0.2194 (95% CI [-0.2443, -0.1948])
- delta(model - baseline_b) spearman: -0.0102

### primary_holdout (n=3804)
- baseline_a: spearman=-0.0273, mse=4.4130
- baseline_b: spearman=0.0508, mse=7.3326
- model: spearman=0.0274, mse=8.2878
- delta(model - baseline_a) spearman: 0.0547 (95% CI [0.0127, 0.0948])
- delta(model - baseline_b) spearman: -0.0234

### external_holdout (n=4256)
- baseline_a: spearman=-0.0490, mse=9.0341
- baseline_b: spearman=-0.0042, mse=4.1856
- model: spearman=-0.0019, mse=3.6672
- delta(model - baseline_a) spearman: 0.0471 (95% CI [0.0098, 0.0843])
- delta(model - baseline_b) spearman: 0.0023

## Rule Evaluation
- Primary delta >= 0.03: 0.0547 -> PASS
- Primary CI excludes 0: [0.0127, 0.0948] -> PASS
- External delta >= 0.01: 0.0471 -> PASS
- No underperform vs baseline_b (primary): -0.0234 -> FAIL
- No underperform vs baseline_b (external): 0.0023 -> PASS
