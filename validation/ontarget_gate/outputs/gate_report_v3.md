# Fast Go/No-Go Validation Report v3 (On-Target)

Generated: 2026-03-09T06:57:36.237913+00:00

## Decision: **INCONCLUSIVE**
- reason: overlap_audit_failed_or_unavailable

## Preconditions
- split_integrity_ok: True
- comparator_self_check_ok: True
- overlap_audit_ok: False
- holdout_min_n_ok: True
- holdout_min_n_required: 200

## Baseline C Comparator
- slot: baseline_c
- name: baseline_c_crisprpred
- version: v0.1-scaffold
- self_check_ok: True
- self_check_message: ok
- provenance_path: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/comparators/crisprpred/PROVENANCE.md
- checksums_path: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/comparators/crisprpred/checksums.sha256

## Overlap Audit
- status: FAIL
- message: overlap_detected
- overlap_count: 1688

## Per-Split Metrics (Spearman / MSE)
### dev_val (doench2016_hg19, n=414)
- baseline_a: spearman=-0.0423, mse=2.7625
- baseline_b: spearman=0.0146, mse=4.6956
- baseline_c: spearman=-0.0474, mse=6.2579
- model: spearman=-0.0402, mse=5.4846
- delta(model - baseline_a) spearman: 0.0020 (95% CI [-0.1074, 0.1091])
- delta(model - baseline_b) spearman: -0.0548
- delta(model - baseline_c) spearman: 0.0072 (95% CI [-0.0809, 0.0870])

### primary_holdout (doench2016_hg19, n=2217)
- baseline_a: spearman=-0.0010, mse=4.3547
- baseline_b: spearman=0.0647, mse=5.7667
- baseline_c: spearman=0.1563, mse=7.1157
- model: spearman=0.0711, mse=6.4274
- delta(model - baseline_a) spearman: 0.0721 (95% CI [0.0184, 0.1295])
- delta(model - baseline_b) spearman: 0.0064
- delta(model - baseline_c) spearman: -0.0851 (95% CI [-0.1268, -0.0392])

### external_holdout (hart2016_hela_lib1_avg, n=4256)
- baseline_a: spearman=-0.0435, mse=9.2447
- baseline_b: spearman=-0.0042, mse=4.1856
- baseline_c: spearman=0.0125, mse=3.2801
- model: spearman=-0.0019, mse=3.6672
- delta(model - baseline_a) spearman: 0.0416 (95% CI [0.0031, 0.0820])
- delta(model - baseline_b) spearman: 0.0023
- delta(model - baseline_c) spearman: -0.0144 (95% CI [-0.0460, 0.0165])

## Rule Evaluation (CI Diagnostic-Only in v3)
- criteria_not_evaluated_due_to_precondition_failure: true
