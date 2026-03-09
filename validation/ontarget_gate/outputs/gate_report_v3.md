# Fast Go/No-Go Validation Report v3 (On-Target)

Generated: 2026-03-09T07:21:47.250070+00:00

## Decision: **NO-GO**
- reason: baseline_c_threshold_failed

## Preconditions
- split_integrity_clean_ok: True
- comparator_self_check_ok: True
- overlap_audit_clean_ok: True
- holdout_min_n_ok: True
- holdout_min_n_required: 200

## Manifest Lineage
- decision_clean_manifest: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/outputs/decision_split_manifest_v3_clean.csv
- exploratory_raw_manifest: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/outputs/exploratory_split_manifest_v3_raw.csv
- dropped_primary_overlap: 1688
- dropped_external_overlap: 0
- dropped_total_overlap: 1688

## Baseline C Comparator
- slot: baseline_c
- name: baseline_c_crisprpred
- version: v0.1-scaffold
- self_check_ok: True
- self_check_message: ok
- provenance_path: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/comparators/crisprpred/PROVENANCE.md
- checksums_path: /Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/comparators/crisprpred/checksums.sha256

## Overlap Audit
- clean_status: PASS
- clean_message: no_overlap
- clean_overlap_count: 0
- raw_status: FAIL
- raw_message: overlap_detected
- raw_overlap_count: 1688

## Decision-Grade Metrics (Clean Manifests)
### dev_val (doench2016_hg19, n=414)
- baseline_a: spearman=-0.0423, mse=2.7625
- baseline_b: spearman=0.0146, mse=4.6956
- baseline_c: spearman=-0.0474, mse=6.2579
- model: spearman=-0.0402, mse=5.4846
- delta(model - baseline_a) spearman: 0.0020 (95% CI [-0.1074, 0.1091])
- delta(model - baseline_b) spearman: -0.0548
- delta(model - baseline_c) spearman: 0.0072 (95% CI [-0.0809, 0.0870])

### primary_holdout (doench2016_hg19, n=529)
- baseline_a: spearman=-0.0240, mse=5.3951
- baseline_b: spearman=0.0429, mse=13.0235
- baseline_c: spearman=-0.0510, mse=16.5125
- model: spearman=-0.0702, mse=14.8312
- delta(model - baseline_a) spearman: -0.0462 (95% CI [-0.1458, 0.0563])
- delta(model - baseline_b) spearman: -0.1131
- delta(model - baseline_c) spearman: -0.0192 (95% CI [-0.1067, 0.0637])

### external_holdout (hart2016_hela_lib1_avg, n=4256)
- baseline_a: spearman=-0.0435, mse=9.2447
- baseline_b: spearman=-0.0042, mse=4.1856
- baseline_c: spearman=0.0125, mse=3.2801
- model: spearman=-0.0019, mse=3.6672
- delta(model - baseline_a) spearman: 0.0416 (95% CI [0.0031, 0.0820])
- delta(model - baseline_b) spearman: 0.0023
- delta(model - baseline_c) spearman: -0.0144 (95% CI [-0.0460, 0.0165])

## Exploratory Metrics (Raw Manifests, Non-Authoritative)
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
- hard_precondition_dev_not_materially_negative: PASS
- primary_delta_vs_baseline_c_min_0p01: FAIL
- external_delta_vs_baseline_c_min_0p01: FAIL
