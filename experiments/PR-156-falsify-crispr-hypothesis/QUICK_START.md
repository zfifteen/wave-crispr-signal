# PR-156 Falsification Experiments: Quick Start Guide

This guide provides examples and quick start commands for running the PR-156 falsification experiments.

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Example Outputs](#example-outputs)
4. [Interpreting Results](#interpreting-results)
5. [Common Use Cases](#common-use-cases)

---

## Overview

The PR-156 experiment tests two hypotheses:

**H1**: θ′(n,k) improves CRISPR guide efficiency prediction  
**H2**: κ(n) creates Z-invariant scoring for variable-length sequences

**Falsification Criteria**:
- H1 is falsified if: ΔROC-AUC ≤ 0 OR p > 0.05 OR 95% CI includes zero
- H2 is falsified if: ANOVA p < 0.05 (not invariant) OR max autocorr < 0.1 (no periodicity)

---

## Quick Start

### 1. Smoke Test (Fastest, <5 seconds)

```bash
make pr156-falsification-smoke
```

This runs both hypotheses with minimal parameters for CI testing.

### 2. Full Run (~5 minutes)

```bash
make run-pr156-falsification
```

This runs both hypotheses with full parameters (100 samples, 1000 bootstrap resamples).

### 3. Individual Hypothesis Tests

**Hypothesis 1 only** (ROC-AUC comparison):
```bash
make run-pr156-h1
```

**Hypothesis 2 only** (Z-invariance):
```bash
make run-pr156-h2
```

### 4. Custom Parameters

**Hypothesis 1 with custom parameters**:
```bash
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --seed 42 \
  --n-samples 200 \
  --n-bootstrap 2000 \
  --n-folds 10 \
  --k-parameter 0.3 \
  --output-dir results/my-h1-test
```

**Hypothesis 2 with custom parameters**:
```bash
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --seed 42 \
  --min-length 15 \
  --max-length 50 \
  --length-step 5 \
  --n-per-length 30 \
  --k-parameter 0.3 \
  --test-motif \
  --output-dir results/my-h2-test
```

---

## Example Outputs

### Hypothesis 1: Console Output

```
======================================================================
Hypothesis 1 Falsification Experiment
======================================================================
Seed: 42
k parameter: 0.3
N samples: 100
N bootstrap: 1000
N folds: 10

Generating synthetic CRISPR guide data...
Generated 100 guide pairs
Label distribution: [52 48]

Computing spectral features with θ′(n,k) phase weighting...
Spectral features shape: (100, 3)
Baseline features shape: (100, 4)

Running 10-fold cross-validation...
Spectral AUC (mean ± std): 0.6234 ± 0.0821
Baseline AUC (mean ± std): 0.5891 ± 0.0742
ΔROC-AUC (CV): 0.0343

Paired t-test: t=2.1456, p=0.0312

Computing bootstrap CI for ΔROC-AUC (1000 resamples)...
Bootstrap ΔROC-AUC: 0.0351 ± 0.0124
95% CI: [0.0112, 0.0589]

======================================================================
FALSIFICATION ANALYSIS
======================================================================
RESULT: Hypothesis 1 is NOT FALSIFIED (supported by data)
  - ΔROC-AUC = 0.0351 > 0
  - p-value = 0.0312 < 0.05
  - 95% CI [0.0112, 0.0589] excludes zero

Results saved to: results/PR-156-falsify-hypothesis1/hypothesis1_results.json
```

### Hypothesis 1: JSON Output

```json
{
  "timestamp": "2026-01-25T17:30:00.000000",
  "parameters": {
    "seed": 42,
    "k_parameter": 0.3,
    "n_samples": 100,
    "n_bootstrap": 1000,
    "n_folds": 10
  },
  "cross_validation": {
    "spectral_auc_mean": 0.6234,
    "spectral_auc_std": 0.0821,
    "baseline_auc_mean": 0.5891,
    "baseline_auc_std": 0.0742,
    "delta_auc": 0.0343,
    "t_statistic": 2.1456,
    "p_value": 0.0312
  },
  "bootstrap": {
    "delta_auc_mean": 0.0351,
    "delta_auc_std": 0.0124,
    "ci_lower": 0.0112,
    "ci_upper": 0.0589
  },
  "falsification": {
    "falsified": false,
    "reasons": []
  }
}
```

### Hypothesis 2: Console Output

```
======================================================================
Hypothesis 2 Falsification Experiment
======================================================================
Seed: 42
k parameter: 0.3
Length range: 20-100 (step 10)
N per length: 20
Test with motif: False

Generating sequences for 9 length groups...
Using random sequences
Generated 180 total sequences

Testing Z-score invariance across sequence lengths...
ANOVA F-statistic: 1.8234
ANOVA p-value: 0.0823
Overall Z-score variance: 0.0123

Mean Z-scores by length:
  Length 20: 0.9812
  Length 30: 0.9834
  Length 40: 0.9851
  Length 50: 0.9869
  Length 60: 0.9881
  Length 70: 0.9892
  Length 80: 0.9901
  Length 90: 0.9909
  Length 100: 0.9916

Testing for emergent periodicity in Z-scores...
Max autocorrelation: 0.2341
Periodicity detected: True

======================================================================
FALSIFICATION ANALYSIS
======================================================================
RESULT: Hypothesis 2 is NOT FALSIFIED (supported by data)
  - Z-scores are invariant (ANOVA p=0.0823 >= 0.05)
  - Emergent periodicity detected (max autocorr=0.2341 > 0.1)

Results saved to: results/PR-156-falsify-hypothesis2/hypothesis2_results.json
```

### Hypothesis 2: JSON Output

```json
{
  "timestamp": "2026-01-25T17:32:00.000000",
  "parameters": {
    "seed": 42,
    "k_parameter": 0.3,
    "min_length": 20,
    "max_length": 100,
    "length_step": 10,
    "n_per_length": 20,
    "test_motif": false
  },
  "invariance_test": {
    "f_statistic": 1.8234,
    "p_value": 0.0823,
    "invariant": true,
    "overall_variance": 0.0123,
    "mean_z_by_length": {
      "20": 0.9812,
      "30": 0.9834,
      "40": 0.9851,
      "50": 0.9869,
      "60": 0.9881,
      "70": 0.9892,
      "80": 0.9901,
      "90": 0.9909,
      "100": 0.9916
    }
  },
  "periodicity_test": {
    "max_autocorrelation": 0.2341,
    "periodicity_detected": true
  },
  "falsification": {
    "falsified": false,
    "reasons": []
  }
}
```

---

## Interpreting Results

### Understanding Falsification

**Hypothesis is FALSIFIED** = Data does not support the claim  
**Hypothesis is NOT FALSIFIED** = Data supports the claim (but doesn't prove it!)

### Hypothesis 1 Interpretation

**Key Metrics**:
- **ΔROC-AUC**: Difference in AUC between spectral and baseline features
  - Positive = spectral features are better
  - Negative = baseline features are better
- **p-value**: Statistical significance of the difference
  - <0.05 = significant difference
  - ≥0.05 = no significant difference
- **95% CI**: Confidence interval for ΔROC-AUC
  - Excludes zero = significant positive/negative difference
  - Includes zero = could be no difference

**Decision Rules**:
1. ✅ NOT FALSIFIED if: ΔROC-AUC > 0 AND p < 0.05 AND CI excludes zero
2. ❌ FALSIFIED if any criterion fails

### Hypothesis 2 Interpretation

**Key Metrics**:
- **ANOVA p-value**: Tests if Z-scores differ significantly across lengths
  - ≥0.05 = invariant (good)
  - <0.05 = not invariant (bad)
- **Max autocorrelation**: Measures periodicity strength
  - >0.1 = periodicity detected (good)
  - ≤0.1 = no periodicity (bad)

**Decision Rules**:
1. ✅ NOT FALSIFIED if: ANOVA p ≥ 0.05 AND max autocorr > 0.1
2. ❌ FALSIFIED if either criterion fails

---

## Common Use Cases

### 1. Testing Different k Values

Compare performance across different k parameters:

```bash
for k in 0.2 0.25 0.3 0.35 0.4; do
  python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
    --k-parameter $k \
    --output-dir results/h1-k-${k}
done
```

### 2. Testing with Preserved Motif

Test if Z-invariance holds when core sequence is preserved:

```bash
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --test-motif \
  --output-dir results/h2-with-motif
```

### 3. Large-Scale Validation

Run with increased sample size and bootstrap resamples:

```bash
python experiments/PR-156-falsify-crispr-hypothesis/run_all_falsification_tests.py \
  --seed 42 \
  --n-samples 500 \
  --n-bootstrap 5000 \
  --n-folds 10 \
  --min-length 15 \
  --max-length 150 \
  --length-step 5 \
  --n-per-length 50 \
  --output-dir results/large-scale-validation
```

### 4. Quick Parameter Sweep

Test multiple configurations quickly:

```bash
# Quick sweep over k values
make pr156-falsification-smoke  # Baseline

# Then test variations
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --n-samples 50 --n-bootstrap 100 --k-parameter 0.25
  
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --n-samples 50 --n-bootstrap 100 --k-parameter 0.35
```

---

## Troubleshooting

### Issue: "TypeError: Object of type bool_ is not JSON serializable"

**Solution**: This has been fixed in the latest version. Update to the latest code.

### Issue: Smoke test fails with non-zero exit code

**Solution**: This is expected behavior when hypotheses are falsified with synthetic data. The experiment still completes successfully. Use `|| true` in scripts to ignore the exit code if needed.

### Issue: "ModuleNotFoundError: No module named 'applications.phase_weighted_scorecard'"

**Solution**: Ensure you're running from the repository root directory:
```bash
cd /path/to/wave-crispr-signal
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py
```

### Issue: Different results each run

**Solution**: Ensure you're using the same `--seed` value for reproducibility:
```bash
--seed 42  # Always use the same seed
```

---

## Next Steps

After running the experiments:

1. **Review results**: Check JSON files in `results/PR-156-*/`
2. **Fill out report**: Use `FALSIFICATION_REPORT_TEMPLATE.md` to document findings
3. **Integrate real data**: Replace synthetic data with Doench 2016, Kim 2025, etc.
4. **Extend analysis**: Add visualization, additional statistical tests, RuleSet3 comparison

---

**Questions or Issues?**  
See [`experiments/PR-156-falsify-crispr-hypothesis/README.md`](README.md) for detailed documentation.
