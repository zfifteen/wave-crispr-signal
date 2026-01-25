# Falsification Experiment: Wave-CRISPR Signal Hypotheses

**Experiment ID**: PR-156-falsify-crispr-hypothesis  
**Created**: 2026-01-25  
**Purpose**: Attempt to falsify two key hypotheses regarding the geodesic resolution function θ′(n,k) and curvature weighting κ(n) in CRISPR guide efficiency prediction.

## Overview

This experiment implements rigorous falsification testing for:

1. **Hypothesis 1**: The geodesic resolution function θ′(n,k) with k≈0.3, embedding golden ratio phasing, quantifies single-nucleotide mutation disruptions via Δentropy and sidelobe metrics, improving CRISPR guide efficiency prediction by ΔROC-AUC +0.047 over baselines like RuleSet3.

2. **Hypothesis 2**: The sequence diversity-dependent curvature weighting κ(n) = d(n) · ln(n+1)/e² creates a Z-invariant scoring abstraction that robustly handles variable-length DNA sequences, revealing emergent periodicity patterns in off-target profiling.

## Falsification Criteria

### Hypothesis 1
- **Falsified if**: ΔROC-AUC ≤ 0 or p > 0.05
- **Method**: Benchmark spectral features (with θ′(n,k) phase weighting) against baseline features using ROC-AUC on CRISPR efficiency prediction
- **Statistical tests**: 
  - 10-fold cross-validation
  - Paired t-test
  - Bootstrap confidence intervals (1,000+ resamples)

### Hypothesis 2
- **Falsified if**: 
  - High variance in Z-scores across lengths (ANOVA p < 0.05), OR
  - No emergent periodicity (autocorrelation < 0.1)
- **Method**: Evaluate Z-scores across variable-length sequences (20-100 nt)
- **Statistical tests**:
  - ANOVA for length dependence
  - Autocorrelation analysis for periodicity

## Scientific Gates

All experiments adhere to the repository's scientific gates:

✅ **Human DNA only**: GRCh38/hg38 reference, A/C/G/T/N validation  
✅ **Fail-fast validation**: Immediate ValueError on invalid nucleotides  
✅ **Z invariants**: Discrete domain Z = A(B / e²)  
✅ **Geometric resolution**: θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3  
✅ **Reproducibility**: Fixed seed, exact version pinning  
✅ **Statistical validity**: Bootstrap CI, permutation tests, FDR correction  

## Usage

### Hypothesis 1: CRISPR Guide Efficiency Prediction

```bash
# Basic run with synthetic data
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --seed 42 \
  --n-samples 100 \
  --n-bootstrap 1000 \
  --n-folds 10 \
  --k-parameter 0.3 \
  --output-dir results/PR-156-falsify-hypothesis1

# Quick test (reduced parameters)
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --seed 42 \
  --n-samples 50 \
  --n-bootstrap 100 \
  --n-folds 5
```

**Expected Runtime**: ~2-5 minutes (100 samples, 1000 bootstrap)

### Hypothesis 2: Z-Invariance Testing

```bash
# Basic run with variable-length sequences
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --seed 42 \
  --min-length 20 \
  --max-length 100 \
  --length-step 10 \
  --n-per-length 20 \
  --k-parameter 0.3 \
  --output-dir results/PR-156-falsify-hypothesis2

# Test with preserved motif
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --seed 42 \
  --min-length 20 \
  --max-length 100 \
  --length-step 10 \
  --n-per-length 20 \
  --test-motif
```

**Expected Runtime**: ~1-2 minutes (9 length groups, 20 sequences each)

### Smoke Test (CI-friendly)

```bash
# Hypothesis 1 smoke test
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --n-samples 20 --n-bootstrap 50 --n-folds 3

# Hypothesis 2 smoke test
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --min-length 20 --max-length 40 --length-step 10 --n-per-length 5
```

**Expected Runtime**: <5 seconds total

## Outputs

Results are saved as JSON to the specified output directory:

### Hypothesis 1 Results

```json
{
  "timestamp": "2026-01-25T...",
  "parameters": { "seed": 42, "k_parameter": 0.3, ... },
  "cross_validation": {
    "spectral_auc_mean": 0.xxx,
    "baseline_auc_mean": 0.xxx,
    "delta_auc": 0.xxx,
    "p_value": 0.xxx
  },
  "bootstrap": {
    "delta_auc_mean": 0.xxx,
    "ci_lower": 0.xxx,
    "ci_upper": 0.xxx
  },
  "falsification": {
    "falsified": false,
    "reasons": []
  }
}
```

### Hypothesis 2 Results

```json
{
  "timestamp": "2026-01-25T...",
  "parameters": { "seed": 42, "k_parameter": 0.3, ... },
  "invariance_test": {
    "f_statistic": 0.xxx,
    "p_value": 0.xxx,
    "invariant": true
  },
  "periodicity_test": {
    "max_autocorrelation": 0.xxx,
    "periodicity_detected": true
  },
  "falsification": {
    "falsified": false,
    "reasons": []
  }
}
```

## Implementation Details

### Hypothesis 1: Spectral Features

The experiment computes the following spectral features using θ′(n,k) phase weighting:

1. **Δentropy**: Change in spectral entropy between wild-type and mutant sequences
2. **Δsidelobes**: Change in number of spectral sidelobes
3. **Frequency shift**: Shift in dominant frequency component

These are compared against baseline features:
- GC content (overall, first 5 nt, last 5 nt)
- Sequence length

### Hypothesis 2: Z-Invariance

The experiment computes Z-scores as:

```
Z = S(entropy / κ)
```

where:
- `entropy` = spectral entropy from phase-weighted FFT
- `κ(n)` = `d(n) · ln(n+1) / e²` (curvature weight)
- `d(n)` = sequence diversity (normalized Shannon entropy of base composition)
- `S(x)` = sigmoid aggregator `1 / (1 + exp(-κ·x))`

## External References

### Datasets
- **Doench 2016**: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4745674/
- **Hsu 2013**: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3773023/
- **Kim 2025**: (if available in repository data)

### Baselines
- **RuleSet3**: https://github.com/MicrosoftResearch/Azimuth

### Related Work
- **Fourier methods in DNA**: https://pubmed.ncbi.nlm.nih.gov/16177379/
- **Golden ratio in biology**: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4345425/

## Limitations

### Current Implementation
- Uses **synthetic data** for demonstration; real datasets (Doench 2016, Kim 2025) should be added
- RuleSet3 baseline not yet integrated; uses simple GC-content features as proxy
- Sample sizes are configurable but default to small values for quick testing

### Future Enhancements
1. Integration with real CRISPR efficiency datasets
2. Direct RuleSet3 comparison (requires Azimuth package)
3. Extended statistical tests (e.g., Benjamini-Hochberg FDR)
4. Visualization of results (ROC curves, Z-score distributions)
5. Off-target profiling dataset integration

## Contributing

To extend this experiment:

1. Add real CRISPR datasets to `data/` directory
2. Implement RuleSet3 baseline in `falsify_hypothesis1.py`
3. Add visualization scripts
4. Create additional falsification tests (e.g., different k values, null models)

## Citation

If you use this experiment in your research, please cite:

```
@software{wave_crispr_signal_2026,
  title={Wave-CRISPR Signal: Falsification Experiment for Geodesic Resolution Hypotheses},
  author={Z Framework Contributors},
  year={2026},
  url={https://github.com/zfifteen/wave-crispr-signal}
}
```

## License

Research use only. See repository LICENSE for details.
