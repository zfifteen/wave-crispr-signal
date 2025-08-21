# Issue #50 Experiment Framework - Implementation Guide

## Overview

This document details the implementation of the Issue #50 experiment specification for validating Z-Framework spectral features in CRISPR guide prediction. The framework provides rigorous statistical validation of 6 falsifiable hypotheses across multiple scales with bio-anchored vs arbitrary encoding bands.

## Experimental Design

### Core Hypotheses Tested

1. **H1 (Lift vs Baseline)**: Z-spectral model improves predictive performance over baseline by ≥15% relative
2. **H2 (k* optimum)**: Performance is maximized near k≈0.3 for geodesic curvature
3. **H3 (Band non-inferiority)**: Arbitrary encoding is non-inferior to bio-anchored with margin δ = 1%
4. **H4 (Phase-bit utility)**: F-phase alternation feature reduces error by ≥5% relative
5. **H5 (Golden-proximity signal)**: |μ_Z − (φ−1)| correlates with guide outcome (partial r ≥ 0.20)
6. **H6 (Non-triviality)**: Under permutation, model lift collapses to 0

### Experimental Configuration

- **Scales**: N∈{10², 10³, 10⁴} sequences (optionally 10⁵ as exploratory)
- **Bands**: Bio-anchored vs Arbitrary complex encoding weights
- **Cross-validation**: Repeated stratified 5×2 CV grouped by gene
- **Statistical validation**: Bootstrap CIs (BCa 95%), permutation tests
- **Random seed**: 42 (for reproducibility)

## Implementation Structure

### Main Components

```
experiments/issue50_experiment.py    # Main experiment framework
tests/test_issue50_experiment.py     # Validation test suite
results/                            # Output directory
├── issue50_summary.json           # Complete results with CIs + p-values
├── issue50_table.csv              # Detailed run-by-run results
└── *.png                          # Visualization plots
logs/                              # Experiment logs
└── issue50_[timestamp].log        # Detailed execution logs
```

### Key Classes

#### `Issue50ExperimentFramework`
Main experimental orchestrator implementing:
- Dataset loading and synthetic data generation
- Z-Framework feature extraction (bio-anchored vs arbitrary)
- Hypothesis testing with statistical validation
- Results visualization and reporting

## Feature Engineering

### Z-Framework Spectral Features

The framework extracts spectral features using complex waveform encoding:

```python
# Complex encoding with geodesic curvature
waveform_component = base_weight * geodesic_factor * exp(1j * position_phase)

# Where:
# - base_weight: Bio-anchored or arbitrary complex weights
# - geodesic_factor: θ'(n,k) = φ·((n mod φ)/φ)^k 
# - position_phase: 2π * position / sequence_length
```

#### Extracted Features:
- **Spectral entropy**: Information content of frequency spectrum
- **Spectral flatness**: Geometric/arithmetic mean ratio
- **Peak analysis**: Density and variance of spectral peaks
- **Z-Framework metrics**: z_mean, z_variance, z_std from discrete mapping
- **Golden proximity**: δφ = |μ_Z − (φ−1)| distance measure
- **Phase bit**: F alternation pattern (period-2)
- **Baseline confounders**: GC content, k-mers, positional features

### Encoding Schemes

#### Bio-anchored Encoding
Uses nucleotide physicochemical properties:
- Polarizability (Å³)
- Hydrogen bonding capacity
- Base stacking energy (kcal/mol)
- Combined into complex weights with normalized scales

#### Arbitrary Control Encoding
Fixed seeded random complex weights with:
- Unit magnitude constraints
- Orthogonal axis distributions
- Reproducible generation for fair comparison

## Statistical Validation

### Cross-Validation Strategy
- **Repeated Stratified 5×2 CV**: 5 repeats of 2-fold CV
- **Group splitting**: By gene to prevent data leakage
- **Stratification**: By efficiency quantiles for balanced folds

### Bootstrap Confidence Intervals
- **BCa 95% CI**: Bias-corrected and accelerated intervals
- **N = 1,000**: Bootstrap resamples on CV fold results
- **Paired resampling**: Maintains fold structure for valid CIs

### Permutation Tests (H6)
- **Label permutation**: Shuffle efficiency within genes
- **Phase scrambling**: Randomize spectral phase features
- **N = 1,000**: Permutation iterations
- **Significance**: p < 0.01 with 99th percentile threshold

## Usage Examples

### Basic Experiment Execution

```python
from experiments.issue50_experiment import Issue50ExperimentFramework

# Initialize experiment
experiment = Issue50ExperimentFramework(random_seed=42)

# Run complete demonstration
results = experiment.run_complete_experiment()

# Results automatically saved to:
# - results/issue50_summary.json
# - results/issue50_table.csv  
# - results/issue50_demo.png
```

### Custom Dataset Testing

```python
# Load custom CRISPR dataset
df_custom = pd.read_csv('custom_crispr_data.csv')  # Must have 'sequence', 'efficiency'

# Run H1 hypothesis test
h1_results = experiment.run_hypothesis_h1_lift_vs_baseline(df_custom, scale=1000)

# Check pass/fail criteria
for band in ['bio_band', 'arbitrary_band']:
    rmr = h1_results[band]['rmr_mean']
    ci_lower = h1_results[band]['rmr_ci_lower']
    passes = h1_results[band]['passes_h1']
    print(f"{band}: RMR = {rmr:.3f}, 95% CI lower = {ci_lower:.3f}, Passes = {passes}")
```

### Reproducibility Validation

```python
# Ensure reproducible results
experiment1 = Issue50ExperimentFramework(random_seed=123)
experiment2 = Issue50ExperimentFramework(random_seed=123)

# Generate identical synthetic datasets
df1 = experiment1._create_synthetic_dataset(n_sequences=100)
df2 = experiment2._create_synthetic_dataset(n_sequences=100)

# Verify identical results
assert df1.equals(df2), "Datasets should be identical with same seed"
```

## Results Interpretation

### Pass/Fail Criteria

| Hypothesis | Pass Criteria | 
|------------|---------------|
| H1 | RMR ≥ 15% AND 95% CI lower bound ≥ 10% (both bands) |
| H2 | k=0.3 is top-1 OR within 1% of top performance |
| H3 | Non-inferiority: lower 95% CI of (Arb − Bio) ≥ −1% |
| H4 | Adding F reduces error by ≥5%; one-sided p < 0.01 |
| H5 | Partial r ≥ 0.20, p < 1e-4 after FDR correction |
| H6 | Permutation p < 0.01 AND observed > 99th percentile |

### Results Structure

```json
{
  "experiment_metadata": {
    "timestamp": "2025-01-20T10:30:00",
    "random_seed": 42,
    "scales": [100, 1000, 10000],
    "n_bootstrap": 1000
  },
  "hypothesis_tests": {
    "h1_scale_100": {
      "bio_band": {
        "rmr_mean": 0.18,
        "rmr_ci_lower": 0.12,
        "rmr_ci_upper": 0.24,
        "p_value": 0.003,
        "passes_h1": true
      },
      "arbitrary_band": { ... }
    }
  }
}
```

## Testing and Validation

### Test Suite
Run comprehensive validation:
```bash
cd /path/to/wave-crispr-signal
python tests/test_issue50_experiment.py
```

### Test Coverage
- **Initialization**: Framework setup and configuration
- **Data generation**: Synthetic dataset creation and validation
- **Feature extraction**: Z-Framework and baseline feature computation
- **Hypothesis testing**: H1 statistical validation framework
- **Reproducibility**: Deterministic behavior with fixed seeds
- **Results storage**: JSON/CSV output validation
- **Integration**: Complete workflow testing

## Performance Considerations

### Computational Complexity
- **Small scale (N=100)**: ~30 seconds for H1 test
- **Medium scale (N=1000)**: ~5 minutes for H1 test  
- **Large scale (N=10000)**: ~30 minutes for H1 test
- **Full experiment (all hypotheses)**: ~2-3 hours estimated

### Memory Requirements
- **Feature matrices**: O(N × F) where F ≈ 30-50 features
- **Bootstrap storage**: O(B × K) where B=1000, K=CV folds
- **Peak memory**: ~500MB for N=10,000 sequences

### Optimization Strategies
- **Batch processing**: Feature extraction in chunks for large datasets
- **Parallel CV**: Multiple fold evaluation (future enhancement)
- **Feature caching**: Store computed features for repeated experiments
- **Early stopping**: Skip expensive tests if preliminary results are clear

## Extensions and Future Work

### Additional Hypotheses
The framework is designed to be extensible for testing additional hypotheses:
- **H7**: Multi-scale consistency across dataset sizes
- **H8**: Guide length sensitivity analysis
- **H9**: PAM sequence context interactions
- **H10**: Cross-species generalization

### Enhanced Statistical Methods
- **Bayesian hypothesis testing**: Posterior probability assessment
- **Multiple comparison correction**: Family-wise error rate control
- **Effect size estimation**: Cohen's d, Glass's delta for practical significance
- **Power analysis**: Sample size determination for future experiments

### Integration with External Tools
- **CRISPRAnalyzeR**: Integration with established CRISPR analysis pipelines
- **DeepCRISPR**: Comparison with deep learning approaches  
- **Benchmarking suite**: Standardized evaluation against published methods

## Troubleshooting

### Common Issues

1. **Import Errors**: Framework includes fallback implementations for missing dependencies
2. **Memory Issues**: Reduce batch sizes or use feature caching for large datasets
3. **Convergence Problems**: Check feature scaling and regularization parameters
4. **Statistical Power**: Ensure sufficient sample sizes for reliable hypothesis testing

### Debug Mode
Enable verbose logging:
```python
experiment = Issue50ExperimentFramework(random_seed=42)
experiment.logger.setLevel(logging.DEBUG)
```

### Contact and Support
- **Repository**: https://github.com/zfifteen/wave-crispr-signal
- **Issues**: Create GitHub issue with experiment configuration and error logs
- **Documentation**: See `docs/` directory for detailed methodology

---

This implementation provides a complete, scientifically rigorous framework for validating the Z-Framework claims in Issue #50, with emphasis on reproducibility, statistical validity, and extensibility for future research directions.