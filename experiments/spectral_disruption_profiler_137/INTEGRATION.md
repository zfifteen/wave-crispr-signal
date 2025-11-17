# Integration Guide: Spectral Disruption Profiler Falsification Experiment

## Overview

This experiment integrates with the `wave-crispr-signal` repository to provide rigorous falsification testing of spectral disruption hypotheses using Z Framework principles.

## Repository Integration Points

### 1. Z Framework Core (`scripts/z_framework.py`)
- **Used by**: `SpectralDisruptionFalsifier.__init__()`
- **Purpose**: High-precision Z Framework calculations (mp.dps=50)
- **Integration**: Direct import and instantiation of `ZFrameworkCalculator`

### 2. FFT CRISPR Disruption (`applications/fft_crispr_disruption.py`)
- **Used by**: Optional integration for advanced FFT features
- **Purpose**: Golden-ratio phase weighting and off-target detection
- **Integration**: Fallback to local implementation if not available

### 3. Topological Analysis (`scripts/topological_analysis.py`)
- **Used by**: Geometric resolution function θ′(n,k)
- **Purpose**: Geodesic curvature analysis with k* ≈ 0.3
- **Integration**: Mathematical constants and formulas

### 4. Synthetic Data Generation (`proof_pack/generate_synthetic_data.py`)
- **Used by**: Pattern reference for sequence generation
- **Purpose**: Generate adversarial test sequences with controlled GC content
- **Integration**: Similar methodology adapted for CRISPR guides

### 5. Validation Framework (`proof_pack/run_validation.py`)
- **Used by**: Statistical testing patterns
- **Purpose**: Bootstrap CI, permutation tests, effect size calculations
- **Integration**: Similar statistical methodology

## Experiment Structure

```
experiments/spectral_disruption_profiler_137/
├── __init__.py                          # Package initialization
├── README.md                            # Technical design document
├── manifest.yml                         # Experiment metadata
├── spectral_disruption_profiler.py      # Main implementation
├── smoke_test.py                        # CI validation (<5s)
└── golden_outputs/                      # Expected results
    ├── expected_structure.json
    └── validation_checksums.txt
```

## CI Integration

### Smoke Test (Fast)
```bash
# Run in CI pipeline for quick validation
python experiments/spectral_disruption_profiler_137/smoke_test.py
# Expected: 7/7 tests pass in <5 seconds
```

### Test Suite Integration
Add to `scripts/run_tests.py`:
```python
{
    "name": "Spectral Disruption Profiler",
    "path": "experiments/spectral_disruption_profiler_137/smoke_test.py"
}
```

## Usage Examples

### Quick Validation (CI/Development)
```bash
# Smoke test only
python experiments/spectral_disruption_profiler_137/smoke_test.py
```

### ⭐ Standard Experiment with Real Human DNA (RECOMMENDED)
```bash
# Use Doench 2016 dataset (real human gRNA sequences with efficiency labels)
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --input data/doench2016.csv \
  --seed 42 \
  --bootstrap 1000 \
  --permutation 1000 \
  --output results/spectral_disruption_profiler_137/doench2016_run.json
```

### Quick Test with Real Data (100 resamples)
```bash
# Faster testing with real data
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --input data/doench2016.csv \
  --seed 42 \
  --bootstrap 100 \
  --permutation 100 \
  --max-sequences 50 \
  --output results/spectral_disruption_profiler_137/quick_test.json
```

### High-Confidence Experiment (Publication)
```bash
# Increased resamples for publication-quality results with real data
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --input data/doench2016.csv \
  --seed 42 \
  --bootstrap 10000 \
  --permutation 10000 \
  --output results/spectral_disruption_profiler_137/publication_run.json
```

### ⚠️ Synthetic Data Mode (NOT RECOMMENDED)
```bash
# For testing only - does NOT provide valid falsification!
python experiments/spectral_disruption_profiler_137/spectral_disruption_profiler.py \
  --use-synthetic \
  --seed 42 \
  --n-sequences 100 \
  --bootstrap 100 \
  --permutation 100 \
  --output results/spectral_disruption_profiler_137/synthetic_test.json
```

**⚠️ WARNING**: Synthetic mode uses random DNA with random labels and does NOT test against real biological signal.

## Dependencies

### Required (from requirements.txt)
- `mpmath==1.3.0` - High-precision mathematics
- `numpy==1.26.4` - Numerical computing
- `scipy==1.16.1` - Scientific computing (FFT, stats)
- `scikit-learn==1.5.1` - ROC-AUC metrics
- `pandas==2.3.1` - CSV data loading (now required)

### Optional
- `matplotlib==3.10.5` - Visualization (for plotting results)

## Scientific Gates Compliance

### Human DNA Only
- ✅ Validates A/C/G/T/N only (case-insensitive)
- ✅ Rejects U (RNA) and IUPAC codes
- ✅ Fail-fast with clear error messages

### Z Framework Invariants
- ✅ Z = A(B/e²) with documented parameters
- ✅ High precision (mp.dps=50)
- ✅ Geometric resolution θ′(n,k) = φ·((n mod φ)/φ)^k

### Statistical Rigor
- ✅ Bootstrap CI (≥1,000 resamples, 95% level)
- ✅ Permutation tests (≥1,000 resamples, two-tailed)
- ✅ Effect size calculations (Cohen's d)
- ✅ Partial correlation (controlling for GC%, length, position)
- ✅ FDR correction (Benjamini-Hochberg)

### Reproducibility
- ✅ Fixed random seed
- ✅ Git commit tracking
- ✅ Metadata persistence (seed, timestamp, environment)
- ✅ Exact dependency versions

### Performance
- ✅ Smoke test: <5s (actual: ~1.3s)
- ✅ Functional test: <30s for 100 sequences
- ✅ Full experiment: <5 minutes for 1000 sequences

## Results Interpretation

### Falsification Status Values

1. **HYPOTHESIS_FALSIFIED** (Expected)
   - Bootstrap CI includes zero
   - Permutation p-value > α (0.05)
   - Conclusion: Phase weighting provides NO significant improvement
   - Interpretation: Null hypothesis cannot be rejected

2. **HYPOTHESIS_SUPPORTED** (Unexpected)
   - Bootstrap CI excludes zero
   - Permutation p-value < α (0.05)
   - Conclusion: Phase weighting provides significant improvement
   - Interpretation: Alternative hypothesis supported

3. **INCONCLUSIVE**
   - Mixed results from different tests
   - Requires further investigation
   - May indicate edge cases or underpowered study

## Troubleshooting

### Import Errors
```bash
# If Z Framework import fails
export PYTHONPATH=/home/runner/work/wave-crispr-signal/wave-crispr-signal:$PYTHONPATH
```

### Memory Issues (Large Experiments)
```bash
# Reduce bootstrap/permutation samples
python spectral_disruption_profiler.py --bootstrap 500 --permutation 500
```

### Slow Performance
```bash
# Use smaller sequence count for testing
python spectral_disruption_profiler.py --n-sequences 20
```

## Future Extensions

### Planned Features
1. Integration with real Kim 2025 dataset
2. RuleSet3 comparison on benchmark data
3. Multi-guide screening analysis
4. k-parameter grid search
5. Conical flow density corrections

### Research Directions
1. Test on high-noise datasets
2. Explore k-parameter sensitivity
3. Compare with alternative phase weighting schemes
4. Integrate QMC (Quasi-Monte Carlo) sampling
5. Add GVA (Geodesic Validation Assault) embedding tests

## References

- **Z Framework**: `scripts/z_framework.py`, `docs/Z_FRAMEWORK.md`
- **Phase Weighting**: `docs/PHASE_WEIGHTED_SCORECARD.md`
- **Repository Policy**: `.github/REPOSITORY_POLICY.md`
- **Main README**: `README.md`

## Contact & Support

For questions about this experiment:
1. Review the technical design: `experiments/spectral_disruption_profiler_137/README.md`
2. Check the manifest: `experiments/spectral_disruption_profiler_137/manifest.yml`
3. Run smoke test: `python experiments/spectral_disruption_profiler_137/smoke_test.py`
4. Check golden outputs: `experiments/spectral_disruption_profiler_137/golden_outputs/`

---

**Last Updated**: 2025-11-17  
**Version**: 1.0  
**Experiment ID**: spectral_disruption_profiler_137
