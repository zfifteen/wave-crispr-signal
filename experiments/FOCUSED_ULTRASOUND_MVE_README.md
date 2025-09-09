# Focused Ultrasound MVE: Z Framework Spatial Targeting Precision

## Overview

This experiment tests the hypothesis that the Z Framework improves targeting precision in simulated focused ultrasound (FUS) by reducing spatial error through nonlinear time-distance transforms and geodesic curvature modeling.

## Hypothesis

**Primary Hypothesis**: The Z Framework improves targeting precision in simulated focused ultrasound (FUS) by reducing spatial error through nonlinear time-distance transforms and geodesic curvature modeling.

**Null Hypothesis**: There is no significant difference in spatial targeting error between the Z Framework and baseline acoustic models.

## Experimental Design

### Simulation Parameters
- **Grid Size**: 100×100 2D tissue grid
- **Acoustic Velocity**: 1540 m/s base with ±10% Gaussian heterogeneity
- **Trials**: 1,000 targeting trials per run
- **Source Position**: Fixed at (5, 5)
- **Target Position**: Random within grid bounds (10-90, 10-90)

### Models Compared

#### Baseline Model
- Euclidean distance calculation
- Constant acoustic velocity assumption
- No heterogeneity compensation
- Error proportional to distance × velocity variance

#### Z Framework Model  
- Discrete domain form: Z = A(B/e²)
- Geodesic resolution: θ'(n,k) = φ·((n mod φ)/φ)^k
- Path-integrated velocity sampling
- Nonlinear error correction through curvature modeling

### Statistical Validation

#### Pre-registered Endpoints
1. **Pearson correlation** with 95% bootstrap CI (≥1,000 resamples)
2. **Permutation test** for difference in means (≥1,000 permutations)
3. **Effect size** (Cohen's d) with confidence intervals
4. **Improvement percentage** with bootstrap CI

#### Null Model
- ≥1,000× permutation testing for empirical p-values
- Two-tailed significance testing (α = 0.05)

#### Reproducibility Controls
- Fixed random seed for reproducibility
- Complete metadata persistence
- Environment state capture

## Usage

### Basic Execution
```bash
# Run with default parameters
cd /home/runner/work/wave-crispr-signal/wave-crispr-signal
python experiments/focused_ultrasound_mve.py --seed 42 --bootstrap 1000 --permutation 1000 --splits single --domain discrete
```

### Advanced Configuration
```bash
# Custom parameters
python experiments/focused_ultrasound_mve.py \
    --seed 42 \
    --bootstrap 2000 \
    --permutation 5000 \
    --splits single \
    --domain discrete \
    --k-parameter 0.3 \
    --grid-size 100 \
    --n-trials 1000 \
    --visualize \
    --output-dir results
```

### Required CLI Flags (Scientific Gates)
- `--seed`: Random seed for reproducibility
- `--bootstrap`: Bootstrap samples (≥1,000)
- `--permutation`: Permutation samples (≥1,000)
- `--splits`: Split strategy
- `--domain`: Z Framework domain form
- `--k-parameter`: Geodesic resolution parameter (default: 0.3)

## Expected Runtime

- **Target**: <5 minutes on consumer hardware
- **Typical**: 2-3 minutes for 1,000 trials
- **Scaling**: O(n) with number of trials

## Output Files

### Results Directory Structure
```
results/focused_ultrasound_mve/run-YYYYMMDD-HHMMSS/
├── results.csv          # Raw trial data
├── analysis.json        # Statistical analysis results
├── metadata.json        # Experiment configuration and environment
├── experiment.log       # Human-readable summary
└── visualization.png    # Plots (if --visualize enabled)
```

### Key Metrics in results.csv
- `trial_id`: Sequential trial identifier
- `target_x`, `target_y`: Target coordinates
- `baseline_error`: Baseline model targeting error
- `z_framework_error`: Z Framework model targeting error
- `baseline_time`, `z_framework_time`: Time-to-target calculations
- `grid_variance`: Local velocity field variance

### Statistical Results in analysis.json
- `improvement_percentage`: Mean improvement with CI
- `cohens_d`: Effect size
- `correlation_r`: Error correlation with CI and p-value
- `permutation_p_value`: Null hypothesis p-value
- `statistically_significant`: Boolean result (p < 0.05)

## Interpretation Guidelines

### Success Criteria
1. **Statistical Significance**: p < 0.05 in permutation test
2. **Effect Size**: |Cohen's d| > 0.2 (small to medium effect)
3. **Consistent Improvement**: Positive improvement percentage
4. **Confidence**: Non-overlapping bootstrap CIs

### Expected Results
Based on Z Framework theory:
- **Improvement**: 5-15% reduction in targeting error
- **Effect Size**: Small to medium (d = 0.2-0.5)
- **Significance**: p < 0.05 if hypothesis is correct
- **Correlation**: Moderate negative correlation between baseline and Z Framework errors

## Scientific Rigor

### Validation Controls
✓ Bootstrap confidence intervals (≥1,000 samples)
✓ Permutation null testing (≥1,000 samples)  
✓ Effect size calculation with CI
✓ Multiple comparison awareness
✓ Reproducible random seeding

### Leakage Prevention
✓ No data leakage between models
✓ Independent error calculations
✓ Controlled random generation
✓ Fixed experimental parameters

### Documentation Standards
✓ Complete methodology documentation
✓ Pre-registered statistical endpoints
✓ Metadata persistence with environment state
✓ Version control integration

## Dependencies

Requires packages from root `requirements.txt`:
- `mpmath==1.3.0` (high-precision arithmetic)
- `numpy==1.26.4` (numerical computing)
- `scipy==1.16.1` (statistical analysis)  
- `matplotlib==3.10.5` (visualization)

## Troubleshooting

### Common Issues
1. **Import Errors**: Ensure all dependencies installed from `requirements.txt`
2. **Memory Issues**: Reduce `--n-trials` for large grids
3. **Runtime Exceeds 5min**: Reduce `--bootstrap` or `--permutation` samples
4. **No Significant Results**: Check `--k-parameter` sensitivity

### Validation Checks
- Grid size must be ≥10 for meaningful heterogeneity
- Bootstrap/permutation samples must be ≥1,000
- K parameter typically in range 0.1-0.9
- Random seed should be consistent across runs

## References

1. Z Framework Core Implementation: `z_framework.py`
2. Repository Policy: `.github/REPOSITORY_POLICY.md` 
3. Scientific Gates: `.github/copilot-instructions.md`
4. Statistical Methods: Standard bootstrap and permutation testing
5. Focused Ultrasound Physics: Simplified acoustic propagation modeling

---

**RESEARCH USE ONLY**: This experiment is for scientific hypothesis testing only, not clinical targeting applications.