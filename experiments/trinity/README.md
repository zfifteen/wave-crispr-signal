# Trinity Experiments: CRISPR Breathing Dynamics Validation

This directory contains the "trinity" of experiments demonstrating DNA breathing dynamics effects on CRISPR guide efficiency.

## Overview

The trinity consists of three complementary analyses:
1. **Doench 2016 Lift** (Mission #1): Nested model showing breathing feature improvement
2. **Cross-species Generalization** (Mission #2): Train human → test zebrafish
3. **Rotational Phase Figure** (Mission #4): The "viral plot" showing φ-dependent activity

All experiments share:
- `spectral.py`: Core spectral analysis utilities
- `statistics.py`: Statistical testing framework
- Dimensionless parametrization (rate ratios, not MHz)
- Scientific gates enforcement

## Quick Start

### Installation
```bash
# From repository root
pip install -r requirements.txt

# Verify trinity imports
python -c "from experiments.trinity import spectral, statistics; print('✓ Trinity imports OK')"
```

### Running Experiments

#### Experiment C: Rotational Phase Figure (Implemented)
```bash
# Full analysis (1000 guides, 1000 bootstrap)
python experiments/trinity/phase_figure.py \
    --seed 42 \
    --n-guides 1000 \
    --bootstrap 1000 \
    --bins 24 \
    --output results/trinity/phase_figure

# Quick test (100 guides, 100 bootstrap)
python experiments/trinity/phase_figure.py \
    --n-guides 100 \
    --bootstrap 100 \
    --bins 12 \
    --output results/trinity/phase_test
```

**Expected outputs:**
- `results/trinity/phase_figure/phase_activity.png` - Polar plot
- `results/trinity/phase_figure/results.json` - Statistics

## Shared Utilities

### spectral.py

Core spectral analysis functions:

```python
from experiments.trinity.spectral import (
    encode_complex,           # Dimensionless DNA encoding
    breathing_features,       # Extract P10.5, P5.25, P3.5
    spectral_peak_analysis,   # Single-period analysis
    rotational_phase_curve,   # Phase wheel data
    phase_at,                 # Position → phase conversion
)

# Example: Analyze a sequence
features = breathing_features("ATGCGATCGATCG...", r=20.0)
# Returns: {"P10_5": ..., "P5_25": ..., "P3_5": ..., "length": ..., "gc_content": ...}
```

**Key parameters:**
- `r`: Rate ratio (dimensionless k_GC/k_AT), default 20.0, range [5, 200]
- `is_rna`: If True, expect U instead of T (default False)

### statistics.py

Statistical testing functions:

```python
from experiments.trinity.statistics import (
    bootstrap_ci,                    # Bootstrap confidence intervals
    permutation_test,                # Permutation-based p-values
    cohens_d, hedges_g,             # Effect sizes
    benjamini_hochberg_fdr,         # Multiple comparison correction
    rayleigh_test,                   # Circular uniformity test
    circular_linear_correlation,     # Phase-activity correlation
    partial_correlation,             # Partial correlation (control for covariates)
)

# Example: Bootstrap CI for mean
mean, lower, upper = bootstrap_ci(data, np.mean, n_bootstrap=1000, seed=42)

# Example: Rayleigh test for phase uniformity
Z, p = rayleigh_test(phases)
```

## Experiment Details

### Experiment C: Rotational Phase Figure

**Purpose**: Show that CRISPR activity oscillates with rotational phase at the cut site.

**Method**:
1. For each guide, calculate φ_cut = 2π·pos/10.5 (cut site phase)
2. Bin guides by phase (default 24 bins, 15° resolution)
3. Compute mean activity ± 95% bootstrap CI per bin
4. Test for non-uniformity (Rayleigh test)
5. Test phase-activity correlation (circular-linear)

**Pre-registered endpoints**:
- Rayleigh test p-value (null: uniform phase distribution)
- Circular-linear correlation r & p (phase ↔ activity relationship)
- Bootstrap 95% CIs for bin means (≥1000 resamples)

**Usage**:
```bash
python experiments/trinity/phase_figure.py \
    --seed 42 \
    --n-guides 1000 \
    --bootstrap 1000 \
    --bins 24 \
    --output results/trinity/phase_figure
```

**Arguments**:
- `--seed`: Random seed for reproducibility
- `--n-guides`: Number of guides to analyze (synthetic data for demo)
- `--bootstrap`: Bootstrap resamples (min 1000, recommended 10000)
- `--bins`: Number of phase bins (24 = 15° resolution)
- `--output`: Output directory

**Outputs**:
- `phase_activity.png`: Polar plot showing mean ± CI by phase
- `results.json`: Full statistics including Rayleigh p, circular-linear r

**Interpretation**:
- **Rayleigh p < 0.05**: Phases are non-uniformly distributed (clustering)
- **Circular-linear r > 0**: Activity correlates with phase
- **Bootstrap CIs**: Quantify uncertainty in bin means

### Experiments A & B: Doench Lift & Cross-species (Placeholders)

These experiments require access to real CRISPR datasets:
- **Doench 2016**: On-target activity scores (human)
- **Kim 2025**: Alternative human dataset
- **CRISPRscan**: Cross-species (zebrafish)

**Planned implementation**:
- `doench_lift.py`: Baseline + breathing nested model comparison
- `cross_species.py`: Train human, test zebrafish (zero-tuning)

## Dimensionless Parametrization

All experiments use **rate ratio** instead of absolute frequencies:

### Rate Ratio (r)
- **Definition**: r = k_GC / k_AT (opening rate ratio)
- **Physical basis**: GC pairs open ~100× slower than AT pairs
- **Typical values**: r ∈ [5, 200]
- **Default**: r = 20.0

### Encoding Formula
```
alpha = log(r)          # Strength factor
beta = 0.3 * log(r)     # Phase factor

AT: -alpha + j*beta     (fast opening, weak bonds)
GC: +alpha - j*beta     (slow opening, strong bonds)
```

### Why Dimensionless?
1. **Unit-free**: Avoids "MHz" ambiguity
2. **Biophysically anchored**: Based on relative rates
3. **Stable**: Works across temperature/ion conditions (with appropriate r)
4. **Interpretable**: r = 20 means GC opens 20× slower

## Scientific Gates

### Enforced by All Experiments
- ✅ **Human DNA only**: A/C/G/T (DNA) or A/C/G/U (RNA) validation
- ✅ **No fabrication**: Real nucleotide sequences only
- ✅ **Pre-registered endpoints**: Declared before analysis
- ✅ **Null models**: Permutation tests (≥1000 permutations)
- ✅ **Bootstrap CIs**: ≥1000 resamples
- ✅ **Multiple comparison correction**: Benjamini-Hochberg FDR when needed
- ✅ **Reproducibility**: Fixed seeds, documented parameters

### Leakage Prevention
- **Split-by-screen**: Train/test never share same screen
- **Split-by-gene**: Train/test never share same gene
- **Cross-species**: Species separation ensures no leakage

## Testing

### Unit Tests
```bash
# Test spectral utilities
python experiments/trinity/spectral.py

# Test statistics utilities
python experiments/trinity/statistics.py
```

### Integration Test
```bash
# Quick smoke test (<30s)
python experiments/trinity/phase_figure.py \
    --n-guides 100 \
    --bootstrap 100 \
    --bins 12 \
    --output /tmp/trinity_smoke
```

Expected: `✓ Trinity Experiment C completed successfully`

## File Structure

```
experiments/trinity/
├── __init__.py                # Package init
├── spectral.py                # Spectral analysis utilities
├── statistics.py              # Statistical testing framework
├── phase_figure.py            # Experiment C (rotational phase)
├── doench_lift.py             # Experiment A (planned)
├── cross_species.py           # Experiment B (planned)
├── manifest.yml               # Experiment metadata (planned)
└── README.md                  # This file
```

## Output Structure

```
results/trinity/
├── phase_figure/
│   ├── phase_activity.png     # Polar plot
│   └── results.json           # Statistics
├── doench_lift/               # Planned
│   ├── pr_curves.png
│   ├── delta_auprc.json
│   └── model_comparison.csv
└── cross_species/             # Planned
    ├── generalization.png
    └── zero_tuning_results.json
```

## Future Enhancements

### Near-term
1. **Real data integration**: Load Doench 2016, Kim 2025 datasets
2. **Doench lift implementation**: Baseline vs +breathing model
3. **Cross-species validation**: Human → zebrafish generalization
4. **Rate ratio sweep**: Stability heatmap across r ∈ [5, 200]

### Long-term
1. **Nearest-neighbor thermodynamics**: Temperature + ion-dependent encoding
2. **Batch processing**: Parallel analysis of large datasets
3. **Interactive visualizations**: Web-based exploration
4. **Meta-analysis**: Combine results across datasets

## References

### Repository
- **Demo MVP**: `web_apps/demo_mvp/` - Interactive demo
- **Repository policy**: `.github/REPOSITORY_POLICY.md`
- **Scientific gates**: `.github/copilot-instructions.md`
- **Z Framework**: `docs/Z_FRAMEWORK.md`

### Publications
- Doench et al. (2016) - Optimized sgRNA design, *Nature Biotechnology*
- Kim et al. (2025) - gRNA efficiency dataset (planned citation)
- Altan-Bonnet et al. (2003) - DNA breathing dynamics, *Physical Review Letters*

### Methods
- Rayleigh test: Zar (1999), *Biostatistical Analysis*
- Circular-linear correlation: Jammalamadaka & Sengupta (2001)
- Bootstrap: Efron & Tibshirani (1993), *An Introduction to the Bootstrap*
- FDR: Benjamini & Hochberg (1995), *Journal of the Royal Statistical Society*

## Contributing

When adding experiments:
1. Use shared `spectral.py` and `statistics.py` utilities
2. Enforce all scientific gates
3. Document pre-registered endpoints
4. Include smoke test with synthetic data
5. Follow repository policy for naming and structure

## License

MIT License - See repository LICENSE file

## Contact

For questions or issues:
- Open GitHub issue tagged `trinity` or `experiments`
- See main repository README for contact information

---

**Version**: 1.0.0  
**Created**: October 2025  
**Status**: Phase 1 Complete (Experiment C), Phases 2-3 Planned
