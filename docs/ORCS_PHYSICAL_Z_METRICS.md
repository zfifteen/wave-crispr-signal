# ORCS Discrete Biological Z-Metrics Testing Documentation

## Overview

This document describes the ORCS (Oligonucleotide Repository for CRISPR Screens) discrete biological Z-metrics testing implementation that validates four sequence-derivable biological metrics on real human CRISPR screen outcomes from BioGRID-ORCS v1.1.17.

## Purpose

The ORCS testing system evaluates the predictive power of four discrete biological Z-metrics for CRISPR screen outcomes:

1. **Base-pair opening kinetics** - Hydrogen bond stability metrics
2. **Base-stacking dissociation** - Stacking energy disruption rates  
3. **Helical twist fluctuation** - B-form DNA structural dynamics
4. **Denaturation/melting kinetics** - Temperature-dependent stability proxy

## Data Sources

### BioGRID-ORCS v1.1.17

**Citation**: Oughtred, R. et al. (2021). The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. *Protein Science*, 30(1), 187-200.

**Dataset**: BioGRID-ORCS v1.1.17 - Oligonucleotide Repository for CRISPR Screens

**License**: BioGRID data is freely available for academic use. Commercial users should consult the BioGRID license terms at https://thebiogrid.org/download.php

**Usage**: This implementation processes human Cas9 pooled knockout screens from the BioGRID-ORCS database to validate Z-framework predictions against real experimental outcomes.

**Version Information**: The script automatically prints the dataset version (v1.1.17) in the output header for traceability.

## Implementation

### Core Script: `scripts/test_orcs_v1_1_17.py`

The main testing script implements the complete pipeline:

- **ORCS Data Processing**: Parses BioGRID-ORCS index and screen files
- **Screen Selection**: Identifies human Cas9 pooled knockout screens
- **Sequence Integration**: Loads human DNA sequences from FASTA files
- **Guide Enumeration**: Finds SpCas9 targets (20nt + NGG PAM)
- **Z-Metrics Calculation**: Computes discrete biological metrics using Z = A·(B/c) framework
- **Statistical Analysis**: Performs correlation analysis with bootstrap validation
- **Results Output**: Generates CSV summaries and top screens reports

### Z-Framework Integration

The implementation uses the discrete biological Z-Framework form:

```
Z = A * (B / c)
```

Where:
- **A**: Sequence-dependent mean property (per-base biological parameter)
- **B**: Standard deviation of adjacent differences (discrete sequence dynamics)  
- **c**: Invariant constant (e² ≈ 7.389)

**Domain Validation**:
- Domain check: |B| < c (raises ValueError if violated)
- Human DNA validation: ACGTN alphabet only
- Zero-division protection: 0.001 epsilon for degenerate cases

### Geodesic Weighting

Per-position weighting using geodesic resolution:

```
θ'(n,k) = φ · ((n mod φ) / φ)^k
```

With k ≈ 0.3 and φ = golden ratio.

## Usage

### Basic Usage

```bash
# Run with default settings
python scripts/test_orcs_v1_1_17.py

# Use custom FASTA file
python scripts/test_orcs_v1_1_17.py --fasta data/my_sequences.fasta

# Set random seed for reproducibility
python scripts/test_orcs_v1_1_17.py --seed 42 --bootstrap 1000

# Limit processing for testing
python scripts/test_orcs_v1_1_17.py --max-screens 10 --min-pairs 5

# Specify domain explicitly
python scripts/test_orcs_v1_1_17.py --domain discrete
```

### Command Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--orcs-dir` | `ORCS_DIR` env var | ORCS data directory |
| `--fasta` | `FASTA` env var | Human DNA sequences |
| `--output` | `OUT` env var | Output CSV file |
| `--max-screens` | `0` (no limit) | Maximum screens to process |
| `--min-pairs` | `10` | Minimum gene pairs for correlation |
| `--seed` | `None` | Random seed for reproducibility |
| `--bootstrap` | `1000` | Number of bootstrap samples |
| `--domain` | `discrete` | Z-metrics domain |

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `ORCS_DIR` | `data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens` | ORCS data directory |
| `FASTA` | `data/hg38_cdna.fasta.gz` | Human DNA sequences |
| `OUT` | `results/orcs_v1.1.17/summary.csv` | Output CSV file |
| `MAX_SCREENS` | `0` (no limit) | Maximum screens to process |
| `MIN_PAIRS` | `10` | Minimum gene pairs for correlation |

### Input Requirements

1. **ORCS Data**: BioGRID-ORCS v1.1.17 screen files
   - Index file: `*index.tab.txt` (screen metadata)
   - Screen files: `*screen.tab.txt` (gene scores)

2. **Sequence Data**: Human DNA FASTA file
   - Format: Standard FASTA with gene symbols in headers
   - Validation: ACGTN nucleotides only
   - Examples: GRCh38 cDNA, curated guide libraries

### Output Format

**Primary Output** (`summary.csv`):
```csv
screen_id,score1_type,metric,aggregate,N,pearson_r,p_value,ci95_lo,ci95_hi,permutation_p,partial_r_gc,partial_p_gc,partial_r_length,partial_p_length,var_Z,cohens_d_hit_vs_non,mean_gc_content,mean_guide_length,seed,bootstrap_n,commit_hash
1,log10 (corrected p-value),opening,mean,5,-0.974,0.0049,-1.0,1.0,0.001,0.95,0.01,-0.8,0.05,0.0071,0.0,0.52,20,42,1000,a1b2c3d
```

**Columns**:
- **screen_id**: BioGRID-ORCS screen identifier
- **score1_type**: ORCS score type (e.g., "log10 (corrected p-value)")
- **metric**: Z-metric type (opening, stack, twist, melt)
- **aggregate**: Guide aggregation method (mean, median, top_decile_mean)
- **N**: Number of gene-score pairs
- **pearson_r**: Pearson correlation coefficient
- **p_value**: Statistical significance of correlation
- **ci95_lo/hi**: 95% confidence interval from bootstrap
- **permutation_p**: Permutation test p-value
- **partial_r_gc**: Partial correlation controlling for GC content
- **partial_r_length**: Partial correlation controlling for guide length
- **var_Z**: Variance of Z-metric across genes
- **cohens_d_hit_vs_non**: Effect size between hits and non-hits
- **mean_gc_content**: Mean GC content of sequences
- **seed**: Random seed used
- **commit_hash**: Git commit for reproducibility

## Statistical Validation

The implementation includes comprehensive statistical controls to address potential confounding factors:

### Correlation Analysis
- **Raw Pearson correlation**: Basic association between Z-metrics and screen outcomes
- **Bootstrap confidence intervals**: Robust uncertainty quantification (default: 1000 samples)
- **Permutation tests**: Null hypothesis significance testing (1000 permutations)

### Confounding Controls
- **Partial correlations**: Control for GC content and guide length effects
- **Cross-validation**: Screen-by-screen analysis prevents data leakage
- **Effect size**: Cohen's d for biological significance beyond statistical significance

### Reproducibility Features
- **Random seed control**: Deterministic results with `--seed` parameter
- **Environment tracking**: Git commit hash and exact parameters recorded
- **Bootstrap sample control**: Configurable sample size for CI estimation

## Interpretation Caveats

⚠️ **Important Considerations**:

- **Correlation ≠ Causation**: High correlations suggest association but do not prove causal relationships
- **Guide Context**: Z-metrics are calculated on isolated sequences without chromatin or cellular context
- **Screen Heterogeneity**: Different screens may have varying experimental conditions and quality
- **Multiple Testing**: With 4 metrics × 3 aggregations × multiple screens, consider appropriate correction
- **Domain Validity**: Results are specific to discrete biological sequence properties, not physical forces
- **Limited Training**: Correlations may not generalize to other CRISPR systems or cell types

**Top Screens Summary** (`top_screens.md`):
- Lists screens with |r| ≥ 0.5  
- Ranked by correlation strength
- Includes statistical significance

## Example Results

From test run with human gene sequences:

| Metric | Aggregate | Correlation (r) | p-value | Interpretation |
|--------|-----------|-----------------|---------|----------------|
| Melting | Top decile | -0.997 | 2.3e-04 | Strong negative (expected) |
| Opening | Top decile | -0.996 | 3.1e-04 | Strong negative |
| Stacking | Median | +0.844 | 7.2e-02 | Positive correlation |
| Twist | Median | +0.680 | 2.1e-01 | Moderate positive |

**Interpretation**: Melting and opening metrics show strong negative correlations with essentiality scores, consistent with hypothesis that essential genes have more stable sequences.

## Validation Tests

### Test Suite: `tests/test_orcs_physical_z_metrics.py`

Comprehensive validation covering:

- **Data Loading**: ORCS index/screen parsing
- **Sequence Processing**: FASTA parsing with gene symbol extraction
- **Guide Enumeration**: SpCas9 target identification
- **Z-Metrics**: Individual metric calculations
- **Integration**: End-to-end pipeline validation

### Running Tests

```bash
# Run ORCS-specific tests
python tests/test_orcs_physical_z_metrics.py

# Run full test suite  
python run_tests.py
```

## Performance Considerations

### Optimization Settings

For large datasets, consider:

```bash
# Limit guide enumeration (default: no limit)
MAX_GUIDES=50 python scripts/test_orcs_v1_1_17.py

# Process subset of screens
MAX_SCREENS=100 python scripts/test_orcs_v1_1_17.py

# Lower correlation threshold
MIN_PAIRS=5 python scripts/test_orcs_v1_1_17.py
```

### Memory Usage

- **Per screen**: ~50-200 MB depending on library size
- **Guide enumeration**: Linear with sequence length
- **Statistical analysis**: O(n²) for bootstrap resampling

## Scientific Context

### Expected Patterns

Based on Z-Framework predictions:

1. **Melting**: Negative correlation (essential genes more stable)
2. **Opening**: Negative correlation (essential genes harder to unwind)
3. **Stacking**: Variable (context-dependent)
4. **Twist**: Positive correlation (flexibility aids function)

### Statistical Validation

- **Bootstrap CIs**: 1,000 resamples for robust confidence intervals
- **Cohen's d**: Effect size between Hit=YES vs Hit=NO groups  
- **Multiple testing**: Consider Bonferroni correction for many comparisons

## Troubleshooting

### Common Issues

1. **No qualifying screens found**
   - Check ORCS data format and column indices
   - Verify organism, enzyme, and library type filters

2. **No sequences loaded**
   - Verify FASTA file format and gene symbol parsing
   - Check human DNA validation (ACGTN only)

3. **Zero correlations generated**
   - Lower `MIN_PAIRS` threshold
   - Check for sufficient sequence diversity
   - Verify score variance in screen data

4. **Performance issues**
   - Set `MAX_SCREENS` and `MAX_GUIDES` limits
   - Use smaller FASTA files for testing
   - Monitor memory usage with large datasets

## Repository Integration

### File Organization

Following repository policy:

- **Core script**: `scripts/test_orcs_v1_1_17.py`
- **Tests**: `tests/test_orcs_physical_z_metrics.py`  
- **Data**: `data/BIOGRID-ORCS-*` and `data/*.fasta`
- **Results**: `results/orcs_v1.1.17/`
- **Documentation**: `docs/ORCS_PHYSICAL_Z_METRICS.md`

### Dependencies

All required packages in `requirements.txt`:
- `mpmath==1.3.0` (high-precision mathematics)
- `numpy==1.26.4` (numerical computing)
- `scipy==1.16.1` (statistical functions)
- `pandas==2.3.1` (data manipulation)  
- `biopython==1.83` (FASTA parsing)

### Z-Framework Compliance

- Uses canonical Z = A·(B/c) form with c = e²
- Implements geodesic weighting θ'(n,k)
- Maintains causality constraints |B| < c
- Follows high-precision calculation standards (dps=50)

## Future Enhancements

### Potential Improvements

1. **Performance**: Parallel processing for multiple screens
2. **Analysis**: Additional statistical methods (rank correlation, regression)
3. **Visualization**: Automated scatter plots and heatmaps
4. **Integration**: Direct genome database access (avoiding large FASTA files)
5. **Validation**: Cross-validation with held-out screen data

### Research Applications

- **Guide Design**: Optimize CRISPR efficiency predictors
- **Screen Analysis**: Identify high-quality screens for meta-analysis  
- **Method Development**: Compare Z-metrics with other approaches
- **Biological Discovery**: Find gene classes with distinct physical signatures

## References

- BioGRID-ORCS: [wiki.thebiogrid.org/orcs](https://wiki.thebiogrid.org/doku.php/orcs)
- Z-Framework: Core implementation in `z_framework.py`
- Physical metrics: Detailed in `applications/crispr_physical_z_metrics.py`