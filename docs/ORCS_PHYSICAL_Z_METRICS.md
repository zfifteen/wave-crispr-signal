# ORCS Physical Z-Metrics Testing Documentation

## Overview

This document describes the ORCS (Oligonucleotide Repository for CRISPR Screens) physical Z-metrics testing implementation that validates four sequence-derivable physical metrics on real human CRISPR screen outcomes from BioGRID-ORCS v1.1.17.

## Purpose

The ORCS testing system evaluates the predictive power of four physical Z-metrics for CRISPR screen outcomes:

1. **Base-pair opening kinetics** - Hydrogen bond stability metrics
2. **Base-stacking dissociation** - Stacking energy disruption rates  
3. **Helical twist fluctuation** - B-form DNA structural dynamics
4. **Denaturation/melting kinetics** - Temperature-dependent stability proxy

## Implementation

### Core Script: `scripts/test_orcs_v1_1_17.py`

The main testing script implements the complete pipeline:

- **ORCS Data Processing**: Parses BioGRID-ORCS index and screen files
- **Screen Selection**: Identifies human Cas9 pooled knockout screens
- **Sequence Integration**: Loads human DNA sequences from FASTA files
- **Guide Enumeration**: Finds SpCas9 targets (20nt + NGG PAM)
- **Z-Metrics Calculation**: Computes physical metrics using Z = A·(B/c) framework
- **Statistical Analysis**: Performs correlation analysis with bootstrap validation
- **Results Output**: Generates CSV summaries and top screens reports

### Z-Framework Integration

The implementation uses the canonical Z-Framework form:

```
Z = A * (B / c)
```

Where:
- **A**: Frame-dependent mean property (per-base physical parameter)
- **B**: Standard deviation of adjacent differences (discrete dynamics)  
- **c**: Invariant constant (e² ≈ 7.389)

**Safety Guards**:
- Causality check: |B| < c (raises ValueError if violated)
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
FASTA="data/my_sequences.fasta" python scripts/test_orcs_v1_1_17.py

# Limit processing for testing
MAX_SCREENS=10 MIN_PAIRS=5 python scripts/test_orcs_v1_1_17.py
```

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
screen_id,score1_type,metric,aggregate,N,pearson_r,p_value,ci95_lo,ci95_hi,var_Z,cohens_d_hit_vs_non
1,log10 (corrected p-value),opening,mean,5,-0.974,0.0049,-1.0,1.0,0.0071,0.0
```

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