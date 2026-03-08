# ORCS v1.1.17 Discrete Biological Z-Metrics Experiment

## Overview

This experiment validates four discrete biological Z-metrics against real human CRISPR screen outcomes using BioGRID-ORCS v1.1.17 data.

## Dataset Information

**Source**: BioGRID-ORCS v1.1.17 - Oligonucleotide Repository for CRISPR Screens  
**URL**: https://thebiogrid.org/download.php  
**Citation**: Oughtred, R. et al. (2021). The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. *Protein Science*, 30(1), 187-200.

**Expected Dataset Structure**:
```
data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens/
├── BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens.index.tab.txt
└── BIOGRID-ORCS-SCREEN_*-1.1.17.screen.tab.txt (multiple files)
```

**SHA256 Checksums** (to be verified):
- Index file: `[to be calculated after download]`
- Total screens: 1,744+ human Cas9 pooled knockout screens

## Environment Requirements

**Python Version**: 3.12+

**Required Packages** (exact versions):
```
numpy==2.3.2
pandas==2.3.2
scipy==1.16.1
biopython==1.85
mpmath==1.3.0
matplotlib==3.10.6
```

Install with:
```bash
pip install numpy==2.3.2 pandas==2.3.2 scipy==1.16.1 biopython==1.85 mpmath==1.3.0 matplotlib==3.10.6
```

## Reproducible Execution

### Step 1: Data Preparation
```bash
# Download BioGRID-ORCS v1.1.17 data
# Place in data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens/

# Prepare human reference sequences (GRCh38 cDNA recommended)
# Place in data/hg38_cdna.fasta.gz
```

### Step 2: Basic Validation Run
```bash
# Quick validation with limited data (5-10 minutes)
python scripts/test_orcs_v1_1_17.py \
    --max-screens 5 \
    --min-pairs 5 \
    --seed 42 \
    --bootstrap 100 \
    --output results/orcs_v1.1.17/validation_summary.csv
```

**Expected Results**:
- Runtime: 5-10 minutes
- Output: ~10-20 correlation results
- Strong correlations (|r| > 0.5): 2-5 results expected

### Step 3: Full Analysis Run
```bash
# Complete analysis with all data (2-6 hours)
python scripts/test_orcs_v1_1_17.py \
    --seed 42 \
    --bootstrap 1000 \
    --min-pairs 10 \
    --output results/orcs_v1.1.17/full_summary.csv
```

**Expected Results**:
- Runtime: 2-6 hours (depending on hardware)
- Output: 500-2000 correlation results
- Processing: ~1744 screens, ~15000+ genes
- Strong correlations (|r| > 0.5): 50-200 results expected

### Step 4: Verification
```bash
# Verify reproducibility with same seed
python scripts/test_orcs_v1_1_17.py \
    --seed 42 \
    --max-screens 3 \
    --output results/orcs_v1.1.17/verification.csv

# Results should be identical to first 3 screens from full run
```

## Expected Performance

### Wall-Clock Time Ranges

| Configuration | Screens | Expected Runtime | Memory Usage |
|---------------|---------|------------------|--------------|
| Validation    | 5       | 5-10 minutes     | <2 GB        |
| Medium        | 50      | 30-60 minutes    | <4 GB        |
| Full          | 1744+   | 2-6 hours        | <8 GB        |

### Hardware Requirements
- **CPU**: 2+ cores recommended
- **RAM**: 8 GB minimum, 16 GB recommended
- **Storage**: 10 GB free space for data and results

## Output Interpretation

### Key Metrics
- **pearson_r**: Primary correlation coefficient
- **permutation_p**: Null hypothesis test p-value
- **partial_r_gc**: Correlation controlling for GC content
- **partial_r_length**: Correlation controlling for sequence length

### Success Criteria
1. **Completion**: All screens process without errors
2. **Statistical Power**: N > 10 gene pairs per correlation
3. **Effect Validation**: |r| > 0.3 for at least 10% of results
4. **Reproducibility**: Identical results with same seed

## Troubleshooting

### Common Issues
1. **Memory Error**: Reduce --max-screens or use smaller FASTA
2. **Missing Data**: Verify ORCS data download and placement
3. **Import Errors**: Check package versions match requirements
4. **Slow Performance**: Use --max-screens to limit scope

### Validation Checks
```bash
# Test human DNA validation
python -c "from scripts.test_orcs_v1_1_17 import validate_seq; validate_seq('ATCG')"

# Check ORCS data presence
ls data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens/*.txt | wc -l
# Should show 1745+ files (1 index + 1744+ screens)
```

## Version Control

**Commit Hash**: Run `git rev-parse HEAD` to get current version  
**Data Version**: BioGRID-ORCS v1.1.17  
**Last Updated**: 2025-01-16

## Contact

For questions or issues with this experiment, refer to the main repository documentation or open an issue.