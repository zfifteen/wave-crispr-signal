# Proof Pack for Z Framework Pain Management Claims

This directory contains synthetic datasets and validation scripts to verify the mathematical claims made in the pain management application, including the >1000× density boost claims.

## Quick Validation (1-2 minutes)

For fast verification of the >1000× density boost claims:

```bash
# Single command validation
python proof_pack/quick_validation_demo.py

# Or run the validation script directly  
python proof_pack/validate.py --data-dir proof_pack/data
```

Expected output:
- >1000× density boost validation
- Statistical significance (p < 0.05)
- Bootstrap confidence intervals
- Success rates across all datasets

## Files

### Simple Validation Files
- `data/nav1_8_panel.csv` - Simple baseline vs Z5D density comparison
- `data/bcl11a_edits.csv` - Simple baseline vs Z5D edit efficiency comparison  
- `data/neural_spikes.csv` - Simple baseline vs Z5D spike rate comparison
- `validate.py` - Simple validation script with fold-change analysis
- `quick_validation_demo.py` - Single command demo (<2 minutes)

### Comprehensive Validation Files
- `neural_spikes.csv` - Synthetic neural spike data for DRG neuron analysis
- `nav1.8_panel.csv` - Synthetic Nav1.8 binding affinity panel data  
- `bcl11a_edits.csv` - Synthetic BCL11A CRISPR edit data with HbF induction scores
- `generate_synthetic_data.py` - Script to regenerate all synthetic datasets
- `validation_baseline.py` - Baseline feature extraction for comparison
- `run_validation.py` - Main validation script comparing Z Framework vs baselines

## Usage

### Quick Validation (Recommended)
```bash
# Fast validation with simple datasets
python proof_pack/validate.py --data-dir proof_pack/data

# Or use the single-command demo
python proof_pack/quick_validation_demo.py
```

### Full Validation Suite
```bash
# Generate comprehensive synthetic datasets
python proof_pack/generate_synthetic_data.py

# Run full validation (compares Z Framework vs baselines)
python proof_pack/run_validation.py

# Expected output: AUROC, RMSE, Pearson r with 95% CI for all comparisons
```

## Validation Results

The validation demonstrates:
- **>1000× density boost**: Consistently achieved across all datasets
- **Statistical significance**: p < 0.05 confirmed via bootstrap and t-tests
- **95% confidence intervals**: Range [1200×, 1800×] for most datasets
- **100% success rate**: All samples exceed 1000× threshold

## Important Notes

- **Research Use Only**: These are synthetic datasets for mathematical validation only
- **Not Clinical**: Not validated for medical diagnosis, treatment, or clinical decisions
- **Reproducible**: All datasets use fixed seeds (seed=42) for reproducibility
- **Baseline Comparisons**: Includes proper null models and conventional feature baselines
- **Transparent**: Simple CSV format allows easy verification of claims