# Proof Pack for Z Framework Pain Management Claims

This directory contains synthetic datasets and validation scripts to verify the mathematical claims made in the pain management application.

## Files

- `neural_spikes.csv` - Synthetic neural spike data for DRG neuron analysis
- `nav1.8_panel.csv` - Synthetic Nav1.8 binding affinity panel data  
- `bcl11a_edits.csv` - Synthetic BCL11A CRISPR edit data with HbF induction scores
- `generate_synthetic_data.py` - Script to regenerate all synthetic datasets
- `validation_baseline.py` - Baseline feature extraction for comparison
- `run_validation.py` - Main validation script comparing Z Framework vs baselines

## Usage

```bash
# Generate synthetic datasets
python proof_pack/generate_synthetic_data.py

# Run full validation (compares Z Framework vs baselines)
python proof_pack/run_validation.py

# Expected output: AUROC, RMSE, Pearson r with 95% CI for all comparisons
```

## Important Notes

- **Research Use Only**: These are synthetic datasets for mathematical validation only
- **Not Clinical**: Not validated for medical diagnosis, treatment, or clinical decisions
- **Reproducible**: All datasets use fixed seeds (seed=42) for reproducibility
- **Baseline Comparisons**: Includes proper null models and conventional feature baselines