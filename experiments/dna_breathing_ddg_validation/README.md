# DNA Breathing Dynamics ΔΔG Validation Experiment

## Overview

This experiment validates the DNA breathing dynamics encoding hypothesis through thermodynamic perturbation analysis (ΔΔG-based mutation testing) following a pre-registered statistical framework.

## Hypothesis

**Primary (H1)**: DNA breathing dynamics encoding shows detectable sensitivity (|Cohen's d| ≥ 0.5) to thermodynamic perturbations (ΔΔG) from single-point mutations in high ΔΔG bin, particularly at the helical band (10.5 bp/turn).

**Secondary (H2)**: Effect size increases monotonically with |ΔΔG| magnitude (low < mid < high bins).

**Falsification (H3)**: Off-band control windows (±20% from f₀) show negligible effects (|d| < 0.2).

## Pre-Registered Parameters

- **Centers**: 10.3, 10.4, 10.5, 10.6, 10.7 bp/turn
- **Widths**: 1%, 2%, 3%, 5%, 6% of center frequency
- **ΔΔG Bins**: Tertiles (low/mid/high) based on absolute |ΔΔG|
- **Features**: Peak magnitude, band power, phase coherence, SNR
- **Effect Size**: Cohen's d with 95% BCa bootstrap CI (1,000 reps minimum)
- **Multiple Comparison Correction**: Benjamini-Hochberg FDR
- **Trend Test**: Jonckheere-Terpstra for dose-response
- **Controls**: Off-band windows, label-shuffle permutations (1,000 reps minimum)

## Data Requirements

- **Source**: Real human CRISPR sequences (no synthetic)
- **Alphabet**: A/C/G/T/N only (DNA), validated for human genome
- **Length**: Uniform guide length (typically 20 nt)
- **Quality**: No ambiguous bases beyond N

## Execution

### Quick Smoke Test

```bash
python ddg_validation.py \
    --input ../../data/doench2016.csv \
    --output results/smoke_test \
    --n-pairs 100 \
    --seed 42 \
    --smoke
```

### Full Validation

```bash
python ddg_validation.py \
    --input ../../data/doench2016.csv \
    --output results/full_validation \
    --n-pairs 6000 \
    --bootstrap 1000 \
    --permutation 1000 \
    --seed 42
```

## Outputs

All artifacts are saved to the specified output directory:

- `config.json` - Full configuration and parameters
- `pairs.csv` - WT-mutant pairs with ΔΔG and bin assignments
- `stats.csv` - Statistical results per (bin, center, width, feature)
- `trend.csv` - Bin-level summaries and trend tests
- `env.txt` - Environment snapshot (pip freeze)
- `FINDINGS.md` - Comprehensive results report

## Acceptance Criteria

### Primary Pass
- At least one (center, width, feature) in ΔΔG-high bin achieves:
  - |Cohen's d| ≥ 0.5
  - FDR-corrected q < 0.05
  - 95% CI excludes 0

### Robustness
- Effect replicates across ≥2 adjacent centers or widths
- Trend test significant across bins (p < 0.05)

### Specificity
- Off-band controls: |d| < 0.2, q ≥ 0.1
- Shuffle controls: Non-significant (q ≥ 0.1)

## Scientific Gates (Mandatory)

Per repository policy (`.github/REPOSITORY_POLICY.md`):

1. **Human DNA Only**: Validates human nucleotide sequences (A/C/G/T/N)
2. **No Fabrication**: Uses real CRISPR data only
3. **Z Invariants**: Applies discrete Z = A(B/e²) with proper guards
4. **Geometric Resolution**: Uses θ′(n,k) with k ≈ 0.3
5. **Statistical Validity**: Pre-registered endpoints, bootstrap CIs, FDR correction
6. **Leakage Control**: No data splitting required (paired design)
7. **Reproducibility**: Fixed seeds, persisted metadata, exact environment

## References

- PR #53: https://github.com/zfifteen/dna-breathing-dynamics-encoding/pull/53
- SantaLucia (1998): Nearest-neighbor thermodynamics
- Helical Period: ~10.5 bp for B-DNA
- BioGRID-ORCS: Dataset v1.1.17 (Homo sapiens)
