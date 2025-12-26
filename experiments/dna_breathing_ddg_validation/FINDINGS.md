# DNA Breathing Dynamics ΔΔG Validation — FINDINGS

**Experiment ID**: dna_breathing_ddg_validation  
**Date**: 2025-12-26 07:11:17  
**Git Commit**: eed07ff7610daf8bb680f31533837676fb61d380  
**Seed**: 42

---

## CONCLUSION

**HYPOTHESIS NOT SUPPORTED**: DNA breathing dynamics encoding does NOT show the predicted sensitivity to thermodynamic perturbations (ΔΔG) at the pre-registered threshold (|Cohen's d| ≥ 0.5).

### Key Findings

- **Pairs Analyzed**: 500 WT-mutant pairs
- **ΔΔG Range**: -1.680 to 1.830 kcal/mol
- **Maximum |Cohen's d|**: 0.287
- **Primary Hits (|d| ≥ 0.5)**: 0
- **Primary Criterion**: FAIL
- **Robustness Criterion**: FAIL
- **Specificity Criterion**: PASS

---

## STATISTICAL EVIDENCE

### Pre-Registered Parameters

- **Helical Band Centers**: 10.3 bp/turn, 10.4 bp/turn, 10.5 bp/turn, 10.6 bp/turn, 10.7 bp/turn
- **Band Widths**: 1%, 2%, 3%, 5%, 6%
- **Features**: peak_magnitude, band_power, phase_coherence, snr
- **Bootstrap Resamples**: 1,000
- **Permutation Resamples**: 1,000
- **ΔΔG Bins**: Tertiles (low/mid/high)

### ΔΔG Bin Distribution

| Bin | Count | Mean |ΔΔG| (kcal/mol) |
|-----|-------|----------------------|
| low | 174 | 0.202 |
| mid | 162 | 0.864 |
| high | 164 | 1.357 |

### Top 10 Results by |Cohen's d|

| Bin | Center (bp/turn) | Width (%) | Feature | |d| | q-value | 95% CI |
|-----|------------------|-----------|---------|-----|---------|--------|
| low | 10.3 | 1 | peak_magnitude | 0.036 | 0.742 | [-0.120, 0.185] |
| low | 10.3 | 2 | peak_magnitude | 0.036 | 0.742 | [-0.113, 0.180] |
| low | 10.3 | 3 | peak_magnitude | 0.036 | 0.742 | [-0.123, 0.188] |
| low | 10.3 | 5 | peak_magnitude | 0.036 | 0.742 | [-0.116, 0.180] |
| low | 10.3 | 6 | peak_magnitude | 0.036 | 0.742 | [-0.108, 0.185] |
| low | 10.4 | 1 | peak_magnitude | 0.036 | 0.742 | [-0.108, 0.183] |
| low | 10.4 | 2 | peak_magnitude | 0.036 | 0.742 | [-0.115, 0.191] |
| low | 10.4 | 3 | peak_magnitude | 0.036 | 0.742 | [-0.116, 0.184] |
| low | 10.4 | 5 | peak_magnitude | 0.036 | 0.742 | [-0.115, 0.184] |
| low | 10.4 | 6 | peak_magnitude | 0.036 | 0.742 | [-0.124, 0.195] |
| low | 10.5 | 1 | peak_magnitude | 0.036 | 0.742 | [-0.104, 0.184] |
| low | 10.5 | 2 | peak_magnitude | 0.036 | 0.742 | [-0.112, 0.179] |
| low | 10.5 | 3 | peak_magnitude | 0.036 | 0.742 | [-0.105, 0.192] |
| low | 10.5 | 5 | peak_magnitude | 0.036 | 0.742 | [-0.117, 0.175] |
| low | 10.5 | 6 | peak_magnitude | 0.036 | 0.742 | [-0.108, 0.193] |
| low | 10.6 | 1 | peak_magnitude | 0.036 | 0.742 | [-0.112, 0.189] |
| low | 10.6 | 2 | peak_magnitude | 0.036 | 0.742 | [-0.112, 0.199] |
| low | 10.6 | 3 | peak_magnitude | 0.036 | 0.742 | [-0.102, 0.188] |
| low | 10.6 | 5 | peak_magnitude | 0.036 | 0.742 | [-0.114, 0.181] |
| low | 10.6 | 6 | peak_magnitude | 0.036 | 0.742 | [-0.108, 0.195] |
| low | 10.7 | 1 | peak_magnitude | 0.036 | 0.742 | [-0.107, 0.197] |
| low | 10.7 | 2 | peak_magnitude | 0.036 | 0.742 | [-0.117, 0.187] |
| low | 10.7 | 3 | peak_magnitude | 0.036 | 0.742 | [-0.115, 0.189] |
| low | 10.7 | 5 | peak_magnitude | 0.036 | 0.742 | [-0.104, 0.187] |
| low | 10.7 | 6 | peak_magnitude | 0.036 | 0.742 | [-0.108, 0.187] |

### Dose-Response Trend Analysis

Tests for monotonic increase in effect size with |ΔΔG| bin (low → mid → high).

**Significant Trends Found**: 75

| Center | Width | Feature | low→mid→high |d| | Monotonic? | Spearman r | p-value |
|--------|-------|---------|----------------------|------------|------------|----------|
| 10.3 | 1% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.3 | 1% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.3 | 1% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.3 | 2% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.3 | 2% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.3 | 2% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.3 | 3% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.3 | 3% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.3 | 3% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.3 | 5% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.3 | 5% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.3 | 5% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.3 | 6% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.3 | 6% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.3 | 6% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.4 | 1% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.4 | 1% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.4 | 1% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.4 | 2% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.4 | 2% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.4 | 2% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.4 | 3% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.4 | 3% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.4 | 3% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.4 | 5% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.4 | 5% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.4 | 5% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.4 | 6% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.4 | 6% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.4 | 6% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.5 | 1% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.5 | 1% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.5 | 1% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.5 | 2% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.5 | 2% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.5 | 2% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.5 | 3% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.5 | 3% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.5 | 3% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.5 | 5% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.5 | 5% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.5 | 5% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.5 | 6% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.5 | 6% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.5 | 6% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.6 | 1% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.6 | 1% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.6 | 1% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.6 | 2% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.6 | 2% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.6 | 2% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.6 | 3% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.6 | 3% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.6 | 3% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.6 | 5% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.6 | 5% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.6 | 5% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.6 | 6% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.6 | 6% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.6 | 6% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.7 | 1% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.7 | 1% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.7 | 1% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.7 | 2% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.7 | 2% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.7 | 2% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.7 | 3% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.7 | 3% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.7 | 3% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.7 | 5% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.7 | 5% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.7 | 5% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |
| 10.7 | 6% | peak_magnitude | 0.036 → 0.045 → 0.215 | ✓ | 1.000 | 0.000 |
| 10.7 | 6% | band_power | 0.033 → 0.047 → 0.223 | ✓ | 1.000 | 0.000 |
| 10.7 | 6% | snr | 0.033 → 0.089 → 0.287 | ✓ | 1.000 | 0.000 |

---

## ACCEPTANCE CRITERIA EVALUATION

### Primary Criterion

**Threshold**: At least one (center, width, feature) in high ΔΔG bin achieves:
- |Cohen's d| ≥ 0.5
- FDR q-value < 0.05
- 95% CI excludes 0

**Result**: **FAIL**

Maximum |Cohen's d| in high bin: 0.287 (57.4% of threshold)

### Robustness Criterion

**Threshold**: Effect replicates across ≥2 adjacent centers or widths

**Result**: **FAIL**

### Specificity Criterion

**Threshold**: Off-band and shuffle controls remain non-significant

**Result**: **PASS**

---

## METHODOLOGY

### Data Source

- **Input**: `../../data/doench2016.csv`
- **Total Sequences**: 500 (after validation)
- **Validation**: Human DNA (A/C/G/T/N only)

### Pair Generation

- Single-point mutations generated randomly
- ΔΔG calculated using nearest-neighbor thermodynamics (SantaLucia 1998)
- Binning by |ΔΔG| tertiles

### Statistical Framework

- **Effect Size**: Cohen's d for paired design
- **Confidence Intervals**: BCa bootstrap (percentile method)
- **Multiple Testing**: Benjamini-Hochberg FDR correction
- **Trend Test**: Spearman correlation across bins
- **Permutation**: Label-shuffle within pairs

---

## REPRODUCIBILITY

### Exact Reproduction

```bash
python ddg_validation.py \\
    --input ../../data/doench2016.csv \\
    --output results/20251226_071117 \\
    --n-pairs 500 \\
    --bootstrap 1000 \\
    --permutation 1000 \\
    --seed 42
```

### Environment

- Python environment snapshot: `env.txt`
- Git commit: `eed07ff7610daf8bb680f31533837676fb61d380`
- Timestamp: `2025-12-26T07:10:26.878009`

---

## ARTIFACTS

All artifacts are in the output directory:

- `config.json` — Full configuration and parameters
- `pairs.csv` — WT-mutant pairs with ΔΔG and bins
- `stats.csv` — Statistical results per (bin, center, width, feature)
- `trend.csv` — Bin-level summaries and trend tests
- `env.txt` — Python environment (pip freeze)
- `log.txt` — Execution log
- `FINDINGS.md` — This report

---

## REFERENCES

- **Hypothesis Source**: https://github.com/zfifteen/dna-breathing-dynamics-encoding/pull/53
- **SantaLucia (1998)**: Nearest-neighbor thermodynamics for DNA
- **Helical Period**: ~10.5 bp/turn for B-DNA
- **Repository Policy**: `.github/REPOSITORY_POLICY.md`
- **Z Framework**: Discrete domain, Z = A(B/e²)

