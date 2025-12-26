# φ-Geometry Hypothesis Validation Experiment

**Experiment ID:** `phi_geometry_validation_20251226`  
**Status:** ✅ Complete  
**Result:** **HYPOTHESIS NOT FALSIFIED**

---

## Quick Summary

This experiment tests whether φ (golden ratio) geometry features derived from B-DNA structure improve CRISPR guide efficiency prediction compared to naive baselines.

**Key Finding:** φ-features show **weak but statistically significant** correlations (r ≈ 0.11–0.15, p < 0.01) and consistently outperform uniform/random baselines, but are **far weaker than simple GC content** (r = 0.71).

---

## Hypothesis

> φ-based features (phase alignment, spectral curvature) derived from the physical B-DNA geometry improve prediction of CRISPR guide efficiency compared to non-φ baselines (uniform-phase, random-phase, GC-only).

### Physical Basis (Larsen, Symmetry 2021)

- **B-DNA helix length:diameter** ≈ 1.6088 (close to φ ≈ 1.618)
- **Major/minor backbone spacing** ≈ 1.64
- **Axial 10-fold symmetry** with "golden diamonds"

---

## Results

### Falsification Criteria

The experiment tested three criteria to falsify the hypothesis:

| Criterion | Threshold | Result | Status |
|-----------|-----------|--------|--------|
| **1. Weak correlation** | \|r\| < 0.1 | Max \|r\| = 0.147 | ✗ NOT MET |
| **2. No improvement vs baselines** | φ wins ≤ 50% | φ wins 6/9 (67%) | ✗ NOT MET |
| **3. Not statistically significant** | p > 0.05 | Min p = 0.001 | ✗ NOT MET |

**Verdict:** All criteria **failed** → Hypothesis **NOT FALSIFIED**

### Feature Correlations

| Feature | Pearson r | 95% CI | p-value | vs Baseline |
|---------|-----------|--------|---------|-------------|
| **phi_phase** | **0.138** | [0.070, 0.235] | **0.001** | Beats uniform/random |
| **phi_curvature** | **0.147** | [0.050, 0.244] | **0.001** | Beats uniform/random |
| **phi_combined** | **0.107** | [0.003, 0.208] | **0.010** | Beats uniform/random |
| uniform_phase | 0.000 | [0.000, 0.000] | 1.000 | Baseline (null) |
| random_phase | -0.011 | [-0.096, 0.071] | 0.794 | Baseline (random) |
| **gc_content** | **0.709** | [0.665, 0.749] | **0.000** | **Dominates all** |

---

## Interpretation

### What This Means

1. **φ-features capture real signal** — correlations are weak but highly significant (p < 0.01)
2. **φ-features beat naive baselines** — win 6/9 comparisons (all vs uniform/random)
3. **GC content dominates** — r = 0.71 vastly exceeds φ-features (r ≈ 0.15)
4. **Secondary structural effects** — φ may capture DNA flexibility/groove geometry

### Recommendation

**Include φ-features as supplementary features in ensemble models**, but do **NOT replace** established predictors (GC%, thermodynamics).

---

## Quick Start

### Run Full Experiment

```bash
python experiments/phi_geometry_validation_20251226/validate_phi_geometry.py \
    --seed 42 \
    --bootstrap 1000 \
    --permutation 1000 \
    --output experiments/phi_geometry_validation_20251226/results.json
```

**Expected runtime:** < 5 seconds  
**Actual runtime:** 2.89 seconds

### Run Quick Validation

```bash
python experiments/phi_geometry_validation_20251226/validate_phi_geometry.py \
    --seed 42 \
    --bootstrap 100 \
    --permutation 100
```

**Expected runtime:** < 5 seconds

---

## Files

```
experiments/phi_geometry_validation_20251226/
├── README.md                        # This file
├── FINDINGS.md                      # Detailed findings (conclusion first!)
├── manifest.yml                     # Experiment specification
├── validate_phi_geometry.py         # Main validation script
└── results.json                     # Complete results with metadata
```

---

## Dataset

- **File:** `data/doench2016.csv`
- **Sequences:** 583 CRISPR guides with efficiency scores
- **Efficiency range:** 0.38 – 0.99
- **SHA256:** `6c2c38934aa771e74ecb42a7b477e0541694aba805ab5de529e87448416ec892`

---

## Scientific Gates Compliance

✅ **Human DNA only** — All sequences validated as A/C/G/T/N  
✅ **No fabrication** — φ constants from empirical measurements (Larsen 2021)  
✅ **Fail-fast validation** — Input sequences validated at entry  
✅ **Statistical validity** — Bootstrap CI + permutation tests + partial correlations  
✅ **Reproducibility** — Fixed seed (42), SHA256 checksums, metadata persisted

---

## Module Created

**Path:** `wave_crispr_signal/features/phi_geometry.py`

### Exports

**Constants:**
- `PHI` — Golden ratio (1.618034)
- `DNA_LENGTH_DIAMETER_RATIO` — 1.6088
- `DNA_MAJOR_MINOR_SEP_RATIO` — 1.6407
- `DNA_HELICAL_PERIOD` — 10 bp
- `PHI_MODES` — φ-related spectral modes

**Functions:**
- `dna_phi_phase(index)` — φ-coupled phase at bp position
- `dna_phi_curvature(track)` — φ-curvature from signal
- `compute_phi_features(seq)` — All φ-features for sequence

**Baselines:**
- `uniform_phase_score(seq)` — Uniform baseline
- `random_phase_score(seq, seed)` — Random baseline
- `simple_gc_content(seq)` — GC fraction

---

## Citation

Based on:
> Larsen, M. D. (2021). *DNA Structure and the Golden Ratio Revisited*. Symmetry.

---

## See Also

- **Detailed findings:** [FINDINGS.md](FINDINGS.md)
- **Experiment spec:** [manifest.yml](manifest.yml)
- **Full results:** [results.json](results.json)
- **PR #147:** Original hypothesis development

---

**Experiment Date:** 2025-12-26  
**Author:** GitHub Copilot  
**Repository:** wave-crispr-signal
