# Implementation Checklist

**Experiment ID**: `spectral_signal_routing_fft_cluster_pr_174_2026_01_14`  
**Status**: Design Complete, Implementation Pending  
**Last Updated**: 2026-01-14

This document tracks the implementation status of the Spectral Signal Routing Hypothesis 1 Falsification Experiment. All work must remain traceable to the H1-strong falsification objectives defined in [TECH_SPEC.md](./TECH_SPEC.md).

---

## 1. Task Wiring

- [ ] **Select EDE task** with clear success/failure definition
  - Candidates: Factor localization, SAT solving (PR #174), other EDE-integrated tasks
  - Must define binary success criterion independent of spectral analysis
  - Document task selection rationale

- [ ] **Implement/reuse runner** for repeated seeded runs
  - Must support per-step snapshot recording
  - Must preserve RNG seed schemes for reproducibility
  - Must track run metadata (configuration hash, task parameters)

---

## 2. Trace Recorder

- [ ] **Implement generic trace recorder**
  - Output format: CSV, JSON, or binary tables (document choice)
  - Per-step cell-level fields:
    - [ ] Cell identifier (index or coordinate)
    - [ ] Strategy/behavior ID (integer-encoded)
    - [ ] Local fitness or error proxy
    - [ ] Other task-relevant scalars
  - Per-step global metadata:
    - [ ] Step index t
    - [ ] Task-level metrics (best score, constraint satisfaction, etc.)
    - [ ] Run identifier
    - [ ] Configuration hash
    - [ ] RNG seed

- [ ] **Document schema**
  - Field names and semantic meanings
  - Data types and units
  - Missing value handling

---

## 3. Series Extraction

- [ ] **Define series types** (pre-declare before success/failure labeling)
  - [ ] Task-level series: x_task(t)
    - E.g., best fitness, residual error
  - [ ] Global aggregation series: x_global(t)
    - E.g., mean/variance of cell-level field
  - [ ] Spatially localized series: x_cell_region(k, t)
    - Front band (first N cells)
    - Cluster cores (from EDE clustering logic)
    - Random cell samples (baseline)

- [ ] **Implement extraction pipeline**
  - Must not reference outcome labels during extraction
  - Document selection rules in configuration file
  - Ensure reproducibility (fixed seeds if randomness involved)

---

## 4. Spectral Analysis Harness

### 4.1 Windowing and Preprocessing

- [ ] **Window segmentation**
  - [ ] Choose T_window (window length)
  - [ ] Choose T_stride (stride/overlap)
  - [ ] Document choices and rationale
  - [ ] Ensure multiple windows fit typical run durations

- [ ] **Detrending** (optional)
  - [ ] Decide: apply or skip
  - [ ] If applied: choose method (mean subtraction, linear fit, etc.)
  - [ ] Apply consistently across all conditions

- [ ] **Tapering**
  - [ ] Choose window function (e.g., Hann, Hamming)
  - [ ] Apply to all windows identically

- [ ] **Normalization**
  - [ ] Choose normalization method (z-score, unit variance, etc.)
  - [ ] Apply per-run/per-channel consistently

### 4.2 FFT Computation

- [ ] **Implement FFT pipeline**
  - [ ] Compute complex coefficients X(f)
  - [ ] Extract magnitude spectrum |X(f)|
  - [ ] Extract phase arg(X(f)) if needed
  - [ ] Document sampling frequency and frequency resolution

### 4.3 Channel Detection

- [ ] **Define operational channel criteria**
  - [ ] Elevated magnitude threshold (e.g., N standard deviations above baseline)
  - [ ] Reproducibility threshold (fraction of windows, number of runs)
  - [ ] Minimum band width
  - [ ] Document all thresholds before success/failure labeling

- [ ] **Implement channel detection algorithm**
  - [ ] Identify contiguous frequency bands [f_low, f_high]
  - [ ] Apply elevation and reproducibility criteria
  - [ ] Return channel metadata per run

### 4.4 Run-Level Feature Aggregation

- [ ] **Channel count metrics**
  - [ ] Number of channels detected
  - [ ] Distribution across frequency ranges (low/mid/high)

- [ ] **Channel strength metrics**
  - [ ] Mean channel magnitude or energy
  - [ ] Maximum channel magnitude

- [ ] **Stability metrics**
  - [ ] Fraction of windows containing channels
  - [ ] Run-level continuity measures

- [ ] **Control metrics**
  - [ ] Total spectral energy
  - [ ] Low-vs-high frequency energy ratios
  - [ ] Spectral centroid and other shape descriptors

---

## 5. Condition Orchestration

### 5.1 Family A: High-Performance vs Low-Performance

- [ ] **A1: High-Performance Runs**
  - [ ] Define success criterion (task-specific, no spectral dependency)
  - [ ] Generate/collect runs meeting criterion within T_max
  - [ ] Record full snapshots and metadata

- [ ] **A2: Low-Performance Runs**
  - [ ] Identical macro-parameters to A1
  - [ ] Collect runs failing success criterion
  - [ ] Record full snapshots and metadata

### 5.2 Family B: Structured vs Randomized Controls

- [ ] **B1: Structured Dynamics**
  - [ ] Reuse A1/A2 data
  - [ ] Apply same FFT pipeline

- [ ] **B2: Phase-Scrambled Controls**
  - [ ] Implement phase randomization (preserve power spectrum)
  - [ ] Document randomization method
  - [ ] Generate controls from B1 data

- [ ] **B3: IID Random Controls**
  - [ ] Implement IID sampling (preserve marginal distributions)
  - [ ] Generate controls from empirical marginals
  - [ ] No temporal/spatial correlation

### 5.3 Family C: Spectral Disruption Interventions

- [ ] **C1: Baseline Replay**
  - [ ] Deterministic replay of A1 runs
  - [ ] Validate determinism assumptions
  - [ ] No perturbations

- [ ] **C2: Spectral Knock-Out**
  - [ ] Identify channels from A1 analysis
  - [ ] Implement inverse FFT with band attenuation/removal
  - [ ] Apply perturbations to cell states or inputs
  - [ ] Preserve amplitude distributions

- [ ] **C3: Non-Spectral Control Perturbations**
  - [ ] Implement broad-spectrum noise (Gaussian, uniform)
  - [ ] Match variance of C2 perturbations
  - [ ] Random cell flips/swaps (if applicable)

---

## 6. Collation and Downstream Hooks

- [ ] **Implement feature aggregation**
  - [ ] Merge per-run features into summary tables
  - [ ] Format: CSV, JSON, or similar for statistical analysis
  - [ ] Do not pre-label "success"/"failure" columns until analysis

- [ ] **Statistical analysis placeholders**
  - [ ] A1 vs A2 comparison (distinctiveness)
  - [ ] B1 vs B2/B3 comparison (specificity)
  - [ ] C1 vs C2 vs C3 comparison (causal necessity)
  - [ ] Document statistical methods when implemented

---

## 7. Documentation Hooks

- [ ] **Link to PR #174**
  - [ ] Reference experiment in PR description
  - [ ] Cross-link to TECH_SPEC.md

- [ ] **Create FINDINGS.md** (when results available)
  - [ ] Link back to TECH_SPEC.md
  - [ ] Document falsification outcomes
  - [ ] Include statistical evidence

- [ ] **Update experiment README**
  - [ ] High-level summary
  - [ ] Usage instructions
  - [ ] References to TECH_SPEC and TODO_IMPLEMENTATION

---

## 8. Reproducibility and Metadata

- [ ] **Environment capture**
  - [ ] Python version (3.12+)
  - [ ] Dependencies snapshot (pip freeze)
  - [ ] Git commit hash
  - [ ] Timestamp

- [ ] **Configuration persistence**
  - [ ] Save all analysis parameters (windowing, thresholds, etc.)
  - [ ] Save RNG seeds
  - [ ] Save dataset checksums (SHA256)

- [ ] **Results structure**
  - [ ] Create `results/<exp_id>/run-YYYYMMDD-HHMMSS/`
  - [ ] Save: results.json, results.csv, env.txt, log.txt

---

## 9. Validation and Testing

- [ ] **Smoke tests**
  - [ ] End-to-end pipeline on small synthetic dataset
  - [ ] Target: < 5 seconds completion
  - [ ] Validate output schema

- [ ] **Unit tests** (if applicable)
  - [ ] FFT computation correctness
  - [ ] Channel detection logic
  - [ ] Phase randomization preserves power spectrum
  - [ ] IID controls preserve marginals

- [ ] **Integration tests**
  - [ ] Full A/B/C condition pipeline
  - [ ] Reproducibility (same seed → same results)

---

## Notes

- **No fabrication**: All metrics, thresholds, and parameters must be justified and non-fabricated.
- **Pre-registration**: Analysis parameters must be fixed before labeling conditions.
- **Scope control**: New ideas must trace back to H1-strong falsification or reusable infrastructure.
- **Traceable work**: All changes reference TECH_SPEC.md or related documentation.

---

**Next Steps**:

1. Select EDE task (Section 1)
2. Implement trace recorder (Section 2)
3. Define and extract time series (Section 3)
4. Build spectral analysis harness (Section 4)

**Status Legend**:
- [ ] Not started
- [x] Completed
- [~] In progress
- [!] Blocked/needs decision
