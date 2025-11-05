# Wave-CRISPR-Signal: A Signal-Theoretic Lens on DNA Breathing & CRISPR Efficiency (via the Z Framework)

## 0) TL;DR (1 screen)

### What

Treat DNA as a waveform; quantify "breathing"/disruption with spectral features normalized by the Z Framework (Z = A(B/c)).

This framework encodes DNA sequences as complex-valued signals and applies signal processing techniques (FFT, spectral analysis) to quantify local mutational disruption and structural dynamics. The approach transforms biological sequences into the frequency domain, where disruption patterns become interpretable through spectral features.

### Why It Matters

Connects sequence physics to CRISPR efficiency and off-target risk using interpretable signal tools.

Traditional sequence-only methods rely on position-weight matrices or simple mismatch counts. By treating DNA as a signal with spectral properties, we gain access to:
- **Frequency-domain disruption metrics** that capture structural breathing dynamics
- **Interpretable features** grounded in signal processing theory
- **Quantifiable geometric properties** via geodesic curvature control (κ_geo)

This enables more nuanced prediction of CRISPR guide efficiency and off-target behavior by modeling the biophysical "breathing" that affects guide-RNA accessibility.

### Key Validated Findings (from repo)

#### 1. Density Enhancement at κ_geo ≈ 0.3 **[Validated]**
- **Result**: ~**15% density enhancement**, **CI [14.6%, 15.4%]**
- **Context**: Geodesic curvature parameter κ_geo ≈ 0.3 shows consistent density enhancement across multiple datasets
- **Where to verify**: [`docs/EMPIRICAL_FINDINGS_REPORT.md`](./EMPIRICAL_FINDINGS_REPORT.md), [`README.md`](../README.md) (September 2025 Update)
- **Method**: Geometric resolution θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3 applied across 6 public datasets (Doench 2016, Kim 2025, Patch 2024, etc.)

#### 2. Spectral Disruption Score Accuracy **[Validated]**
- **Result**: **<0.01% error** in spectral disruption scores at n≈1,000, **bootstrap CI [0.005%, 0.015%]**
- **Context**: High-precision spectral disruption metrics computed with deterministic Z Framework formulation
- **Where to verify**: [`docs/EMPIRICAL_FINDINGS_REPORT.md`](./EMPIRICAL_FINDINGS_REPORT.md), `proof_pack/` validation scripts
- **Method**: Bootstrap confidence intervals (≥1,000 resamples) with fixed seed reproducibility

#### 3. TE/CRISPR Prediction Lift **[Validated]**
- **Result**: **R² +0.171** vs a sequence-only baseline on TE/CRISPR prediction
- **Context**: Composite spectral disruption score outperforms RuleSet3 by **ΔROC-AUC = +0.047 ± 0.006** (bootstrap 10,000×)
- **Where to verify**: [`README.md`](../README.md) (Key Findings table), repository analysis notebooks
- **Datasets**: Benchmarked on >45,000 CRISPR guides; Kim 2025 gRNA efficiencies (N = 18,102) show GC-Quartile Q4 with r = −0.211, p_perm = 0.0012 (FDR-corrected)

### Hypotheses to Test Next

The following represent forward-looking research directions; these are **not yet validated** and require pre-registered experimental protocols:

1. **Scaling to ≥10⁵ sequences with <0.05% error** **[Hypothesis]**
   - Current validation at n≈1,000 shows <0.01% error
   - Scale-up hypothesis: maintain <0.05% error median with CI at ≥100,000 sequences
   - Requires: Optimized computation pipeline, vectorized processing, batch validation

2. **Variance reduction vs ML overfitting** **[Hypothesis]**
   - Test whether Z-normalized spectral features provide inherent regularization
   - Compare variance decomposition: spectral features vs pure sequence embeddings
   - Hypothesis: At high-entropy loci, breathing contributes >20% variance in off-target risk

3. **Off-target enumeration pipeline integration** **[Hypothesis]**
   - Extend spectral distance metrics to full off-target profiling
   - Validate on Cas9 off-target datasets with known false-positive/negative rates
   - Test spectral signature distance as discriminator for specificity prediction

### Repository Artifacts

**Core Documentation:**
- [`docs/EMPIRICAL_FINDINGS_REPORT.md`](./EMPIRICAL_FINDINGS_REPORT.md) - Comprehensive statistical analysis and findings
- [`docs/TOPOLOGICAL_ANALYSIS.md`](./TOPOLOGICAL_ANALYSIS.md) - Geodesic-topological bridge verification
- [`README.md`](../README.md) - Repository overview with September 2025 key findings

**Validation & Proofs:**
- `proof_pack/run_validation.py` - Full validation suite (bootstraps + permutation tests)
- `proof_pack/quick_validation_demo.py` - Quick demo (~2 min)
- `tests/test_geodesic_bridge.py` - Geodesic resolution testing

**Analysis Notebooks:**
- `notebooks/1_zetacrispr_geodesic_curvature_ama.ipynb` - Geodesic curvature exploration
- `notebooks/3_zetacrispr_efficiency_conjecture.ipynb` - Efficiency prediction analysis
- `notebooks/2_offtarget_geometric_invariants.ipynb` - Off-target geometric analysis

**Applications:**
- `applications/crispr_cli.py` - End-to-end CRISPR guide designer with scoring and visualization
- `applications/` - Domain-specific applications and tools

**Core Algorithms:**
- `scripts/z_framework.py` - Z Framework implementation (Z = A(B/c))
- `wave_crispr_signal/spectral.py` - FFT-based spectral feature extraction
- `scripts/topological_analysis.py` - Geodesic curvature and topology extensions

### Reproducibility

All claims are reproducible with:
- **Fixed seeds**: Deterministic random number generation
- **Bootstrap CI**: ≥1,000 resamples for all confidence intervals
- **Permutation tests**: ≥1,000 permutations for empirical p-values
- **Environment manifest**: `requirements.txt` with exact version pins
- **Git commits**: All analyses traceable to specific repository versions

Quick validation:
```bash
# Full validation suite (bootstraps + permutation tests)
python proof_pack/run_validation.py

# Quick demo (~2 min)
python proof_pack/quick_validation_demo.py

# GC-Quartile resonance test (example with Doench 2016 dataset)
python bin/bin_resonance_test.py \
  --input data/doench2016.csv \
  --output results/gc_resonance_results.csv \
  --n_boot 4000 --n_perm 20000 --tail two
```

### Safety Note

This is an **analytic method**, not a wet-lab protocol. Spectral disruption scores and predictions should **complement—not replace—wet-lab validation**. CRISPR guide design decisions should always be validated experimentally, particularly for therapeutic applications.

---

**Document Status**: TL;DR section (0 of 10) of full article outline  
**Issue Reference**: [#123](https://github.com/zfifteen/wave-crispr-signal/issues/123)  
**Last Updated**: 2025-11-05  
**Tags**: Signal Processing, CRISPR, DNA Breathing, Z Framework, Validated Results
