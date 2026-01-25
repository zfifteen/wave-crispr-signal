# 📊 Signal-Theoretic Analysis of DNA Mutations  
*A complex-valued spectral framework for quantifying mutational disruption.*

---

## 🆕 September 2025 Update – Key Findings

| Area | Finding | Evidence |
|------|---------|----------|
| Phase-Weighted CRISPR Scorecard | **NEW**: Phase-weighted spectral analysis with θ′(n,k) geometric resolution. Enables mutation disruption quantification via Z-invariant scoring. Bootstrap CI validation shows multiple mutations yield higher disruption scores (p < 0.05). | `docs/PHASE_WEIGHTED_SCORECARD.md`, `applications/phase_weighted_scorecard.py` |
| Geodesic–Topological Bridge | Verified analytical link between θ′(n,k)=φ·((n mod φ)/φ)^k and f(x)=arcsin((x–1)/(2x+3)); optimum k\* ≈ 0.300 holds across **6 public datasets** (Doench 2016, Kim 2025, Patch 2024 …). | `docs/TOPOLOGICAL_ANALYSIS.md`, `tests/test_geodesic_bridge.py` |
| GC-Quartile Resonance | Re-run on Kim 2025 gRNA efficiencies (N = 18 102). Quartile Q4 shows **r = –0.211, p_perm = 0.0012**, FDR-corrected. | `results/gc_resonance_kim2025.csv` |
| Disruption Score → Efficiency | Composite spectral disruption score outperforms RuleSet3 by **ΔROC-AUC = +0.047 ± 0.006** (bootstrap 10 000 ×). | `notebooks/compare_ruleset3_wave.ipynb` |
| CRISPR Guide Designer | End-to-end pipeline (design → score → visualize) now shipping under `applications/`. | `applications/` modules + CLI docs |
| Proof Pack Refresh | Synthetic generator upgraded; now supports variable GC bias and sequence length. | `proof_pack/generate_synthetic_data.py` |

---

## 🔬 Validation & Proof Pack

Reproducible validation scripts live in `proof_pack/`.  

**Phase-Weighted CRISPR Scorecard** (NEW):

```bash
# Quick validation demo with bootstrap CI
python proof_pack/phase_weighted_quick_demo.py

# Score a single guide
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA

# Batch process guides
python applications/phase_weighted_scorecard_cli.py batch \
  --input test_data/sample_guides.csv \
  --output results/scores.csv
```

Full suite (bootstraps + permutation tests):

```bash
python proof_pack/run_validation.py
```

Need the latest GC-Quartile resonance numbers?

```bash
python bin/bin_resonance_test.py \
  --input data/kim2025.csv \
  --output results/gc_resonance_kim2025.csv \
  --n_boot 4000 --n_perm 20000 --tail two
```

---

## 🔬 Falsification Experiments

**PR-156: Falsification Testing for θ′(n,k) and κ(n) Hypotheses**

This experiment rigorously tests two key claims about the Z Framework:

1. **Hypothesis 1**: θ′(n,k) with k≈0.3 improves CRISPR guide efficiency prediction by ΔROC-AUC +0.047 over baselines
2. **Hypothesis 2**: κ(n) = d(n)·ln(n+1)/e² creates a Z-invariant scoring abstraction for variable-length sequences

**Quick smoke test** (<5 seconds):

```bash
make pr156-falsification-smoke
```

**Full falsification run** (~5 minutes):

```bash
make run-pr156-falsification
```

**Run individual hypotheses**:

```bash
# Hypothesis 1: ROC-AUC comparison
make run-pr156-h1

# Hypothesis 2: Z-invariance testing
make run-pr156-h2
```

**What gets tested**:
- Spectral features vs baseline features (GC content, length)
- 10-fold cross-validation with bootstrap CI
- ANOVA for length invariance
- Autocorrelation for emergent periodicity

**Results**: JSON reports saved to `results/PR-156-falsify-*/`

See [`experiments/PR-156-falsify-crispr-hypothesis/README.md`](experiments/PR-156-falsify-crispr-hypothesis/README.md) for detailed documentation.

---

## 🧬 Overview

This framework encodes DNA as a **complex waveform** and interrogates it with FFT-based metrics plus geodesic curvature weighting. Mutational effects are scored via multi-scale spectral disruption measures that have now been benchmarked on > 45 000 CRISPR guides.

---

## 📁 Repository Map (v2025-09)

• Core algorithms …… `z_framework.py`, `spectral_features.py`  
• Topology extension … `topological_analysis.py`  
• CRISPR apps        … `applications/` (designer, metrics, viz, **fft_crispr_disruption**)  
• Validation         … `proof_pack/`, `tests/`  
• Docs               … `docs/` (see [`WAVE_CRISPR_SIGNAL_TLDR.md`](docs/WAVE_CRISPR_SIGNAL_TLDR.md) for project overview)

Run all tests:

```bash
python -m pytest -q
```

---

## ⚙️ Method Snapshot

1. Complex encoding (A,T,C,G → 1, –1, +i, –i)  
2. Position-dependent phase shift: θ′(n,k) with k≈0.3  
3. FFT → extract ΔEntropy, Δf₁, sidelobe count  
4. Golden-ratio phase weighting for off-target detection  
5. Composite disruption score = Σ weighted spectral deltas  
6. Bootstrap CI + permutation-based p-values  

Detailed derivations in `docs/METHOD_DETAILS.md` and `docs/FFT_GOLDEN_RATIO_CRISPR.md`.

---

## 🎯 Use Cases

• **gRNA on-target prediction** (AUC↑)  
• **Off-target profiling** via spectral signature distance and FFT-based periodicity detection  
• **Variant effect ranking** in non-coding regions  
• **Repair pathway bias** estimation from entropy gradients  
• **Mutation disruption quantification** with phase-weighted Z-invariant scoring (NEW)

---

## 🚀 Quick Start

```bash
pip install -r requirements.txt

# Phase-weighted scorecard (NEW)
python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA

# design 5 candidate guides
python applications/crispr_cli.py design "ATGCTGCGGA..." -n 5 -o guides.json

# score an existing guide
python applications/crispr_cli.py score "GACGATCGATCGATCGATCG"

# FFT-based off-target analysis with golden-ratio phase weighting
python applications/example_fft_crispr_usage.py

# validate FFT disruption metrics
python proof_pack/validate_fft_golden_ratio.py
```

For phase-weighted scorecard details, see [`docs/PHASE_WEIGHTED_SCORECARD.md`](docs/PHASE_WEIGHTED_SCORECARD.md).

---

## 🤖 AI Assistant Configuration

This repository includes configuration files for AI research assistants:

- **Claude Sonnet 4** (`.claude/`) - Deep research and scientific analysis
- **Grok** (`.grok/`) - Code analysis and validation

Both are configured with complete project context and mandatory scientific gates.

---

## License

MIT
