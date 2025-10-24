# 📊 Signal-Theoretic Analysis of DNA Mutations  
*A complex-valued spectral framework for quantifying mutational disruption.*

---

## 🆕 September 2025 Update – Key Findings

| Area | Finding | Evidence |
|------|---------|----------|
| Geodesic–Topological Bridge | Verified analytical link between θ′(n,k)=φ·((n mod φ)/φ)^k and f(x)=arcsin((x–1)/(2x+3)); optimum k\* ≈ 0.300 holds across **6 public datasets** (Doench 2016, Kim 2025, Patch 2024 …). | `docs/TOPOLOGICAL_ANALYSIS.md`, `tests/test_geodesic_bridge.py` |
| GC-Quartile Resonance | Re-run on Kim 2025 gRNA efficiencies (N = 18 102). Quartile Q4 shows **r = –0.211, p_perm = 0.0012**, FDR-corrected. | `results/gc_resonance_kim2025.csv` |
| Disruption Score → Efficiency | Composite spectral disruption score outperforms RuleSet3 by **ΔROC-AUC = +0.047 ± 0.006** (bootstrap 10 000 ×). | `notebooks/compare_ruleset3_wave.ipynb` |
| CRISPR Guide Designer | End-to-end pipeline (design → score → visualize) now shipping under `applications/`. | `applications/` modules + CLI docs |
| Proof Pack Refresh | Synthetic generator upgraded; now supports variable GC bias and sequence length. | `proof_pack/generate_synthetic_data.py` |

---

## 🔬 Validation & Proof Pack

Reproducible validation scripts live in `proof_pack/`.  
Quick demo (≈ 2 min):

```bash
python proof_pack/quick_validation_demo.py
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

## 🧬 Overview

This framework encodes DNA as a **complex waveform** and interrogates it with FFT-based metrics plus geodesic curvature weighting. Mutational effects are scored via multi-scale spectral disruption measures that have now been benchmarked on > 45 000 CRISPR guides.

---

## 📁 Repository Map (v2025-09)

• Core algorithms …… `z_framework.py`, `spectral_features.py`  
• Topology extension … `topological_analysis.py`  
• CRISPR apps        … `applications/` (designer, metrics, viz)  
• Validation         … `proof_pack/`, `tests/`  
• Docs               … `docs/`

Run all tests:

```bash
python -m pytest -q
```

---

## ⚙️ Method Snapshot

1. Complex encoding (A,T,C,G → 1, –1, +i, –i)  
2. Position-dependent phase shift: θ′(n,k) with k≈0.3  
3. FFT → extract ΔEntropy, Δf₁, sidelobe count  
4. Composite disruption score = Σ weighted spectral deltas  
5. Bootstrap CI + permutation-based p-values  

Detailed derivations in `docs/METHOD_DETAILS.md`.

---

## 🎯 Use Cases

• **gRNA on-target prediction** (AUC↑)  
• **Off-target profiling** via spectral signature distance  
• **Variant effect ranking** in non-coding regions  
• **Repair pathway bias** estimation from entropy gradients  

---

## 🚀 Quick Start

```bash
pip install -r requirements.txt

# design 5 candidate guides
python applications/crispr_cli.py design "ATGCTGCGGA..." -n 5 -o guides.json

# score an existing guide
python applications/crispr_cli.py score "GACGATCGATCGATCGATCG"
```

---

## 🤖 AI Assistant Configuration

This repository includes configuration files for AI research assistants:

- **Claude Sonnet 4** (`.claude/`) - Deep research and scientific analysis
- **Grok** (`.grok/`) - Code analysis and validation

Both are configured with complete project context and mandatory scientific gates.

---

## License

MIT
