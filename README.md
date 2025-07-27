# 📊 Signal-Theoretic Analysis of DNA Mutations

### A computational method for encoding and quantifying mutational disruptions using complex-valued spectral analysis of nucleotide sequences.

---

## 🧬 Overview

This project introduces a **novel computational framework** that encodes DNA sequences as **complex-valued waveforms**, enabling mutation analysis through **signal processing techniques**. Using this representation, we define **spectral disruption scores** to quantify how single-nucleotide variants alter the mathematical structure of DNA.

> ⚠️ This is a purely computational method — it does **not** model physical DNA vibrations or molecular dynamics.

---

## 🎯 Purpose

* Provide a **new feature space** for variant analysis and machine learning models
* Quantify mutational effects using **sequence-encoded spectral properties**
* Explore non-biological representations of DNA that may correlate with biological function

---

## ⚙️ Method Summary

### 1. **Sequence Encoding**

* Each nucleotide is mapped to a complex value:

  ```
  A → 1 + 0j
  T → -1 + 0j
  C → 0 + 1j
  G → 0 - 1j
  ```
* A synthetic waveform is generated using position-based phase modulation:

  $$
  Ψ_n = w_n \cdot e^{2πi s_n}
  $$

  where $s_n$ is the cumulative position using uniform or mutation-scaled spacing.

### 2. **Spectral Disruption from Mutation**

* For a given point mutation:

  * The waveform is rebuilt with local positional scaling (Z-tuning)
  * FFT is applied to extract spectral features
  * Differences from baseline include:

    * Δf₁: Frequency magnitude shift at selected harmonic
    * ΔEntropy: Spectral entropy change
    * ΔPeaks: Side-lobe count increase
* A **composite disruption score** is computed:

  $$
  \text{Score} = Z_n \cdot |\Delta f_1| + \text{ΔPeaks} + \text{ΔEntropy}
  $$

---

## 📦 Repository Structure

```
.
├── wave_crispr_signal.py         # Main script for mutation scoring
├── README.md                     # Documentation (this file)
├── notebooks/
│   └── validation.ipynb          # Benchmarking with public datasets (WIP)
├── data/                         # Sample sequences or variant datasets
└── examples/                     # Demonstration outputs
```

---

## ✅ What This Method **Is**

* A **mathematical model** for encoding DNA as a symbolic waveform
* A tool to quantify **mutational disruption in the signal domain**
* A generator of **novel numerical features** for machine learning

---

## ❌ What This Method **Is Not**

* ❌ A model of DNA vibrational physics or THz spectroscopy
* ❌ A predictor of gene function or expression on its own
* ❌ A substitute for biochemical or chromatin-based CRISPR scoring

---

## 📊 Experimental Validation Plan

### 🔬 Phase 1: Correlation Studies

* **ClinVar** pathogenicity comparison (ROC/AUC)
* CRISPR guide efficiency datasets (correlation with activity)
* PhyloP conservation alignment (signal vs. conservation)

### 🔬 Phase 2: Functional Overlays

* eQTL effect size correlation (GTEx)
* TF binding and chromatin accessibility (ENCODE)

### 🔬 Phase 3: Predictive Modeling

* Integrate disruption scores with:

  * CADD / DeepSEA / Basenji features
  * CRISPR base editing prediction tools
* Evaluate predictive performance improvements

---

## 📈 Success Criteria

| Level      | Outcome                                                                                |
| ---------- | -------------------------------------------------------------------------------------- |
| ✅ Minimal  | Statistically significant correlation (p < 0.001) with at least one biological dataset |
| ✅ Moderate | Comparable or superior to baseline conservation metrics                                |
| ✅ High     | Improvement in existing predictive models or discovery of new interpretable patterns   |

---

## 📚 Usage

```bash
python wave_crispr_signal.py
```

Outputs top mutational "hotspots" in terms of spectral disruption score.

*NOTE:*
This tool currently uses a hardcoded 150bp mock sequence. Swap in real sequences for production use.

---

## 🧠 License & Attribution

MIT License.
Original concept developed under the reframed idea:

> “Spectral Disruption Profiling (SDP) for DNA Sequence Analysis”
