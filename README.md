# 📊 Signal-Theoretic Analysis of DNA Mutations

### A computational method for encoding and quantifying mutational disruptions using complex-valued spectral analysis of nucleotide sequences.

---

## 🧬 Overview

This project introduces a **novel computational framework** that encodes DNA sequences as **complex-valued waveforms**,
enabling mutation analysis through **signal processing techniques**. Using this representation, we define **spectral
disruption scores** to quantify how single-nucleotide variants alter the mathematical structure of DNA.

> ⚠️ This is a purely computational method — it does **not** model physical DNA vibrations or molecular dynamics.

---

## 🎯 Purpose

* Provide a **new feature space** for variant analysis and machine learning models
* Quantify mutational effects using **sequence-encoded spectral properties**
* Explore non-biological representations of DNA that may correlate with biological function

---

# Advantages

Complex-valued waveform encoding offers several notable advantages over traditional methods for analyzing DNA sequences
and mutations:

- **Captures Phase and Amplitude Information**
  Traditional DNA encodings (e.g., one-hot vectors) are purely real-valued and typically binary, recording only the
  presence or absence of each nucleotide at a position. In contrast, complex-valued encoding assigns both real and
  imaginary components to each nucleotide, allowing representation of both amplitude and phase. This richer
  representation enables the capture of more nuanced sequential and structural information in the DNA.
- **Enables Frequency-Domain Analysis**
  By encoding DNA as a synthetic complex waveform, it becomes possible to apply Fourier/spectral analysis directly to
  genomic data. This facilitates the detection of periodicities, motifs, harmonics, and global sequence features that
  are otherwise difficult to summarize quantitatively using classical metrics.
- **Quantifies Mutational Disruption Mathematically**
  Changes caused by single-nucleotide variants can be precisely described in the frequency domain (via magnitude shifts
  at specific harmonics, entropy changes, and alterations in spectral peaks). Such mathematical descriptors provide a
  novel, interpretable, and continuous feature space for mutation effect prediction and machine learning
  pipelines—enabling distinction between subtle and dramatic mutational impacts in a single framework.
- **Retains Spatial Relationships**
  The position-based phase modulation in the encoding scheme means that changes at different sequence positions yield
  distinct spectral fingerprints. This property is valuable for modeling spatial context and mutation locality, which
  can be important in both regulatory and coding regions.
- **Promotes Generalization and Integrability**
  Since spectral features are independent of sequence length and can be combined or compared across sequences, this
  approach is easily integrated into scalable machine learning models for genome-wide analyses.

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

## ✅ What This Method **Is**

* A **mathematical model** for encoding DNA as a symbolic waveform
* A tool to quantify **mutational disruption in the signal domain**
* A generator of **novel numerical features** for machine learning

---
## Why Some CRISPR Guides Fail
The signal-theoretic techniques demonstrated in the script offer innovative approaches to solving key CRISPR research problems through **mathematical modeling of DNA sequence properties**. Here's how they address specific challenges:

### 1. **gRNA On-Target Efficiency Prediction**
- **Problem**: Existing tools (e.g., RuleSet3, DeepHF) rely on sequence motifs but ignore global structural properties.
- **Solution**:  
  Spectral entropy and sidelobe metrics capture **sequence harmonic stability** – low-entropy regions with dominant frequencies (like the baseline plot shows) indicate structurally stable DNA sites. High `hotspot_score` mutations disrupt these harmonics, potentially identifying:
  - Optimal gRNA binding sites (high disruption → easy cleavage)
  - Fragile genomic contexts vulnerable to off-targeting

### 2. **Off-Target Effect Identification**
- **Problem**: Off-target sites often share local sequence similarity but differ in global structure.
- **Solution**:  
  The **complex waveform encoding** (phase-modulated cumulative positions) detects subtle structural differences:
  ```python
  # Position-dependent phase scaling
  s_n = d * (1 + zn_map.get(i, 0))  # zn = position-dependent scaling
  Ψ_n = w_n * np.exp(2j * np.pi * s_n) 
  ```
  - Compare target vs. off-target spectral signatures using KL divergence of FFT magnitudes
  - Sites with similar sequences but different spectral entropy/peak profiles would be flagged

### 3. **Predicting Functional Consequences**
- **Problem**: Non-coding variants are hard to interpret; most tools focus on amino acid changes.
- **Solution**:  
  The composite disruption score quantifies **regulatory impact**:
  - High `hotspot_score` at promoter/enhancer regions → likely disruptive regulatory variants
  - Correlation studies could link entropy changes (ΔEntropy) with epigenetic modifications

### 4. **CRISPR Repair Outcome Bias**
- **Problem**: Indel profiles vary by genomic context (microhomology, sequence stability).
- **Solution**:  
  **Spectral stability metrics** predict repair tendencies:
  - Low-entropy sites → template-independent NHEJ bias
  - High sidelobe regions → microhomology-mediated repair (MMEJ)

### 5. **Multiplexed gRNA Design**
- **Problem**: Simultaneous cuts may cause chromosomal rearrangements.
- **Solution**:  
  Use **spectral coherence analysis**:
  ```python
  # Cross-correlation of gRNA target spectra
  corr = np.correlate(spec_gRNA1, spec_gRNA2, mode='same')
  ```
  - High spectral correlation → risk of genomic instability
  - Prioritize gRNA pairs with orthogonal spectral signatures

### Validation & Integration Pathway
1. **Benchmarking**:  
   Train ML models using spectral features against datasets like:
   - CRISPR-Cas9 efficiency (e.g., Doench 2016)
   - Off-target cleavage data (e.g., GUIDE-seq)

2. **Biological interpretability**:  
   Correlate spectral metrics with:
   ```python
   # Epigenetic feature integration
   d_adjusted = d * (1 + dnase_signal[pos])  # Scale spacing by openness
   ```
   - Chromatin accessibility (ATAC-seq)
   - DNA shape parameters (minor groove width)

3. **Tool integration**:  
   Embed in existing pipelines:
   ```bash
   crisprScore --spectral-features input.fa
   ```

### Key Advantages Over Traditional Methods
- **Positional context sensitivity**: `zn` scaling weights central mutations more heavily
- **Phase awareness**: Complex encoding captures base transition dynamics
- **Length invariance**: FFT features enable cross-sequence comparison
- **Noise robustness**: Sidelobe thresholds ignore low-magnitude variations

### Limitations to Address
- **Parameter calibration**: Harmonic index (f1_index=10) needs optimization
- **Physical basis**: Incorporate DNA biophysical models (e.g., *ab initio* charge distributions)
- **Runtime**: Optimize FFT computation for genome-scale analysis

This approach provides a **mathematically rigorous framework** to model DNA as an information-carrying waveform – bridging digital sequence analysis and analog structural biology. By quantifying mutational disruptions in spectral space, it offers new dimensions for predicting CRISPR behavior beyond sequence-level patterns.
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
