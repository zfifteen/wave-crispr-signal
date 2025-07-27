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
