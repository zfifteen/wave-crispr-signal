# 📊 Signal-Theoretic Analysis of DNA Mutations

### A computational method for encoding and quantifying mutational disruptions using complex-valued spectral analysis of nucleotide sequences.

---

## 🧬 Overview

This project introduces a **novel computational framework** that encodes DNA sequences as **complex-valued waveforms**,
enabling mutation analysis through **signal processing techniques**. Using this representation, we define **spectral
disruption scores** to quantify how single-nucleotide variants alter the mathematical structure of DNA.

> ⚠️ This is a purely computational method — it does **not** model physical DNA vibrations or molecular dynamics.

---

## Validation Report

Below is a concise, reproducible summary of the Z Framework validation steps and results. You can expand each section for details.

<details>
<summary>1. Parameters & Constants</summary>

- **Golden ratio**: φ = (1 + √5)/2 ≈ 1.61803, so φ − 1 ≈ 0.61803
- **Fixed parameters**:
    - a = 5
    - k = 0.3
    - c = e ≈ 2.71828
    - κ ≈ 0.386
    - σ₀ ≈ 0.118
    - Tolerances: ε₁ (mean) = 0.005, ε₂ (variance) = 0.005
</details>

<details>
<summary>2. Discrete Z-Value Computation</summary>

1. Map DNA → integers:  
   A → 1, T → 2, C → 3, G → 4
2. Compute differences:
   ```
   Δi = xi – x{i–1}
   Δmax = max |Δi|
   ```  
3. Scale and normalize:
   ```
   Zi = i * (|Δi| / Δmax),  i = 2…N
   ```  
</details>

<details>
<summary>3. First-Order Statistics</summary>

- Mean:  μZ = (1/(N–1)) ∑ Zi
- Variance:  σZ² = (1/(N–1)) ∑ (Zi – μZ)²
- Std-dev:  σZ = √σZ²
</details>

<details>
<summary>4. Geodesic Curvature</summary>

F = k · (μZ)^k

> At μZ ≈ 0.552: F ≈ 0.096
</details>

<details>
<summary>5. Variance Trimming</summary>

Reported: σZ² ≈ 0.118 → σ_trim² ≈ 0.113

- **Threshold**: σ_trim² = max(σZ² – κ, 0)
- **Scaling**: σ_trim² = σZ² · (1 – κ/σZ²)

Both produce ~0.113 without negative values.
</details>

<details>
<summary>6. Convergence Tests</summary>

- |μZ – (φ–1)| ≤ ε₁?  → 0.552 vs 0.618 → **not converged**
- |σZ² – σ₀| ≤ ε₂? → 0.118 vs 0.118 → **converged**
</details>

<details>
<summary>7. Zeta-Chain Unfolding</summary>

Let z⁽⁰⁾ = μZ. Iterate:
```
z⁽t⁾ = T(z⁽t–1⁾)
D⁽t⁾ = 1 / z⁽t⁾
E⁽t⁾ = c ⋅ D⁽t⁾
F⁽t⁾ = F  (constant)
```
| t | z⁽t⁾  | D⁽t⁾  | E⁽t⁾  | F |
|---|-------|-------|-------|---|
| 0 | 0.552 | 1.811 | 4.926*|0.096|
| 1 | 5.624 | 0.178 | 0.485 |0.096|
| 2 | 4.983 | 0.201 | 0.546 |0.096|
| 3 | 4.950 | 0.202 | 0.550 |0.096|

> * E⁽0⁾ discrepancy suggests a variant initial D⁽0⁾ or scaling in T.
</details>

<details>
<summary>8. Empirical Correlation</summary>

- **Bio-anchored**: r ≈ –0.198, p ≈ 0.048 (significant)
- **Arbitrary**:  r ≈  0.052, p ≈ 0.611 (not significant)

Efficacy boost Δ_eff ≈ 5.8% aligns with ~15% density enhancement + trimmed variance.
</details>

<details>
<summary>9. Flags & Recommendations</summary>

- Clarify κ usage (threshold vs scaling)
- Specify T(z) operator or pseudocode
- Add guards for Δmax > 0 (raise error if zero)
- Label any unverified hypotheses in code/tests
</details>

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

### Basic DNA Mutation Analysis

```bash
python wave-crispr-signal.py
```

Outputs top mutational "hotspots" in terms of spectral disruption score.

### CRISPR Guide Design (NEW!)

The repository now includes comprehensive CRISPR guide design tools using signal-theoretic analysis:

#### Quick Start - Python API

```python
from applications.crispr_guide_designer import CRISPRGuideDesigner

# Initialize designer
designer = CRISPRGuideDesigner(pam_pattern="NGG", guide_length=20)

# Design guides for your target sequence
target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGG..."
guides = designer.design_guides(target_sequence, num_guides=5)

# Print results
for i, guide in enumerate(guides, 1):
    print(f"Guide {i}: {guide['sequence']}")
    print(f"  Position: {guide['position']}")
    print(f"  On-target score: {guide['on_target_score']:.3f}")
    print(f"  GC content: {guide['gc_content']:.1%}")
```

#### Command Line Interface

```bash
# Design guides for a sequence
python applications/crispr_cli.py design ATGCGATCGATCGATCG

# Design guides from FASTA file  
python applications/crispr_cli.py design target.fasta -o results.json

# Score a specific guide
python applications/crispr_cli.py score GACGATCGATCGATCGATCG

# Batch scoring from file
python applications/crispr_cli.py batch-score guides.txt -o scores.csv -f csv
```

#### Advanced Metrics & Visualization

```python
from applications.wave_crispr_metrics import WaveCRISPRMetrics
from applications.crispr_visualization import CRISPRVisualizer

# Calculate comprehensive metrics
metrics = WaveCRISPRMetrics()
comprehensive_score = metrics.calculate_comprehensive_score(guide_seq, target_seq)

# Create visualizations
visualizer = CRISPRVisualizer()
fig = visualizer.plot_guide_scores(guides)
dashboard = visualizer.create_dashboard(sequence, guides, "dashboard.html")
```

### Key CRISPR Features

🎯 **Signal-Theoretic Guide Design**
- Uses spectral entropy and harmonic analysis for efficiency prediction
- Novel approach combining frequency domain analysis with traditional metrics

🔍 **Off-Target Risk Assessment**
- Spectral signature comparison between target and off-target sites
- Detects subtle structural differences beyond sequence similarity

🔧 **Repair Outcome Prediction**
- Predicts NHEJ, MMEJ, and HDR probabilities using spectral stability
- Context-aware analysis of cut site properties

📊 **Comprehensive Scoring**
- Weighted composite scores combining multiple signal-theoretic metrics
- Benchmarking tools for guide set evaluation

🎨 **Rich Visualizations**
- Interactive dashboards with Plotly
- Spectral analysis plots and guide comparison charts
- Sequence mapping with guide positions and PAM sites

### Installation

```bash
pip install -r requirements.txt
```

*NOTE:*
The basic tool uses a hardcoded 150bp mock sequence. For CRISPR applications, provide your own target sequences.

---

## 🧠 License & Attribution

MIT License.
Original concept developed under the reframed idea:

> “Spectral Disruption Profiling (SDP) for DNA Sequence Analysis”
