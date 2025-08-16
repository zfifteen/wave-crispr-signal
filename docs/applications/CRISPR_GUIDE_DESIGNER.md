# CRISPR Guide Designer Documentation

## Overview

The CRISPR Guide Designer uses signal-theoretic analysis of DNA sequences to design and evaluate CRISPR guide RNAs. This approach combines traditional bioinformatics methods with novel spectral features derived from complex-valued DNA encodings.

## Key Features

- **Signal-Theoretic Analysis**: DNA sequences are encoded as complex waveforms and analyzed in frequency domain
- **On-Target Efficiency Prediction**: Uses spectral entropy and harmonic content to predict guide efficiency
- **Off-Target Risk Assessment**: Compares spectral signatures between target and off-target sites
- **Repair Outcome Prediction**: Predicts NHEJ, MMEJ, and HDR probabilities using spectral stability metrics
- **Comprehensive Scoring**: Combines multiple metrics for robust guide evaluation

## Installation

```bash
pip install -r requirements.txt
```

Required dependencies:
- numpy >= 1.21.0
- matplotlib >= 3.5.0  
- scipy >= 1.7.0
- scikit-learn >= 1.0.0
- plotly >= 5.0.0
- seaborn >= 0.11.0
- mpmath >= 1.2.0
- pandas >= 1.3.0

## Quick Start

### Basic Guide Design

```python
from applications.crispr_guide_designer import CRISPRGuideDesigner

# Initialize designer
designer = CRISPRGuideDesigner(pam_pattern="NGG", guide_length=20)

# Design guides for target sequence
sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGG..."
guides = designer.design_guides(sequence, num_guides=5)

# Print results
for i, guide in enumerate(guides, 1):
    print(f"Guide {i}: {guide['sequence']}")
    print(f"  Position: {guide['position']}")
    print(f"  Score: {guide['on_target_score']:.3f}")
    print(f"  GC Content: {guide['gc_content']:.1%}")
```

### Command Line Interface

```bash
# Design guides for a sequence
python applications/crispr_cli.py design ATGCGATCGATCGATCG

# Design guides from FASTA file
python applications/crispr_cli.py design target.fasta -o results.json

# Score a specific guide
python applications/crispr_cli.py score GACGATCGATCGATCGATCG

# Batch scoring
python applications/crispr_cli.py batch-score guides.txt -o scores.csv -f csv
```

## Core Methodology

### 1. DNA Sequence Encoding

Each nucleotide is mapped to a complex value:
- A → 1 + 0j
- T → -1 + 0j  
- C → 0 + 1j
- G → 0 - 1j

The sequence is then converted to a synthetic waveform using position-based phase modulation:

```
Ψₙ = wₙ × e^(2πi × sₙ)
```

where `sₙ` is the cumulative position with optional mutation-dependent scaling.

### 2. Spectral Analysis

The complex waveform is transformed using Fast Fourier Transform (FFT) to obtain:
- **Frequency spectrum**: Magnitude of FFT components
- **Spectral entropy**: Information content measure
- **Sidelobe count**: Number of significant spectral peaks
- **Harmonic content**: Power distribution across harmonics

### 3. Scoring Algorithms

#### On-Target Efficiency Score

Combines spectral features with traditional metrics:

```python
entropy_score = 1.0 - (spectral_entropy / 10.0)    # Lower entropy → higher efficiency
sidelobe_score = 1.0 - (sidelobes / spectrum_length) # Fewer sidelobes → higher efficiency  
gc_score = 1.0 - abs(gc_content - 0.5) * 2         # Optimal GC ~50%

on_target_score = entropy_score * 0.4 + sidelobe_score * 0.4 + gc_score * 0.2
```

#### Off-Target Risk Score

Uses spectral similarity between target and off-target sites:

```python
spectral_correlation = correlate(target_spectrum, off_target_spectrum)
sequence_similarity = hamming_distance(target_seq, off_target_seq)

risk_score = spectral_correlation * 0.6 + sequence_similarity * 0.4
```

#### Repair Outcome Prediction

Predicts repair pathway probabilities based on spectral stability:

- **Low entropy** → NHEJ bias (template-independent repair)
- **High sidelobes** → MMEJ bias (microhomology-mediated repair)  
- **High stability** → HDR efficiency (homology-directed repair)

## API Reference

### CRISPRGuideDesigner Class

#### Constructor
```python
CRISPRGuideDesigner(pam_pattern="NGG", guide_length=20)
```

#### Key Methods

**design_guides(sequence, num_guides=5)**
- Designs CRISPR guides for target sequence
- Returns list of guide dictionaries with scores and positions

**calculate_on_target_score(guide_seq)**  
- Calculates on-target efficiency score (0-1)
- Higher scores indicate better predicted efficiency

**calculate_off_target_risk(guide_seq, target_seq, off_target_seq)**
- Calculates off-target cleavage risk (0-1)
- Higher scores indicate higher risk

**predict_repair_outcomes(guide_seq, target_context)**
- Predicts repair pathway probabilities
- Returns NHEJ, MMEJ, and HDR efficiency estimates

**find_pam_sites(sequence)**
- Finds PAM sites in sequence using regex pattern
- Returns list of (position, PAM_sequence) tuples

### Advanced Metrics (WaveCRISPRMetrics)

#### Spectral Complexity Analysis
```python
from applications.wave_crispr_metrics import WaveCRISPRMetrics

metrics = WaveCRISPRMetrics()
complexity = metrics.calculate_spectral_complexity(guide_seq)

# Returns:
# - spectral_entropy: Information content
# - spectral_flatness: Wiener entropy  
# - spectral_centroid: Center of mass frequency
# - spectral_rolloff: 95% energy frequency
# - zero_crossing_rate: Spectrum variation
# - peak_count: Number of significant peaks
```

#### Harmonic Content Analysis
```python
harmonics = metrics.calculate_harmonic_content(guide_seq)

# Returns:
# - fundamental_frequency: Strongest frequency component
# - harmonic_power: Total harmonic energy
# - total_harmonic_distortion: Harmonic to fundamental ratio
# - harmonic_to_noise_ratio: Signal quality metric
```

#### Mutational Sensitivity
```python
sensitivity = metrics.calculate_mutational_sensitivity(guide_seq)

# Returns:
# - average_sensitivity: Mean disruption across positions
# - hotspot_positions: Most sensitive positions
# - position_effects: Per-position sensitivity scores
```

#### Comprehensive Scoring
```python
comprehensive = metrics.calculate_comprehensive_score(guide_seq, target_seq)

# Returns weighted composite score combining:
# - On-target efficiency (35%)
# - Spectral quality (25%) 
# - Mutational stability (20%)
# - Off-target safety (15%)
# - Repair confidence (5%)
```

## Visualization Tools

### Basic Plotting
```python
from applications.crispr_visualization import CRISPRVisualizer

visualizer = CRISPRVisualizer()

# Plot sequence spectrum
fig = visualizer.plot_sequence_spectrum(sequence)

# Plot guide scores  
fig = visualizer.plot_guide_scores(guides)

# Compare guides
fig = visualizer.plot_guide_comparison(guides)
```

### Interactive Dashboard
```python
# Create interactive Plotly dashboard
dashboard = visualizer.create_dashboard(sequence, guides, "dashboard.html")
dashboard.show()

# Interactive guide map
guide_map = visualizer.plot_interactive_guide_map(sequence, guides)
guide_map.show()
```

## Output Formats

### JSON Output (Default)
```json
{
  "sequence_name": "target_1",
  "sequence_length": 150,
  "guides": [
    {
      "sequence": "GACGATCGATCGATCGATCG",
      "position": 45,
      "pam_position": 65,
      "pam_sequence": "TGG",
      "on_target_score": 0.823,
      "gc_content": 0.6,
      "length": 20,
      "repair_outcomes": {
        "nhej_probability": 0.65,
        "mmej_probability": 0.25,
        "hdr_efficiency": 0.10
      }
    }
  ]
}
```

### CSV Output
```csv
sequence_name,guide_sequence,position,pam_sequence,on_target_score,gc_content,nhej_prob,mmej_prob,hdr_eff
target_1,GACGATCGATCGATCGATCG,45,TGG,0.823,0.600,0.650,0.250,0.100
```

## Best Practices

### Guide Design Strategy

1. **Target Multiple Sites**: Design 5-10 guides per target for redundancy
2. **Check GC Content**: Aim for 40-60% GC content
3. **Avoid Repetitive Sequences**: Screen for poly-T runs and repetitive elements
4. **Consider Context**: Include flanking sequences for off-target analysis
5. **Validate Experimentally**: Computational predictions should be validated

### Quality Thresholds

- **Excellent guides**: On-target score > 0.8, Off-target risk < 0.3
- **Good guides**: On-target score > 0.6, Off-target risk < 0.5  
- **Acceptable guides**: On-target score > 0.4, Off-target risk < 0.7

### Troubleshooting

#### Low Scores
- Check sequence quality and length
- Verify PAM sites are present
- Consider alternative PAM patterns (e.g., "NRG" for SpRY-Cas9)

#### High Off-Target Risk
- Use more stringent scoring thresholds
- Consider base editors or prime editors for precision
- Validate with orthogonal off-target detection methods

#### Poor Repair Outcomes
- Optimize target site selection for desired repair pathway
- Consider delivery timing and template design for HDR
- Use repair pathway modulators (e.g., small molecules)

## Advanced Usage

### Custom PAM Patterns
```python
# SpRY-Cas9 (relaxed PAM)
designer_spry = CRISPRGuideDesigner(pam_pattern="NRG", guide_length=20)

# Cas12a (TTTV PAM)  
designer_cas12a = CRISPRGuideDesigner(pam_pattern="TTT[ACG]", guide_length=23)

# xCas9 (relaxed NGG)
designer_xcas9 = CRISPRGuideDesigner(pam_pattern="NG", guide_length=20)
```

### Batch Processing
```python
# Process multiple targets
targets = {"gene1": seq1, "gene2": seq2, "gene3": seq3}
all_results = {}

for name, sequence in targets.items():
    guides = designer.design_guides(sequence)
    benchmark = metrics.benchmark_guide_set(guides, sequence)
    all_results[name] = benchmark
```

### Integration with Other Tools
```python
# Export for other CRISPR tools
def export_to_chopchop(guides):
    """Export guides for CHOPCHOP analysis."""
    return [guide['sequence'] for guide in guides]

def export_to_crispor(guides, genome="hg38"):
    """Export guides for CRISPOR scoring."""
    return [(guide['sequence'], guide['position']) for guide in guides]
```

## Scientific Background

### Signal Theory Rationale

The signal-theoretic approach is based on the hypothesis that DNA structural properties encoded in spectral features correlate with CRISPR efficiency and specificity:

1. **Spectral Entropy**: Reflects sequence complexity and accessibility
2. **Harmonic Content**: Indicates periodic patterns that may affect protein-DNA interactions  
3. **Spectral Stability**: Measures sensitivity to mutations and off-target binding

### Validation Studies

This methodology builds on several key observations:

- Chromatin accessibility correlates with spectral flatness
- Efficient guides often have lower spectral entropy
- Off-target sites with similar spectral signatures pose higher risk
- Repair outcomes correlate with local sequence stability measures

### Comparison to Traditional Methods

| Feature | Traditional | Signal-Theoretic | Advantage |
|---------|------------|------------------|-----------|
| Sequence features | Position-specific | Global patterns | Captures long-range effects |
| Efficiency prediction | Rule-based | Entropy-based | Adapts to sequence context |
| Off-target detection | Alignment-based | Spectral similarity | Detects subtle structural differences |
| Repair prediction | Context rules | Stability analysis | Quantitative probability estimates |

## References

1. Doench, J.G. *et al.* Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. *Nat Biotechnol* **34**, 184–191 (2016).

2. Hsu, P.D. *et al.* DNA targeting specificity of RNA-guided Cas9 nucleases. *Nat Biotechnol* **31**, 827–832 (2013).

3. Concordet, J.P. & Haeussler, M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing. *Nucleic Acids Res* **46**, W242–W245 (2018).

4. Listgarten, J. *et al.* Prediction of off-target activities for the end-to-end design of CRISPR guide RNAs. *Nat Biomed Eng* **2**, 38–47 (2018).

## License

MIT License. See LICENSE file for details.

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Check the documentation for common solutions
- Refer to the examples in the applications/ directory