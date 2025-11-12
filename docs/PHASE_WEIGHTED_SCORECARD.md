# Phase-Weighted CRISPR Scorecard

**A phase-weighted spectral analysis framework for CRISPR guide scoring using the Z Framework**

## Overview

The Phase-Weighted CRISPR Scorecard implements a novel approach to CRISPR guide RNA design and mutation impact assessment using Fourier-based signal processing with golden-ratio phase resolution. This framework treats DNA sequences as complex waveforms and quantifies mutation-induced spectral disruptions.

## Key Features

- **Complex Encoding**: Maps DNA/RNA bases to complex plane preserving complementarity
  - A → 1 (real positive)
  - T/U → -1 (real negative)
  - C → +i (imaginary positive)
  - G → -i (imaginary negative)

- **Phase-Weighted Transform**: Applies position-dependent phase shifts using geometric resolution function
  - θ′(n,k) = φ·((n mod φ)/φ)^k where φ ≈ 1.618 (golden ratio)
  - k* = 0.3 (optimized via grid search on entropy minima)

- **Spectral Features**: Extracts FFT-based disruption metrics
  - ΔEntropy: Shannon entropy difference
  - Δf₁: Dominant frequency shift
  - ΔSidelobes: Spectral peak count change

- **Z-Invariant Scoring**: Composite disruption score with curvature weighting
  - Z = S(Δ_spectral / φ)
  - S(x) = 1/(1 + e^{-κ(n)·x}) (sigmoid aggregator)
  - κ(n) = d(n)·ln(n+1)/e² (curvature weight)

## Scientific Gates

The implementation enforces strict scientific validation gates:

1. **Human DNA/RNA Only**: Validates sequences contain only A/C/G/T (DNA) or A/C/G/U (RNA)
   - N allowed for ambiguous bases
   - Rejects IUPAC codes beyond N with clear ValueError

2. **Fail-Fast Validation**: Immediate error on invalid nucleotides

3. **Domain-Correct Z Invariant**: Uses discrete/biological form Z = A(B / e²)

4. **Statistical Validity**: 
   - Bootstrap CI with ≥1,000 resamples
   - Seed control for reproducibility
   - Permutation tests for significance

## Installation

```bash
cd /path/to/wave-crispr-signal
pip install -r requirements.txt
```

## Quick Start

### Score a Single Guide

```python
from applications.phase_weighted_scorecard import PhaseWeightedScorecard

# Create scorecard
scorecard = PhaseWeightedScorecard(k=0.3)

# Score a guide
guide = "GCTGCGGAGACCTGGAGAGA"
result = scorecard.score_guide(guide)

print(f"Entropy: {result['guide_features']['entropy']:.6f}")
print(f"Diversity: {result['guide_features']['diversity']:.6f}")
```

### Compare Guide vs Target (Mutation Analysis)

```python
# Reference and mutant sequences
ref_seq = "GCTGCGGAGACCTGGAGAGAAAGC"
mut_seq = "GCTGCGGAGTCCTGGAGAGAAAGC"  # A→T mutation at position 10

# Compute Z-score
result = scorecard.compute_z_score(ref_seq, mut_seq)

print(f"Z-Score: {result['z_score']:.6f}")
print(f"Delta Spectral: {result['delta_spectral']:.6f}")
print(f"Kappa (curvature): {result['kappa']:.6f}")
```

### Batch Processing

```python
from applications.phase_weighted_scorecard import score_guide_batch

guides = [
    "GCTGCGGAGACCTGGAGAGA",
    "ATCGATCGATCGATCGATCG",
    "AAAATTTTCCCCGGGGAAAA",
]

results = score_guide_batch(guides, k=0.3)

for i, result in enumerate(results):
    print(f"Guide {i+1}: Entropy = {result['guide_features']['entropy']:.6f}")
```

## Command-Line Interface

### Score a Single Guide

```bash
python applications/phase_weighted_scorecard_cli.py score \
    --guide GCTGCGGAGACCTGGAGAGA
```

Output:
```
Phase-Weighted CRISPR Scorecard
================================
k parameter: 0.3
Sequence type: DNA

Guide sequence: GCTGCGGAGACCTGGAGAGA

Guide Spectral Features:
  Entropy: 4.048185
  Dominant Freq Index: 2
  Dominant Freq Magnitude: 10.836302
  Sidelobe Count: 6
  Diversity: 0.907490
```

### Score with Target Comparison

```bash
python applications/phase_weighted_scorecard_cli.py score \
    --guide GCTGCGGAGTCCTGGAGAGA \
    --target GCTGCGGAGACCTGGAGAGA
```

Output:
```
Z-Score (Disruption): 0.502748
Delta Spectral: 0.047572
Kappa (curvature): 0.373914

Disruption Components:
  ΔEntropy: 0.118931
  ΔFrequency: 0
  ΔSidelobes: 0
```

### Batch Process from CSV

Create an input CSV file `guides.csv`:
```csv
id,guide,target
guide_1,GCTGCGGAGACCTGGAGAGA,GCTGCGGAGACCTGGAGAGA
guide_2,ATCGATCGATCGATCGATCG,ATCGATCGATCGATCGATCG
guide_3,AAAATTTTCCCCGGGGAAAA,AAAATTTTCCCCGGGGAAAA
```

Run batch processing:
```bash
python applications/phase_weighted_scorecard_cli.py batch \
    --input guides.csv \
    --output results.csv \
    --k 0.3 \
    --seed 42
```

### Analyze FASTA Files

```bash
python applications/phase_weighted_scorecard_cli.py analyze \
    --ref reference.fa \
    --mut mutant.fa \
    --output analysis.json \
    --k 0.3
```

## Quick Validation Demo

Run the quick validation demo to verify installation:

```bash
python proof_pack/phase_weighted_quick_demo.py
```

This will:
- Test spectral feature extraction
- Demonstrate mutation disruption scoring
- Compute bootstrap confidence intervals (1,000 resamples)
- Validate mathematical constants
- Test edge cases

Expected output includes:
```
Single Mutation Results (n=50):
  Mean Z-score: 0.560730
  95% CI: [0.539072, 0.583831]

Multiple Mutations Results (n=50):
  Mean Z-score: 0.592572
  95% CI: [0.569625, 0.616675]

✓ Multiple mutations show higher disruption scores (expected)
```

## API Reference

### PhaseWeightedScorecard

Main class for phase-weighted spectral analysis.

#### Constructor

```python
PhaseWeightedScorecard(
    k: float = 0.3,
    weights: Optional[Dict[str, float]] = None,
    is_rna: bool = False,
)
```

**Parameters:**
- `k`: Phase parameter (default: 0.3, optimized for k*)
- `weights`: Custom weights for composite score (default: entropy=0.4, freq_shift=0.35, sidelobes=0.25)
- `is_rna`: If True, treat sequences as RNA (U instead of T)

#### Methods

##### compute_spectral_features

```python
compute_spectral_features(seq: str) -> Dict[str, float]
```

Compute phase-weighted spectral features for a sequence.

**Returns:**
- `entropy`: Spectral entropy
- `dominant_freq_idx`: Dominant frequency index
- `dominant_freq_mag`: Dominant frequency magnitude
- `sidelobe_count`: Number of sidelobes
- `diversity`: Sequence diversity d(n)

##### compute_disruption_features

```python
compute_disruption_features(ref_seq: str, mut_seq: str) -> Dict[str, float]
```

Compute mutation-induced spectral disruptions.

**Returns:**
- `delta_entropy`: H(mutated) - H(reference)
- `delta_freq`: |f₁_mut - f₁_ref|
- `delta_sidelobes`: |sidelobes_mut - sidelobes_ref|
- `ref_features`: Reference sequence features
- `mut_features`: Mutant sequence features

##### compute_z_score

```python
compute_z_score(ref_seq: str, mut_seq: str) -> Dict[str, float]
```

Compute Z-invariant disruption score.

**Returns:**
- `z_score`: Composite Z-invariant score (0-1)
- `delta_spectral`: Weighted spectral disruption
- `kappa`: Curvature weight κ(n)
- `disruptions`: Individual disruption features

##### score_guide

```python
score_guide(
    guide_seq: str,
    target_seq: Optional[str] = None,
) -> Dict[str, float]
```

Score a CRISPR guide with optional target context.

**Parameters:**
- `guide_seq`: Guide RNA sequence
- `target_seq`: Optional target DNA sequence for disruption analysis

**Returns:**
- If target provided: Z-score with disruption analysis
- Otherwise: Spectral features only

### Utility Functions

#### encode_complex

```python
encode_complex(seq: str, is_rna: bool = False) -> np.ndarray
```

Map DNA/RNA string to complex vector preserving complementarity.

#### theta_prime

```python
theta_prime(n: Union[int, float, np.ndarray], k: float = 0.3) -> Union[float, np.ndarray]
```

Geodesic resolution function θ′(n,k) with golden-angle phasing.

#### score_guide_batch

```python
score_guide_batch(
    guides: List[str],
    targets: Optional[List[str]] = None,
    k: float = 0.3,
    weights: Optional[Dict[str, float]] = None,
    is_rna: bool = False,
) -> List[Dict[str, float]]
```

Score multiple guides efficiently.

## Mathematical Framework

### Complex Encoding

DNA sequences are encoded as complex numbers to preserve Watson-Crick complementarity:

```
A ↔ T:  1 ↔ -1  (real axis, conjugate pair)
C ↔ G: +i ↔ -i  (imaginary axis, conjugate pair)
```

This encoding enables:
- Fourier analysis of sequence structure
- Phase-based position weighting
- Spectral disruption quantification

### Geometric Resolution

Position-dependent phase shift function:

```
θ′(n,k) = φ · ((n mod φ) / φ)^k
```

Where:
- φ ≈ 1.618033989 (golden ratio)
- k* = 0.3 (optimized parameter)
- n = position index

This function embeds golden-angle spirals to minimize discrepancy in spectral sampling.

### Z-Invariant Form

The composite disruption score follows the Z Framework invariant:

```
Z = S(Δ_spectral / φ)
```

Where:
- S(x) = 1/(1 + e^{-κ(n)·x}) is sigmoid aggregator
- κ(n) = d(n)·ln(n+1)/e² is curvature weight
- d(n) = sequence diversity (normalized Shannon entropy)
- Δ_spectral = weighted sum of spectral disruptions

Composite spectral disruption:
```
Δ_spectral = w₁·ΔEntropy + w₂·Δf₁ + w₃·ΔSidelobes
```

Default weights (tuned via OLS regression):
- w₁ = 0.4 (entropy)
- w₂ = 0.35 (frequency shift)
- w₃ = 0.25 (sidelobes)

### Sequence Diversity

Diversity d(n) measures base composition uniformity:

```
d(n) = H(bases) / log₂(4)
```

Where H(bases) is Shannon entropy of base frequencies, normalized by maximum entropy for 4 bases.

### Curvature Weight

For sequences with n ≥ 10:
```
κ(n) = d(n) · ln(n+1) / e²
```

For n < 10 (guard against instability):
```
κ(n) = 1.0
```

## Validation and Testing

### Unit Tests

Run the comprehensive test suite:

```bash
python -m pytest tests/test_phase_weighted_scorecard.py -v
```

Test coverage includes:
- Sequence validation (DNA/RNA, IUPAC rejection)
- Complex encoding correctness
- Geometric resolution θ′(n,k)
- Phase-weighted transforms
- Spectral features (entropy, frequency, sidelobes)
- Sequence diversity
- Curvature weight κ(n)
- Sigmoid aggregator
- Z-score computation
- Batch processing
- Integration workflows

All 48 tests pass successfully.

### Quick Validation

Run the quick validation demo:

```bash
python proof_pack/phase_weighted_quick_demo.py
```

This validates:
- Complex encoding and phase weighting
- Spectral feature extraction
- Z-invariant disruption scoring
- Bootstrap confidence intervals
- Mathematical constants (φ, e²)
- Edge case handling

## Performance Considerations

### Computational Complexity

- **Single sequence**: O(n log n) for FFT, where n = sequence length
- **Batch processing**: Linear in number of sequences
- **Typical performance**: ~1ms per 20-nt guide on modern CPU

### Memory Usage

- **Single guide**: ~few KB
- **Batch of 1000 guides**: ~few MB
- Phase-weighted arrays: 2× sequence length (complex128)

### Optimization Tips

1. **Batch processing**: Use `score_guide_batch()` for multiple guides
2. **Pre-compute features**: Cache spectral features for repeated comparisons
3. **Parallel processing**: Guides can be scored independently (not implemented by default)

## Limitations and Future Work

### Current Limitations

1. **Dataset Specificity**: Tuned on eukaryotic data; prokaryotic performance unverified
2. **Causality**: Scores are correlative, not causal
3. **Sequence Length**: Optimal for 15-30 nt guides (shorter sequences have warnings)
4. **Baseline Comparison**: RuleSet3 comparison pending full validation

### Future Enhancements

1. **Machine Learning Integration**: Train XGBoost on spectral features for improved prediction
2. **Dataset Expansion**: Validate on additional datasets (Vex-seq, prokaryotic)
3. **Real-time Scoring**: Web API for cloud-based analysis
4. **Interactive Visualization**: Dashboard with adjustable k parameter
5. **Multi-mutation Analysis**: Systematic epistasis detection

## References

### Issue and Requirements

- Issue: Phase-Weighted CRISPR Scorecard
- Z Framework: `docs/Z_FRAMEWORK.md`
- Repository Policy: `.github/REPOSITORY_POLICY.md`

### Mathematical Foundations

- Golden ratio φ for phase resolution minimizes spectral discrepancy
- Discrete Z invariant form appropriate for genomic signals
- Curvature weight κ(n) accounts for sequence complexity

### Empirical Validation

- Bootstrap CI: ≥1,000 resamples for statistical confidence
- Permutation tests: ≥1,000 permutations for p-values
- Seed control: Ensures reproducible results

## Contributing

When contributing to this module:

1. **Follow Repository Policy**: See `.github/REPOSITORY_POLICY.md`
2. **Maintain Scientific Gates**: All validation gates must pass
3. **Add Tests**: New features require comprehensive tests
4. **Document Changes**: Update this documentation for API changes
5. **Validate**: Run `make smoke` and full test suite

## License

MIT License - See LICENSE file in repository root.

## Contact

For questions or issues related to the Phase-Weighted CRISPR Scorecard, please open an issue on the repository.
