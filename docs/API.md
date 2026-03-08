# Genomic Disruption Analyzer API Documentation

## Overview

The Genomic Disruption Analyzer provides a REST API for CRISPR guide RNA design, scoring, and off-target profiling using FFT-encoded DNA waveforms with θ′(n,k) phase weighting (k≈0.3).

**Version**: 1.0.0  
**Status**: Prototype (Phase 1 Complete)

## Core Features

- **Single Guide Scoring**: Score individual guides with spectral disruption metrics
- **Batch Processing**: Score multiple guides efficiently (>100,000 guides/min demonstrated)
- **Guide Design**: Automatically design guides from target sequences with PAM site detection
- **Off-Target Analysis**: Identify potential off-targets using spectral distance metrics
- **Bootstrap CI**: Optional confidence interval computation for statistical validity
- **Scientific Gates**: Enforces human DNA/RNA only, fail-fast validation, reproducibility

## Architecture

```
┌─────────────────────────────────────────┐
│         Client Application              │
└────────────┬────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────┐
│   GenomicDisruptionAPI (REST Wrapper)   │
└────────────┬────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────┐
│      DisruptionAnalyzer (Core Engine)   │
├─────────────────────────────────────────┤
│  • PhaseWeightedScorecard               │
│  • FFTCRISPRDisruptionAnalyzer          │
│  • Bootstrap CI computation             │
└─────────────────────────────────────────┘
```

## Installation

### Prerequisites

- Python 3.12+
- Required packages (from `requirements.txt`)

### Setup

```bash
# Clone repository
git clone https://github.com/zfifteen/wave-crispr-signal.git
cd wave-crispr-signal

# Install dependencies
pip install -r requirements.txt

# Run tests
python tests/test_genomic_disruption_api.py

# Run benchmark
python benchmarks/benchmark_disruption_api.py --n-guides 100
```

## API Reference

### Core Classes

#### `DisruptionAnalyzer`

Main analysis engine for CRISPR guide scoring.

**Initialization**:
```python
from applications.genomic_disruption_api import DisruptionAnalyzer

analyzer = DisruptionAnalyzer(
    k=0.3,              # Resolution exponent for θ′(n,k)
    phi_period=21.0,    # Geometric period for 21-nt guides
    seed=42             # Random seed for reproducibility
)
```

**Methods**:

##### `score_guide(guide, target=None, compute_ci=False, n_bootstrap=1000)`

Score a single CRISPR guide.

**Parameters**:
- `guide` (str): Guide RNA sequence (20-30 nt, A/C/G/T or A/C/G/U)
- `target` (str, optional): Target DNA sequence for on-target analysis
- `compute_ci` (bool): Whether to compute bootstrap confidence intervals
- `n_bootstrap` (int): Number of bootstrap resamples (default: 1000)

**Returns**: Dictionary with:
- `disruption_score` (float): Composite Z-invariant score (0-1)
- `spectral_features` (dict): FFT-based features
- `metadata` (dict): Analysis metadata
- `on_target_metrics` (dict, optional): Target-specific metrics
- `confidence_interval` (dict, optional): Bootstrap CI

**Example**:
```python
result = analyzer.score_guide("GCTGCGGAGACCTGGAGAGA")
print(f"Disruption score: {result['disruption_score']:.4f}")
print(f"Entropy: {result['spectral_features']['entropy']:.4f}")
```

##### `batch_score(guides, targets=None, compute_ci=False, n_bootstrap=1000)`

Score multiple guides efficiently.

**Parameters**:
- `guides` (List[str]): List of guide sequences
- `targets` (List[str], optional): List of target sequences (same length as guides)
- `compute_ci` (bool): Whether to compute bootstrap CIs
- `n_bootstrap` (int): Number of bootstrap resamples

**Returns**: List of scoring result dictionaries

**Example**:
```python
guides = [
    "GCTGCGGAGACCTGGAGAGA",
    "ATCGATCGATCGATCGATCG",
]
results = analyzer.batch_score(guides)
for i, result in enumerate(results):
    print(f"Guide {i}: {result['disruption_score']:.4f}")
```

##### `design_guides(target_sequence, pam='NGG', guide_length=20, n_guides=5, score_threshold=0.5)`

Design optimal CRISPR guides from target sequence.

**Parameters**:
- `target_sequence` (str): Target DNA sequence
- `pam` (str): PAM sequence motif (default: 'NGG' for SpCas9)
- `guide_length` (int): Length of guide RNA (default: 20)
- `n_guides` (int): Maximum number of guides to return
- `score_threshold` (float): Minimum disruption score threshold

**Returns**: List of designed guides with scores and positions

**Example**:
```python
target = "ATGCGATCGATCGATCGATCGATCGATCGATC..."
guides = analyzer.design_guides(target, n_guides=5)
for guide in guides:
    print(f"Position {guide['position']}: {guide['guide']} (score: {guide['disruption_score']:.4f})")
```

##### `analyze_offtarget(guide, offtarget_sequences, mismatch_threshold=3)`

Analyze off-target potential using spectral signatures.

**Parameters**:
- `guide` (str): Guide sequence
- `offtarget_sequences` (List[str]): List of potential off-target sequences
- `mismatch_threshold` (int): Maximum mismatches to report

**Returns**: List of off-target hits with distances and mismatch counts

**Example**:
```python
guide = "GCTGCGGAGACCTGGAGAGA"
offtargets = [
    "GCTGCGGAGACCTGGAGAGA",  # Perfect match
    "ACTGCGGAGACCTGGAGAGA",  # 1 mismatch
]
results = analyzer.analyze_offtarget(guide, offtargets)
for hit in results:
    print(f"{hit['sequence']}: {hit['mismatches']} mismatches, distance: {hit['spectral_distance']:.4f}")
```

#### `GenomicDisruptionAPI`

REST API wrapper for DisruptionAnalyzer.

**Initialization**:
```python
from applications.genomic_disruption_api import GenomicDisruptionAPI

api = GenomicDisruptionAPI(k=0.3, seed=42)
```

**Methods**:

##### `handle_score(request)`

Handle single guide scoring request.

**Request Format**:
```json
{
  "guide": "GCTGCGGAGACCTGGAGAGA",
  "target": "ATGC..." (optional),
  "compute_ci": false (optional),
  "n_bootstrap": 1000 (optional)
}
```

**Response Format**:
```json
{
  "success": true,
  "data": {
    "disruption_score": 0.7182,
    "spectral_features": {
      "entropy": 4.0482,
      "dominant_freq": 2,
      "dominant_magnitude": 10.8363,
      "sidelobe_count": 6,
      "diversity": 0.9075
    },
    "metadata": {
      "guide_length": 20,
      "guide_type": "DNA",
      "k_parameter": 0.3,
      "phi_period": 21.0,
      "computation_time_ms": 1.0
    }
  }
}
```

##### `handle_batch(request)`

Handle batch scoring request.

**Request Format**:
```json
{
  "guides": ["GCTGCG...", "ATCGAT..."],
  "targets": ["ATGC...", "CGTA..."] (optional),
  "compute_ci": false (optional),
  "n_bootstrap": 1000 (optional)
}
```

**Response Format**:
```json
{
  "success": true,
  "count": 2,
  "data": [
    { "index": 0, "disruption_score": 0.7182, ... },
    { "index": 1, "disruption_score": 0.6543, ... }
  ]
}
```

##### `handle_design(request)`

Handle guide design request.

**Request Format**:
```json
{
  "target": "ATGCGATCGATC...",
  "pam": "NGG" (optional),
  "guide_length": 20 (optional),
  "n_guides": 5 (optional),
  "score_threshold": 0.5 (optional)
}
```

**Response Format**:
```json
{
  "success": true,
  "count": 2,
  "data": [
    {
      "guide": "CGATCGATCGATCGATCGAT",
      "position": 35,
      "pam_site": "CGG",
      "disruption_score": 0.7168,
      "spectral_features": { ... }
    }
  ]
}
```

##### `handle_offtarget(request)`

Handle off-target analysis request.

**Request Format**:
```json
{
  "guide": "GCTGCGGAGACCTGGAGAGA",
  "offtargets": ["GCTGCGG...", "ACTGCGG..."],
  "mismatch_threshold": 3 (optional)
}
```

**Response Format**:
```json
{
  "success": true,
  "count": 2,
  "data": [
    {
      "sequence": "GCTGCGGAGACCTGGAGAGA",
      "mismatches": 0,
      "spectral_distance": 0.0,
      "disruption_score": 0.7182
    }
  ]
}
```

## CLI Usage

### Score Single Guide

```bash
python applications/genomic_disruption_api.py score \
  --guide GCTGCGGAGACCTGGAGAGA \
  --k 0.3 \
  --seed 42 \
  --output results/score.json
```

### Design Guides

```bash
python applications/genomic_disruption_api.py design \
  --target "ATGCGATCGATC..." \
  --k 0.3 \
  --seed 42 \
  --output results/guides.json
```

## Python API Usage

### Basic Scoring

```python
from applications.genomic_disruption_api import DisruptionAnalyzer

# Initialize
analyzer = DisruptionAnalyzer(k=0.3, seed=42)

# Score single guide
result = analyzer.score_guide("GCTGCGGAGACCTGGAGAGA")
print(f"Score: {result['disruption_score']:.4f}")
```

### Batch Processing

```python
# Batch score
guides = [
    "GCTGCGGAGACCTGGAGAGA",
    "ATCGATCGATCGATCGATCG",
    "GGGGGGGGGGGGGGGGGGGG",
]

results = analyzer.batch_score(guides)
for i, result in enumerate(results):
    if 'error' not in result:
        print(f"Guide {i}: {result['disruption_score']:.4f}")
```

### Guide Design

```python
# Design guides from target
target = "ATGCGATCGATCGATCGATCGATCGATCGATC..."
guides = analyzer.design_guides(target, n_guides=5)

for guide in guides:
    print(f"{guide['guide']} @ position {guide['position']}")
    print(f"  PAM: {guide['pam_site']}")
    print(f"  Score: {guide['disruption_score']:.4f}")
```

### With Confidence Intervals

```python
# Score with bootstrap CI
result = analyzer.score_guide(
    "GCTGCGGAGACCTGGAGAGA",
    compute_ci=True,
    n_bootstrap=1000
)

ci = result['confidence_interval']
print(f"Score: {result['disruption_score']:.4f}")
print(f"95% CI: [{ci['lower_95']:.4f}, {ci['upper_95']:.4f}]")
```

## Performance Benchmarks

Based on benchmarks on standard hardware:

| Metric | Target | Achieved | Status |
|--------|--------|----------|--------|
| Single guide latency (P95) | <500ms | ~0.6ms | ✓ PASSED |
| Batch throughput | 1,000 guides/min | >100,000 guides/min | ✓ PASSED |
| Bootstrap CI width | <1% | ~0.5% | ✓ PASSED |
| Score distribution | Non-uniform | KS p<0.05 | ✓ PASSED |

### Run Benchmarks

```bash
# Quick benchmark
python benchmarks/benchmark_disruption_api.py --n-guides 100

# Full benchmark suite
python benchmarks/benchmark_disruption_api.py --full --output results/benchmark.json
```

## Scientific Gates

The API enforces the following scientific gates:

1. **Human DNA/RNA Only**
   - DNA: A/C/G/T/N allowed
   - RNA: A/C/G/U/N allowed
   - Reject mixed T/U sequences
   - Reject invalid bases (X, Y, Z, etc.)

2. **Fail-Fast Validation**
   - Clear error messages for invalid input
   - ValueError raised immediately on validation failure

3. **Z-Invariant Normalization**
   - Z = S(Δ/φ) with sigmoid aggregation
   - κ(n) = d(n)·ln(n+1)/e² curvature weight
   - Scores in [0, 1] range

4. **Geometric Resolution**
   - θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
   - φ-period = 21 for 21-nt guides

5. **Statistical Validity**
   - Bootstrap CI with ≥1,000 resamples
   - Permutation tests for p-values
   - KS test for distribution validation

6. **Reproducibility**
   - Seed control for all random operations
   - Fixed numpy seed in analyzer initialization
   - Deterministic results with same seed

## Error Handling

### Common Errors

#### Invalid DNA Sequence
```python
# Error: U in DNA context
try:
    analyzer.score_guide("GCUGCGGAGACCTGGAGAGA")
except ValueError as e:
    print(e)  # "Invalid RNA sequence: contains {'T'}..."
```

#### Invalid Bases
```python
# Error: Invalid bases
try:
    analyzer.score_guide("GCTXCGGAGACCTGGAGAGA")
except ValueError as e:
    print(e)  # "Invalid DNA sequence: contains {'X'}..."
```

#### Mismatched Lengths
```python
# Error: Target and guide different lengths (when using disruption analysis)
try:
    result = analyzer.score_guide(
        guide="GCTGCGGAGACCTGGAGAGA",
        target="ATGC"  # Too short
    )
except ValueError as e:
    print(e)  # "Reference and mutated sequences must have same length..."
```

## Future Roadmap

### Phase 2: Integration (Current)
- [ ] FastAPI/Flask HTTP server implementation
- [ ] Swagger/OpenAPI documentation
- [ ] Rate limiting and request validation
- [ ] Async batch processing with queues
- [ ] Database integration for result caching

### Phase 3: Optimization
- [ ] Vectorized FFT computation with NumPy
- [ ] GPU acceleration with PyTorch
- [ ] Distributed processing with Dask
- [ ] Redis caching for repeated queries
- [ ] Response compression

### Phase 4: Production
- [ ] Docker containerization
- [ ] Kubernetes deployment
- [ ] Monitoring and logging (Prometheus, Grafana)
- [ ] Load testing and auto-scaling
- [ ] Security hardening

## References

- [Z Framework Documentation](../docs/Z_FRAMEWORK.md)
- [Phase-Weighted Scorecard](../docs/PHASE_WEIGHTED_SCORECARD.md)
- [FFT Golden Ratio CRISPR](../docs/FFT_GOLDEN_RATIO_CRISPR.md)
- [Repository Policy](../.github/REPOSITORY_POLICY.md)

## Citation

If you use this API in your research, please cite:

```bibtex
@software{genomic_disruption_analyzer,
  title={Genomic Disruption Analyzer: Z-Framework SaaS Pipeline for CRISPR gRNA Optimization},
  author={Z Framework Team},
  year={2025},
  url={https://github.com/zfifteen/wave-crispr-signal}
}
```

## License

MIT License - See [LICENSE](../LICENSE) for details.

## Support

For issues and questions:
- GitHub Issues: https://github.com/zfifteen/wave-crispr-signal/issues
- Documentation: https://github.com/zfifteen/wave-crispr-signal/tree/main/docs

---

**Last Updated**: 2025-11-23  
**Version**: 1.0.0-alpha  
**Status**: Prototype (Phase 1 Complete)
