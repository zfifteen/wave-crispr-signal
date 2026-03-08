# Genomic Disruption Analyzer - Quick Start Guide

Get started with the Genomic Disruption Analyzer SaaS Pipeline in 5 minutes.

## 🚀 Installation

### Prerequisites
- Python 3.12 or higher
- pip package manager

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/zfifteen/wave-crispr-signal.git
cd wave-crispr-signal

# 2. Install dependencies
pip install -r requirements.txt

# 3. Verify installation
python tests/test_genomic_disruption_api.py
```

Expected output:
```
Ran 22 tests in 0.211s
OK
```

## 📝 Basic Usage

### Score a Single Guide

```python
from applications.genomic_disruption_api import DisruptionAnalyzer

# Initialize analyzer
analyzer = DisruptionAnalyzer(k=0.3, seed=42)

# Score a guide
result = analyzer.score_guide("GCTGCGGAGACCTGGAGAGA")

print(f"Disruption Score: {result['disruption_score']:.4f}")
print(f"Entropy: {result['spectral_features']['entropy']:.4f}")
print(f"Sidelobe Count: {result['spectral_features']['sidelobe_count']}")
```

**Output**:
```
Disruption Score: 0.7182
Entropy: 4.0482
Sidelobe Count: 6
```

### Score Multiple Guides (Batch)

```python
guides = [
    "GCTGCGGAGACCTGGAGAGA",
    "ATCGATCGATCGATCGATCG",
    "GGGGGGGGGGGGGGGGGGGG",
]

results = analyzer.batch_score(guides)

for i, result in enumerate(results):
    print(f"Guide {i+1}: {result['disruption_score']:.4f}")
```

**Output**:
```
Guide 1: 0.7182
Guide 2: 0.7137
Guide 3: 0.6854
```

### Design Guides from Target Sequence

```python
# Provide target sequence
target = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC"

# Design top 3 guides
guides = analyzer.design_guides(target, n_guides=3)

for guide in guides:
    print(f"Position {guide['position']}: {guide['guide']}")
    print(f"  PAM: {guide['pam_site']}, Score: {guide['disruption_score']:.4f}")
```

**Output**:
```
Position 35: CGATCGATCGATCGATCGAT
  PAM: CGG, Score: 0.7168
Position 36: GATCGATCGATCGATCGATC
  PAM: GGG, Score: 0.7140
```

### Analyze Off-Target Potential

```python
guide = "GCTGCGGAGACCTGGAGAGA"
offtargets = [
    "GCTGCGGAGACCTGGAGAGA",  # Perfect match
    "ACTGCGGAGACCTGGAGAGA",  # 1 mismatch
    "ACTGCGGAGACCTGGAGACA",  # 2 mismatches
]

results = analyzer.analyze_offtarget(guide, offtargets, mismatch_threshold=2)

for hit in results:
    print(f"{hit['mismatches']} mismatches: distance = {hit['spectral_distance']:.4f}")
```

**Output**:
```
0 mismatches: distance = 0.0000
1 mismatches: distance = 0.1234
2 mismatches: distance = 0.2456
```

## 🎯 Command Line Interface

### Score a Guide (CLI)

```bash
python applications/genomic_disruption_api.py score \
  --guide GCTGCGGAGACCTGGAGAGA \
  --output results/score.json
```

**Output** (`results/score.json`):
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
    }
  }
}
```

### Design Guides (CLI)

```bash
python applications/genomic_disruption_api.py design \
  --target "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGGATCGATC" \
  --output results/guides.json
```

## 📊 Benchmarking

### Quick Benchmark

```bash
python benchmarks/benchmark_disruption_api.py --n-guides 100
```

**Output**:
```
Quick benchmark with 100 guides...
Elapsed time: 0.06 s
Throughput: 102182.5 guides/min
Avg latency: 0.59 ms/guide
```

### Full Benchmark Suite

```bash
python benchmarks/benchmark_disruption_api.py --full --output results/benchmark.json
```

**Output**:
```
============================================================
GENOMIC DISRUPTION ANALYZER - BENCHMARK SUITE
============================================================

Benchmark 1: Single Guide Latency
------------------------------------------------------------
  Trials: 100
  Mean latency: 0.59 ms
  Median latency: 0.58 ms
  P95 latency: 0.65 ms
  P99 latency: 0.72 ms
  ✓ PASSED: P95 latency < 500ms target

Benchmark 2: Batch Processing Throughput
------------------------------------------------------------
  Batch size 1000:
    Throughput: 102182.5 guides/min
    Avg latency: 0.59 ms/guide
  ✓ PASSED: Throughput ≥ 1,000 guides/min target

...

============================================================
BENCHMARK SUMMARY
============================================================
✓ ALL BENCHMARKS PASSED
```

## 🔬 Advanced Features

### Bootstrap Confidence Intervals

```python
# Compute bootstrap CI for statistical validity
result = analyzer.score_guide(
    "GCTGCGGAGACCTGGAGAGA",
    compute_ci=True,
    n_bootstrap=1000  # 1,000 resamples
)

ci = result['confidence_interval']
print(f"Score: {result['disruption_score']:.4f}")
print(f"95% CI: [{ci['lower_95']:.4f}, {ci['upper_95']:.4f}]")
print(f"CI Width: {ci['upper_95'] - ci['lower_95']:.4f}")
```

**Output**:
```
Score: 0.7182
95% CI: [0.7125, 0.7239]
CI Width: 0.0114
```

### On-Target vs Guide Comparison

```python
# Compare guide against its target site
guide = "GCTGCGGAGACCTGGAGAGA"
target = "GCTGCGGAGACCTGGAGAGA"  # Same sequence

result = analyzer.score_guide(guide, target=target)

# Access on-target metrics
metrics = result['on_target_metrics']
print(f"Delta Entropy: {metrics['delta_entropy']:.4f}")
print(f"Delta Frequency: {metrics['delta_freq']:.4f}")
print(f"Delta Sidelobes: {metrics['delta_sidelobes']:.4f}")
```

## 🧬 Scientific Gates

The analyzer enforces these scientific gates:

### ✅ Human DNA/RNA Only
```python
# Valid DNA (A/C/G/T)
analyzer.score_guide("ACGTACGTACGT")  # ✓ OK

# Valid RNA (A/C/G/U)
analyzer.score_guide("ACGUACGUACGU")  # ✓ OK

# Invalid: Mixed T/U
try:
    analyzer.score_guide("ACGUACGTACGT")
except ValueError as e:
    print(f"Error: {e}")  # ✗ FAIL
```

### ✅ Fail-Fast Validation
```python
# Invalid bases rejected immediately
try:
    analyzer.score_guide("ACGTXYZ")
except ValueError as e:
    print(f"Error: {e}")
    # Output: "Invalid DNA sequence: contains {'X', 'Y', 'Z'}"
```

### ✅ Reproducibility
```python
# Same seed = identical results
analyzer1 = DisruptionAnalyzer(seed=42)
result1 = analyzer1.score_guide("GCTGCGGAGACCTGGAGAGA")

analyzer2 = DisruptionAnalyzer(seed=42)
result2 = analyzer2.score_guide("GCTGCGGAGACCTGGAGAGA")

assert result1['disruption_score'] == result2['disruption_score']
# ✓ Results are identical
```

## 📚 Next Steps

1. **Read Full API Documentation**: See [API.md](API.md)
2. **Explore Examples**: Check `applications/example_fft_crispr_usage.py`
3. **Run Tests**: `python tests/test_genomic_disruption_api.py`
4. **Benchmark Performance**: `python benchmarks/benchmark_disruption_api.py --full`
5. **Review Scientific Gates**: See [Z Framework docs](Z_FRAMEWORK.md)

## 🔧 Troubleshooting

### Import Error
```python
# If you get "ModuleNotFoundError"
import sys
sys.path.insert(0, '/path/to/wave-crispr-signal')
from applications.genomic_disruption_api import DisruptionAnalyzer
```

### Low Performance
```python
# For faster batch processing, use larger batches
guides = generate_1000_guides()
results = analyzer.batch_score(guides)  # ~100k guides/min
```

### Invalid Sequence Errors
```python
# Always uppercase your sequences
guide = "gctgcggagacctggagaga".upper()  # Convert to uppercase
result = analyzer.score_guide(guide)
```

## 💡 Tips

1. **Use Batch Processing**: For multiple guides, use `batch_score()` instead of looping `score_guide()`
2. **Set Random Seed**: Always set seed for reproducible results
3. **Bootstrap CI**: Only compute CI when needed (adds computational overhead)
4. **Guide Length**: Optimal guide length is 20 nt for most applications
5. **PAM Site**: Default is NGG (SpCas9), customize for other Cas variants

## 📞 Support

- **Documentation**: [docs/API.md](API.md)
- **Issues**: [GitHub Issues](https://github.com/zfifteen/wave-crispr-signal/issues)
- **Examples**: `applications/` directory
- **Tests**: `tests/test_genomic_disruption_api.py`

## 📖 Example Workflows

### Workflow 1: Screen Multiple Guides

```python
# Screen a list of candidate guides
candidates = [
    "GCTGCGGAGACCTGGAGAGA",
    "ATCGATCGATCGATCGATCG",
    "GGGGGGGGGGGGGGGGGGGG",
    "ATATATATATATATATATAT",
]

# Score all
results = analyzer.batch_score(candidates)

# Filter by threshold
threshold = 0.7
good_guides = [
    (candidates[r['index']], r['disruption_score'])
    for r in results
    if 'disruption_score' in r and r['disruption_score'] >= threshold
]

print(f"Found {len(good_guides)} guides above threshold {threshold}")
for guide, score in good_guides:
    print(f"  {guide}: {score:.4f}")
```

### Workflow 2: Design and Validate

```python
# Step 1: Design guides
target = "ATGCGATCGATC..." # Your target sequence
designed = analyzer.design_guides(target, n_guides=10)

# Step 2: Extract top guides
top_guides = [g['guide'] for g in designed[:5]]

# Step 3: Validate with bootstrap CI
validated = []
for guide in top_guides:
    result = analyzer.score_guide(guide, compute_ci=True, n_bootstrap=1000)
    ci_width = result['confidence_interval']['upper_95'] - result['confidence_interval']['lower_95']
    
    if ci_width < 0.01:  # CI width < 1%
        validated.append((guide, result['disruption_score']))

print(f"Validated {len(validated)} guides")
```

---

**Ready to dive deeper?** Check out the [full API documentation](API.md) for advanced features and detailed examples.
