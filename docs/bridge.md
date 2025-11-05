# Bridge API Documentation

## Overview

The bridge module provides implementations of geometric resolution (θ′) and curvature (κ) features for CRISPR guide scoring. These features emerge from the geodesic-topological bridge connecting discrete sequence space to topological structures.

---

## Quick Start

```python
from wave_crispr_signal.features import theta_prime, kappa, ThetaPrimeConfig, KappaConfig
from wave_crispr_signal.features.curvature import compute_coupled_features
import numpy as np

# Calculate phase resolution for position 21 with k=0.3
result = theta_prime(21, k=0.3)
print(f"θ′(21, 0.3) = {float(result):.6f}")

# Calculate curvature for position 21
kappa_val = kappa(21, mode='discrete')
print(f"κ(21) = {float(kappa_val):.6f}")

# Compute coupled features (multiplicative)
coupled = compute_coupled_features(21, k=0.3, coupling='multiplicative')
print(f"κ(21) · θ′(21, 0.3) = {float(coupled):.6f}")

# Process array of positions
positions = np.arange(1, 22)
theta_vals = theta_prime(positions, k=0.3)
kappa_vals = kappa(positions, mode='discrete')
```

---

## API Reference

### `theta_prime(n, config=None, k=None)`

Calculate geometric resolution function θ′(n,k) = φ·((n mod φ)/φ)^k

**Parameters:**
- `n` (int, array, mpf): Position index (must be positive)
- `config` (ThetaPrimeConfig, optional): Configuration object
- `k` (float, optional): Override k parameter (default: 0.3)

**Returns:**
- mpf (scalar) or ndarray of mpf (array): Geometric resolution values

**Raises:**
- `ValueError`: If n ≤ 0 or invalid parameters

**Example:**
```python
# Scalar with default k=0.3
result = theta_prime(21)

# Array with custom k
positions = np.arange(1, 22)
results = theta_prime(positions, k=0.35)

# With configuration
config = ThetaPrimeConfig(k=0.28, precision_dps=30)
result = theta_prime(21, config=config)
```

---

### `ThetaPrimeConfig`

Configuration for theta_prime calculation.

**Fields:**
- `k` (float): Geometric resolution exponent, range [0, 1], default 0.3
- `phi` (float, optional): Geometric period, defaults to golden ratio φ
- `precision_dps` (int): Decimal precision for mpmath, range [15, 100], default 50

**Example:**
```python
from wave_crispr_signal.features import ThetaPrimeConfig

config = ThetaPrimeConfig(
    k=0.30,              # Optimal for CRISPR
    phi=None,            # Use golden ratio (default)
    precision_dps=50     # High precision
)
```

---

### `kappa(n, config=None, mode=None)`

Calculate geodesic curvature κ(n) in discrete or continuous mode.

**Modes:**
- **Discrete** (default, biological): κ(n) = (1/φ) / (1 + (n mod φ))
- **Continuous** (analytical): κ(n) = (1/φ) / (1 + n/φ)

**Parameters:**
- `n` (int, array, mpf): Position index (must be positive)
- `config` (KappaConfig, optional): Configuration object
- `mode` (str, optional): 'discrete' or 'continuous' (overrides config)

**Returns:**
- mpf (scalar) or ndarray of mpf (array): Curvature values

**Raises:**
- `ValueError`: If n ≤ 0 or invalid parameters

**Example:**
```python
# Discrete mode (default)
kappa_discrete = kappa(21, mode='discrete')

# Continuous mode
kappa_continuous = kappa(21, mode='continuous')

# Array processing
positions = np.arange(1, 22)
kappa_vals = kappa(positions, mode='discrete')

# With configuration
config = KappaConfig(mode='discrete', scale=0.5)
result = kappa(21, config=config)
```

---

### `KappaConfig`

Configuration for kappa calculation.

**Fields:**
- `mode` (str): Calculation mode, 'discrete' or 'continuous', default 'discrete'
- `precision_dps` (int): Decimal precision for mpmath, range [15, 100], default 50
- `scale` (float, optional): Scaling factor, defaults to 1/φ

**Example:**
```python
from wave_crispr_signal.features import KappaConfig

config = KappaConfig(
    mode='discrete',     # For biological sequences
    scale=None,          # Use default 1/φ
    precision_dps=50
)
```

---

### `compute_coupled_features(n, k=0.3, coupling='multiplicative', kappa_config=None)`

Compute coupled κ(n) + θ′(n,k) features with specified coupling mode.

**Coupling Modes:**
- **Multiplicative** (recommended): F(n) = κ(n) · θ′(n,k)
- **Additive**: F(n) = κ(n) + θ′(n,k)

**Parameters:**
- `n` (int, array): Position index (must be positive)
- `k` (float): Phase resolution parameter, default 0.3
- `coupling` (str): 'additive' or 'multiplicative', default 'multiplicative'
- `kappa_config` (KappaConfig, optional): Configuration for kappa

**Returns:**
- mpf (scalar) or ndarray of mpf (array): Coupled feature values

**Example:**
```python
from wave_crispr_signal.features.curvature import compute_coupled_features

# Multiplicative coupling (default, best performance)
coupled = compute_coupled_features(21, k=0.3, coupling='multiplicative')

# Additive coupling (for ablation studies)
coupled_add = compute_coupled_features(21, k=0.3, coupling='additive')

# Array processing
positions = np.arange(1, 22)
coupled_vals = compute_coupled_features(positions, k=0.3)
```

---

## Vectorized Functions (Performance)

For large-scale calculations (>10,000 positions), use the vectorized versions:

### `theta_prime_vectorized(n, k=0.3, phi=None)`

Fast numpy-based version (float64 precision).

**Warning:** Uses float64 instead of mpmath. Suitable for production pipelines but not for scientific validation requiring high precision.

### `kappa_vectorized(n, mode='discrete', scale=None)`

Fast numpy-based version (float64 precision).

**Warning:** Uses float64 instead of mpmath. Suitable for production pipelines but not for scientific validation requiring high precision.

---

## Configuration Keys

### Default Parameters

```yaml
# Default configuration
theta_prime:
  k: 0.3                # Optimal for CRISPR (k*)
  phi: 1.618033988749895  # Golden ratio
  precision_dps: 50

kappa:
  mode: discrete        # Biological sequences
  scale: 0.618033988749895  # 1/φ
  precision_dps: 50

coupling:
  mode: multiplicative  # Consistent with derivation
```

### Config File Format (YAML)

Create `configs/bridge.yaml`:

```yaml
# Bridge Configuration
experiment:
  name: "bridge_validation"
  seed: [0, 1, 2, 3, 4]
  
features:
  theta_prime:
    k_range: [0.1, 0.6]
    k_step: 0.02
    k_star: 0.3
    precision_dps: 50
  
  kappa:
    mode: discrete
    scale: null  # Use default 1/φ
  
  coupling:
    mode: multiplicative
    ablations: [additive, multiplicative]

validation:
  bootstrap_resamples: 5000
  permutation_tests: 1000
  confidence_level: 0.95
  
datasets:
  - name: doench_2016
    path: data/doench_2016.csv
    sha256: "expected_hash_here"
  
  - name: kim_2025
    path: data/kim_2025.csv
    sha256: "expected_hash_here"

output:
  run_dir: runs/
  metrics_file: metrics.json
  figures_dir: figures/
```

---

## Exact Commands

### Run Stability Study

```bash
# Full validation on all datasets
python -m wave_crispr_signal.experiments.bridge \
  --config configs/bridge.yaml \
  --mode stability \
  --seeds 0 1 2 3 4

# Quick smoke test (micro shard, 1 seed)
python -m wave_crispr_signal.experiments.bridge \
  --config configs/bridge.yaml \
  --mode stability \
  --quick
```

### Run Benchmark vs RuleSet3

```bash
python -m wave_crispr_signal.experiments.bridge \
  --config configs/bridge.yaml \
  --mode benchmark \
  --baseline ruleset3 \
  --seeds 0 1 2 3 4
```

### Run Ablation Study

```bash
python -m wave_crispr_signal.experiments.bridge \
  --config configs/bridge.yaml \
  --mode ablation \
  --ablate coupling
```

---

## Interpreting k* Deviations

### Normal Operation

k* should be within [0.28, 0.32] for most CRISPR datasets:
- **k* = 0.300 ± 0.006**: Nominal optimal value
- **CI overlaps 0.30**: Dataset consistent with theory
- **15% variance reduction**: Expected density enhancement

### Deviation Diagnostics

If k* is outside [0.25, 0.35]:

1. **Check dataset quality**:
   - Look for outliers or data entry errors
   - Verify guide lengths are 20-23 nt
   - Check for technical artifacts

2. **Control for confounders**:
   - GC% distribution anomalies
   - Position bias (5' vs 3' end effects)
   - PAM sequence variations

3. **Map to κ(n) anomaly**:
   - Plot κ(n) and θ′(n,k*) vs position
   - Look for phase shifts or amplitude changes
   - Check if deviation corresponds to known biological feature

4. **Test arcsin bridge**:
   - Use bidirectional oracle notebook
   - Map k* deviation to f(x) domain
   - Check if it corresponds to critical point proximity

### Failure Modes

**Numeric underflow:**
- Symptom: κ or θ′ returns 0 or NaN
- Remedy: Increase precision_dps or use mpmath throughout

**Data leakage:**
- Symptom: k* estimate unstable across splits
- Remedy: Use split_by_guide or split_by_gene, never split_by_sample

**Overfitting:**
- Symptom: k* varies >0.05 between train/val
- Remedy: Increase dataset size or use regularization

---

## Integration with Existing Pipeline

### Add to Existing Feature Extraction

```python
from wave_crispr_signal.features import theta_prime, kappa

def extract_spectral_features(guide_sequence, position):
    """Extract all spectral features including bridge features."""
    # Existing features (FFT, etc.)
    fft_features = compute_fft_features(guide_sequence)
    
    # Bridge features
    theta = float(theta_prime(position, k=0.3))
    kappa_val = float(kappa(position, mode='discrete'))
    coupled = theta * kappa_val  # Multiplicative coupling
    
    return {
        **fft_features,
        'theta_prime': theta,
        'kappa': kappa_val,
        'coupled': coupled
    }
```

### Add to Machine Learning Pipeline

```python
from sklearn.ensemble import RandomForestRegressor
from wave_crispr_signal.features import compute_coupled_features
import numpy as np

# Extract features for all guides
def featurize_guides(guides_df):
    positions = np.arange(1, len(guides_df) + 1)
    
    # Compute bridge features
    coupled_features = compute_coupled_features(
        positions, 
        k=0.3, 
        coupling='multiplicative'
    )
    
    # Convert to float for sklearn
    guides_df['coupled_feature'] = [float(f) for f in coupled_features]
    
    return guides_df

# Train model
features = featurize_guides(train_df)
X = features[['gc_content', 'coupled_feature', ...]]
y = features['efficiency']

model = RandomForestRegressor(random_state=42)
model.fit(X, y)
```

---

## Reproducibility Checklist

Before running experiments:

- [ ] Set PYTHONHASHSEED=0
- [ ] Use fixed seeds: [0, 1, 2, 3, 4]
- [ ] Verify exact package versions match requirements.txt
- [ ] Compute and verify dataset SHA-256 checksums
- [ ] Document git commit hash
- [ ] Ensure no guide/gene leakage in splits
- [ ] Set mpmath precision to 50 (or document if different)
- [ ] Save config.yaml to run directory
- [ ] Log pip freeze to run directory

---

## References

- **Mathematical derivation**: `docs/bridge_kappa_to_theta_prime.md`
- **Z Framework**: `docs/Z_FRAMEWORK.md`
- **Topological analysis**: `docs/TOPOLOGICAL_ANALYSIS.md`
- **Repository policy**: `.github/REPOSITORY_POLICY.md`

---

**API Version**: 1.0  
**Last Updated**: 2025-11-05  
**Compatibility**: Python 3.12+
