# Topological Analysis: Linking Geodesic Curvature to f(x) Properties

## Overview

This module implements the mathematical hypothesis connecting the topological properties of f(x) = arcsin((x-1)/(2x+3)) with the Z Framework's geodesic curvature mapping θ'(n, k) = φ · ((n mod φ)/φ)^k.

## Mathematical Foundation

### The f(x) Function

The function f(x) = arcsin((x-1)/(2x+3)) exhibits critical topological properties:

- **Pole Singularity**: x = -3/2 (where denominator 2x+3 = 0)
- **Branch Point**: x = -2/3 (critical for domain analysis)
- **Domain Constraints**: -1 ≤ (x-1)/(2x+3) ≤ 1 for valid arcsin argument

### Admissible Intervals

The function's domain creates "admissible universes":
- **Interval 1**: (-∞, -3/2) - approaching the pole from negative infinity
- **Interval 2**: (-4, -2/3] - bounded interval with specific endpoints

### Connection to Z Framework

The hypothesis establishes that:

1. **Singularities** in f(x) mirror "prime-like ruptures" in discrete Z space
2. **Domain constraints** align with the invariant c = e² ≈ 7.389
3. **Geodesic mapping** provides a mathematical bridge between spaces
4. **Optimal parameter** k* ≈ 0.3 emerges naturally in both systems

## Implementation

### TopologicalAnalyzer Class

```python
from topological_analysis import TopologicalAnalyzer

# Initialize with high precision
analyzer = TopologicalAnalyzer(precision_dps=50)

# Calculate f(x) at valid points
fx_value = analyzer.f_x(0.5)

# Analyze domain constraints
domain_info = analyzer.analyze_domain_constraints()

# Map to geodesic space
x_values = [0, 1, -0.6]
mapping = analyzer.map_fx_to_geodesic(x_values)

# Comprehensive analysis
results = analyzer.comprehensive_analysis()
```

### Key Methods

#### `f_x(x)`
Calculates f(x) = arcsin((x-1)/(2x+3)) with:
- Pole detection at x = -3/2
- Domain validation for arcsin argument
- High-precision arithmetic

#### `analyze_domain_constraints()`
Returns:
- Domain bounds where arg = ±1
- Pole locations
- Admissible intervals

#### `map_fx_to_geodesic(x_values, k=0.3)`
Establishes correspondence between:
- f(x) values in original space
- θ'(n, k) values in geodesic space
- Resonance measures

#### `analyze_invariant_alignment()`
Compares:
- Main interval span vs e² invariant
- Golden ratio φ relationships
- Topological resonance measures

#### `demonstrate_density_enhancement()`
Validates k* ≈ 0.3 optimality:
- Baseline variance (k=0)
- Enhanced variance (k=0.3) 
- Achievement vs 15% target

## Usage Examples

### Basic Function Evaluation

```python
analyzer = TopologicalAnalyzer()

# Valid domain points
print(analyzer.f_x(0))      # -0.339837...
print(analyzer.f_x(1))      # 0.0
print(analyzer.f_x(-0.6))   # -1.094914...

# Pole singularity (raises ValueError)
try:
    analyzer.f_x(-1.5)
except ValueError as e:
    print(f"Pole detected: {e}")
```

### Domain Analysis

```python
domain = analyzer.analyze_domain_constraints()

print(f"Lower bound: {domain['domain_lower_bound']}")  # -2/3
print(f"Upper bound: {domain['domain_upper_bound']}")  # -4  
print(f"Pole: {domain['pole']}")                      # -3/2

for name, start, end in domain['admissible_intervals']:
    print(f"Interval {name}: ({start}, {end}]")
```

### Geodesic Mapping

```python
x_points = [0, 1, -0.6]
mapping = analyzer.map_fx_to_geodesic(x_points, k=0.3)

for corr in mapping['correspondences']:
    if corr['fx'] is not None:
        print(f"x={corr['x']:.3f} -> f(x)={corr['fx']:.3f}, θ'={corr['geodesic']:.3f}")
```

### Comprehensive Analysis

```python
results = analyzer.comprehensive_analysis()

# Check hypothesis validation
validation = results['hypothesis_validation']
print(f"Domain bounded correctly: {validation['domain_bounded_correctly']}")
print(f"Bridge established: {validation['topological_bridge_established']}")

# Examine invariant alignment  
alignment = results['invariant_alignment']
print(f"e² resonance: {alignment['e_squared_resonance']}")
print(f"φ resonance: {alignment['phi_resonance']}")
```

## Testing

Run the test suite to validate implementation:

```bash
python test_topological_analysis.py
```

The test suite includes:
- ✅ Basic f(x) function calculation
- ✅ Domain violation detection  
- ✅ Pole singularity handling
- ✅ Geodesic resolution validation
- ✅ Mapping correspondence verification
- ✅ Mathematical constant accuracy
- ✅ Comprehensive analysis integration

## Demonstration

Run the interactive demonstration:

```bash
python demo_topological_analysis.py
```

This showcases:
- Complete topological analysis workflow
- Mathematical connections visualization
- Hypothesis validation results
- Theoretical implications summary

## Mathematical Insights

### Key Findings

1. **Topological Correspondence**: Singularities in f(x) space correspond to "ruptures" in Z Framework discrete space

2. **Universal Scaling**: Golden ratio φ provides consistent scaling across both mathematical domains

3. **Optimal Curvature**: k* ≈ 0.3 emerges as optimal parameter for both geodesic resolution and density enhancement

4. **Invariant Alignment**: Arcsine domain constraints create bounded "admissible universes" similar to c = e² invariant

### Theoretical Implications

- **Cross-Domain Mathematics**: Demonstrates universal principles spanning different mathematical frameworks
- **Topological Invariants**: Shows how singularities encode structural information across domains  
- **Optimization Convergence**: k* ≈ 0.3 appears as natural optimum in multiple contexts
- **Geometric Resonance**: φ and e² emerge as fundamental constants bridging spaces

## References

- Original Z Framework implementation: `z_framework.py`
- Geodesic mapping validation: `bio_v_arbitrary.py`
- Experimental setup documentation: `EXPERIMENTAL_SETUP.md`
- Issue specification: GitHub Issue #26

---

*Developed by Dionisio A. Lopez ("Big D"), Z Framework Developer*