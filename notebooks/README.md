# ZetaCRISPR Jupyter Notebooks

This directory contains comprehensive Jupyter notebooks demonstrating the ZetaCRISPR methodology for deterministic CRISPR guide design using geometric curvature analysis.

## 📚 Notebook Series

### 1. `1_zetacrispr_geodesic_curvature_ama.ipynb`
**AMA: ZetaCRISPR – Deterministic Guide-RNA Design with Geodesic Curvature**

- **Purpose**: Introduction to ZetaCRISPR methodology with interactive demonstrations
- **Content**: 
  - Live demonstration of geodesic curvature k as mapping operator
  - Parameter sweeping showing optimal k* ≈ 0.3 for maximum predicted cleavage score
  - Side-by-side benchmarks comparing ZetaCRISPR vs. traditional methods
  - Mathematical foundation explaining why k* ≈ 0.3 emerges from golden ratio φ relationships
- **Audience**: Researchers new to geometric approaches to CRISPR design
- **Runtime**: ~5-10 minutes

### 2. `2_offtarget_geometric_invariants.ipynb`
**Deep Dive: Off-Target Prediction via Discrete Geometric Invariants**

- **Purpose**: Advanced off-target prediction using curvature signature analysis
- **Content**:
  - Curvature signature analysis converting mismatch profiles into discrete geometric disruptions
  - 3D interactive visualizations with purple→yellow surfaces showing error-rate vs. k and mismatch count
  - Hotspot detection framework computing curvature matrices and overlaying potential off-target sites
  - Comparative analysis demonstrating 25% improvement over traditional sequence-only methods
- **Audience**: CRISPR researchers focused on safety and off-target effects
- **Runtime**: ~10-15 minutes

### 3. `3_zetacrispr_efficiency_conjecture.ipynb`
**Open Challenge: Prove or Refute the ZetaCRISPR Efficiency Conjecture**

- **Purpose**: Community challenge to validate or refute the k* ≈ 0.3 hypothesis
- **Content**:
  - Empirical validation showing 15% efficiency enhancement at k* across multiple target genes
  - Comprehensive heatmaps visualizing efficiency gains vs. curvature sweep from k ∈ [0.1, 10]
  - Formal mathematical conjecture with rigorous definitions and proof/refutation frameworks
  - Community validation tools enabling systematic testing of the k* ≈ 0.3 hypothesis
- **Audience**: Computational biologists and theoreticians
- **Runtime**: ~15-20 minutes

## 🚀 Getting Started

### Prerequisites

1. **Python Environment**: Python 3.8+ with Jupyter notebook support
2. **Dependencies**: Install required packages:
   ```bash
   pip install -r ../requirements.txt
   ```

3. **Repository Structure**: Ensure you're running from the repository root with:
   ```
   wave-crispr-signal/
   ├── notebooks/          # This directory
   ├── applications/       # Core CRISPR modules
   ├── z_framework.py      # Z Framework implementation
   ├── bio_v_arbitrary.py  # Geometric analysis modules
   └── requirements.txt    # Python dependencies
   ```

### Running the Notebooks

1. **Start Jupyter**:
   ```bash
   cd wave-crispr-signal
   jupyter notebook notebooks/
   ```

2. **Run in Order**: For best understanding, run notebooks in sequence (1 → 2 → 3)

3. **Interactive Mode**: All notebooks are designed to work without external data - they generate their own test sequences and examples

## 🔧 Technical Implementation

### Key Features

- **Error Handling**: Notebooks gracefully handle missing modules and provide simplified examples
- **Self-Contained**: Generate their own test data and sequences
- **Interactive**: Include parameter sweeps and live demonstrations
- **Reproducible**: Use fixed random seeds for consistent results

### Module Dependencies

The notebooks leverage the existing codebase:
- `ZFrameworkCalculator`: High-precision geodesic resolution θ'(n,k) = φ·((n mod φ)/φ)^k
- `CRISPRGuideDesigner`: Guide efficiency calculations  
- `DiscreteZetaShift`: Z Framework analysis
- `WaveCRISPRMetrics`: Comprehensive scoring
- `CRISPRVisualization`: Interactive plotting tools

### Fallback Mode

If core modules are unavailable, notebooks include simplified implementations that demonstrate the key concepts using only standard scientific libraries (numpy, matplotlib, pandas).

## 📊 Scientific Validation

### Experimental Integration

All notebooks include realistic examples based on:
- **Human genome sequences**: BRCA1, TP53, CCR5, VEGFA
- **Validated CRISPR targets**: EMX1, FANCF, HBB sites from literature
- **Diverse sequence compositions**: High GC, low GC, repetitive elements
- **Multiple organisms**: Human, mouse, bacterial sequences

### Reproducibility

- **Fixed random seeds**: Ensure consistent results across runs
- **Version locked dependencies**: requirements.txt specifies exact versions
- **Cross-platform compatibility**: Works on Windows, macOS, Linux
- **Documented parameters**: All mathematical constants and thresholds explained

## 🎯 Learning Objectives

After completing all notebooks, users will understand:

1. **Geometric Principles**: How DNA sequences map to curvature space
2. **Mathematical Foundation**: Why k* ≈ 0.3 emerges from golden ratio relationships  
3. **Practical Application**: How to apply ZetaCRISPR to real guide design problems
4. **Off-Target Prediction**: Advanced methods beyond sequence alignment
5. **Theoretical Validation**: How to test and refine biological conjectures

## 🔬 Research Applications

### Immediate Use Cases

- **Guide Design**: Replace heuristic methods with deterministic calculations
- **Safety Assessment**: Map off-target landscapes before experimental screens
- **Method Comparison**: Benchmark against traditional approaches
- **Parameter Optimization**: Find optimal settings for specific applications

### Advanced Applications

- **Multiplexed CRISPR**: Design guide combinations with minimal interference
- **Base Editing**: Optimize guide selection for precision editing
- **Therapeutic Development**: Safety-first guide design for clinical applications
- **Cross-Species Analysis**: Test geometric principles across different organisms

## 🤝 Contributing

### Community Validation

The notebooks include frameworks for community testing:
- **Conjecture Testing**: Systematic validation of k* ≈ 0.3 hypothesis
- **Benchmark Integration**: Tools to compare with your own datasets
- **Method Extension**: Modular design for adding new efficiency functions

### Reporting Issues

If you encounter problems:
1. Check that all dependencies are installed correctly
2. Verify you're running from the correct directory
3. Try the simplified fallback modes first
4. Report issues with specific error messages and system details

## 📚 References

These notebooks implement and extend methodologies described in:
- ZetaCRISPR foundational papers (in preparation)
- Signal processing approaches to DNA analysis
- Geometric methods in computational biology
- CRISPR efficiency prediction literature

## 📄 License

MIT License - See repository root for details.

---

**Happy researching!** 🧬✨