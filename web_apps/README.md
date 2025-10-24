# DNA Breathing Dynamics Interactive Web Application

## Overview

This interactive web application demonstrates the breakthrough findings from **PR #103** showing that DNA breathing dynamics (base pair opening frequencies) significantly outperform arbitrary encodings in spectral DNA analysis.

## Key Findings Demonstrated

### The Breakthrough
- **First positive result**: DNA breathing dynamics is the first biological property to outperform arbitrary weights
- **Effect size**: Cohen's d = +4.130 (very large effect)
- **Statistical significance**: p < 0.000001

### DNA Breathing Frequencies (Experimental)
- **AT pairs**: ~10 MHz (10⁷ Hz) - Fast opening, 2 hydrogen bonds
- **GC pairs**: ~1 GHz (10⁹ Hz) - Slow opening, 3 hydrogen bonds
- **Separation**: 100× frequency difference

### Biological Selectivity
The encoding shows intelligent selectivity:
- **GC-affecting mutations** (AT↔GC): Breathing wins with large effect
- **Within-class mutations** (A↔T or G↔C): Shows ~0 signal (correct behavior)
- **Random mutations**: Mixed results (expected dilution)

## Scientific Gates Enforced

This application enforces strict scientific gates as per repository policy:

1. **Human DNA Only**: Only accepts A, C, G, T nucleotides (no U or IUPAC codes)
2. **No Fabrication**: Uses real human DNA sequences from test data
3. **Validation**: Fail-fast validation for all inputs
4. **Reproducibility**: Fixed random seeds, documented parameters

## Features

### Interactive Analysis
- Real-time DNA sequence analysis
- Comparison of breathing dynamics vs arbitrary encoding
- Mutation effect visualization
- Frequency spectrum plots

### Visualizations
1. **Sequence Composition**: Base counts and GC content
2. **Mutation Analysis**: Z-scores for different mutation types
3. **Power Spectrum**: Frequency domain comparison
4. **Encoding Weights**: Side-by-side comparison of weight schemes

### Example Sequences
Pre-loaded real human DNA sequences from test data:
- Human TOP2A (GC-rich)
- Human CDK6 (balanced GC content)
- Human MYB (repeat-rich)

## Installation & Usage

### Prerequisites
```bash
# Python 3.12+
python --version

# Install required packages
pip install -r requirements.txt
```

### Running the Application
```bash
# From the repository root
cd web_apps

# Run the Flask application
python breathing_dynamics_demo.py
```

The application will start on `http://127.0.0.1:5000`

### Using the Web Interface

1. **Open your browser** to `http://127.0.0.1:5000`
2. **Enter a DNA sequence** (10-500 bases, A/C/G/T only)
   - Or click an example button to load real human DNA
3. **Click "Analyze Sequence"**
4. **View results**:
   - Sequence composition and GC content
   - Mutation analysis showing breathing vs arbitrary performance
   - Frequency spectrum visualization
   - Encoding weights comparison

## Technical Implementation

### Encoding Strategy

#### Breathing Dynamics Encoder
Based on experimental base pair opening frequencies:
```python
# Real part: log-normalized frequency using formula (log10(freq) - 8.0) * 10.0
# This centers the encoding around 10^8 Hz and scales to ±10 range
# 
# AT: 10^7 Hz → (log10(10^7) - 8.0) * 10.0 = (7 - 8) * 10 = -10.00 (fast opening)
# GC: 10^9 Hz → (log10(10^9) - 8.0) * 10.0 = (9 - 8) * 10 = +10.00 (slow opening)

# Imaginary part: phase from kinetics
# AT: +3.0j (weak bonds, positive phase)
# GC: -3.0j (strong bonds, negative phase)

weights = {
    'A': -10.00 + 3.00j,
    'T': -10.00 + 3.00j,
    'C': +10.00 - 3.00j,
    'G': +10.00 - 3.00j
}
```

#### Phase Modulation
Both encoders use identical phase modulation:
- **Helical periodicity**: 2π × i / 10.5 (DNA double helix)
- **Positional context**: 2π × i / seq_length × 0.3
- **Combined**: Fair comparison between encoders

### Z Framework
Uses discrete domain form: **Z = A(B/c)**
- **A**: Frame-dependent sequence entropy
- **B**: Spectral shift from mutation (normalized)
- **c**: e² ≈ 7.389 (discrete invariant)

### Mutation Types Tested
1. **GC-affecting**: Changes frequency class (AT→GC or GC→AT)
2. **Within-class**: Stays in same class (A↔T or G↔C)
3. **Random**: Any base substitution

## File Structure

```
web_apps/
├── breathing_dynamics_demo.py     # Flask application
└── README.md                       # This file

templates/
└── breathing_dynamics.html         # Web interface

data/
└── test_human_cdna.fasta          # Human DNA test sequences
```

## API Endpoints

### POST /analyze
Analyzes a DNA sequence
```json
Request:
{
  "sequence": "ATGCGATCGATCG..."
}

Response:
{
  "success": true,
  "results": {
    "sequence_length": 100,
    "gc_content": 55.0,
    "base_composition": {...},
    "breathing_spectrum": {...},
    "arbitrary_spectrum": {...},
    "mutation_analysis": {...}
  }
}
```

### GET /example_sequences
Returns pre-loaded human DNA examples
```json
Response:
{
  "success": true,
  "sequences": [
    {
      "name": "Human TOP2A",
      "sequence": "ATGAGTCCG...",
      "length": 100
    }
  ]
}
```

## References

### Original Research
- **PR #103**: Breathing Dynamics Encoding Breakthrough
- **Documentation**: `docs/PR-breathing-dynamics/`
- **Experiment**: `experiments/test_breathing_dynamics_encoding.py`
- **Tests**: `tests/test_breathing_dynamics.py`

### Scientific Background
- Altan-Bonnet et al. (2003) - Base pair opening dynamics, PRL
- Doench et al. (2016) - CRISPR guide efficiency, Nature Biotechnology
- Z Framework: `docs/Z_FRAMEWORK.md`

### Repository Policy
- `.github/REPOSITORY_POLICY.md` - File structure and naming
- `.github/copilot-instructions.md` - Scientific gates

## Validation

### Testing the Application
```bash
# Test with valid human DNA
curl -X POST http://127.0.0.1:5000/analyze \
  -H "Content-Type: application/json" \
  -d '{"sequence":"ATGCGATCGATCGATCG"}'

# Should return analysis results

# Test with invalid sequence
curl -X POST http://127.0.0.1:5000/analyze \
  -H "Content-Type: application/json" \
  -d '{"sequence":"ATGCXYZ"}'

# Should return error about invalid bases
```

### Expected Behavior
- ✅ Validates DNA sequences (A/C/G/T only)
- ✅ Rejects sequences with U, N, or other IUPAC codes
- ✅ Shows higher Z-scores for breathing encoding on GC-affecting mutations
- ✅ Shows ~0 signal for within-class mutations (correct selectivity)
- ✅ Displays frequency spectrum with 100× separation (AT vs GC)

## Limitations

### Current Scope
- **Synthetic testing**: Uses controlled mutation tests, not full CRISPR datasets
- **Single property**: Only breathing dynamics (not combined with other properties)
- **Web performance**: Limited to 500 bases for responsive UI

### Future Enhancements
1. **Real CRISPR data**: Integration with Doench 2016 / Kim 2025 datasets
2. **Combined encoding**: Breathing + electronic transitions + torsional
3. **Batch analysis**: Upload FASTA files for multiple sequences
4. **3D visualization**: Interactive molecular structure view
5. **Export results**: Download analysis as JSON/CSV

## Troubleshooting

### Common Issues

**Problem**: Application won't start
```bash
# Check Python version
python --version  # Should be 3.12+

# Check dependencies
pip install flask numpy scipy biopython
```

**Problem**: No example sequences loading
```bash
# Verify test data exists
ls -la data/test_human_cdna.fasta

# If missing, examples will use hardcoded fallbacks
```

**Problem**: "Invalid bases" error
- Ensure sequence contains only A, C, G, T (uppercase or lowercase)
- Remove any whitespace, numbers, or other characters
- Do NOT use U (RNA) or IUPAC ambiguity codes

## Contributing

When modifying this application:
1. **Maintain scientific gates**: Only human DNA, no fabrication
2. **Test thoroughly**: Verify with real human sequences
3. **Document changes**: Update README and inline comments
4. **Follow repository policy**: See `.github/REPOSITORY_POLICY.md`

## License

MIT License - See repository LICENSE file

## Contact

For questions about this application or the breathing dynamics findings:
- **Issues**: Open a GitHub issue tagged `web-app` or `breathing-dynamics`
- **Documentation**: See `docs/PR-breathing-dynamics/` for full details
- **Repository**: https://github.com/zfifteen/wave-crispr-signal

---

**Created**: October 2025
**Version**: 1.0
**Status**: Production Ready
**Impact**: 🔥 Demonstrates Breakthrough Finding (Cohen's d = +4.130)
