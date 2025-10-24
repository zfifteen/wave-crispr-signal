# Demo MVP: Interactive CRISPR Breathing Dynamics Analyzer

This is the public-facing interactive demo showing DNA breathing dynamics and the 10.5-bp spectral peak phenomenon.

## Features

### 🎯 What It Shows
- **10.5-bp spectral peak**: Fundamental helical period of B-DNA
- **Harmonics**: 5.25 bp (2nd) and 3.5 bp (3rd) harmonics
- **Rotational phase wheel**: Polar visualization of phase-dependent activity
- **Breathing Δscore**: Proxy for model improvement from adding breathing features

### 🔒 Scientific Gates Enforced
- ✅ Human DNA only (A/C/G/T for DNA, A/C/G/U for RNA)
- ✅ No fabrication (real nucleotides only)
- ✅ Dimensionless parametrization (rate ratios r, not MHz)
- ✅ Privacy-safe logging (opt-in, SHA256 + features only, no raw sequences)

## Quick Start

### Installation
```bash
# Install dependencies (from repository root)
pip install fastapi uvicorn pydantic numpy scipy

# Or use repository requirements.txt
pip install -r requirements.txt
```

### Running the Demo
```bash
# Navigate to demo directory
cd web_apps/demo_mvp

# Start the server
python app.py

# Or use uvicorn directly
uvicorn app:app --reload
```

The demo will be available at:
- **Main UI**: http://127.0.0.1:8000
- **API docs**: http://127.0.0.1:8000/docs
- **Health check**: http://127.0.0.1:8000/health

## Usage

### Web Interface
1. Open http://127.0.0.1:8000 in your browser
2. Enter:
   - **Guide RNA**: 20-30 nucleotides (A/C/G/U)
   - **Target DNA**: 30-100 base pairs (A/C/G/T)
3. (Optional) Adjust advanced parameters:
   - **Rate Ratio**: Dimensionless k_GC/k_AT (default: 20.0)
   - **Weight Mode**: Currently rate_ratio only
   - **Logging consent**: Opt-in for privacy-safe analytics
4. Click "Analyze Sequence"
5. View results:
   - Spectral peaks at 10.5, 5.25, 3.5 bp
   - Rotational phase wheel (polar plot)
   - Bar chart of harmonic magnitudes
   - Sequence properties

### API Usage

#### POST /analyze
Analyze a guide RNA + target DNA pair.

**Request:**
```bash
curl -X POST http://127.0.0.1:8000/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "guide": "GACGAUCGAUCGAUCGAUCG",
    "target": "ATGCGATCGATCGATCGATCGCTAGCTAGCTA",
    "rate_ratio": 20.0,
    "weight_mode": "rate_ratio",
    "log_consent": false
  }'
```

**Response:**
```json
{
  "success": true,
  "peak": {
    "p10": 0.5876,
    "p5_25": 1.6384,
    "p3_5": 1.8724,
    "angle": 0.2699
  },
  "rotation": {
    "phi": [0.2618, 0.7854, ...],
    "curve": [3.1276, 3.1276, ...]
  },
  "breathing_features": {
    "P10_5": 3.7517,
    "P5_25": 8.3302,
    "P3_5": 14.5463,
    "length": 32,
    "gc_content": 0.5
  },
  "meta": {
    "guide_length": 20,
    "target_length": 32,
    "weight_mode": "rate_ratio",
    "rate_ratio": 20.0,
    "breathing_score": 0.0375,
    "temperature_c": null,
    "mg_mM": null
  },
  "message": "Analysis completed successfully"
}
```

#### GET /health
Check service health.

```bash
curl http://127.0.0.1:8000/health
```

#### GET /info
Get information about available parameters.

```bash
curl http://127.0.0.1:8000/info
```

## Technical Details

### Dimensionless Parametrization

The demo uses **rate ratios** instead of absolute frequencies to avoid unit ambiguity:

- **Rate ratio**: r = k_GC / k_AT (dimensionless)
- **Default**: r = 20.0 (AT opens ~20× faster than GC)
- **Typical range**: 5.0 to 200.0

**Encoding formula:**
```
alpha = log(r)          # Strength factor
beta = 0.3 * log(r)     # Phase factor

AT: -alpha + j*beta     (weak bonds, fast opening)
GC: +alpha - j*beta     (strong bonds, slow opening)
```

### Spectral Analysis

Uses **Chirp Z-Transform (CZT)** approach for precise frequency evaluation:
1. Remove DC component (mean subtraction)
2. Apply Hann window (reduce spectral leakage)
3. Compute Fourier coefficient at f = 1/period
4. Extract magnitude (power) and phase

**Key periods:**
- **10.5 bp**: Helical period (1 turn of B-DNA)
- **5.25 bp**: First harmonic (2 turns/wavelength)
- **3.5 bp**: Second harmonic (3 turns/wavelength)

### Privacy & Logging

When logging is enabled (opt-in only):
- **Stored**: SHA256 hash + features (length, GC%, peaks)
- **NOT stored**: Raw sequences
- **Purpose**: Aggregate analytics, no individual identification

## File Structure

```
web_apps/demo_mvp/
├── app.py              # FastAPI backend
├── index.html          # Web interface (HTMX + Plotly)
├── test_api.py         # API tests
└── README.md           # This file
```

## Testing

Run the test suite:
```bash
cd web_apps/demo_mvp
python test_api.py
```

Expected output:
```
Testing Demo MVP API...
✓ Health check passed
✓ Info endpoint passed
✓ Analyze endpoint passed
✓ Invalid guide rejection passed
✓ Invalid target rejection passed
✓ All tests passed!
```

## Scientific Background

### DNA Breathing Dynamics
- **AT pairs**: 2 hydrogen bonds → fast opening (~10⁷ Hz experimental)
- **GC pairs**: 3 hydrogen bonds → slow opening (~10⁹ Hz experimental)
- **Separation**: ~100× frequency difference

### Helical Periodicity
- **B-DNA**: ~10.5 base pairs per helical turn
- **Effect**: Rotational phasing modulates accessibility
- **CRISPR**: Guide RNA binding affected by target accessibility

### Z Framework Connection
- **Discrete domain**: Z = A(B / e²)
- **Geometric resolution**: θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
- **Period**: φ = 10.5 bp for helical structure

## References

### Repository
- **Spectral utilities**: `experiments/trinity/spectral.py`
- **Repository policy**: `.github/REPOSITORY_POLICY.md`
- **Scientific gates**: `.github/copilot-instructions.md`

### Related Work
- Doench et al. (2016) - On-target CRISPR prediction
- Altan-Bonnet et al. (2003) - DNA breathing dynamics measurements
- Z Framework documentation: `docs/Z_FRAMEWORK.md`

## Future Enhancements

### Planned Features
1. **Nearest-neighbor thermodynamics**: Temperature + ion-dependent weights
2. **Batch analysis**: Upload FASTA files
3. **3D visualization**: Interactive molecular structure
4. **Export**: Download results as JSON/CSV
5. **Rate ratio sweep**: Stability heatmap across r ∈ [5, 200]

### Integration with Trinity Experiments
This demo uses the same spectral utilities as the trinity experiments:
- **Experiment A**: Doench-2016 lift (nested model)
- **Experiment B**: Cross-species generalization
- **Experiment C**: Rotational phase figure

See `experiments/trinity/` for full research implementations.

## Troubleshooting

### Server won't start
```bash
# Check Python version (needs 3.12+)
python --version

# Install missing dependencies
pip install fastapi uvicorn pydantic numpy scipy
```

### Invalid sequence errors
- **Guide RNA**: Must use A/C/G/U (not T)
- **Target DNA**: Must use A/C/G/T (not U)
- Remove whitespace, numbers, and special characters
- No IUPAC ambiguity codes (N, R, Y, etc.)

### Plots not showing
- Ensure internet connection (Plotly loads from CDN)
- Check browser console for JavaScript errors
- Try a different browser (Chrome/Firefox recommended)

## License

MIT License - See repository LICENSE file

## Contact

For questions or issues:
- Open a GitHub issue tagged `demo-mvp`
- See `experiments/trinity/README.md` for research details

---

**Version**: 1.0.0  
**Created**: October 2025  
**Status**: Production Ready
