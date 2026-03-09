# Phase-Weighted Scorecard (Authoritative Spec)

This document describes the behavior currently implemented in:
- `applications/phase_weighted_scorecard.py`
- `applications/phase_weighted_scorecard_cli.py`

## Scope

The scorecard computes phase-weighted spectral features for guide sequences and, when a target is provided, computes disruption-based Z scores.

## Implemented Computation Pipeline

1. Sequence validation:
- `validate_dna_sequence(seq, allow_rna=...)` is delegated to `wave_crispr_signal.sequence_utils`.

2. Encoding and phase weighting:
- Complex encoding is delegated to `wave_crispr_signal.sequence_utils.encode_complex_phase_weighted`.
- Phase weighting uses `theta_prime_zero_based` via wrapper `theta_prime(n, k)`.
- Default phase parameter: `K_STAR = 0.3`.

3. Spectral feature extraction (`compute_spectral_features`):
- FFT magnitude spectrum.
- Shannon entropy of normalized spectrum.
- Dominant frequency index/magnitude (excluding DC).
- Sidelobe count from `scipy.signal.find_peaks` with default threshold ratio `0.1`.
- Sequence diversity from normalized base-composition entropy.

4. Disruption features (`compute_disruption_features`):
- `delta_entropy = H(mut) - H(ref)`
- `delta_freq = |f_mut - f_ref|`
- `delta_sidelobes = |sidelobes_mut - sidelobes_ref|`

5. Composite Z score (`compute_z_score`):
- `delta_spectral = 0.4*delta_entropy + 0.35*delta_freq + 0.25*delta_sidelobes`
- `kappa = diversity * ln(n+1) / e^2` for `n >= 10` else `1.0`
- `z_score = sigmoid((delta_spectral / PHI), kappa)` where `PHI = 1.6180339887498949`

6. Guide scoring (`score_guide`):
- Without target: returns `guide_features`, and `z_score`/`delta_spectral` as `None`.
- With target: returns disruption components, `z_score`, and both guide/target features.

7. Batch scoring (`score_guide_batch`):
- Scores each guide, optionally with aligned targets.
- Returns per-item error objects on failure instead of aborting entire batch.

## CLI Contract (`phase_weighted_scorecard_cli.py`)

Commands:
- `score --guide <SEQ> [--target <SEQ>] [--k <float>] [--rna] [--output <json>]`
- `batch --input <csv> --output <csv> [--k <float>] [--rna] [--seed <int>]`
- `analyze --ref <fasta> --mut <fasta> --output <json> [--k <float>] [--rna] [--seed <int>]`

Input details:
- `batch` expects CSV columns `guide` and optional `id`, `target`.
- `analyze` compares matching FASTA record IDs between reference and mutant files.

Outputs:
- `score` prints spectral features or disruption metrics.
- `batch` writes CSV feature/disruption summaries.
- `analyze` writes per-sequence comparison plus summary stats JSON.
