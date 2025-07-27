import numpy as np
from scipy.fft import fft
from scipy.stats import entropy
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import random

# =====================
# WAVECRISPR-SIGNAL CORE
# =====================
weights = {'A': 1 + 0j, 'T': -1 + 0j, 'C': 0 + 1j, 'G': 0 - 1j}


def build_waveform(seq):
    """Convert DNA sequence to complex waveform"""
    s = np.cumsum([0.34] * len(seq))
    return np.array([weights[base] * np.exp(2j * np.pi * s[i]) for i, base in enumerate(seq)])


def compute_spectrum(waveform):
    """Compute FFT magnitude spectrum"""
    return np.abs(fft(waveform))


def spectral_features(seq):
    """Extract key spectral properties from sequence"""
    wave = build_waveform(seq)
    spec = compute_spectrum(wave)
    ps = spec / np.sum(spec)  # Normalized power spectrum

    # Skip DC component (index 0)
    return {
        'entropy': entropy(ps[1:], base=2),
        'sidelobes': np.sum(spec[1:] > (0.25 * np.max(spec[1:]))),
        'dominant_freq': np.argmax(spec[1:]) + 1,
        'harmonic_power': np.mean(spec[10:30])  # Mid-frequency harmonics
    }


# ==============================
# RIBONN-INSPIRED TE PREDICTION
# ==============================
def codon_adaptation_index(seq):
    """Calculate codon adaptation index (simplified)"""
    optimal_codons = ['GCT', 'GCC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
                      'AAT', 'AAC', 'GAT', 'GAC', 'TGT', 'TGC', 'CAA', 'CAG',
                      'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC',
                      'ATT', 'ATC', 'ATA', 'AAA', 'AAG', 'TTG', 'CTT', 'CTC',
                      'CTA', 'CTG', 'TTA', 'ATG', 'TTC', 'CCT', 'CCC', 'CCA',
                      'CCG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'ACT',
                      'ACC', 'ACA', 'ACG', 'TGG', 'TAT', 'TAC', 'GTT', 'GTC',
                      'GTA', 'GTG', 'TAA', 'TAG', 'TGA']
    codon_count = 0
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i + 3] in optimal_codons:
            codon_count += 1
    return codon_count / (len(seq) // 3) if len(seq) > 2 else 0


def extract_features(seq):
    """Extract combined features for TE prediction"""
    # Traditional sequence features
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    utr_length = len(seq)
    cai = codon_adaptation_index(seq)

    # Spectral features
    spec_feats = spectral_features(seq)

    return [
        gc_content,
        utr_length,
        cai,
        spec_feats['entropy'],
        spec_feats['sidelobes'],
        spec_feats['harmonic_power']
    ]


# ========================
# SYNTHETIC DATA GENERATOR
# ========================
def generate_sequence(length=50):
    """Generate random DNA sequence"""
    return ''.join(random.choices('ATCG', k=length))


def synthetic_te(features):
    """Synthetic translation efficiency model"""
    # Weights: [gc, length, cai, entropy, sidelobes, harmonic]
    weights = [0.3, -0.002, 0.4, -0.25, -0.1, 0.15]
    te = sum(w * f for w, f in zip(weights, features))

    # Apply non-linear transformation
    te = 1 / (1 + np.exp(-te))  # Sigmoid

    # Add noise
    return te + random.gauss(0, 0.1)


def generate_dataset(n_samples=200, min_len=30, max_len=120):
    """Generate synthetic dataset with TE labels"""
    sequences = []
    features = []
    te_values = []

    for _ in range(n_samples):
        length = random.randint(min_len, max_len)
        seq = generate_sequence(length)
        seq_features = extract_features(seq)
        te = synthetic_te(seq_features)

        sequences.append(seq)
        features.append(seq_features)
        te_values.append(te)

    return sequences, np.array(features), np.array(te_values)


# ===================
# DEMONSTRATION & VIS
# ===================
def run_proof_of_concept():
    # Generate synthetic dataset
    sequences, X, y = generate_dataset(n_samples=300)

    # Split features: traditional vs traditional+spectral
    X_traditional = X[:, :3]  # GC, Length, CAI
    X_combined = X[:, :]  # All features

    # Train models
    model_trad = LinearRegression().fit(X_traditional, y)
    model_combined = LinearRegression().fit(X_combined, y)

    # Predictions
    y_pred_trad = model_trad.predict(X_traditional)
    y_pred_combined = model_combined.predict(X_combined)

    # Evaluate
    r2_trad = r2_score(y, y_pred_trad)
    r2_combined = r2_score(y, y_pred_combined)

    # Visualization
    plt.figure(figsize=(10, 6))
    plt.scatter(y, y_pred_trad, alpha=0.5, label=f'Traditional Features (R²={r2_trad:.3f})')
    plt.scatter(y, y_pred_combined, alpha=0.5, label=f'Combined Features (R²={r2_combined:.3f})')
    plt.plot([min(y), max(y)], [min(y), max(y)], 'k--', lw=1)
    plt.xlabel('True TE (Synthetic)')
    plt.ylabel('Predicted TE')
    plt.title('TE Prediction: Traditional vs Spectral-Enhanced Features')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Feature importance
    plt.figure(figsize=(10, 4))
    features = ['GC Content', 'UTR Length', 'CAI', 'Spectral Entropy', 'Sidelobes', 'Harmonic Power']
    importance = np.abs(model_combined.coef_)
    plt.bar(features, importance / np.sum(importance))
    plt.title('Feature Importance in Combined Model')
    plt.ylabel('Normalized Importance')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Sequence-spectrum example
    plt.figure(figsize=(10, 4))
    sample_seq = sequences[np.argmax(y)]  # Highest TE sequence
    wave = build_waveform(sample_seq)
    spec = compute_spectrum(wave)
    plt.plot(spec[1:50])
    plt.title(f'Frequency Spectrum: High-TE Sequence\n"{sample_seq[:15]}..."')
    plt.xlabel('Frequency Index')
    plt.ylabel('Magnitude')
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.show()

    # Return performance metrics
    return {
        'traditional_r2': r2_trad,
        'combined_r2': r2_combined,
        'improvement': r2_combined - r2_trad
    }


# ========
# EXECUTION
# ========
if __name__ == "__main__":
    results = run_proof_of_concept()
    print(f"\nTraditional Model R²: {results['traditional_r2']:.3f}")
    print(f"Combined Model R²:    {results['combined_r2']:.3f}")
    print(f"Improvement:          +{results['improvement']:.3f}")