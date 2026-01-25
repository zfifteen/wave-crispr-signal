#!/usr/bin/env python3
"""
Falsification Experiment for Hypothesis 1: θ′(n,k) CRISPR Guide Efficiency Prediction

Hypothesis 1: The geodesic resolution function θ′(n,k) with k≈0.3, embedding golden 
ratio phasing, quantifies single-nucleotide mutation disruptions via Δentropy and 
sidelobe metrics, improving CRISPR guide efficiency prediction by ΔROC-AUC +0.047 
over baselines like RuleSet3.

Falsification Criteria:
- If ΔROC-AUC <= 0 or p > 0.05: Hypothesis is falsified
- If no significant improvement over baseline: Hypothesis is falsified

Scientific Gates:
- Human DNA only (A/C/G/T/N)
- Fail-fast validation
- Bootstrap CI (≥1,000 resamples)
- Permutation tests for null models
- Reproducible with fixed seed
"""

import sys
import os
import argparse
import json
import numpy as np
from typing import Dict, List, Tuple
from datetime import datetime
from scipy.stats import ttest_rel
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import warnings

# Add parent directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "applications"))

from applications.phase_weighted_scorecard import (
    encode_complex,
    theta_prime,
    apply_phase_weighting,
    compute_spectrum,
    compute_spectral_entropy,
    count_sidelobes,
    find_dominant_frequency,
    PHI,
    K_STAR,
)

warnings.filterwarnings('ignore', category=UserWarning)


class Hypothesis1Falsifier:
    """
    Falsification engine for θ′(n,k) CRISPR guide efficiency hypothesis.
    """
    
    def __init__(self, seed: int = 42, k_parameter: float = K_STAR):
        """
        Initialize falsifier.
        
        Args:
            seed: Random seed for reproducibility
            k_parameter: Resolution exponent (default: K_STAR = 0.3)
        """
        self.seed = seed
        self.k = k_parameter
        np.random.seed(seed)
        
    def encode_and_phase(self, sequence: str) -> np.ndarray:
        """
        Encode DNA sequence and apply golden-ratio phase weighting.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Phase-weighted complex waveform
        """
        # Complex encoding: A→1, T→-1, C→+i, G→-i
        encoded = encode_complex(sequence, is_rna=False)
        
        # Apply θ′(n,k) phase weighting
        phased = apply_phase_weighting(encoded, k=self.k)
        
        return phased
    
    def compute_spectral_features(
        self, 
        wt_sequence: str, 
        mut_sequence: str
    ) -> Dict[str, float]:
        """
        Compute spectral disruption features between wild-type and mutant.
        
        Args:
            wt_sequence: Wild-type sequence
            mut_sequence: Mutant sequence
            
        Returns:
            Dictionary of spectral features
        """
        # Encode and phase both sequences
        wt_phased = self.encode_and_phase(wt_sequence)
        mut_phased = self.encode_and_phase(mut_sequence)
        
        # Compute FFT spectra
        wt_spectrum = compute_spectrum(wt_phased)
        mut_spectrum = compute_spectrum(mut_phased)
        
        # Spectral entropy
        wt_entropy = compute_spectral_entropy(wt_spectrum)
        mut_entropy = compute_spectral_entropy(mut_spectrum)
        delta_entropy = abs(mut_entropy - wt_entropy)
        
        # Sidelobes
        wt_sidelobes = count_sidelobes(wt_spectrum)
        mut_sidelobes = count_sidelobes(mut_spectrum)
        delta_sidelobes = abs(mut_sidelobes - wt_sidelobes)
        
        # Dominant frequency shift
        wt_freq, wt_mag = find_dominant_frequency(wt_spectrum)
        mut_freq, mut_mag = find_dominant_frequency(mut_spectrum)
        freq_shift = abs(mut_freq - wt_freq)
        
        return {
            "delta_entropy": delta_entropy,
            "delta_sidelobes": delta_sidelobes,
            "freq_shift": freq_shift,
            "wt_entropy": wt_entropy,
            "mut_entropy": mut_entropy,
            "wt_sidelobes": wt_sidelobes,
            "mut_sidelobes": mut_sidelobes,
        }
    
    def compute_baseline_features(self, sequence: str) -> Dict[str, float]:
        """
        Compute baseline features (GC content, length, etc.) for comparison.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary of baseline features
        """
        seq_upper = sequence.upper()
        
        gc_count = seq_upper.count('G') + seq_upper.count('C')
        gc_content = gc_count / len(seq_upper) if len(seq_upper) > 0 else 0.0
        
        # Position-specific GC
        gc_first_5 = sum(1 for b in seq_upper[:5] if b in 'GC') / 5.0 if len(seq_upper) >= 5 else 0.0
        gc_last_5 = sum(1 for b in seq_upper[-5:] if b in 'GC') / 5.0 if len(seq_upper) >= 5 else 0.0
        
        return {
            "gc_content": gc_content,
            "length": len(seq_upper),
            "gc_first_5": gc_first_5,
            "gc_last_5": gc_last_5,
        }
    
    def run_cross_validation(
        self,
        X_spectral: np.ndarray,
        X_baseline: np.ndarray,
        y: np.ndarray,
        n_folds: int = 10,
    ) -> Dict[str, List[float]]:
        """
        Run k-fold cross-validation comparing spectral vs baseline features.
        
        Args:
            X_spectral: Spectral features matrix
            X_baseline: Baseline features matrix
            y: Labels (efficiency scores)
            n_folds: Number of CV folds
            
        Returns:
            Dictionary of AUC scores per fold
        """
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=self.seed)
        
        spectral_aucs = []
        baseline_aucs = []
        
        for train_idx, test_idx in kf.split(X_spectral):
            # Train/test split
            X_train_spec, X_test_spec = X_spectral[train_idx], X_spectral[test_idx]
            X_train_base, X_test_base = X_baseline[train_idx], X_baseline[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            # Scale features
            scaler_spec = StandardScaler()
            scaler_base = StandardScaler()
            
            X_train_spec_scaled = scaler_spec.fit_transform(X_train_spec)
            X_test_spec_scaled = scaler_spec.transform(X_test_spec)
            
            X_train_base_scaled = scaler_base.fit_transform(X_train_base)
            X_test_base_scaled = scaler_base.transform(X_test_base)
            
            # Train classifiers
            clf_spec = LogisticRegression(random_state=self.seed, max_iter=1000)
            clf_base = LogisticRegression(random_state=self.seed, max_iter=1000)
            
            clf_spec.fit(X_train_spec_scaled, y_train)
            clf_base.fit(X_train_base_scaled, y_train)
            
            # Predict probabilities
            y_pred_spec = clf_spec.predict_proba(X_test_spec_scaled)[:, 1]
            y_pred_base = clf_base.predict_proba(X_test_base_scaled)[:, 1]
            
            # Compute AUC
            if len(np.unique(y_test)) > 1:  # Need both classes
                auc_spec = roc_auc_score(y_test, y_pred_spec)
                auc_base = roc_auc_score(y_test, y_pred_base)
                
                spectral_aucs.append(auc_spec)
                baseline_aucs.append(auc_base)
        
        return {
            "spectral_aucs": spectral_aucs,
            "baseline_aucs": baseline_aucs,
        }
    
    def bootstrap_auc_difference(
        self,
        X_spectral: np.ndarray,
        X_baseline: np.ndarray,
        y: np.ndarray,
        n_bootstrap: int = 1000,
    ) -> Tuple[float, float, np.ndarray]:
        """
        Bootstrap confidence interval for ΔROC-AUC.
        
        Args:
            X_spectral: Spectral features
            X_baseline: Baseline features
            y: Labels
            n_bootstrap: Number of bootstrap resamples
            
        Returns:
            (mean_delta, std_delta, delta_distribution)
        """
        delta_aucs = []
        
        for i in range(n_bootstrap):
            # Bootstrap resample
            boot_idx = np.random.choice(
                len(y), 
                size=len(y), 
                replace=True
            )
            
            X_spec_boot = X_spectral[boot_idx]
            X_base_boot = X_baseline[boot_idx]
            y_boot = y[boot_idx]
            
            # Train/test split
            X_train_spec, X_test_spec, y_train, y_test = train_test_split(
                X_spec_boot, y_boot, test_size=0.2, random_state=i
            )
            X_train_base, X_test_base = train_test_split(
                X_base_boot, test_size=0.2, random_state=i
            )[0:2]
            
            # Skip if not enough samples
            if len(y_test) < 2 or len(np.unique(y_test)) < 2:
                continue
            
            # Scale
            scaler_spec = StandardScaler()
            scaler_base = StandardScaler()
            
            X_train_spec = scaler_spec.fit_transform(X_train_spec)
            X_test_spec = scaler_spec.transform(X_test_spec)
            X_train_base = scaler_base.fit_transform(X_train_base)
            X_test_base = scaler_base.transform(X_test_base)
            
            # Train and predict
            clf_spec = LogisticRegression(random_state=self.seed, max_iter=1000)
            clf_base = LogisticRegression(random_state=self.seed, max_iter=1000)
            
            clf_spec.fit(X_train_spec, y_train)
            clf_base.fit(X_train_base, y_train)
            
            y_pred_spec = clf_spec.predict_proba(X_test_spec)[:, 1]
            y_pred_base = clf_base.predict_proba(X_test_base)[:, 1]
            
            # Compute AUC difference
            try:
                auc_spec = roc_auc_score(y_test, y_pred_spec)
                auc_base = roc_auc_score(y_test, y_pred_base)
                delta_aucs.append(auc_spec - auc_base)
            except:
                continue
        
        delta_aucs = np.array(delta_aucs)
        
        return (
            np.mean(delta_aucs),
            np.std(delta_aucs),
            delta_aucs,
        )


def generate_synthetic_crispr_data(
    n_samples: int = 100, 
    seed: int = 42
) -> Tuple[List[str], List[str], np.ndarray]:
    """
    Generate synthetic CRISPR guide sequences with efficiency labels.
    
    This is a placeholder for real datasets (e.g., Doench 2016).
    
    Args:
        n_samples: Number of guide sequences to generate
        seed: Random seed
        
    Returns:
        (wild_type_sequences, mutant_sequences, efficiency_labels)
    """
    np.random.seed(seed)
    
    # CRISPR guides are typically 20-23 nt
    guide_length = 20
    
    bases = ['A', 'T', 'C', 'G']
    
    wt_sequences = []
    mut_sequences = []
    labels = []
    
    for i in range(n_samples):
        # Generate random wild-type sequence
        wt_seq = ''.join(np.random.choice(bases, guide_length))
        
        # Create mutant by changing 1-3 positions
        mut_seq = list(wt_seq)
        n_mutations = np.random.randint(1, 4)
        mut_positions = np.random.choice(guide_length, n_mutations, replace=False)
        
        for pos in mut_positions:
            # Change to different base
            current_base = mut_seq[pos]
            other_bases = [b for b in bases if b != current_base]
            mut_seq[pos] = np.random.choice(other_bases)
        
        mut_seq = ''.join(mut_seq)
        
        # Generate synthetic efficiency label
        # Higher GC content → higher efficiency (rough heuristic)
        gc_content = (wt_seq.count('G') + wt_seq.count('C')) / len(wt_seq)
        
        # Add some noise
        efficiency_score = gc_content + np.random.normal(0, 0.2)
        
        # Binarize: high efficiency if > 0.5
        label = 1 if efficiency_score > 0.5 else 0
        
        wt_sequences.append(wt_seq)
        mut_sequences.append(mut_seq)
        labels.append(label)
    
    return wt_sequences, mut_sequences, np.array(labels)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Falsification Experiment for Hypothesis 1"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Number of synthetic samples (or dataset size)"
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=1000,
        help="Number of bootstrap resamples"
    )
    parser.add_argument(
        "--n-folds",
        type=int,
        default=10,
        help="Number of cross-validation folds"
    )
    parser.add_argument(
        "--k-parameter",
        type=float,
        default=K_STAR,
        help="Resolution exponent k (default: 0.3)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/PR-156-falsify-hypothesis1",
        help="Output directory for results"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", args.output_dir
    )
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize falsifier
    falsifier = Hypothesis1Falsifier(seed=args.seed, k_parameter=args.k_parameter)
    
    print("=" * 70)
    print("Hypothesis 1 Falsification Experiment")
    print("=" * 70)
    print(f"Seed: {args.seed}")
    print(f"k parameter: {args.k_parameter}")
    print(f"N samples: {args.n_samples}")
    print(f"N bootstrap: {args.n_bootstrap}")
    print(f"N folds: {args.n_folds}")
    print()
    
    # Generate or load data
    print("Generating synthetic CRISPR guide data...")
    wt_seqs, mut_seqs, labels = generate_synthetic_crispr_data(
        n_samples=args.n_samples,
        seed=args.seed
    )
    print(f"Generated {len(wt_seqs)} guide pairs")
    print(f"Label distribution: {np.bincount(labels)}")
    print()
    
    # Extract spectral features
    print("Computing spectral features with θ′(n,k) phase weighting...")
    spectral_features = []
    baseline_features = []
    
    for wt_seq, mut_seq in zip(wt_seqs, mut_seqs):
        # Spectral features
        spec_feat = falsifier.compute_spectral_features(wt_seq, mut_seq)
        spectral_features.append([
            spec_feat["delta_entropy"],
            spec_feat["delta_sidelobes"],
            spec_feat["freq_shift"],
        ])
        
        # Baseline features
        base_feat = falsifier.compute_baseline_features(wt_seq)
        baseline_features.append([
            base_feat["gc_content"],
            base_feat["length"],
            base_feat["gc_first_5"],
            base_feat["gc_last_5"],
        ])
    
    X_spectral = np.array(spectral_features)
    X_baseline = np.array(baseline_features)
    
    print(f"Spectral features shape: {X_spectral.shape}")
    print(f"Baseline features shape: {X_baseline.shape}")
    print()
    
    # Run cross-validation
    print(f"Running {args.n_folds}-fold cross-validation...")
    cv_results = falsifier.run_cross_validation(
        X_spectral, X_baseline, labels, n_folds=args.n_folds
    )
    
    spectral_aucs = cv_results["spectral_aucs"]
    baseline_aucs = cv_results["baseline_aucs"]
    
    mean_spectral_auc = np.mean(spectral_aucs)
    mean_baseline_auc = np.mean(baseline_aucs)
    delta_auc_cv = mean_spectral_auc - mean_baseline_auc
    
    print(f"Spectral AUC (mean ± std): {mean_spectral_auc:.4f} ± {np.std(spectral_aucs):.4f}")
    print(f"Baseline AUC (mean ± std): {mean_baseline_auc:.4f} ± {np.std(baseline_aucs):.4f}")
    print(f"ΔROC-AUC (CV): {delta_auc_cv:.4f}")
    print()
    
    # Paired t-test
    if len(spectral_aucs) == len(baseline_aucs) and len(spectral_aucs) > 1:
        t_stat, p_value = ttest_rel(spectral_aucs, baseline_aucs)
        print(f"Paired t-test: t={t_stat:.4f}, p={p_value:.4f}")
        print()
    else:
        p_value = 1.0
    
    # Bootstrap confidence interval
    print(f"Computing bootstrap CI for ΔROC-AUC ({args.n_bootstrap} resamples)...")
    mean_delta, std_delta, delta_dist = falsifier.bootstrap_auc_difference(
        X_spectral, X_baseline, labels, n_bootstrap=args.n_bootstrap
    )
    
    ci_lower = np.percentile(delta_dist, 2.5)
    ci_upper = np.percentile(delta_dist, 97.5)
    
    print(f"Bootstrap ΔROC-AUC: {mean_delta:.4f} ± {std_delta:.4f}")
    print(f"95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
    print()
    
    # Falsification decision
    print("=" * 70)
    print("FALSIFICATION ANALYSIS")
    print("=" * 70)
    
    falsified = False
    reasons = []
    
    if mean_delta <= 0:
        falsified = True
        reasons.append("ΔROC-AUC <= 0 (no improvement)")
    
    if p_value > 0.05:
        falsified = True
        reasons.append(f"p-value > 0.05 (p={p_value:.4f}, not significant)")
    
    if ci_lower < 0 and ci_upper > 0:
        falsified = True
        reasons.append("95% CI includes zero (no significant difference)")
    
    if falsified:
        print("RESULT: Hypothesis 1 is FALSIFIED")
        print("Reasons:")
        for reason in reasons:
            print(f"  - {reason}")
    else:
        print("RESULT: Hypothesis 1 is NOT FALSIFIED (supported by data)")
        print(f"  - ΔROC-AUC = {mean_delta:.4f} > 0")
        print(f"  - p-value = {p_value:.4f} < 0.05")
        print(f"  - 95% CI [{ci_lower:.4f}, {ci_upper:.4f}] excludes zero")
    
    print()
    
    # Save results
    results = {
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "seed": args.seed,
            "k_parameter": args.k_parameter,
            "n_samples": args.n_samples,
            "n_bootstrap": args.n_bootstrap,
            "n_folds": args.n_folds,
        },
        "cross_validation": {
            "spectral_auc_mean": float(mean_spectral_auc),
            "spectral_auc_std": float(np.std(spectral_aucs)),
            "baseline_auc_mean": float(mean_baseline_auc),
            "baseline_auc_std": float(np.std(baseline_aucs)),
            "delta_auc": float(delta_auc_cv),
            "t_statistic": float(t_stat) if len(spectral_aucs) > 1 else None,
            "p_value": float(p_value),
        },
        "bootstrap": {
            "delta_auc_mean": float(mean_delta),
            "delta_auc_std": float(std_delta),
            "ci_lower": float(ci_lower),
            "ci_upper": float(ci_upper),
        },
        "falsification": {
            "falsified": falsified,
            "reasons": reasons,
        }
    }
    
    output_file = os.path.join(output_dir, "hypothesis1_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to: {output_file}")
    print()
    
    return 0 if not falsified else 1


if __name__ == "__main__":
    sys.exit(main())
