#!/usr/bin/env python3
"""
Falsification Experiment for Hypothesis 2: Z-Invariant Scoring Abstraction

Hypothesis 2: The sequence diversity-dependent curvature weighting κ(n) = d(n) · ln(n+1)/e² 
creates a Z-invariant scoring abstraction that robustly handles variable-length DNA sequences, 
revealing emergent periodicity patterns in off-target profiling.

Falsification Criteria:
- If variance in Z-scores across variable lengths is high (ANOVA p < 0.05): Non-invariant
- If no emergent periodicity (autocorrelation < 0.1): No periodicity detected
- If either criterion fails: Hypothesis is falsified

Scientific Gates:
- Human DNA only (A/C/G/T/N)
- Fail-fast validation
- Statistical validity with ANOVA and autocorrelation analysis
- Reproducible with fixed seed
"""

import sys
import os
import argparse
import json
import numpy as np
from typing import Dict, List, Tuple
from datetime import datetime
from scipy.stats import f_oneway, pearsonr
from scipy.signal import correlate
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
    compute_sequence_diversity,
    kappa_curvature,
    sigmoid_aggregator,
    PHI,
    K_STAR,
    E_SQUARED,
)

warnings.filterwarnings('ignore', category=UserWarning)


class Hypothesis2Falsifier:
    """
    Falsification engine for κ(n) Z-invariant hypothesis.
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
    
    def compute_z_score(self, sequence: str) -> float:
        """
        Compute Z-invariant score for a DNA sequence.
        
        Z-score incorporates:
        - Spectral entropy (from phase-weighted FFT)
        - Curvature weighting κ(n) = d(n) · ln(n+1) / e²
        - Sigmoid aggregation
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Z-score (0-1 range via sigmoid)
        """
        # Encode and phase weight
        encoded = encode_complex(sequence, is_rna=False)
        phased = apply_phase_weighting(encoded, k=self.k)
        
        # Compute spectrum
        spectrum = compute_spectrum(phased)
        
        # Spectral entropy
        entropy = compute_spectral_entropy(spectrum)
        
        # Sequence diversity
        diversity = compute_sequence_diversity(sequence)
        
        # Curvature weight
        n = len(sequence)
        kappa = kappa_curvature(n, diversity)
        
        # Z-score via sigmoid aggregator
        # Z = S(entropy / κ) where S is sigmoid
        z_score = sigmoid_aggregator(entropy / (kappa + 1e-10), kappa=kappa)
        
        return z_score
    
    def test_length_invariance(
        self,
        sequences_by_length: Dict[int, List[str]]
    ) -> Dict[str, any]:
        """
        Test invariance of Z-scores across different sequence lengths.
        
        Uses ANOVA to test if Z-score means differ significantly by length.
        
        Args:
            sequences_by_length: Dictionary mapping length -> list of sequences
            
        Returns:
            Dictionary with ANOVA results
        """
        z_scores_by_length = {}
        
        for length, seqs in sequences_by_length.items():
            z_scores = [self.compute_z_score(seq) for seq in seqs]
            z_scores_by_length[length] = z_scores
        
        # Prepare for ANOVA
        groups = list(z_scores_by_length.values())
        
        # Check if we have enough data
        if len(groups) < 2:
            return {
                "f_statistic": None,
                "p_value": None,
                "invariant": None,
                "reason": "Need at least 2 length groups"
            }
        
        # ANOVA test
        f_stat, p_value = f_oneway(*groups)
        
        # Compute variance across length groups
        all_z_scores = []
        for z_list in groups:
            all_z_scores.extend(z_list)
        
        overall_variance = np.var(all_z_scores)
        
        # Mean Z-scores per length
        mean_z_by_length = {
            length: np.mean(scores) 
            for length, scores in z_scores_by_length.items()
        }
        
        # Invariance criterion: p >= 0.05 (no significant length dependence)
        invariant = p_value >= 0.05
        
        return {
            "f_statistic": float(f_stat),
            "p_value": float(p_value),
            "invariant": bool(invariant),
            "overall_variance": float(overall_variance),
            "mean_z_by_length": mean_z_by_length,
            "z_scores_by_length": z_scores_by_length,
        }
    
    def test_periodicity(self, sequences: List[str]) -> Dict[str, any]:
        """
        Test for emergent periodicity in Z-scores.
        
        Computes autocorrelation of Z-scores to detect periodic patterns.
        
        Args:
            sequences: List of DNA sequences (should be sorted or ordered)
            
        Returns:
            Dictionary with periodicity analysis results
        """
        # Compute Z-scores
        z_scores = np.array([self.compute_z_score(seq) for seq in sequences])
        
        # Autocorrelation
        z_centered = z_scores - np.mean(z_scores)
        autocorr = correlate(z_centered, z_centered, mode='full')
        autocorr = autocorr[len(autocorr) // 2:]  # Keep only positive lags
        
        # Normalize
        autocorr = autocorr / autocorr[0] if autocorr[0] != 0 else autocorr
        
        # Find peaks in autocorrelation (excluding lag 0)
        if len(autocorr) > 1:
            # Look for first significant peak
            max_autocorr = np.max(autocorr[1:]) if len(autocorr) > 1 else 0.0
            
            # Detect periodicity: max autocorr > threshold
            periodicity_detected = max_autocorr > 0.1
        else:
            max_autocorr = 0.0
            periodicity_detected = False
        
        return {
            "max_autocorrelation": float(max_autocorr),
            "periodicity_detected": bool(periodicity_detected),
            "autocorrelation_values": autocorr[:20].tolist() if len(autocorr) >= 20 else autocorr.tolist(),
            "z_scores": z_scores.tolist(),
        }


def generate_variable_length_sequences(
    lengths: List[int],
    n_per_length: int = 20,
    seed: int = 42,
) -> Dict[int, List[str]]:
    """
    Generate random DNA sequences of variable lengths.
    
    Args:
        lengths: List of sequence lengths to generate
        n_per_length: Number of sequences per length
        seed: Random seed
        
    Returns:
        Dictionary mapping length -> list of sequences
    """
    np.random.seed(seed)
    
    bases = ['A', 'T', 'C', 'G']
    sequences_by_length = {}
    
    for length in lengths:
        sequences = []
        for _ in range(n_per_length):
            seq = ''.join(np.random.choice(bases, length))
            sequences.append(seq)
        sequences_by_length[length] = sequences
    
    return sequences_by_length


def generate_sequences_with_motif(
    base_motif: str,
    lengths: List[int],
    n_per_length: int = 20,
    seed: int = 42,
) -> Dict[int, List[str]]:
    """
    Generate sequences of different lengths while preserving a core motif.
    
    This tests if Z-scores remain invariant when sequence length changes
    but core sequence content is preserved.
    
    Args:
        base_motif: Core motif to preserve
        lengths: List of target lengths
        n_per_length: Number of sequences per length
        seed: Random seed
        
    Returns:
        Dictionary mapping length -> list of sequences
    """
    np.random.seed(seed)
    
    bases = ['A', 'T', 'C', 'G']
    sequences_by_length = {}
    
    motif_len = len(base_motif)
    
    for length in lengths:
        sequences = []
        
        if length < motif_len:
            # If target shorter than motif, truncate
            for _ in range(n_per_length):
                seq = base_motif[:length]
                sequences.append(seq)
        else:
            # Extend motif with random bases
            for _ in range(n_per_length):
                # Random prefix/suffix
                extra_bases = length - motif_len
                prefix_len = np.random.randint(0, extra_bases + 1)
                suffix_len = extra_bases - prefix_len
                
                prefix = ''.join(np.random.choice(bases, prefix_len))
                suffix = ''.join(np.random.choice(bases, suffix_len))
                
                seq = prefix + base_motif + suffix
                sequences.append(seq)
        
        sequences_by_length[length] = sequences
    
    return sequences_by_length


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Falsification Experiment for Hypothesis 2"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=20,
        help="Minimum sequence length"
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=100,
        help="Maximum sequence length"
    )
    parser.add_argument(
        "--length-step",
        type=int,
        default=10,
        help="Step size between lengths"
    )
    parser.add_argument(
        "--n-per-length",
        type=int,
        default=20,
        help="Number of sequences per length"
    )
    parser.add_argument(
        "--k-parameter",
        type=float,
        default=K_STAR,
        help="Resolution exponent k (default: 0.3)"
    )
    parser.add_argument(
        "--test-motif",
        action="store_true",
        help="Test with preserved motif across lengths"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/PR-156-falsify-hypothesis2",
        help="Output directory for results"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", args.output_dir
    )
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize falsifier
    falsifier = Hypothesis2Falsifier(seed=args.seed, k_parameter=args.k_parameter)
    
    print("=" * 70)
    print("Hypothesis 2 Falsification Experiment")
    print("=" * 70)
    print(f"Seed: {args.seed}")
    print(f"k parameter: {args.k_parameter}")
    print(f"Length range: {args.min_length}-{args.max_length} (step {args.length_step})")
    print(f"N per length: {args.n_per_length}")
    print(f"Test with motif: {args.test_motif}")
    print()
    
    # Generate length list
    lengths = list(range(args.min_length, args.max_length + 1, args.length_step))
    
    # Generate sequences
    print(f"Generating sequences for {len(lengths)} length groups...")
    
    if args.test_motif:
        # Use a conserved motif
        base_motif = "ATCGATCGATCGATCG"  # 16-bp motif
        sequences_by_length = generate_sequences_with_motif(
            base_motif, lengths, args.n_per_length, args.seed
        )
        print(f"Using preserved motif: {base_motif}")
    else:
        # Random sequences
        sequences_by_length = generate_variable_length_sequences(
            lengths, args.n_per_length, args.seed
        )
        print("Using random sequences")
    
    total_sequences = sum(len(seqs) for seqs in sequences_by_length.values())
    print(f"Generated {total_sequences} total sequences")
    print()
    
    # Test 1: Length invariance
    print("Testing Z-score invariance across sequence lengths...")
    invariance_results = falsifier.test_length_invariance(sequences_by_length)
    
    print(f"ANOVA F-statistic: {invariance_results['f_statistic']:.4f}")
    print(f"ANOVA p-value: {invariance_results['p_value']:.4f}")
    print(f"Overall Z-score variance: {invariance_results['overall_variance']:.4f}")
    print()
    
    print("Mean Z-scores by length:")
    for length, mean_z in sorted(invariance_results['mean_z_by_length'].items()):
        print(f"  Length {length}: {mean_z:.4f}")
    print()
    
    # Test 2: Periodicity
    print("Testing for emergent periodicity in Z-scores...")
    
    # Flatten sequences for periodicity test
    all_sequences = []
    for length in sorted(sequences_by_length.keys()):
        all_sequences.extend(sequences_by_length[length])
    
    periodicity_results = falsifier.test_periodicity(all_sequences)
    
    print(f"Max autocorrelation: {periodicity_results['max_autocorrelation']:.4f}")
    print(f"Periodicity detected: {periodicity_results['periodicity_detected']}")
    print()
    
    # Falsification decision
    print("=" * 70)
    print("FALSIFICATION ANALYSIS")
    print("=" * 70)
    
    falsified = False
    reasons = []
    
    # Criterion 1: Invariance
    if not invariance_results['invariant']:
        falsified = True
        reasons.append(
            f"Z-scores NOT invariant across lengths (ANOVA p={invariance_results['p_value']:.4f} < 0.05)"
        )
    
    # Criterion 2: Periodicity
    if not periodicity_results['periodicity_detected']:
        falsified = True
        reasons.append(
            f"No emergent periodicity detected (max autocorr={periodicity_results['max_autocorrelation']:.4f} < 0.1)"
        )
    
    if falsified:
        print("RESULT: Hypothesis 2 is FALSIFIED")
        print("Reasons:")
        for reason in reasons:
            print(f"  - {reason}")
    else:
        print("RESULT: Hypothesis 2 is NOT FALSIFIED (supported by data)")
        print(f"  - Z-scores are invariant (ANOVA p={invariance_results['p_value']:.4f} >= 0.05)")
        print(f"  - Emergent periodicity detected (max autocorr={periodicity_results['max_autocorrelation']:.4f} > 0.1)")
    
    print()
    
    # Save results
    results = {
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "seed": args.seed,
            "k_parameter": args.k_parameter,
            "min_length": args.min_length,
            "max_length": args.max_length,
            "length_step": args.length_step,
            "n_per_length": args.n_per_length,
            "test_motif": args.test_motif,
        },
        "invariance_test": {
            "f_statistic": invariance_results['f_statistic'],
            "p_value": invariance_results['p_value'],
            "invariant": invariance_results['invariant'],
            "overall_variance": invariance_results['overall_variance'],
            "mean_z_by_length": invariance_results['mean_z_by_length'],
        },
        "periodicity_test": {
            "max_autocorrelation": periodicity_results['max_autocorrelation'],
            "periodicity_detected": periodicity_results['periodicity_detected'],
        },
        "falsification": {
            "falsified": falsified,
            "reasons": reasons,
        }
    }
    
    output_file = os.path.join(output_dir, "hypothesis2_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to: {output_file}")
    print()
    
    return 0 if not falsified else 1


if __name__ == "__main__":
    sys.exit(main())
