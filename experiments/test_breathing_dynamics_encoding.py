#!/usr/bin/env python3
"""
Test: DNA Breathing Dynamics vs Arbitrary Encoding

Hypothesis: Electronic polarizability (breathing dynamics) provide better
spectral encoding than arbitrary weights because they reflect real oscillatory
phenomena that affect CRISPR accessibility.

Breathing Frequencies (experimental values):
- AT pairs: ~10^7 Hz (fast opening, weaker bonds)
- GC pairs: ~10^9 Hz (slow opening, stronger bonds)

Prediction: Breathing-based encoding will show significantly higher Z-scores
than arbitrary encodings, especially for mutations affecting GC content.
"""

import numpy as np
import random
from scipy import stats
from scipy.fft import fft
from scipy import signal
from typing import Dict, List, Tuple
import sys
import os
from Bio import SeqIO

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# Import DiscreteZetaShift implementation inline (to avoid complex dependencies)
E_SQUARED = np.e**2  # ≈ 7.389

class DiscreteZetaShift:
    """
    Simple implementation of Z Framework's universal equation Z = A(B/c)
    """
    def __init__(self, sequence_length: int):
        self.length = sequence_length
        self.c = E_SQUARED

    def compute_frame_entropy(self, sequence: str) -> float:
        """Calculate frame-dependent sequence entropy"""
        base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for base in sequence:
            if base in base_counts:
                base_counts[base] += 1

        total = sum(base_counts.values())
        if total == 0:
            return 0.0

        entropy = 0.0
        for count in base_counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        return entropy

    def compute_spectral_shift(self, original_spectrum: np.ndarray,
                               mutated_spectrum: np.ndarray) -> float:
        """Calculate spectral mutation shift"""
        diff = np.abs(mutated_spectrum - original_spectrum)
        return np.sum(diff)

    def compute_z_score(self, sequence: str, original_spectrum: np.ndarray,
                       mutated_spectrum: np.ndarray) -> float:
        """Compute Z = A(B/c) with normalization"""
        A = self.compute_frame_entropy(sequence)
        B_raw = self.compute_spectral_shift(original_spectrum, mutated_spectrum)

        # Normalize B by spectrum magnitude
        spectrum_magnitude = np.sum(np.abs(original_spectrum))
        B = B_raw / (spectrum_magnitude + 1e-10)

        Z = A * (B / self.c)
        return Z

# Experimental breathing frequencies (Hz)
# Dimensionless breathing rate ratios - AT opens 100x faster than GC (μs-ms timescales)
BREATHING_RATE_RATIO = {
    'A': 1,     # AT pair: fast opening (μs scale)
    'T': 1,     # AT pair: fast opening
    'C': 100,   # GC pair: slow opening (ms scale)
    'G': 100    # GC pair: slow opening
}

# For comparison: helical periodicity (structural)
HELICAL_PERIOD = 10.5  # bp per turn


class BreathingDynamicsEncoder:
    """Encode DNA using dimensionless breathing rate ratios""" 
    

    def __init__(self):
        """Initialize with experimentally-derived breathing frequencies"""
        # Normalize frequencies to suitable complex weight range
        # Log scale because frequencies span orders of magnitude
        self.weights = {}

        for base in 'ATCG':
            freq = BREATHING_RATE_RATIO[base]

            # Real part: log-normalized frequency (captures magnitude)
            # Map 10^7 to 10^9 Hz -> -10 to +10 range
            real_part = (np.log10(freq) - 8.0) * 10.0

            # Imaginary part: phase derived from opening kinetics
            # AT (weak) = positive phase, GC (strong) = negative phase
            if base in 'AT':
                imag_part = 3.0  # Fast opening = positive phase
            else:
                imag_part = -3.0  # Slow opening = negative phase

            self.weights[base] = real_part + 1j * imag_part

        print("Breathing Dynamics Encoding Weights:")
        for base in 'ATCG':
            w = self.weights[base]
            freq = BREATHING_RATE_RATIO[base]
            print(f"  {base}: {w.real:+.2f}{w.imag:+.2f}j (freq: {freq:.0e} Hz)")

    def encode_sequence(self, sequence: str) -> np.ndarray:
        """Encode sequence with breathing dynamics"""
        encoded = []

        for i, base in enumerate(sequence):
            if base not in self.weights:
                encoded.append(0 + 0j)
                continue

            # Base weight from breathing frequency
            base_weight = self.weights[base]

            # Add helical periodicity phase (DNA wraps every 10.5 bp)
            helical_phase = 2 * np.pi * i / HELICAL_PERIOD

            # Add positional phase (for context)
            positional_phase = 2 * np.pi * i / len(sequence)

            # Combined phase: helical structure + positional context
            total_phase = helical_phase + positional_phase * 0.3

            # Apply phase modulation
            encoded_base = base_weight * np.exp(1j * total_phase)
            encoded.append(encoded_base)

        return np.array(encoded)


class ArbitraryEncoder:
    """Arbitrary random encoding for control"""

    def __init__(self, seed: int = None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        # Random weights in similar magnitude range as breathing encoder
        self.weights = {
            'A': random.uniform(-10, 10) + random.uniform(-3, 3) * 1j,
            'T': random.uniform(-10, 10) + random.uniform(-3, 3) * 1j,
            'C': random.uniform(-10, 10) + random.uniform(-3, 3) * 1j,
            'G': random.uniform(-10, 10) + random.uniform(-3, 3) * 1j,
        }

    def encode_sequence(self, sequence: str) -> np.ndarray:
        """Encode with arbitrary weights"""
        encoded = []

        for i, base in enumerate(sequence):
            if base not in self.weights:
                encoded.append(0 + 0j)
                continue

            # Same phase modulation as breathing encoder (fair comparison)
            helical_phase = 2 * np.pi * i / HELICAL_PERIOD
            positional_phase = 2 * np.pi * i / len(sequence)
            total_phase = helical_phase + positional_phase * 0.3

            encoded_base = self.weights[base] * np.exp(1j * total_phase)
            encoded.append(encoded_base)

        return np.array(encoded)


class BreathingDynamicsValidator:
    """Test breathing dynamics encoding against arbitrary controls"""
    def __init__(self):
        self.encoder = BreathingDynamicsEncoder()

    def generate_crispr_sequences(self, n_sequences: int, seq_length: int = 20) -> List[str]:
        """Generate CRISPR guide sequences from real human cDNA"""
        records = list(SeqIO.parse('data/test_human_cdna.fasta', 'fasta'))
        sequences = []
       
        while len(sequences) < n_sequences:
            record = random.choice(records)
            full_seq = str(record.seq)
            if len(full_seq) < seq_length:
                continue
            start = random.randint(0, len(full_seq) - seq_length)
            seq = full_seq[start:start + seq_length].upper()
            seq = ''.join(c for c in seq if c in 'ATCG')
            if len(seq) == seq_length:
                sequences.append(seq)
        return sequences

    def compute_spectrum(self, encoded_seq: np.ndarray) -> np.ndarray:
        """Compute FFT spectrum"""
        return np.abs(fft(encoded_seq))

    def analyze_mutation_sensitivity(self, sequence: str, encoder,
                                      mutation_type: str = 'gc_affecting') -> float:
        """
        Analyze how mutations affect spectral signature

        mutation_type:
            'gc_affecting': Mutations that change GC content (most relevant for breathing)
            'at_affecting': Mutations that change AT content
            'random': Random single-point mutations
        """
        zeta_shift = DiscreteZetaShift(len(sequence))

        # Original spectrum
        encoded_orig = encoder.encode_sequence(sequence)
        spectrum_orig = self.compute_spectrum(encoded_orig)

        z_scores = []

        # Create targeted mutations
        for pos in range(len(sequence)):
            original_base = sequence[pos]

            # Choose mutation based on type
            if mutation_type == 'gc_affecting':
                # AT -> GC or GC -> AT (changes breathing dynamics significantly)
                if original_base in 'AT':
                    new_bases = ['G', 'C']
                else:
                    new_bases = ['A', 'T']
            elif mutation_type == 'at_affecting':
                # Within AT or within GC (smaller effect on breathing)
                if original_base in 'AT':
                    new_bases = ['A' if original_base == 'T' else 'T']
                else:
                    new_bases = ['C' if original_base == 'G' else 'G']
            else:  # random
                new_bases = [b for b in 'ATCG' if b != original_base]

            for new_base in new_bases[:1]:  # Test one mutation per position
                mutated_seq = sequence[:pos] + new_base + sequence[pos+1:]

                # Mutated spectrum
                encoded_mut = encoder.encode_sequence(mutated_seq)
                spectrum_mut = self.compute_spectrum(encoded_mut)

                # Compute Z score
                z_score = zeta_shift.compute_z_score(sequence, spectrum_orig, spectrum_mut)
                z_scores.append(z_score)

        return np.mean(z_scores) if z_scores else 0.0

    def run_comparative_test(self, n_sequences: int = 100,
                            n_arbitrary_trials: int = 10) -> Dict:
        """Run breathing dynamics vs arbitrary encoding test"""

        print("\n" + "="*70)
        print("BREATHING DYNAMICS ENCODING TEST")
        print("="*70)

        print(f"\nGenerating {n_sequences} CRISPR-like test sequences...")
        sequences = self.generate_crispr_sequences(n_sequences)

        print("\n--- Testing Breathing Dynamics Encoding ---")
        breathing_scores_gc = []
        breathing_scores_at = []
        breathing_scores_random = []

        for seq in sequences:
            z_gc = self.analyze_mutation_sensitivity(seq, self.encoder,
                                                     'gc_affecting')
            z_at = self.analyze_mutation_sensitivity(seq, self.encoder,
                                                     'at_affecting')
            z_rand = self.analyze_mutation_sensitivity(seq, self.encoder,
                                                       'random')

            breathing_scores_gc.append(z_gc)
            breathing_scores_at.append(z_at)
            breathing_scores_random.append(z_rand)

        print(f"  GC-affecting mutations: Mean Z = {np.mean(breathing_scores_gc):.4f}")
        print(f"  AT-affecting mutations: Mean Z = {np.mean(breathing_scores_at):.4f}")
        print(f"  Random mutations: Mean Z = {np.mean(breathing_scores_random):.4f}")

        print(f"\n--- Testing Arbitrary Encodings ({n_arbitrary_trials} trials) ---")
        arbitrary_results = []

        for trial in range(n_arbitrary_trials):
            arb_encoder = ArbitraryEncoder(seed=trial)

            arb_scores_gc = []
            arb_scores_at = []
            arb_scores_random = []

            for seq in sequences:
                z_gc = self.analyze_mutation_sensitivity(seq, arb_encoder, 'gc_affecting')
                z_at = self.analyze_mutation_sensitivity(seq, arb_encoder, 'at_affecting')
                z_rand = self.analyze_mutation_sensitivity(seq, arb_encoder, 'random')

                arb_scores_gc.append(z_gc)
                arb_scores_at.append(z_at)
                arb_scores_random.append(z_rand)

            arbitrary_results.append({
                'gc_mean': np.mean(arb_scores_gc),
                'at_mean': np.mean(arb_scores_at),
                'random_mean': np.mean(arb_scores_random),
                'gc_scores': arb_scores_gc,
                'at_scores': arb_scores_at,
                'random_scores': arb_scores_random
            })

        # Aggregate arbitrary results
        arb_gc_means = [r['gc_mean'] for r in arbitrary_results]
        arb_at_means = [r['at_mean'] for r in arbitrary_results]
        arb_random_means = [r['random_mean'] for r in arbitrary_results]

        print(f"  GC-affecting: Mean Z = {np.mean(arb_gc_means):.4f} ± {np.std(arb_gc_means):.4f}")
        print(f"  AT-affecting: Mean Z = {np.mean(arb_at_means):.4f} ± {np.std(arb_at_means):.4f}")
        print(f"  Random: Mean Z = {np.mean(arb_random_means):.4f} ± {np.std(arb_random_means):.4f}")

        # Statistical comparison
        print("\n" + "="*70)
        print("STATISTICAL ANALYSIS")
        print("="*70)

        results = {}

        for mut_type, breath_scores, arb_means_list in [
            ('GC-affecting', breathing_scores_gc, arb_gc_means),
            ('AT-affecting', breathing_scores_at, arb_at_means),
            ('Random', breathing_scores_random, arb_random_means)
        ]:
            # t-test
            t_stat, p_value = stats.ttest_ind(breath_scores, arb_means_list)

            # Effect size (Cohen's d)
            pooled_std = max(np.sqrt(
                ((len(breath_scores) - 1) * np.var(breath_scores, ddof=1) +
                 (len(arb_means_list) - 1) * np.var(arb_means_list, ddof=1)) /
                (len(breath_scores) + len(arb_means_list) - 2)
            ), 1e-8)
            cohens_d = (np.mean(breath_scores) - np.mean(arb_means_list)) / pooled_std

            breath_better = np.mean(breath_scores) > np.mean(arb_means_list)

            print(f"\n{mut_type} Mutations:")
            print(f"  Breathing Mean Z: {np.mean(breath_scores):.4f}")
            print(f"  Arbitrary Mean Z: {np.mean(arb_means_list):.4f}")
            print(f"  Difference: {np.mean(breath_scores) - np.mean(arb_means_list):+.4f}")
            print(f"  t-statistic: {t_stat:.4f}")
            print(f"  p-value: {p_value:.6f}")
            print(f"  Cohen's d: {cohens_d:+.4f}")
            print(f"  Significant (p<0.05): {'YES' if p_value < 0.05 else 'NO'}")
            print(f"  Winner: {'BREATHING' if breath_better else 'ARBITRARY'}")

            results[mut_type] = {
                'breathing_mean': np.mean(breath_scores),
                'arbitrary_mean': np.mean(arb_means_list),
                'difference': np.mean(breath_scores) - np.mean(arb_means_list),
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'significant': p_value < 0.05,
                'breathing_wins': breath_better
            }

        return results

    def interpret_results(self, results: Dict) -> None:
        """Interpret and print conclusions"""
        print("\n" + "="*70)
        print("INTERPRETATION")
        print("="*70)

        # Check if breathing wins on GC-affecting mutations (key prediction)
        gc_result = results['GC-affecting']

        if gc_result['breathing_wins'] and gc_result['significant']:
            print("\n✓ HYPOTHESIS CONFIRMED")
            print("  Breathing dynamics encoding OUTPERFORMS arbitrary encoding")
            print("  for GC-affecting mutations (as predicted).")
            print(f"  Effect size: Cohen's d = {gc_result['cohens_d']:+.3f}")

            if abs(gc_result['cohens_d']) > 0.8:
                print("  This is a LARGE effect size!")
            elif abs(gc_result['cohens_d']) > 0.5:
                print("  This is a MEDIUM effect size.")

            print("\n  IMPLICATION: DNA breathing frequencies (real oscillatory")
            print("  phenomena) translate better to spectral encoding than")
            print("  arbitrary weights. Frequency-native properties work!")

        elif not gc_result['breathing_wins'] and gc_result['significant']:
            print("\n✗ HYPOTHESIS REJECTED")
            print("  Arbitrary encoding OUTPERFORMS breathing dynamics")
            print("  (unexpected result).")
            print(f"  Effect size: Cohen's d = {gc_result['cohens_d']:+.3f}")

            print("\n  IMPLICATION: Even frequency-native properties like")
            print("  breathing dynamics don't improve spectral encoding.")
            print("  May need different frequency properties or different")
            print("  mathematical framework.")

        else:
            print("\n⚠ NO SIGNIFICANT DIFFERENCE")
            print("  Breathing dynamics encoding shows no significant")
            print("  advantage over arbitrary encoding.")
            print(f"  p-value: {gc_result['p_value']:.4f} (not < 0.05)")

            print("\n  IMPLICATION: Either breathing frequencies don't affect")
            print("  CRISPR function enough, or the Z Framework doesn't")
            print("  capture their contribution effectively.")

        # Summary table
        print("\n" + "-"*70)
        print("SUMMARY TABLE")
        print("-"*70)
        print(f"{'Mutation Type':<20} {'Winner':<12} {'p-value':<10} {'Effect Size'}")
        print("-"*70)

        for mut_type, data in results.items():
            winner = 'Breathing' if data['breathing_wins'] else 'Arbitrary'
            significance = '*' if data['significant'] else ''
            print(f"{mut_type:<20} {winner:<12} {data['p_value']:<10.6f} {data['cohens_d']:+.3f}{significance}")

        print("-"*70)
        print("* = statistically significant (p < 0.05)")


def main():
    """Run breathing dynamics encoding test"""

    # Set random seed for reproducibility
    random.seed(42)
    np.random.seed(42)

    print("\n" + "="*70)
    print("DNA BREATHING DYNAMICS vs ARBITRARY ENCODING")
    print("="*70)
    print("\nHypothesis: Base pair opening rates provide better spectral")
    print("encoding than arbitrary weights because they reflect real oscillatory")
    print("phenomena affecting CRISPR accessibility.")
    print("\nExperimental basis:")
    print("  - AT pairs open at ~10 MHz (2 H-bonds)")
    print("  - GC pairs open at ~1 GHz (3 H-bonds)")
    print("  - 100× frequency difference")
    print("\nPrediction: Breathing encoding will show significantly higher Z-scores")
    print("than arbitrary encodings, especially for GC-content-affecting mutations.")

    # Run test
    validator = BreathingDynamicsValidator()
    results = validator.run_comparative_test(n_sequences=10, n_arbitrary_trials=2)

    # Interpret
    validator.interpret_results(results)

    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70)

    return results


if __name__ == '__main__':
    results = main()
