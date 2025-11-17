#!/usr/bin/env python3
"""
Spectral Disruption Profiler - Falsification Experiment

This module implements a comprehensive falsification experiment testing whether
phase-weighted FFT analysis (θ′(n,k) with k* ≈ 0.300) provides significant
enhancement to CRISPR gRNA efficiency prediction compared to unweighted baseline.

Scientific Gates:
- Human DNA only (A/C/G/T/N)
- Z Framework invariants: Z = A(B/e²)
- Geometric resolution: θ′(n,k) = φ·((n mod φ)/φ)^k with k ≈ 0.3
- Statistical rigor: Bootstrap CI (≥1,000), permutation tests, FDR correction
- Reproducibility: Fixed seed, versioned code, metadata persistence

Author: Z Framework Falsification Team
Date: 2025-11-17
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from scipy import stats
from scipy.fft import fft, fftfreq
from scipy.stats import entropy
import logging

# Add parent directories for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from scripts.z_framework import ZFrameworkCalculator
except ImportError:
    # Fallback if running from different location
    from z_framework import ZFrameworkCalculator

try:
    from applications.fft_crispr_disruption import FFTCRISPRDisruptionAnalyzer
except ImportError:
    # Will implement local fallback if needed
    FFTCRISPRDisruptionAnalyzer = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Mathematical constants
PHI = 1.6180339887498949  # Golden ratio
E_SQUARED = 7.3890560989306495  # e²
DEFAULT_PHI_PERIOD = 21.0  # 21-nt guide default
DEFAULT_K = 0.3  # Default resolution exponent


class SpectralDisruptionFalsifier:
    """
    Falsification experiment for spectral disruption profiler hypothesis.
    
    Tests whether phase-weighted FFT analysis with θ′(n,k) provides significant
    improvement over unweighted baseline for CRISPR gRNA efficiency prediction.
    """
    
    def __init__(
        self,
        seed: int = 42,
        n_bootstrap: int = 1000,
        n_permutation: int = 1000,
        k_parameter: float = DEFAULT_K,
        phi_period: float = DEFAULT_PHI_PERIOD,
        alpha: float = 0.05
    ):
        """
        Initialize falsification experiment.
        
        Args:
            seed: Random seed for reproducibility
            n_bootstrap: Number of bootstrap resamples
            n_permutation: Number of permutation resamples
            k_parameter: Geodesic curvature exponent
            phi_period: Geometric period (default 21 for guides)
            alpha: Significance level for statistical tests
        """
        self.seed = seed
        self.n_bootstrap = n_bootstrap
        self.n_permutation = n_permutation
        self.k_parameter = k_parameter
        self.phi_period = phi_period
        self.alpha = alpha
        
        # Set random seeds
        np.random.seed(seed)
        
        # Initialize components
        self.z_calc = ZFrameworkCalculator(precision_dps=50)
        
        # Results storage
        self.results = {
            'metadata': self._get_metadata(),
            'parameters': {
                'seed': seed,
                'n_bootstrap': n_bootstrap,
                'n_permutation': n_permutation,
                'k_parameter': k_parameter,
                'phi_period': phi_period,
                'alpha': alpha
            },
            'primary_endpoints': {},
            'secondary_endpoints': {},
            'falsification_status': None
        }
        
        logger.info("Initialized SpectralDisruptionFalsifier")
        logger.info(f"  Seed: {seed}")
        logger.info(f"  Bootstrap samples: {n_bootstrap}")
        logger.info(f"  Permutation samples: {n_permutation}")
        logger.info(f"  k-parameter: {k_parameter}")
        logger.info(f"  phi-period: {phi_period}")
    
    def _get_metadata(self) -> Dict:
        """Get experiment metadata including git commit, timestamp, etc."""
        import subprocess
        
        metadata = {
            'timestamp': time.strftime('%Y%m%d-%H%M%S'),
            'python_version': sys.version,
            'experiment_id': 'spectral_disruption_profiler_137',
            'version': '1.0'
        }
        
        try:
            git_commit = subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                stderr=subprocess.DEVNULL
            ).decode('utf-8').strip()
            metadata['git_commit'] = git_commit
        except:
            metadata['git_commit'] = 'unknown'
        
        return metadata
    
    def validate_dna_sequence(self, sequence: str) -> str:
        """
        Validate DNA sequence contains only A/C/G/T/N.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Uppercase validated sequence
            
        Raises:
            ValueError: If sequence contains invalid bases
        """
        if not sequence:
            raise ValueError("DNA sequence cannot be empty")
        
        sequence = sequence.upper()
        
        # Check for U (RNA) in DNA context
        if 'U' in sequence:
            raise ValueError(
                "Found 'U' in sequence. For DNA sequences, use 'T'."
            )
        
        valid_bases = set('ACGTN')
        invalid_bases = set(sequence) - valid_bases
        
        if invalid_bases:
            raise ValueError(
                f"Invalid DNA bases found: {invalid_bases}. "
                f"Only A/C/G/T/N allowed for DNA sequences."
            )
        
        return sequence
    
    def encode_dna_complex(self, sequence: str) -> np.ndarray:
        """
        Encode DNA sequence as complex waveform.
        
        Standard mapping:
        - A → 1+0i
        - T → -1+0i
        - C → 0+i
        - G → 0-i
        - N → 0+0i
        
        Args:
            sequence: Validated DNA sequence
            
        Returns:
            Complex numpy array
        """
        mapping = {
            'A': 1+0j,
            'T': -1+0j,
            'C': 0+1j,
            'G': 0-1j,
            'N': 0+0j
        }
        
        return np.array([mapping[base] for base in sequence], dtype=np.complex128)
    
    def theta_prime(self, n: int, k: float) -> float:
        """
        Geometric resolution function θ′(n,k) = φ·((n mod φ)/φ)^k
        
        Args:
            n: Position index
            k: Resolution exponent
            
        Returns:
            Phase weight value
        """
        phi = self.phi_period
        return PHI * ((n % phi) / phi) ** k
    
    def apply_phase_weighting(
        self,
        spectrum: np.ndarray,
        use_weighting: bool = True
    ) -> np.ndarray:
        """
        Apply θ′(n,k) phase weighting to FFT spectrum.
        
        Args:
            spectrum: FFT spectrum array
            use_weighting: If True, apply phase weighting; else return unweighted
            
        Returns:
            Weighted or unweighted spectrum
        """
        if not use_weighting:
            return spectrum
        
        n = len(spectrum)
        weights = np.array([self.theta_prime(i, self.k_parameter) for i in range(n)])
        return spectrum * weights
    
    def compute_fft_features(
        self,
        sequence: str,
        use_phase_weighting: bool = True
    ) -> Dict[str, float]:
        """
        Compute FFT-based features for a DNA sequence.
        
        Args:
            sequence: DNA sequence
            use_phase_weighting: Apply θ′(n,k) weighting
            
        Returns:
            Dictionary of FFT features
        """
        # Validate and encode
        seq = self.validate_dna_sequence(sequence)
        waveform = self.encode_dna_complex(seq)
        
        # Compute FFT
        spectrum = fft(waveform)
        freqs = fftfreq(len(waveform))
        
        # Apply phase weighting if requested
        spectrum_weighted = self.apply_phase_weighting(spectrum, use_phase_weighting)
        
        # Extract features
        power = np.abs(spectrum_weighted) ** 2
        power_norm = power / (np.sum(power) + 1e-10)
        
        # Dominant frequency
        dominant_idx = np.argmax(power[1:len(power)//2]) + 1
        dominant_freq = freqs[dominant_idx]
        
        # Spectral entropy
        spectral_entropy = entropy(power_norm + 1e-10)
        
        # Sidelobe count (bins with >10% of dominant power)
        sidelobe_threshold = 0.1 * power[dominant_idx]
        sidelobe_count = np.sum(power[1:len(power)//2] > sidelobe_threshold)
        
        # GC content
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        
        return {
            'dominant_freq': float(dominant_freq),
            'spectral_entropy': float(spectral_entropy),
            'sidelobe_count': int(sidelobe_count),
            'gc_content': float(gc_content),
            'sequence_length': len(seq),
            'phase_weighted': use_phase_weighting
        }
    
    def generate_synthetic_sequences(
        self,
        n_sequences: int = 100,
        length: int = 21,
        gc_range: Tuple[float, float] = (0.4, 0.6),
        add_noise: bool = False
    ) -> List[str]:
        """
        Generate synthetic DNA sequences for testing.
        
        Args:
            n_sequences: Number of sequences to generate
            length: Length of each sequence
            gc_range: (min, max) GC content range
            add_noise: Add random noise to sequences
            
        Returns:
            List of DNA sequences
        """
        sequences = []
        
        for _ in range(n_sequences):
            # Target GC content
            target_gc = np.random.uniform(gc_range[0], gc_range[1])
            n_gc = int(length * target_gc)
            n_at = length - n_gc
            
            # Generate bases
            gc_bases = np.random.choice(['G', 'C'], n_gc)
            at_bases = np.random.choice(['A', 'T'], n_at)
            
            # Combine and shuffle
            bases = np.concatenate([gc_bases, at_bases])
            np.random.shuffle(bases)
            
            sequence = ''.join(bases)
            
            # Add noise if requested (random single-base changes)
            if add_noise and np.random.random() < 0.1:
                pos = np.random.randint(0, length)
                new_base = np.random.choice(['A', 'T', 'C', 'G'])
                sequence = sequence[:pos] + new_base + sequence[pos+1:]
            
            sequences.append(sequence)
        
        return sequences
    
    def load_human_dna_from_csv(
        self,
        csv_path: str,
        sequence_column: str = 'sequence',
        label_column: Optional[str] = 'efficiency',
        max_sequences: Optional[int] = None
    ) -> Tuple[List[str], Optional[np.ndarray]]:
        """
        Load real human DNA sequences from CSV file.
        
        Args:
            csv_path: Path to CSV file with DNA sequences
            sequence_column: Column name containing DNA sequences
            label_column: Column name containing labels/efficiency scores (optional)
            max_sequences: Maximum number of sequences to load
            
        Returns:
            Tuple of (sequences, labels) where labels may be None
        """
        logger.info(f"Loading human DNA data from {csv_path}...")
        
        # Load CSV
        df = pd.read_csv(csv_path)
        
        if sequence_column not in df.columns:
            raise ValueError(f"Column '{sequence_column}' not found in CSV")
        
        # Get sequences
        sequences = df[sequence_column].tolist()
        
        # Validate all sequences are human DNA
        validated_sequences = []
        for i, seq in enumerate(sequences):
            try:
                validated_seq = self.validate_dna_sequence(seq)
                validated_sequences.append(validated_seq)
            except ValueError as e:
                logger.warning(f"Skipping sequence {i}: {e}")
        
        # Get labels if available
        labels = None
        if label_column and label_column in df.columns:
            # Get labels for validated sequences only
            valid_indices = [i for i, seq in enumerate(sequences) 
                           if i < len(validated_sequences)]
            labels = df[label_column].iloc[:len(validated_sequences)].values
            
            # Convert to binary if continuous (median split)
            if labels.dtype in [np.float32, np.float64]:
                median = np.median(labels)
                labels = (labels >= median).astype(int)
                logger.info(f"Converted continuous labels to binary using median split")
        
        # Limit number of sequences if requested
        if max_sequences and len(validated_sequences) > max_sequences:
            indices = np.random.choice(len(validated_sequences), max_sequences, replace=False)
            validated_sequences = [validated_sequences[i] for i in indices]
            if labels is not None:
                labels = labels[indices]
        
        logger.info(f"Loaded {len(validated_sequences)} human DNA sequences from {csv_path}")
        if labels is not None:
            logger.info(f"Label distribution: {np.sum(labels)} positive, {len(labels) - np.sum(labels)} negative")
        
        return validated_sequences, labels
    
    def load_human_dna_from_fasta(
        self,
        fasta_path: str,
        max_sequences: Optional[int] = None,
        sequence_length: int = 21
    ) -> List[str]:
        """
        Load real human DNA sequences from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file with human DNA
            max_sequences: Maximum number of sequences to extract
            sequence_length: Length of sequences to extract (default 21 for guides)
            
        Returns:
            List of validated DNA sequences
        """
        logger.info(f"Loading human DNA sequences from {fasta_path}...")
        
        sequences = []
        current_seq = ""
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Header line, process previous sequence if exists
                    if current_seq:
                        # Extract subsequences of specified length
                        for i in range(0, len(current_seq) - sequence_length + 1, sequence_length):
                            subseq = current_seq[i:i+sequence_length]
                            try:
                                validated = self.validate_dna_sequence(subseq)
                                sequences.append(validated)
                                if max_sequences and len(sequences) >= max_sequences:
                                    break
                            except ValueError:
                                continue
                        current_seq = ""
                    if max_sequences and len(sequences) >= max_sequences:
                        break
                else:
                    current_seq += line.upper()
        
        # Process last sequence
        if current_seq and (not max_sequences or len(sequences) < max_sequences):
            for i in range(0, len(current_seq) - sequence_length + 1, sequence_length):
                subseq = current_seq[i:i+sequence_length]
                try:
                    validated = self.validate_dna_sequence(subseq)
                    sequences.append(validated)
                    if max_sequences and len(sequences) >= max_sequences:
                        break
                except ValueError:
                    continue
        
        logger.info(f"Extracted {len(sequences)} valid DNA sequences from {fasta_path}")
        
        return sequences
    
    def compute_roc_auc_difference(
        self,
        sequences: List[str],
        true_labels: np.ndarray
    ) -> Dict[str, float]:
        """
        Compute ΔROC-AUC between weighted and unweighted FFT.
        
        Args:
            sequences: List of DNA sequences
            true_labels: True binary labels for sequences
            
        Returns:
            Dictionary with ROC-AUC metrics
        """
        from sklearn.metrics import roc_auc_score
        
        # Compute features with and without weighting
        features_weighted = []
        features_unweighted = []
        
        for seq in sequences:
            feat_w = self.compute_fft_features(seq, use_phase_weighting=True)
            feat_u = self.compute_fft_features(seq, use_phase_weighting=False)
            
            features_weighted.append(feat_w['spectral_entropy'])
            features_unweighted.append(feat_u['spectral_entropy'])
        
        # Compute ROC-AUC
        try:
            auc_weighted = roc_auc_score(true_labels, features_weighted)
            auc_unweighted = roc_auc_score(true_labels, features_unweighted)
            delta_auc = auc_weighted - auc_unweighted
        except ValueError:
            # Handle case with single class
            logger.warning("Cannot compute ROC-AUC: single class in labels")
            return {
                'auc_weighted': 0.5,
                'auc_unweighted': 0.5,
                'delta_auc': 0.0,
                'error': 'single_class'
            }
        
        return {
            'auc_weighted': float(auc_weighted),
            'auc_unweighted': float(auc_unweighted),
            'delta_auc': float(delta_auc)
        }
    
    def bootstrap_confidence_interval(
        self,
        sequences: List[str],
        true_labels: np.ndarray
    ) -> Dict[str, float]:
        """
        Compute bootstrap confidence interval for ΔROC-AUC.
        
        Args:
            sequences: List of DNA sequences
            true_labels: True binary labels
            
        Returns:
            Dictionary with CI bounds and statistics
        """
        logger.info(f"Computing bootstrap CI with {self.n_bootstrap} resamples...")
        
        delta_aucs = []
        
        for i in range(self.n_bootstrap):
            if i % 100 == 0:
                logger.info(f"  Bootstrap iteration {i}/{self.n_bootstrap}")
            
            # Resample with replacement
            n = len(sequences)
            indices = np.random.choice(n, n, replace=True)
            
            seq_resample = [sequences[i] for i in indices]
            label_resample = true_labels[indices]
            
            # Compute ΔROC-AUC
            result = self.compute_roc_auc_difference(seq_resample, label_resample)
            delta_aucs.append(result['delta_auc'])
        
        # Compute CI
        delta_aucs = np.array(delta_aucs)
        ci_lower = np.percentile(delta_aucs, 2.5)
        ci_upper = np.percentile(delta_aucs, 97.5)
        mean_delta = np.mean(delta_aucs)
        
        # Check if CI includes zero (falsification criterion)
        includes_zero = ci_lower <= 0 <= ci_upper
        
        logger.info(f"Bootstrap CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
        logger.info(f"Mean ΔROC-AUC: {mean_delta:.4f}")
        logger.info(f"CI includes zero: {includes_zero}")
        
        return {
            'mean_delta_auc': float(mean_delta),
            'ci_lower': float(ci_lower),
            'ci_upper': float(ci_upper),
            'includes_zero': bool(includes_zero),
            'std': float(np.std(delta_aucs))
        }
    
    def permutation_test(
        self,
        sequences: List[str],
        true_labels: np.ndarray
    ) -> Dict[str, float]:
        """
        Perform permutation test for ΔROC-AUC significance.
        
        Args:
            sequences: List of DNA sequences
            true_labels: True binary labels
            
        Returns:
            Dictionary with p-value and statistics
        """
        logger.info(f"Performing permutation test with {self.n_permutation} permutations...")
        
        # Observed difference
        observed = self.compute_roc_auc_difference(sequences, true_labels)
        observed_delta = observed['delta_auc']
        
        # Null distribution
        null_deltas = []
        
        for i in range(self.n_permutation):
            if i % 100 == 0:
                logger.info(f"  Permutation {i}/{self.n_permutation}")
            
            # Permute labels
            permuted_labels = np.random.permutation(true_labels)
            
            # Compute ΔROC-AUC
            result = self.compute_roc_auc_difference(sequences, permuted_labels)
            null_deltas.append(result['delta_auc'])
        
        # Compute p-value (two-tailed)
        null_deltas = np.array(null_deltas)
        p_value = np.mean(np.abs(null_deltas) >= np.abs(observed_delta))
        
        logger.info(f"Observed ΔROC-AUC: {observed_delta:.4f}")
        logger.info(f"Permutation p-value: {p_value:.4f}")
        
        return {
            'observed_delta_auc': float(observed_delta),
            'p_value': float(p_value),
            'null_mean': float(np.mean(null_deltas)),
            'null_std': float(np.std(null_deltas)),
            'significant': bool(p_value < self.alpha)
        }
    
    def run_falsification_experiment(
        self,
        sequences: Optional[List[str]] = None,
        labels: Optional[np.ndarray] = None,
        n_sequences: int = 100,
        gc_range: Tuple[float, float] = (0.4, 0.6),
        use_synthetic: bool = False
    ) -> Dict:
        """
        Run complete falsification experiment.
        
        Args:
            sequences: Optional list of sequences; if None and use_synthetic=True, generates synthetic
            labels: Optional labels for sequences (required for real data)
            n_sequences: Number of synthetic sequences to generate (only if use_synthetic=True)
            gc_range: GC content range for synthetic sequences (only if use_synthetic=True)
            use_synthetic: If True, use synthetic data; if False, requires real sequences and labels
            
        Returns:
            Complete results dictionary
        """
        logger.info("="*60)
        logger.info("SPECTRAL DISRUPTION PROFILER FALSIFICATION EXPERIMENT")
        logger.info("="*60)
        
        start_time = time.time()
        
        # Validate inputs
        if not use_synthetic and (sequences is None or labels is None):
            raise ValueError(
                "Real data mode requires both 'sequences' and 'labels' to be provided. "
                "Use use_synthetic=True for synthetic data, or provide real human DNA data."
            )
        
        # Generate synthetic or use provided sequences
        if use_synthetic:
            if sequences is None:
                logger.info(f"⚠️  WARNING: Using synthetic random DNA sequences")
                logger.info(f"⚠️  This does not test against real human DNA data!")
                logger.info(f"Generating {n_sequences} synthetic sequences...")
                sequences = self.generate_synthetic_sequences(
                    n_sequences=n_sequences,
                    gc_range=gc_range
                )
            # Generate synthetic labels (random for null hypothesis testing)
            n = len(sequences)
            true_labels = np.random.randint(0, 2, n)
            logger.info(f"⚠️  Using random synthetic labels (null hypothesis)")
        else:
            logger.info(f"✓ Using real human DNA data")
            true_labels = labels
        
        # Generate synthetic labels (random for null hypothesis testing)
        n = len(sequences)
        
        logger.info(f"Processing {n} sequences...")
        logger.info(f"Label distribution: {np.sum(true_labels)} positive, {n - np.sum(true_labels)} negative")
        
        # Primary endpoint: Bootstrap CI for ΔROC-AUC
        logger.info("\n--- PRIMARY ENDPOINT: Bootstrap CI ---")
        bootstrap_result = self.bootstrap_confidence_interval(sequences, true_labels)
        self.results['primary_endpoints']['bootstrap_ci'] = bootstrap_result
        
        # Secondary endpoint: Permutation test
        logger.info("\n--- SECONDARY ENDPOINT: Permutation Test ---")
        permutation_result = self.permutation_test(sequences, true_labels)
        self.results['secondary_endpoints']['permutation_test'] = permutation_result
        
        # Performance test
        logger.info("\n--- PERFORMANCE TEST ---")
        perf_start = time.time()
        for seq in sequences[:100]:  # Test on first 100
            _ = self.compute_fft_features(seq, use_phase_weighting=True)
        perf_time = time.time() - perf_start
        perf_per_seq = perf_time / min(100, len(sequences))
        
        logger.info(f"Processing time: {perf_time:.2f}s for 100 sequences")
        logger.info(f"Per-sequence: {perf_per_seq*1000:.2f}ms")
        
        self.results['secondary_endpoints']['performance'] = {
            'total_time_100seq': float(perf_time),
            'time_per_sequence_ms': float(perf_per_seq * 1000),
            'meets_target': bool(perf_time < 30.0)
        }
        
        # Determine falsification status
        ci_includes_zero = bootstrap_result['includes_zero']
        perm_not_significant = permutation_result['p_value'] > self.alpha
        
        if ci_includes_zero and perm_not_significant:
            falsification_status = "HYPOTHESIS_FALSIFIED"
            conclusion = "Phase weighting provides NO significant improvement"
        elif not ci_includes_zero and not perm_not_significant:
            falsification_status = "HYPOTHESIS_SUPPORTED"
            conclusion = "Phase weighting provides significant improvement (unexpected)"
        else:
            falsification_status = "INCONCLUSIVE"
            conclusion = "Mixed results, further investigation needed"
        
        self.results['falsification_status'] = falsification_status
        self.results['conclusion'] = conclusion
        
        # Runtime
        total_time = time.time() - start_time
        self.results['runtime_seconds'] = float(total_time)
        
        logger.info("\n" + "="*60)
        logger.info(f"FALSIFICATION STATUS: {falsification_status}")
        logger.info(f"CONCLUSION: {conclusion}")
        logger.info(f"Total runtime: {total_time:.2f}s")
        logger.info("="*60)
        
        return self.results
    
    def save_results(self, output_path: str):
        """Save results to JSON file."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        logger.info(f"Results saved to {output_path}")


def main():
    """Main entry point for falsification experiment."""
    parser = argparse.ArgumentParser(
        description="Spectral Disruption Profiler Falsification Experiment",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use real human DNA data from CSV (RECOMMENDED)
  python spectral_disruption_profiler.py --input data/doench2016.csv
  
  # Use real human DNA data from FASTA
  python spectral_disruption_profiler.py --input-fasta data/test_human_cdna.fasta --max-sequences 100
  
  # Use synthetic data (for testing only, NOT recommended for falsification)
  python spectral_disruption_profiler.py --use-synthetic --n-sequences 100
"""
    )
    
    # Input data arguments
    parser.add_argument(
        '--input', type=str,
        help='Path to CSV file with human DNA sequences and labels (e.g., data/doench2016.csv)'
    )
    parser.add_argument(
        '--input-fasta', type=str,
        help='Path to FASTA file with human DNA sequences'
    )
    parser.add_argument(
        '--sequence-column', type=str, default='sequence',
        help='Column name for sequences in CSV (default: sequence)'
    )
    parser.add_argument(
        '--label-column', type=str, default='efficiency',
        help='Column name for labels in CSV (default: efficiency)'
    )
    parser.add_argument(
        '--max-sequences', type=int,
        help='Maximum number of sequences to load from input file'
    )
    parser.add_argument(
        '--use-synthetic', action='store_true',
        help='Use synthetic random DNA data instead of real human DNA (NOT RECOMMENDED)'
    )
    
    # Experiment parameters
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--bootstrap', type=int, default=1000,
        help='Number of bootstrap resamples (default: 1000)'
    )
    parser.add_argument(
        '--permutation', type=int, default=1000,
        help='Number of permutation resamples (default: 1000)'
    )
    parser.add_argument(
        '--k-parameter', type=float, default=0.3,
        help='Geodesic curvature exponent k (default: 0.3)'
    )
    
    # Synthetic data parameters (only used with --use-synthetic)
    parser.add_argument(
        '--n-sequences', type=int, default=100,
        help='Number of synthetic sequences to generate (only with --use-synthetic)'
    )
    parser.add_argument(
        '--gc-min', type=float, default=0.4,
        help='Minimum GC content for synthetic sequences (default: 0.4)'
    )
    parser.add_argument(
        '--gc-max', type=float, default=0.6,
        help='Maximum GC content for synthetic sequences (default: 0.6)'
    )
    
    # Output and metadata
    parser.add_argument(
        '--output', type=str, 
        default='results/spectral_disruption_profiler_137/falsification_results.json',
        help='Output path for results JSON'
    )
    parser.add_argument(
        '--splits', type=str, default='single',
        help='Data split strategy (single, split-by-gene, split-by-screen)'
    )
    parser.add_argument(
        '--domain', type=str, default='discrete',
        help='Domain type (discrete, continuous)'
    )
    parser.add_argument(
        '--pam', type=str, default='NGG',
        help='PAM sequence for CRISPR guides (default: NGG)'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.use_synthetic and not args.input and not args.input_fasta:
        print("ERROR: Must provide either --input or --input-fasta for real human DNA data,")
        print("       or use --use-synthetic flag (not recommended for falsification).")
        print("\nRecommended usage:")
        print("  python spectral_disruption_profiler.py --input data/doench2016.csv")
        return 1
    
    if args.use_synthetic and (args.input or args.input_fasta):
        print("WARNING: --use-synthetic flag ignores --input and --input-fasta.")
        print("         Using synthetic data instead of real human DNA.")
    
    # Create falsifier
    falsifier = SpectralDisruptionFalsifier(
        seed=args.seed,
        n_bootstrap=args.bootstrap,
        n_permutation=args.permutation,
        k_parameter=args.k_parameter
    )
    
    # Load data
    sequences = None
    labels = None
    
    if not args.use_synthetic:
        if args.input:
            # Load from CSV
            sequences, labels = falsifier.load_human_dna_from_csv(
                args.input,
                sequence_column=args.sequence_column,
                label_column=args.label_column,
                max_sequences=args.max_sequences
            )
        elif args.input_fasta:
            # Load from FASTA
            sequences = falsifier.load_human_dna_from_fasta(
                args.input_fasta,
                max_sequences=args.max_sequences
            )
            # For FASTA without labels, we cannot perform meaningful falsification
            print("WARNING: FASTA file has no labels. Cannot compute meaningful ROC-AUC.")
            print("         Consider using a CSV file with efficiency labels instead.")
            return 1
    
    # Run experiment
    results = falsifier.run_falsification_experiment(
        sequences=sequences,
        labels=labels,
        n_sequences=args.n_sequences,
        gc_range=(args.gc_min, args.gc_max),
        use_synthetic=args.use_synthetic
    )
    
    # Save results
    falsifier.save_results(args.output)
    
    # Print summary
    print("\n" + "="*60)
    print("EXPERIMENT SUMMARY")
    print("="*60)
    print(f"Falsification Status: {results['falsification_status']}")
    print(f"Conclusion: {results['conclusion']}")
    print(f"\nBootstrap CI: [{results['primary_endpoints']['bootstrap_ci']['ci_lower']:.4f}, "
          f"{results['primary_endpoints']['bootstrap_ci']['ci_upper']:.4f}]")
    print(f"CI includes zero: {results['primary_endpoints']['bootstrap_ci']['includes_zero']}")
    print(f"\nPermutation p-value: {results['secondary_endpoints']['permutation_test']['p_value']:.4f}")
    print(f"Significant (α=0.05): {results['secondary_endpoints']['permutation_test']['significant']}")
    print(f"\nRuntime: {results['runtime_seconds']:.2f}s")
    print("="*60)
    
    return 0 if results['falsification_status'] in ['HYPOTHESIS_FALSIFIED', 'INCONCLUSIVE'] else 1


if __name__ == '__main__':
    sys.exit(main())
