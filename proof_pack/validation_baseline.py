#!/usr/bin/env python3
"""
Baseline Feature Extraction for Z Framework Validation

Provides conventional feature extraction methods for comparison against
Z Framework geodesic mapping and Z5D predictor approaches.

Usage:
    from validation_baseline import BaselineFeatureExtractor
    extractor = BaselineFeatureExtractor()
    features = extractor.extract_features(sequences)
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Union
from collections import Counter
import logging

# Set reproducible seed
np.random.seed(42)

class BaselineFeatureExtractor:
    """Extract conventional sequence and numerical features for baseline comparison."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def extract_sequence_features(self, sequences: List[str]) -> np.ndarray:
        """Extract basic sequence composition features."""
        features = []
        
        for seq in sequences:
            seq_features = {}
            
            # Length
            seq_features['length'] = len(seq)
            
            # Base composition
            base_counts = Counter(seq.upper())
            total_bases = len(seq)
            
            for base in ['A', 'T', 'G', 'C']:
                seq_features[f'{base}_count'] = base_counts.get(base, 0)
                seq_features[f'{base}_freq'] = base_counts.get(base, 0) / total_bases if total_bases > 0 else 0
            
            # GC content
            gc_count = base_counts.get('G', 0) + base_counts.get('C', 0)
            seq_features['gc_content'] = gc_count / total_bases if total_bases > 0 else 0
            
            # Dinucleotide frequencies
            dinucs = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 
                     'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
            
            for dinuc in dinucs:
                count = seq.upper().count(dinuc)
                seq_features[f'dinuc_{dinuc}'] = count / (total_bases - 1) if total_bases > 1 else 0
            
            # Simple complexity measures
            seq_features['entropy'] = self._calculate_entropy(seq)
            seq_features['complexity'] = len(set(seq)) / len(seq) if len(seq) > 0 else 0
            
            features.append(list(seq_features.values()))
        
        return np.array(features)
    
    def extract_numeric_features(self, numeric_sequences: List[List[int]]) -> np.ndarray:
        """Extract statistical features from numeric sequences."""
        features = []
        
        for seq in numeric_sequences:
            seq_array = np.array(seq)
            
            if len(seq_array) == 0:
                # Handle empty sequences
                seq_features = [0] * 15
            else:
                seq_features = [
                    np.mean(seq_array),
                    np.std(seq_array),
                    np.var(seq_array),
                    np.min(seq_array),
                    np.max(seq_array),
                    np.median(seq_array),
                    np.sum(seq_array),
                    len(seq_array),
                    np.percentile(seq_array, 25),
                    np.percentile(seq_array, 75),
                    np.ptp(seq_array),  # Range
                    np.mean(np.diff(seq_array)) if len(seq_array) > 1 else 0,  # Mean difference
                    np.std(np.diff(seq_array)) if len(seq_array) > 1 else 0,   # Diff std
                    len(np.unique(seq_array)),  # Unique values
                    np.sum(seq_array > np.mean(seq_array)) / len(seq_array)  # Above mean ratio
                ]
            
            features.append(seq_features)
        
        return np.array(features)
    
    def extract_histogram_features(self, numeric_sequences: List[List[int]], bins: int = 20) -> np.ndarray:
        """Extract histogram-based density features for comparison with Z Framework."""
        features = []
        
        for seq in numeric_sequences:
            if len(seq) == 0:
                seq_features = [0] * bins
            else:
                hist, _ = np.histogram(seq, bins=bins, density=True)
                seq_features = hist.tolist()
            
            features.append(seq_features)
        
        return np.array(features)
    
    def extract_fourier_features(self, numeric_sequences: List[List[int]], n_components: int = 5) -> np.ndarray:
        """Extract Fourier transform features for frequency analysis."""
        features = []
        
        for seq in numeric_sequences:
            if len(seq) < 2:
                seq_features = [0] * (n_components * 2)  # Real + imaginary
            else:
                # Pad or truncate to standard length for FFT
                standard_length = 64
                if len(seq) < standard_length:
                    seq_padded = np.pad(seq, (0, standard_length - len(seq)), 'wrap')
                else:
                    seq_padded = seq[:standard_length]
                
                # Compute FFT
                fft_vals = np.fft.fft(seq_padded)[:n_components]
                
                # Extract real and imaginary parts
                seq_features = []
                for val in fft_vals:
                    seq_features.append(val.real)
                    seq_features.append(val.imag)
            
            features.append(seq_features)
        
        return np.array(features)
    
    def _calculate_entropy(self, sequence: str) -> float:
        """Calculate Shannon entropy of sequence."""
        if len(sequence) == 0:
            return 0
        
        counts = Counter(sequence.upper())
        probs = [count / len(sequence) for count in counts.values()]
        entropy = -sum(p * np.log2(p) for p in probs if p > 0)
        return entropy
    
    def extract_all_features(self, sequences: List[str], numeric_sequences: List[List[int]]) -> Dict[str, np.ndarray]:
        """Extract all baseline feature types."""
        return {
            'sequence': self.extract_sequence_features(sequences),
            'numeric': self.extract_numeric_features(numeric_sequences),
            'histogram': self.extract_histogram_features(numeric_sequences),
            'fourier': self.extract_fourier_features(numeric_sequences)
        }

class PoissonnNullModel:
    """Generate Poisson-based null model for comparison."""
    
    def __init__(self, random_state: int = 42):
        self.random_state = random_state
        np.random.seed(random_state)
    
    def generate_null_sequences(self, n_samples: int, mean_length: int = 100) -> List[str]:
        """Generate random sequences with Poisson-distributed lengths."""
        sequences = []
        bases = ['A', 'T', 'G', 'C']
        
        for _ in range(n_samples):
            length = max(1, np.random.poisson(mean_length))
            sequence = ''.join(np.random.choice(bases, length))
            sequences.append(sequence)
        
        return sequences
    
    def generate_null_spikes(self, n_samples: int, mean_rate: float = 10.0) -> List[Dict]:
        """Generate Poisson spike trains as null model."""
        spikes = []
        
        for i in range(n_samples):
            rate = max(0.1, np.random.exponential(mean_rate))
            burst_ratio = np.random.beta(1, 3)  # Skewed toward low bursting
            pain_score = np.random.uniform(0, 10)
            
            spikes.append({
                'sample_id': f'NULL_{i:04d}',
                'condition': 'null',
                'spike_rate_hz': rate,
                'burst_ratio': burst_ratio,
                'pain_score': pain_score,
                'therapeutic_response': np.random.randint(0, 2)
            })
        
        return spikes