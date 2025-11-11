#!/usr/bin/env python3
"""
Arm A: Batch-Scan Baseline Implementation

This module implements the "Spring Batch" style baseline processing pipeline:
Reader → Processor → Writer

Features extracted:
- GC content and composition
- Homopolymer runs
- Seed region mismatches
- Basic thermodynamic surrogates
- PAM context analysis

No spectral features are included in this baseline arm.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import re
from collections import Counter
import logging
from sklearn.linear_model import LinearRegression, Ridge, LogisticRegression
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.model_selection import cross_val_score
import warnings

warnings.filterwarnings('ignore')


class PAMScanner:
    """Detect and analyze PAM sequences in DNA."""
    
    def __init__(self, pam_pattern: str = "NGG"):
        """
        Initialize PAM scanner.
        
        Args:
            pam_pattern: PAM sequence pattern (e.g., "NGG" for Cas9)
        """
        self.pam_pattern = pam_pattern
        self.pam_regex = self._create_pam_regex(pam_pattern)
        
    def _create_pam_regex(self, pattern: str) -> re.Pattern:
        """Convert PAM pattern to regex."""
        # Convert IUPAC codes to regex
        iupac_map = {
            'N': '[ATCG]',
            'R': '[AG]',
            'Y': '[CT]',
            'S': '[GC]',
            'W': '[AT]',
            'K': '[GT]',
            'M': '[AC]',
            'B': '[CGT]',
            'D': '[AGT]',
            'H': '[ACT]',
            'V': '[ACG]'
        }
        
        regex_pattern = pattern
        for code, replacement in iupac_map.items():
            regex_pattern = regex_pattern.replace(code, replacement)
            
        return re.compile(regex_pattern)
    
    def find_pam_sites(self, sequence: str) -> List[Tuple[int, str]]:
        """
        Find all PAM sites in sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            List of (position, pam_sequence) tuples
        """
        sites = []
        for match in self.pam_regex.finditer(sequence.upper()):
            sites.append((match.start(), match.group()))
        return sites
    
    def extract_guide_candidates(self, sequence: str, guide_length: int = 20) -> List[Dict[str, Any]]:
        """
        Extract guide RNA candidates upstream of PAM sites.
        
        Args:
            sequence: DNA sequence
            guide_length: Length of guide sequence
            
        Returns:
            List of guide candidate dictionaries
        """
        candidates = []
        pam_sites = self.find_pam_sites(sequence)
        
        for pam_pos, pam_seq in pam_sites:
            # Extract guide sequence upstream of PAM
            guide_start = pam_pos - guide_length
            if guide_start >= 0:
                guide_seq = sequence[guide_start:pam_pos].upper()
                if len(guide_seq) == guide_length:
                    candidates.append({
                        'guide_sequence': guide_seq,
                        'pam_sequence': pam_seq,
                        'guide_start': guide_start,
                        'pam_start': pam_pos,
                        'cut_site': pam_pos - 3  # Approximate cut site
                    })
        
        return candidates


class BaselineFeatureExtractor:
    """Extract simple baseline features from CRISPR guide sequences."""
    
    def __init__(self, pam_scanner: PAMScanner):
        """
        Initialize feature extractor.
        
        Args:
            pam_scanner: PAM scanner instance
        """
        self.pam_scanner = pam_scanner
        self.logger = logging.getLogger(__name__)
        
    def gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence) if len(sequence) > 0 else 0.0
    
    def gc_composition(self, sequence: str) -> Dict[str, float]:
        """Calculate detailed nucleotide composition."""
        sequence = sequence.upper()
        total = len(sequence)
        
        if total == 0:
            return {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        
        composition = {}
        for base in 'ATCG':
            composition[f'{base}_content'] = sequence.count(base) / total
            
        return composition
    
    def homopolymer_analysis(self, sequence: str) -> Dict[str, float]:
        """Analyze homopolymer runs in sequence."""
        sequence = sequence.upper()
        
        # Find runs of identical nucleotides
        runs = {'A': [], 'T': [], 'C': [], 'G': []}
        current_base = None
        current_length = 0
        
        for base in sequence:
            if base == current_base:
                current_length += 1
            else:
                if current_base and current_length > 1:
                    runs[current_base].append(current_length)
                current_base = base
                current_length = 1
        
        # Handle final run
        if current_base and current_length > 1:
            runs[current_base].append(current_length)
        
        # Calculate statistics
        features = {
            'max_homopolymer': max([max(runs[base]) if runs[base] else 0 for base in 'ATCG']),
            'total_homopolymer_bases': sum([sum(runs[base]) for base in 'ATCG']),
            'num_homopolymer_runs': sum([len(runs[base]) for base in 'ATCG'])
        }
        
        return features
    
    def seed_region_analysis(self, guide_sequence: str, seed_positions: Tuple[int, int] = (1, 8)) -> Dict[str, float]:
        """
        Analyze seed region properties.
        
        Args:
            guide_sequence: Guide RNA sequence
            seed_positions: (start, end) positions of seed region (1-indexed)
            
        Returns:
            Seed region features
        """
        if len(guide_sequence) < max(seed_positions):
            return {'seed_gc': 0.0, 'seed_mismatches': 0.0}
        
        # Extract seed region (convert to 0-indexed)
        seed_start = seed_positions[0] - 1
        seed_end = seed_positions[1]
        seed_region = guide_sequence[seed_start:seed_end].upper()
        
        features = {
            'seed_gc': self.gc_content(seed_region),
            'seed_length': len(seed_region),
            'seed_a_content': seed_region.count('A') / len(seed_region) if seed_region else 0,
            'seed_t_content': seed_region.count('T') / len(seed_region) if seed_region else 0
        }
        
        return features
    
    def thermodynamic_surrogates(self, sequence: str) -> Dict[str, float]:
        """
        Calculate simple thermodynamic surrogates.
        
        These are rough approximations, not true thermodynamic calculations.
        """
        sequence = sequence.upper()
        
        # Simple nearest-neighbor approximation weights
        dinucleotide_weights = {
            'AA': -1.0, 'AT': -0.88, 'AC': -1.45, 'AG': -1.29,
            'TA': -0.58, 'TT': -1.0, 'TC': -1.29, 'TG': -1.45,
            'CA': -1.45, 'CT': -1.29, 'CC': -1.84, 'CG': -2.27,
            'GA': -1.29, 'GT': -1.45, 'GC': -2.27, 'GG': -1.84
        }
        
        # Calculate rough stability score
        stability_score = 0.0
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            stability_score += dinucleotide_weights.get(dinuc, 0.0)
        
        # Calculate melting temperature surrogate
        # Simple GC counting approximation
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        tm_surrogate = (gc_count * 4) + (at_count * 2)  # Very rough approximation
        
        return {
            'stability_surrogate': stability_score,
            'tm_surrogate': tm_surrogate,
            'gc_tm_ratio': gc_count / len(sequence) if len(sequence) > 0 else 0.0
        }
    
    def extract_features(self, guide_data: Dict[str, Any]) -> Dict[str, float]:
        """
        Extract all baseline features for a guide.
        
        Args:
            guide_data: Dictionary with guide information
            
        Returns:
            Feature dictionary
        """
        guide_seq = guide_data.get('guide_sequence', '')
        context_seq = guide_data.get('context_sequence', guide_seq)  # 201-bp window if available
        
        features = {}
        
        # Basic composition features
        features.update(self.gc_composition(guide_seq))
        features['gc_content'] = self.gc_content(guide_seq)
        
        # Homopolymer features
        features.update(self.homopolymer_analysis(guide_seq))
        
        # Seed region features
        features.update(self.seed_region_analysis(guide_seq))
        
        # Thermodynamic surrogates
        features.update(self.thermodynamic_surrogates(guide_seq))
        
        # PAM context features if available
        if 'pam_sequence' in guide_data:
            pam_seq = guide_data['pam_sequence']
            features['pam_gc'] = self.gc_content(pam_seq)
        
        # Context window features if available
        if len(context_seq) > len(guide_seq):
            features['context_gc'] = self.gc_content(context_seq)
            features['context_length'] = len(context_seq)
        
        # Guide length and basic properties
        features['guide_length'] = len(guide_seq)
        features['purine_content'] = (guide_seq.count('A') + guide_seq.count('G')) / len(guide_seq) if guide_seq else 0
        features['pyrimidine_content'] = (guide_seq.count('C') + guide_seq.count('T')) / len(guide_seq) if guide_seq else 0
        
        return features


class BaselinePipeline:
    """Spring Batch style baseline processing pipeline."""
    
    def __init__(self, pam_pattern: str = "NGG", seed: int = 42):
        """
        Initialize baseline pipeline.
        
        Args:
            pam_pattern: PAM sequence pattern
            seed: Random seed for reproducibility
        """
        self.pam_scanner = PAMScanner(pam_pattern)
        self.feature_extractor = BaselineFeatureExtractor(self.pam_scanner)
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Set random seeds
        np.random.seed(seed)
        
    def read_data(self, data_path: str) -> pd.DataFrame:
        """
        Reader: Load and validate input data.
        
        Args:
            data_path: Path to input data file
            
        Returns:
            Loaded DataFrame
        """
        self.logger.info(f"Reading data from {data_path}")
        
        if data_path.endswith('.csv'):
            df = pd.read_csv(data_path)
        else:
            raise ValueError(f"Unsupported file format: {data_path}")
        
        # Basic validation
        required_columns = ['sequence']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        self.logger.info(f"Loaded {len(df)} records")
        return df
    
    def process_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Processor: Extract baseline features from sequences.
        
        Args:
            df: Input DataFrame with sequences
            
        Returns:
            DataFrame with extracted features
        """
        self.logger.info("Processing sequences to extract baseline features")
        
        processed_records = []
        
        for idx, row in df.iterrows():
            try:
                # Create guide data structure
                guide_data = {
                    'guide_sequence': row['sequence'],
                    'original_index': idx
                }
                
                # Add efficiency if available
                if 'efficiency' in row:
                    guide_data['efficiency'] = row['efficiency']
                
                # Extract features
                features = self.feature_extractor.extract_features(guide_data)
                
                # Combine with original data
                processed_record = {**guide_data, **features}
                processed_records.append(processed_record)
                
            except Exception as e:
                self.logger.warning(f"Failed to process sequence {idx}: {e}")
                continue
        
        result_df = pd.DataFrame(processed_records)
        self.logger.info(f"Successfully processed {len(result_df)} sequences")
        
        return result_df
    
    def write_results(self, df: pd.DataFrame, output_path: str) -> None:
        """
        Writer: Save processed results.
        
        Args:
            df: Processed DataFrame
            output_path: Output file path
        """
        self.logger.info(f"Writing results to {output_path}")
        df.to_csv(output_path, index=False)
    
    def run_pipeline(self, input_path: str, output_path: str) -> pd.DataFrame:
        """
        Run complete baseline processing pipeline.
        
        Args:
            input_path: Input data path
            output_path: Output results path
            
        Returns:
            Processed DataFrame
        """
        # Reader
        df = self.read_data(input_path)
        
        # Processor
        processed_df = self.process_data(df)
        
        # Writer
        self.write_results(processed_df, output_path)
        
        return processed_df


class BaselineModels:
    """Baseline models for on-target and off-target prediction."""
    
    def __init__(self, seed: int = 42):
        """
        Initialize baseline models.
        
        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Initialize models
        self.on_target_model = Ridge(alpha=1.0, random_state=seed)
        self.off_target_model = LogisticRegression(random_state=seed, max_iter=1000)
        
        self.feature_columns = None
        self.is_fitted = False
    
    def prepare_features(self, df: pd.DataFrame) -> np.ndarray:
        """
        Prepare feature matrix from processed DataFrame.
        
        Args:
            df: Processed DataFrame with features
            
        Returns:
            Feature matrix
        """
        # Select numeric feature columns (exclude metadata)
        exclude_columns = [
            'guide_sequence', 'original_index', 'efficiency', 
            'pam_sequence', 'context_sequence'
        ]
        
        numeric_columns = [col for col in df.columns 
                          if col not in exclude_columns and pd.api.types.is_numeric_dtype(df[col])]
        
        if self.feature_columns is None:
            self.feature_columns = numeric_columns
            self.logger.info(f"Using {len(self.feature_columns)} baseline features")
        
        # Ensure consistent feature ordering
        X = df[self.feature_columns].fillna(0).values
        return X
    
    def train_on_target(self, df: pd.DataFrame) -> Dict[str, float]:
        """
        Train on-target efficiency regression model.
        
        Args:
            df: Training data with efficiency labels
            
        Returns:
            Training metrics
        """
        if 'efficiency' not in df.columns:
            raise ValueError("No efficiency labels found for on-target training")
        
        X = self.prepare_features(df)
        y = df['efficiency'].values
        
        # Remove any invalid labels
        valid_mask = ~np.isnan(y)
        X = X[valid_mask]
        y = y[valid_mask]
        
        # Train model
        self.on_target_model.fit(X, y)
        
        # Calculate training metrics
        y_pred = self.on_target_model.predict(X)
        r2 = r2_score(y, y_pred)
        mae = mean_absolute_error(y, y_pred)
        
        # Cross-validation score
        cv_scores = cross_val_score(self.on_target_model, X, y, cv=5, scoring='r2')
        
        metrics = {
            'train_r2': r2,
            'train_mae': mae,
            'cv_r2_mean': np.mean(cv_scores),
            'cv_r2_std': np.std(cv_scores),
            'n_samples': len(y),
            'n_features': X.shape[1]
        }
        
        self.logger.info(f"On-target training complete: R² = {r2:.3f}, MAE = {mae:.3f}")
        return metrics
    
    def predict_on_target(self, df: pd.DataFrame) -> np.ndarray:
        """
        Predict on-target efficiency.
        
        Args:
            df: Data to predict
            
        Returns:
            Predicted efficiency scores
        """
        X = self.prepare_features(df)
        return self.on_target_model.predict(X)
    
    def train_off_target(self, df: pd.DataFrame) -> Dict[str, float]:
        """
        Train off-target classification model.
        
        Args:
            df: Training data with off-target labels
            
        Returns:
            Training metrics
        """
        if 'off_target' not in df.columns:
            raise ValueError("No off-target labels found for classification training")
        
        X = self.prepare_features(df)
        y = df['off_target'].values
        
        # Train model
        self.off_target_model.fit(X, y)
        
        # Calculate training metrics
        train_score = self.off_target_model.score(X, y)
        cv_scores = cross_val_score(self.off_target_model, X, y, cv=5, scoring='accuracy')
        
        metrics = {
            'train_accuracy': train_score,
            'cv_accuracy_mean': np.mean(cv_scores),
            'cv_accuracy_std': np.std(cv_scores),
            'n_samples': len(y),
            'n_features': X.shape[1]
        }
        
        self.logger.info(f"Off-target training complete: Accuracy = {train_score:.3f}")
        return metrics
    
    def predict_off_target(self, df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict off-target probability.
        
        Args:
            df: Data to predict
            
        Returns:
            (predictions, probabilities)
        """
        X = self.prepare_features(df)
        predictions = self.off_target_model.predict(X)
        probabilities = self.off_target_model.predict_proba(X)[:, 1]  # Positive class probability
        return predictions, probabilities


if __name__ == "__main__":
    # Quick test of baseline pipeline
    import tempfile
    import os
    
    # Create test data
    test_data = pd.DataFrame({
        'sequence': [
            'ATCGATCGATCGATCGAT',
            'GGGGGGGGGGGGGGGGGG',
            'AAAAAAAAAAAAAAAAAAA',
            'CCCCCCCCCCCCCCCCCC'
        ],
        'efficiency': [0.72, 0.85, 0.45, 0.91]
    })
    
    # Test pipeline
    pipeline = BaselinePipeline()
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_input:
        test_data.to_csv(tmp_input.name, index=False)
        tmp_input_path = tmp_input.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_output:
        tmp_output_path = tmp_output.name
    
    try:
        # Run pipeline
        result_df = pipeline.run_pipeline(tmp_input_path, tmp_output_path)
        print(f"✓ Baseline pipeline test passed: {len(result_df)} sequences processed")
        print(f"Features extracted: {[col for col in result_df.columns if col.endswith('_content') or col.startswith('seed_')]}")
        
        # Test model training
        models = BaselineModels()
        metrics = models.train_on_target(result_df)
        print(f"✓ Model training test passed: R² = {metrics['train_r2']:.3f}")
        
    finally:
        # Cleanup
        os.unlink(tmp_input_path)
        os.unlink(tmp_output_path)
    
    print("Baseline module self-test complete")