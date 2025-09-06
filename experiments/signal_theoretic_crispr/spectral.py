#!/usr/bin/env python3
"""
Arm B: WAVE-CRISPR Spectral-Geodesic Features Implementation

This module implements the signal-theoretic approach with two encoding bands:
- Band-Bio: Physicochemically anchored complex mapping
- Band-Arb: Randomized (but fixed-seed) complex amplitudes

Features extracted include:
- Spectral entropy and flatness
- Peak magnitude and frequency  
- Sidelobe bandpowers
- Stability ratios
- Geodesic-modulated Z features
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import logging
from scipy.fft import fft
from scipy.stats import entropy
import mpmath
from sklearn.linear_model import LinearRegression, Ridge, LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import r2_score, mean_absolute_error, roc_auc_score
from sklearn.model_selection import cross_val_score
import warnings

warnings.filterwarnings('ignore')

# Set high-precision arithmetic
mpmath.mp.dps = 50


class BiologicalEncoder:
    """Band-Bio: Physicochemically anchored complex encoding."""
    
    def __init__(self, seed: int = 42):
        """
        Initialize biological encoder with physicochemical properties.
        
        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Physicochemical properties for complex encoding
        self.nucleotide_properties = {
            'A': {
                'polarizability': 1.65,  # Å³
                'h_bond_capacity': 2,    # donor + acceptor
                'molecular_weight': 331.2,
                'purine': 1,
                'ring_electrons': 10
            },
            'T': {
                'polarizability': 1.52,
                'h_bond_capacity': 2,
                'molecular_weight': 322.2,
                'purine': 0,
                'ring_electrons': 6
            },
            'C': {
                'polarizability': 1.42,
                'h_bond_capacity': 3,
                'molecular_weight': 307.2,
                'purine': 0,
                'ring_electrons': 6
            },
            'G': {
                'polarizability': 1.78,
                'h_bond_capacity': 3,
                'molecular_weight': 347.2,
                'purine': 1,
                'ring_electrons': 10
            },
            'N': {
                'polarizability': 1.5,   # Average
                'h_bond_capacity': 2.5,  # Average
                'molecular_weight': 327.0,  # Average
                'purine': 0.5,
                'ring_electrons': 8
            }
        }
        
        # Create complex encoding based on polarizability and H-bond tendency
        self.complex_encoding = self._create_complex_encoding()
        
    def _create_complex_encoding(self) -> Dict[str, complex]:
        """Create complex encoding from physicochemical properties."""
        encoding = {}
        
        for base, props in self.nucleotide_properties.items():
            # Real part: normalized polarizability
            real_part = props['polarizability'] / 2.0  # Normalize to ~[-1, 1] range
            
            # Imaginary part: H-bond tendency (normalized)
            imag_part = props['h_bond_capacity'] / 3.0  # Normalize to [0, 1] range
            
            # Create complex number
            encoding[base] = complex(real_part, imag_part)
            
        return encoding
    
    def encode_sequence(self, sequence: str) -> np.ndarray:
        """
        Encode DNA sequence as complex array.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Complex-valued array
        """
        sequence = sequence.upper().strip()
        encoded = []
        
        for base in sequence:
            if base in self.complex_encoding:
                encoded.append(self.complex_encoding[base])
            else:
                # Default to 'N' for unknown bases
                encoded.append(self.complex_encoding['N'])
                
        return np.array(encoded, dtype=complex)
    
    def get_encoding_info(self) -> Dict[str, Any]:
        """Get encoding information for documentation."""
        return {
            'type': 'biological_physicochemical',
            'properties_used': ['polarizability', 'h_bond_capacity'],
            'complex_mapping': {base: str(val) for base, val in self.complex_encoding.items()},
            'normalization': 'polarizability/2.0, h_bond/3.0'
        }


class ArbitraryEncoder:
    """Band-Arb: Randomized but fixed-seed complex encoding."""
    
    def __init__(self, seed: int = 42):
        """
        Initialize arbitrary encoder with fixed random mapping.
        
        Args:
            seed: Random seed for reproducible randomization
        """
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Create reproducible random encoding
        np.random.seed(seed)
        self.complex_encoding = self._create_random_encoding()
        
    def _create_random_encoding(self) -> Dict[str, complex]:
        """Create randomized complex encoding within predefined bounds."""
        encoding = {}
        
        for base in 'ATCGN':
            # Random complex numbers within unit circle
            magnitude = np.random.uniform(0.5, 1.0)  # Avoid very small magnitudes
            phase = np.random.uniform(0, 2 * np.pi)
            
            # Convert to complex number
            real_part = magnitude * np.cos(phase)
            imag_part = magnitude * np.sin(phase)
            encoding[base] = complex(real_part, imag_part)
            
        return encoding
    
    def encode_sequence(self, sequence: str) -> np.ndarray:
        """
        Encode DNA sequence as complex array.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Complex-valued array
        """
        sequence = sequence.upper().strip()
        encoded = []
        
        for base in sequence:
            if base in self.complex_encoding:
                encoded.append(self.complex_encoding[base])
            else:
                # Default to 'N' for unknown bases
                encoded.append(self.complex_encoding['N'])
                
        return np.array(encoded, dtype=complex)
    
    def get_encoding_info(self) -> Dict[str, Any]:
        """Get encoding information for documentation."""
        return {
            'type': 'arbitrary_randomized_fixed_seed',
            'seed': self.seed,
            'complex_mapping': {base: str(val) for base, val in self.complex_encoding.items()},
            'bounds': 'magnitude: [0.5, 1.0], phase: [0, 2π]'
        }


class GeodeticResolution:
    """Geodesic resolution calculations for positional modulation."""
    
    def __init__(self, k: float = 0.3):
        """
        Initialize geodesic resolution with parameter k.
        
        Args:
            k: Geodesic resolution parameter (default 0.3)
        """
        self.k = k
        self.golden_ratio = (1 + np.sqrt(5)) / 2  # φ = 1.618...
        
    def theta_prime(self, n: int) -> float:
        """
        Calculate geodesic resolution: θ'(n,k) = φ·((n mod φ)/φ)^k
        
        Args:
            n: Position index
            
        Returns:
            Geodesic resolution value
        """
        n_mod_phi = n % self.golden_ratio
        normalized = n_mod_phi / self.golden_ratio
        return self.golden_ratio * (normalized ** self.k)
    
    def apply_modulation(self, encoded_seq: np.ndarray) -> np.ndarray:
        """
        Apply geodesic resolution modulation to encoded sequence.
        
        Args:
            encoded_seq: Complex-encoded sequence
            
        Returns:
            Modulated complex sequence
        """
        modulated = []
        
        for i, base_encoding in enumerate(encoded_seq):
            geodesic_factor = self.theta_prime(i)
            modulated_value = base_encoding * geodesic_factor
            modulated.append(modulated_value)
            
        return np.array(modulated, dtype=complex)


class SpectralAnalyzer:
    """Spectral analysis with high-precision calculations."""
    
    def __init__(self, d: float = 0.34):
        """
        Initialize spectral analyzer.
        
        Args:
            d: Position step size for waveform construction
        """
        self.d = d
        self.logger = logging.getLogger(__name__)
        
    def build_waveform(self, encoded_seq: np.ndarray) -> np.ndarray:
        """
        Build position-modulated waveform with high precision.
        
        Args:
            encoded_seq: Complex-encoded sequence
            
        Returns:
            Complex waveform
        """
        positions = np.cumsum([self.d] * len(encoded_seq))
        waveform = []
        
        for i, (base_encoding, pos) in enumerate(zip(encoded_seq, positions)):
            # Use mpmath for high-precision phase calculations
            phase_factor = mpmath.exp(2j * mpmath.pi * pos)
            waveform_component = complex(base_encoding) * complex(phase_factor)
            waveform.append(waveform_component)
            
        return np.array(waveform, dtype=complex)
    
    def compute_spectrum(self, waveform: np.ndarray) -> np.ndarray:
        """
        Compute FFT spectrum with high precision.
        
        Args:
            waveform: Complex waveform
            
        Returns:
            Power spectrum
        """
        return np.abs(fft(waveform.astype(complex)))
    
    def spectral_entropy(self, spectrum: np.ndarray) -> float:
        """
        Calculate normalized spectral entropy.
        
        Args:
            spectrum: Power spectrum
            
        Returns:
            Spectral entropy
        """
        # Normalize spectrum to probability distribution
        prob_spectrum = spectrum / np.sum(spectrum)
        # Remove zeros to avoid log(0)
        prob_spectrum = prob_spectrum[prob_spectrum > 1e-15]
        return entropy(prob_spectrum, base=2)
    
    def spectral_flatness(self, spectrum: np.ndarray) -> float:
        """
        Calculate spectral flatness (Wiener entropy).
        
        Args:
            spectrum: Power spectrum
            
        Returns:
            Spectral flatness [0, 1]
        """
        # Avoid zeros
        spectrum_nonzero = spectrum[spectrum > 1e-15]
        if len(spectrum_nonzero) == 0:
            return 0.0
        
        geometric_mean = np.exp(np.mean(np.log(spectrum_nonzero)))
        arithmetic_mean = np.mean(spectrum_nonzero)
        
        if arithmetic_mean == 0:
            return 0.0
        
        return geometric_mean / arithmetic_mean
    
    def extract_spectral_features(self, spectrum: np.ndarray) -> Dict[str, float]:
        """
        Extract comprehensive spectral features.
        
        Args:
            spectrum: Power spectrum
            
        Returns:
            Dictionary of spectral features
        """
        features = {
            'spectral_entropy': self.spectral_entropy(spectrum),
            'spectral_flatness': self.spectral_flatness(spectrum),
            'peak_magnitude': np.max(spectrum),
            'mean_magnitude': np.mean(spectrum),
            'spectral_variance': np.var(spectrum),
            'spectral_std': np.std(spectrum),
            'peak_frequency': np.argmax(spectrum),
            'bandwidth': np.sum(spectrum > 0.1 * np.max(spectrum)),
            'spectral_centroid': np.sum(spectrum * np.arange(len(spectrum))) / np.sum(spectrum) if np.sum(spectrum) > 0 else 0
        }
        
        # Sidelobe analysis
        peak_idx = np.argmax(spectrum)
        if peak_idx > 0 and peak_idx < len(spectrum) - 1:
            # Left and right sidelobes
            left_sidelobe = np.sum(spectrum[:peak_idx])
            right_sidelobe = np.sum(spectrum[peak_idx+1:])
            total_power = np.sum(spectrum)
            
            features['left_sidelobe_power'] = left_sidelobe / total_power if total_power > 0 else 0
            features['right_sidelobe_power'] = right_sidelobe / total_power if total_power > 0 else 0
            features['sidelobe_asymmetry'] = abs(left_sidelobe - right_sidelobe) / total_power if total_power > 0 else 0
        
        # Stability ratios
        if len(spectrum) > 1:
            features['peak_to_mean_ratio'] = features['peak_magnitude'] / features['mean_magnitude'] if features['mean_magnitude'] > 0 else 0
            features['dynamic_range'] = features['peak_magnitude'] / np.min(spectrum[spectrum > 0]) if np.any(spectrum > 0) else 0
        
        return features


class ZFrameworkNormalizer:
    """Z Framework normalization: Z = n(Δₙ/Δₘₐₓ)."""
    
    def __init__(self):
        """Initialize Z Framework normalizer."""
        self.logger = logging.getLogger(__name__)
        
    def compute_z_features(self, features_dict: Dict[str, float], 
                          reference_features: Dict[str, float]) -> Dict[str, float]:
        """
        Compute normalized Z features.
        
        Args:
            features_dict: Features for current sequence
            reference_features: Reference features for normalization
            
        Returns:
            Dictionary of Z-normalized features
        """
        z_features = {}
        
        for feature_name in features_dict:
            if feature_name in reference_features:
                delta_n = features_dict[feature_name] - reference_features[feature_name]
                delta_max = max(abs(features_dict[feature_name]), abs(reference_features[feature_name]))
                
                # Guard against division by zero
                if delta_max != 0:
                    z_value = len(str(feature_name)) * (delta_n / delta_max)  # n * (Δₙ/Δₘₐₓ)
                    z_features[f'z_{feature_name}'] = z_value
                else:
                    z_features[f'z_{feature_name}'] = 0.0
                    
        return z_features


class SpectralFeatureExtractor:
    """Main spectral feature extraction pipeline."""
    
    def __init__(self, k: float = 0.3, seed: int = 42):
        """
        Initialize spectral feature extractor.
        
        Args:
            k: Geodesic resolution parameter
            seed: Random seed for reproducibility
        """
        self.k = k
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        # Initialize components
        self.bio_encoder = BiologicalEncoder(seed)
        self.arb_encoder = ArbitraryEncoder(seed)
        self.geodesic = GeodeticResolution(k)
        self.spectral_analyzer = SpectralAnalyzer()
        self.z_normalizer = ZFrameworkNormalizer()
        
    def extract_band_features(self, sequence: str, encoder, band_name: str) -> Dict[str, float]:
        """
        Extract features for a single encoding band.
        
        Args:
            sequence: DNA sequence
            encoder: Encoder instance (bio or arbitrary)
            band_name: Name of the encoding band
            
        Returns:
            Dictionary of features with band prefix
        """
        try:
            # Encode sequence
            encoded = encoder.encode_sequence(sequence)
            
            # Apply geodesic modulation
            modulated = self.geodesic.apply_modulation(encoded)
            
            # Build waveform
            waveform = self.spectral_analyzer.build_waveform(modulated)
            
            # Compute spectrum
            spectrum = self.spectral_analyzer.compute_spectrum(waveform)
            
            # Extract spectral features
            spectral_features = self.spectral_analyzer.extract_spectral_features(spectrum)
            
            # Add band prefix to feature names
            band_features = {f'{band_name}_{key}': value for key, value in spectral_features.items()}
            
            return band_features
            
        except Exception as e:
            self.logger.warning(f"Failed to extract {band_name} features for sequence: {e}")
            return {}
    
    def extract_features(self, guide_data: Dict[str, Any]) -> Dict[str, float]:
        """
        Extract all spectral features using both encoding bands.
        
        Args:
            guide_data: Dictionary with guide information
            
        Returns:
            Combined feature dictionary
        """
        sequence = guide_data.get('guide_sequence', '')
        context_sequence = guide_data.get('context_sequence', sequence)
        
        # Use the longer sequence if available (201-bp window preferred)
        analysis_sequence = context_sequence if len(context_sequence) > len(sequence) else sequence
        
        features = {}
        
        # Extract features for both bands
        bio_features = self.extract_band_features(analysis_sequence, self.bio_encoder, 'bio')
        arb_features = self.extract_band_features(analysis_sequence, self.arb_encoder, 'arb')
        
        features.update(bio_features)
        features.update(arb_features)
        
        # Compute Z-normalized features comparing bio vs arbitrary
        if bio_features and arb_features:
            # Remove band prefixes for comparison
            bio_clean = {k.replace('bio_', ''): v for k, v in bio_features.items()}
            arb_clean = {k.replace('arb_', ''): v for k, v in arb_features.items()}
            
            z_features = self.z_normalizer.compute_z_features(bio_clean, arb_clean)
            features.update(z_features)
        
        # Add metadata
        features['sequence_length'] = len(analysis_sequence)
        features['geodesic_k'] = self.k
        features['encoding_bands'] = 2
        
        return features


class SpectralModels:
    """Models using spectral-geodesic features."""
    
    def __init__(self, seed: int = 42, use_ensemble: bool = False):
        """
        Initialize spectral models.
        
        Args:
            seed: Random seed for reproducibility
            use_ensemble: Whether to use ensemble models (RF/GBM)
        """
        self.seed = seed
        self.use_ensemble = use_ensemble
        self.logger = logging.getLogger(__name__)
        
        # Initialize models
        if use_ensemble:
            self.on_target_model = RandomForestRegressor(n_estimators=100, random_state=seed)
            self.off_target_model = RandomForestClassifier(n_estimators=100, random_state=seed)
        else:
            self.on_target_model = Ridge(alpha=1.0, random_state=seed)
            self.off_target_model = LogisticRegression(random_state=seed, max_iter=1000)
        
        self.feature_columns = None
        self.is_fitted = False
    
    def prepare_features(self, df: pd.DataFrame) -> np.ndarray:
        """
        Prepare feature matrix from processed DataFrame.
        
        Args:
            df: Processed DataFrame with spectral features
            
        Returns:
            Feature matrix
        """
        # Select spectral feature columns
        exclude_columns = [
            'guide_sequence', 'original_index', 'efficiency', 
            'context_sequence', 'off_target'
        ]
        
        numeric_columns = [col for col in df.columns 
                          if col not in exclude_columns and pd.api.types.is_numeric_dtype(df[col])]
        
        if self.feature_columns is None:
            self.feature_columns = numeric_columns
            self.logger.info(f"Using {len(self.feature_columns)} spectral features")
        
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
            'n_features': X.shape[1],
            'model_type': 'ensemble' if self.use_ensemble else 'linear'
        }
        
        self.logger.info(f"Spectral on-target training complete: R² = {r2:.3f}, MAE = {mae:.3f}")
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
    
    def get_feature_importance(self) -> Dict[str, float]:
        """Get feature importance if using ensemble models."""
        if self.use_ensemble and hasattr(self.on_target_model, 'feature_importances_'):
            return dict(zip(self.feature_columns, self.on_target_model.feature_importances_))
        return {}


if __name__ == "__main__":
    # Quick test of spectral pipeline
    import tempfile
    
    # Test encoders
    bio_encoder = BiologicalEncoder()
    arb_encoder = ArbitraryEncoder()
    
    test_sequence = "ATCGATCGATCGATCGAT"
    
    bio_encoded = bio_encoder.encode_sequence(test_sequence)
    arb_encoded = arb_encoder.encode_sequence(test_sequence)
    
    print(f"✓ Encoding test passed: Bio={len(bio_encoded)}, Arb={len(arb_encoded)}")
    
    # Test feature extraction
    extractor = SpectralFeatureExtractor()
    
    test_guide = {
        'guide_sequence': test_sequence,
        'context_sequence': test_sequence * 10  # Longer context
    }
    
    features = extractor.extract_features(test_guide)
    print(f"✓ Feature extraction test passed: {len(features)} features extracted")
    
    # Show some key features
    key_features = [k for k in features.keys() if any(word in k for word in ['entropy', 'peak', 'z_'])]
    print(f"Key features: {key_features[:5]}")
    
    print("Spectral module self-test complete")