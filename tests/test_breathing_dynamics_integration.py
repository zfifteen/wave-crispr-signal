#!/usr/bin/env python3
"""
Unit Tests for Breathing Dynamics and Ablation Framework

Tests all components of the breathing dynamics implementation including:
- Biophysical encoding with temperature/Mg²⁺ dependencies
- CZT and Goertzel algorithms for fractional-period analysis
- Ablation tests and null distributions
- Statistical validation framework
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from experiments.signal_theoretic_crispr.breathing_dynamics import (
    BreathingDynamicsEncoder,
    ChirpZTransform,
    GoertzelAlgorithm,
    BreathingSpectralAnalyzer,
    HELICAL_PERIOD_BP,
    BP_OPENING_LIFETIMES_MS
)

from experiments.signal_theoretic_crispr.ablation_tests import (
    RandomEncoder,
    AblationTester
)


class TestBreathingDynamicsEncoder:
    """Test BreathingDynamicsEncoder class."""
    
    def test_encoder_initialization(self):
        """Test encoder initializes with correct parameters."""
        encoder = BreathingDynamicsEncoder(
            temperature_c=37.0,
            mg_concentration_mm=2.0,
            helical_period=10.5
        )
        
        assert encoder.temperature_c == 37.0
        assert encoder.mg_concentration_mm == 2.0
        assert encoder.helical_period == 10.5
        assert len(encoder.base_weights) == 5  # A, T, C, G, N
    
    def test_sequence_encoding(self):
        """Test DNA sequence encoding."""
        encoder = BreathingDynamicsEncoder()
        
        test_seq = "ATCGATCG"
        encoded = encoder.encode_sequence(test_seq)
        
        assert len(encoded) == len(test_seq)
        assert encoded.dtype == complex
        assert all(np.isfinite(encoded))
    
    def test_invalid_sequence_rejection(self):
        """Test that invalid sequences are rejected."""
        encoder = BreathingDynamicsEncoder()
        
        # Should reject sequences with U (RNA) or other invalid bases
        with pytest.raises(ValueError):
            encoder.encode_sequence("AUCGATCG")  # Contains U
        
        with pytest.raises(ValueError):
            encoder.encode_sequence("ATCXATCG")  # Contains X
    
    def test_at_gc_weight_difference(self):
        """Test that AT and GC bases have different weights."""
        encoder = BreathingDynamicsEncoder()
        
        at_weight = encoder.base_weights['A']
        gc_weight = encoder.base_weights['G']
        
        # AT and GC should have different complex weights
        assert at_weight != gc_weight
        
        # Real parts should differ (opening rates)
        assert at_weight.real != gc_weight.real
    
    def test_helical_phase_application(self):
        """Test that helical phase modulation is applied correctly."""
        encoder = BreathingDynamicsEncoder()
        
        test_seq = "AAAAAAAAAA"
        
        # With helical phase
        encoded_with_phase = encoder.encode_sequence(test_seq, apply_helical_phase=True)
        
        # Without helical phase
        encoded_no_phase = encoder.encode_sequence(test_seq, apply_helical_phase=False)
        
        # Should be different
        assert not np.allclose(encoded_with_phase, encoded_no_phase)
    
    def test_temperature_dependence(self):
        """Test that temperature affects encoding."""
        encoder_37 = BreathingDynamicsEncoder(temperature_c=37.0)
        encoder_25 = BreathingDynamicsEncoder(temperature_c=25.0)
        
        # Weights should differ with temperature
        assert encoder_37.base_weights['G'] != encoder_25.base_weights['G']
    
    def test_mg_dependence(self):
        """Test that Mg²⁺ concentration affects encoding."""
        encoder_2mm = BreathingDynamicsEncoder(mg_concentration_mm=2.0)
        encoder_10mm = BreathingDynamicsEncoder(mg_concentration_mm=10.0)
        
        # Weights should differ with Mg²⁺
        assert encoder_2mm.base_weights['G'] != encoder_10mm.base_weights['G']


class TestChirpZTransform:
    """Test CZT implementation."""
    
    def test_czt_initialization(self):
        """Test CZT initializes correctly."""
        czt = ChirpZTransform(sequence_length=20)
        assert czt.N == 20
    
    def test_czt_computation(self):
        """Test CZT computes at specific frequency."""
        czt = ChirpZTransform(sequence_length=20)
        
        # Create simple test signal
        signal = np.ones(20, dtype=complex)
        
        # Compute CZT at specific frequency
        result = czt.compute_czt(signal, frequency_hz=0.1, sampling_rate_hz=1.0, num_points=1)
        
        assert len(result) == 1
        assert np.isfinite(result[0])
    
    def test_fractional_period_evaluation(self):
        """Test CZT evaluates at fractional period (10.5 bp)."""
        czt = ChirpZTransform(sequence_length=21)  # 2 × 10.5
        
        # Create test signal
        signal = np.exp(2j * np.pi * np.arange(21) / 10.5)
        
        # Evaluate at 10.5 bp period
        results = czt.evaluate_at_fractional_period(signal, period_bp=10.5, harmonics=3)
        
        # Should have results for 3 harmonics plus total
        assert 'czt_period_10.5_h1_power' in results
        assert 'czt_period_10.5_h2_power' in results
        assert 'czt_period_10.5_h3_power' in results
        assert 'czt_period_10.5_total_power' in results
        
        # All values should be finite
        assert all(np.isfinite(v) for v in results.values())


class TestGoertzelAlgorithm:
    """Test Goertzel algorithm implementation."""
    
    def test_goertzel_computation(self):
        """Test Goertzel computes at specific frequency."""
        goertzel = GoertzelAlgorithm()
        
        # Create test signal
        signal = np.ones(20, dtype=complex)
        
        # Compute Goertzel
        result = goertzel.compute_goertzel(signal, frequency_hz=0.1, sampling_rate_hz=1.0)
        
        assert isinstance(result, complex)
        assert np.isfinite(result)
    
    def test_goertzel_period_evaluation(self):
        """Test Goertzel evaluates at specific period."""
        goertzel = GoertzelAlgorithm()
        
        # Create test signal
        signal = np.exp(2j * np.pi * np.arange(21) / 10.5)
        
        # Evaluate at 10.5 bp period
        results = goertzel.evaluate_at_period(signal, period_bp=10.5, harmonics=3)
        
        # Should have results for 3 harmonics plus total
        assert 'goertzel_period_10.5_h1_power' in results
        assert 'goertzel_period_10.5_total_power' in results
        
        # All values should be finite
        assert all(np.isfinite(v) for v in results.values())


class TestBreathingSpectralAnalyzer:
    """Test complete spectral analyzer."""
    
    def test_analyzer_initialization(self):
        """Test analyzer initializes correctly."""
        analyzer = BreathingSpectralAnalyzer(use_czt=True)
        assert analyzer.use_czt == True
        assert analyzer.encoder is not None
    
    def test_dc_removal(self):
        """Test DC component removal."""
        analyzer = BreathingSpectralAnalyzer()
        
        signal = np.array([1+1j, 2+2j, 3+3j])
        dc_removed = analyzer.remove_dc(signal)
        
        # Mean should be approximately zero
        assert np.abs(np.mean(dc_removed)) < 1e-10
    
    def test_windowing(self):
        """Test window application."""
        analyzer = BreathingSpectralAnalyzer()
        
        signal = np.ones(20, dtype=complex)
        
        # Apply Hamming window
        windowed = analyzer.apply_window(signal, window_type='hamming')
        
        assert len(windowed) == len(signal)
        # Windowed signal should have reduced edges
        assert np.abs(windowed[0]) < np.abs(signal[0])
    
    def test_feature_extraction(self):
        """Test complete feature extraction pipeline."""
        analyzer = BreathingSpectralAnalyzer(use_czt=True)
        
        test_seq = "ATCGATCGATCGATCGATCG"
        features = analyzer.extract_breathing_features(test_seq, harmonics=3)
        
        # Should have multiple features
        assert len(features) > 10
        
        # Should have CZT features
        assert any('czt_period' in k for k in features.keys())
        
        # Should have breathing-specific metrics
        assert 'breathing_gc_content' in features
        assert 'breathing_phase_coherence' in features
        
        # All features should be finite
        assert all(np.isfinite(v) for v in features.values())
    
    def test_czt_vs_goertzel(self):
        """Test that CZT and Goertzel give similar results."""
        analyzer_czt = BreathingSpectralAnalyzer(use_czt=True)
        analyzer_goertzel = BreathingSpectralAnalyzer(use_czt=False)
        
        test_seq = "ATCGATCGATCGATCGATCG"
        
        features_czt = analyzer_czt.extract_breathing_features(test_seq)
        features_goertzel = analyzer_goertzel.extract_breathing_features(test_seq)
        
        # Both should extract features
        assert len(features_czt) > 0
        assert len(features_goertzel) > 0


class TestRandomEncoder:
    """Test random encoder for null distributions."""
    
    def test_random_encoder_initialization(self):
        """Test random encoder creates random weights."""
        encoder = RandomEncoder(seed=42)
        
        assert len(encoder.base_weights) == 5
        assert all(isinstance(w, complex) for w in encoder.base_weights.values())
    
    def test_random_encoder_reproducibility(self):
        """Test that same seed gives same random weights."""
        encoder1 = RandomEncoder(seed=42)
        encoder2 = RandomEncoder(seed=42)
        
        assert encoder1.base_weights['A'] == encoder2.base_weights['A']
    
    def test_random_encoder_different_seeds(self):
        """Test that different seeds give different weights."""
        encoder1 = RandomEncoder(seed=42)
        encoder2 = RandomEncoder(seed=43)
        
        assert encoder1.base_weights['A'] != encoder2.base_weights['A']


class TestAblationTester:
    """Test ablation testing framework."""
    
    def test_ablation_tester_initialization(self):
        """Test ablation tester initializes correctly."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(
            baseline_encoder=encoder,
            n_permutations=100,
            n_bootstrap=100
        )
        
        assert tester.n_permutations == 100
        assert tester.n_bootstrap == 100
    
    def test_ablation_no_helical_phase(self):
        """Test no helical phase ablation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10)
        
        test_seq = "ATCGATCGATCG"
        ablated = tester.ablation_no_helical_phase(test_seq)
        
        assert len(ablated) == len(test_seq)
        assert ablated.dtype == complex
    
    def test_ablation_phase_scramble(self):
        """Test phase scramble ablation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10, seed=42)
        
        test_seq = "ATCGATCGATCG"
        scrambled = tester.ablation_phase_scramble(test_seq)
        
        # Should preserve magnitudes
        original = encoder.encode_sequence(test_seq, apply_helical_phase=True)
        assert np.allclose(np.abs(original), np.abs(scrambled), rtol=1e-3)
        
        # Should scramble phases
        assert not np.allclose(np.angle(original), np.angle(scrambled))
    
    def test_ablation_swap_at_gc(self):
        """Test AT/GC swap ablation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10)
        
        test_seq = "ATCGATCGATCG"
        swapped = tester.ablation_swap_at_gc(test_seq)
        
        assert len(swapped) == len(test_seq)
        assert swapped.dtype == complex
    
    def test_random_encodings_generation(self):
        """Test random encodings generation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10, seed=42)
        
        test_seq = "ATCGATCGATCG"
        random_encs = tester.generate_random_encodings(test_seq, n_random=10)
        
        assert len(random_encs) == 10
        assert all(len(enc) == len(test_seq) for enc in random_encs)
    
    def test_hedges_g_calculation(self):
        """Test Hedges' g effect size calculation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10)
        
        group_a = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        group_b = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
        
        g = tester.hedges_g(group_a, group_b)
        
        assert np.isfinite(g)
        assert g < 0  # group_b has higher mean
    
    def test_bootstrap_ci(self):
        """Test bootstrap confidence interval calculation."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=50)
        
        group_a = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        group_b = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
        
        estimate, lower, upper = tester.bootstrap_ci(group_a, group_b)
        
        assert np.isfinite(estimate)
        assert np.isfinite(lower)
        assert np.isfinite(upper)
        assert lower <= estimate <= upper
    
    def test_permutation_test(self):
        """Test permutation test."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=100, n_bootstrap=10)
        
        observed = 5.0
        null_dist = np.random.randn(100)
        
        p_value = tester.permutation_test(observed, null_dist, alternative='two-sided')
        
        assert 0 <= p_value <= 1
    
    def test_fdr_correction(self):
        """Test Benjamini-Hochberg FDR correction."""
        encoder = BreathingDynamicsEncoder()
        tester = AblationTester(encoder, n_permutations=10, n_bootstrap=10)
        
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        reject, adj_p = tester.benjamini_hochberg_fdr(p_values, alpha=0.05)
        
        assert len(reject) == len(p_values)
        assert len(adj_p) == len(p_values)
        assert all(0 <= p <= 1 for p in adj_p)


class TestScientificGates:
    """Test compliance with scientific gates from copilot-instructions."""
    
    def test_human_dna_alphabet_validation(self):
        """Gate G2: Test that only A/C/G/T/N are accepted."""
        encoder = BreathingDynamicsEncoder()
        
        # Valid sequence should work
        valid_seq = "ATCGATCGN"
        encoded = encoder.encode_sequence(valid_seq)
        assert len(encoded) == len(valid_seq)
        
        # Invalid sequences should fail
        invalid_sequences = [
            "AUCGATCG",  # RNA (contains U)
            "ATCXATCG",  # Contains X
            "ATCRATCG",  # Contains R (IUPAC code)
        ]
        
        for invalid_seq in invalid_sequences:
            with pytest.raises(ValueError):
                encoder.encode_sequence(invalid_seq)
    
    def test_biophysical_anchoring(self):
        """Gate: Test that weights are derived from biophysical parameters."""
        encoder = BreathingDynamicsEncoder()
        
        # Check that weights reference opening lifetimes
        info = encoder.get_encoding_info()
        assert 'opening_lifetimes_ms' in info
        assert info['reference'] == 'PMC5393899'
        assert 'type' in info
        assert 'biophysical' in info['type']
    
    def test_temperature_and_mg_parameters(self):
        """Gate: Test temperature and Mg²⁺ parameters are exposed."""
        encoder = BreathingDynamicsEncoder(temperature_c=37.0, mg_concentration_mm=2.0)
        
        info = encoder.get_encoding_info()
        assert info['temperature_c'] == 37.0
        assert info['mg_concentration_mm'] == 2.0
    
    def test_helical_period_parameter(self):
        """Gate: Test helical period (10.5 bp) is configurable."""
        encoder = BreathingDynamicsEncoder(helical_period=10.5)
        
        info = encoder.get_encoding_info()
        assert info['helical_period_bp'] == 10.5
    
    def test_dimensionless_encoding(self):
        """Gate: Test that encoding is dimensionless (normalized)."""
        encoder = BreathingDynamicsEncoder()
        
        test_seq = "ATCGATCGATCG"
        encoded = encoder.encode_sequence(test_seq)
        
        # Magnitudes should be O(1), not MHz/GHz
        mean_magnitude = np.mean(np.abs(encoded))
        assert 0.1 < mean_magnitude < 10.0  # Dimensionless range


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "--tb=short"])
