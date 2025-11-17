"""
Unit Tests for Spectral Disruption Profiler

Tests all core modules with scientific gate validation.
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from experiments.spectral_disruption_profiler.encoding import (
    validate_dna_sequence,
    encode_sequence,
    phase_weighted_encoding,
    geometric_phase_weight,
    batch_encode
)
from experiments.spectral_disruption_profiler.analysis import (
    compute_spectral_features,
    analyze_disruption,
    compute_fft_spectrum,
    compute_fundamental_frequency,
    compute_spectral_entropy,
    count_sidelobes
)
from experiments.spectral_disruption_profiler.scoring import (
    compute_z_score,
    compute_composite_score,
    bootstrap_confidence_interval,
    score_with_confidence,
    auto_optimize_k
)
from experiments.spectral_disruption_profiler.detection import (
    detect_off_targets,
    compute_gc_resonance,
    benjamini_hochberg_fdr,
    effect_size_cohens_d
)


class TestEncoding:
    """Test encoding module."""
    
    def test_validate_dna_sequence_valid(self):
        """Test validation accepts valid DNA sequences."""
        validate_dna_sequence("ATCG")
        validate_dna_sequence("ATCGATCGATCG")
        validate_dna_sequence("NNNATCGNN")  # N allowed
    
    def test_validate_dna_sequence_invalid(self):
        """Test validation rejects invalid DNA sequences."""
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_dna_sequence("ATCGXYZ")
        
        with pytest.raises(ValueError, match="Invalid DNA sequence"):
            validate_dna_sequence("AUCG")  # U not allowed for DNA
    
    def test_validate_rna_sequence_valid(self):
        """Test validation accepts valid RNA sequences."""
        validate_dna_sequence("AUCG", is_rna=True)
        validate_dna_sequence("AUCGAUCGAUCG", is_rna=True)
    
    def test_validate_rna_sequence_invalid(self):
        """Test validation rejects invalid RNA sequences."""
        with pytest.raises(ValueError, match="Invalid RNA sequence"):
            validate_dna_sequence("ATCG", is_rna=True)  # T not allowed for RNA
    
    def test_encode_sequence_dna(self):
        """Test DNA encoding produces correct waveform."""
        waveform = encode_sequence("ATCG")
        
        assert len(waveform) == 4
        assert waveform.dtype == np.complex128
        
        # Check encoding: A→1+0i, T→-1+0i, C→0+1i, G→0-1i
        np.testing.assert_almost_equal(waveform[0], 1.0 + 0.0j)
        np.testing.assert_almost_equal(waveform[1], -1.0 + 0.0j)
        np.testing.assert_almost_equal(waveform[2], 0.0 + 1.0j)
        np.testing.assert_almost_equal(waveform[3], 0.0 - 1.0j)
    
    def test_encode_sequence_rna(self):
        """Test RNA encoding produces correct waveform."""
        waveform = encode_sequence("AUCG", is_rna=True)
        
        assert len(waveform) == 4
        # A→1+0i, U→-1+0i, C→0+1i, G→0-1i
        np.testing.assert_almost_equal(waveform[0], 1.0 + 0.0j)
        np.testing.assert_almost_equal(waveform[1], -1.0 + 0.0j)
    
    def test_geometric_phase_weight(self):
        """Test geometric phase weighting function."""
        # Test k=0.3 (optimal)
        weight = geometric_phase_weight(0, phi=21.0, k=0.3)
        assert weight >= 0
        
        weight = geometric_phase_weight(10, phi=21.0, k=0.3)
        assert weight >= 0
    
    def test_phase_weighted_encoding(self):
        """Test phase-weighted encoding."""
        seq = "ATCGATCGATCG"
        waveform = phase_weighted_encoding(seq, phi=21.0, k=0.3)
        
        assert len(waveform) == len(seq)
        assert waveform.dtype == np.complex128
        
        # Should be different from standard encoding
        standard = encode_sequence(seq)
        assert not np.allclose(waveform, standard)
    
    def test_batch_encode(self):
        """Test batch encoding."""
        sequences = ["ATCG", "GCTA", "TTAA"]
        waveforms = batch_encode(sequences, use_phase_weighting=False)
        
        assert len(waveforms) == 3
        assert all(len(w) == 4 for w in waveforms)


class TestAnalysis:
    """Test analysis module."""
    
    def test_compute_fft_spectrum(self):
        """Test FFT spectrum computation."""
        waveform = np.array([1+0j, -1+0j, 1+0j, -1+0j])
        freqs, spectrum = compute_fft_spectrum(waveform)
        
        assert len(freqs) == len(spectrum)
        assert len(spectrum) == 4
        assert np.all(spectrum >= 0)  # Magnitude is non-negative
    
    def test_compute_fundamental_frequency(self):
        """Test fundamental frequency extraction."""
        # Simple periodic signal
        t = np.linspace(0, 1, 100, endpoint=False)
        signal = np.exp(2j * np.pi * 5 * t)  # 5 Hz signal
        
        freqs, spectrum = compute_fft_spectrum(signal)
        f1 = compute_fundamental_frequency(spectrum, freqs)
        
        # Should detect peak near 5 Hz (normalized frequency ≈ 0.05)
        assert f1 > 0
    
    def test_compute_spectral_entropy(self):
        """Test spectral entropy computation."""
        # Uniform spectrum → high entropy
        uniform_spectrum = np.ones(100)
        entropy_uniform = compute_spectral_entropy(uniform_spectrum)
        
        # Single peak → low entropy
        peaked_spectrum = np.zeros(100)
        peaked_spectrum[50] = 1.0
        entropy_peaked = compute_spectral_entropy(peaked_spectrum)
        
        assert entropy_uniform > entropy_peaked
    
    def test_count_sidelobes(self):
        """Test sidelobe counting."""
        # Create spectrum with known peaks
        spectrum = np.zeros(100)
        spectrum[10] = 1.0
        spectrum[30] = 0.8
        spectrum[50] = 0.6
        
        freqs = np.linspace(-0.5, 0.5, 100)
        n_lobes = count_sidelobes(spectrum, freqs, threshold_percentile=50)
        
        assert n_lobes >= 0
    
    def test_compute_spectral_features(self):
        """Test spectral feature extraction."""
        waveform = encode_sequence("ATCGATCGATCG")
        features = compute_spectral_features(waveform)
        
        assert 'f1' in features
        assert 'entropy' in features
        assert 'sidelobes' in features
        assert 'length' in features
        
        assert features['length'] == 12
    
    def test_analyze_disruption(self):
        """Test disruption analysis."""
        ref_seq = "ATCGATCGATCG"
        mut_seq = "ATCGATCGATCG"
        
        features = analyze_disruption(mut_seq, ref_seq)
        
        assert 'delta_f1' in features
        assert 'delta_entropy' in features
        assert 'delta_sidelobes' in features
        assert 'gc_content' in features


class TestScoring:
    """Test scoring module."""
    
    def test_compute_z_score(self):
        """Test Z-invariant score computation."""
        z = compute_z_score(A=1.0, B=7.389)  # B ≈ e²
        assert np.isclose(z, 1.0, rtol=0.01)
        
        z = compute_z_score(A=2.0, B=14.778)
        assert np.isclose(z, 4.0, rtol=0.01)  # 2 * (14.778 / 7.389) ≈ 4
    
    def test_compute_composite_score(self):
        """Test composite score computation."""
        features = {
            'delta_entropy': 0.5,
            'delta_f1': 0.1,
            'delta_sidelobes': 2,
            'gc_content': 0.6
        }
        
        score = compute_composite_score(features)
        
        assert isinstance(score, float)
        assert score >= 0
    
    def test_bootstrap_confidence_interval(self):
        """Test bootstrap CI computation."""
        data = np.random.randn(100)
        mean, lower, upper = bootstrap_confidence_interval(
            data,
            np.mean,
            n_bootstrap=100,
            seed=42
        )
        
        assert lower <= mean <= upper
        assert isinstance(mean, float)
    
    def test_score_with_confidence(self):
        """Test scoring with confidence intervals."""
        features_list = [
            {'delta_entropy': 0.5, 'delta_f1': 0.1, 'delta_sidelobes': 2, 'gc_content': 0.6}
            for _ in range(10)
        ]
        
        results = score_with_confidence(
            features_list,
            n_bootstrap=100,
            seed=42
        )
        
        assert 'mean_score' in results
        assert 'ci_lower' in results
        assert 'ci_upper' in results
        assert results['ci_lower'] <= results['mean_score'] <= results['ci_upper']
    
    def test_auto_optimize_k(self):
        """Test k parameter auto-optimization."""
        seq = "ATCGATCGATCGATCGATCG"
        optimal_k = auto_optimize_k(seq)
        
        assert 0.0 <= optimal_k <= 1.0
        assert isinstance(optimal_k, float)


class TestDetection:
    """Test detection module."""
    
    def test_detect_off_targets(self):
        """Test off-target detection."""
        features_list = [
            {'delta_entropy': 0.01, 'gc_content': 0.5},  # Normal
            {'delta_entropy': 0.8, 'gc_content': 0.5},   # High entropy → flag
            {'delta_entropy': 0.01, 'gc_content': 0.8},  # High GC → flag
        ]
        
        flags = detect_off_targets(features_list, entropy_threshold=0.05)
        
        assert len(flags) == 3
        assert flags[0] == False  # Normal - low entropy and normal GC
        assert flags[1] == True   # High entropy
        assert flags[2] == True   # High GC
    
    def test_compute_gc_resonance(self):
        """Test GC-resonance computation."""
        features_list = [
            {'gc_content': 0.3, 'delta_entropy': 0.5},
            {'gc_content': 0.5, 'delta_entropy': 0.3},
            {'gc_content': 0.7, 'delta_entropy': 0.1},
        ]
        
        results = compute_gc_resonance(
            features_list,
            n_permutations=100,
            seed=42
        )
        
        assert 'r' in results
        assert 'p_permutation' in results
        assert -1 <= results['r'] <= 1
    
    def test_benjamini_hochberg_fdr(self):
        """Test BH-FDR correction."""
        p_values = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        adjusted_p, significant = benjamini_hochberg_fdr(p_values, alpha=0.05)
        
        assert len(adjusted_p) == len(p_values)
        assert len(significant) == len(p_values)
        assert np.all(adjusted_p >= p_values)  # Adjusted p should be >= original
    
    def test_effect_size_cohens_d(self):
        """Test Cohen's d effect size."""
        group1 = np.array([1, 2, 3, 4, 5])
        group2 = np.array([3, 4, 5, 6, 7])
        
        d = effect_size_cohens_d(group1, group2)
        
        assert isinstance(d, float)
        # Group2 has higher mean, so d should be negative
        assert d < 0


class TestIntegration:
    """Integration tests for end-to-end workflow."""
    
    def test_full_pipeline_dna(self):
        """Test complete pipeline with DNA sequences."""
        ref_seq = "ATCGATCGATCGATCGATCG"
        mut_seq = "ATCGATCGATCGATCGATCG"
        
        # 1. Analyze disruption
        features = analyze_disruption(mut_seq, ref_seq)
        
        # 2. Compute score
        score = compute_composite_score(features)
        
        # 3. Detect off-targets
        flags = detect_off_targets([features])
        
        assert isinstance(score, float)
        assert len(flags) == 1
        assert isinstance(flags[0], bool)
    
    def test_full_pipeline_rna(self):
        """Test complete pipeline with RNA sequences."""
        ref_seq = "AUCGAUCGAUCGAUCGAUCG"
        mut_seq = "AUCGAUCGAUCGAUCGAUCG"
        
        features = analyze_disruption(mut_seq, ref_seq, is_rna=True)
        score = compute_composite_score(features)
        
        assert isinstance(score, float)
    
    def test_batch_workflow(self):
        """Test batch processing workflow."""
        ref_seqs = ["ATCGATCGATCGATCG" for _ in range(5)]
        mut_seqs = ["GCTAGCTAGCTAGCTA" for _ in range(5)]
        
        # Batch analyze
        from experiments.spectral_disruption_profiler.analysis import batch_analyze_disruption
        features_list = batch_analyze_disruption(mut_seqs, ref_seqs)
        
        # Score with CI
        results = score_with_confidence(features_list, n_bootstrap=100, seed=42)
        
        # GC resonance
        gc_results = compute_gc_resonance(features_list, n_permutations=100, seed=42)
        
        assert len(features_list) == 5
        assert 'mean_score' in results
        assert 'r' in gc_results


# Smoke test marker
@pytest.mark.smoke
class TestSmoke:
    """Quick smoke tests for CI (< 5s)."""
    
    def test_import_modules(self):
        """Test all modules import successfully."""
        from experiments.spectral_disruption_profiler import (
            encode_sequence,
            phase_weighted_encoding,
            compute_spectral_features,
            analyze_disruption,
            compute_z_score,
            compute_composite_score,
            detect_off_targets,
            compute_gc_resonance
        )
    
    def test_encode_simple(self):
        """Test simple encoding."""
        waveform = encode_sequence("ATCG")
        assert len(waveform) == 4
    
    def test_analyze_simple(self):
        """Test simple analysis."""
        features = analyze_disruption("ATCGATCG", "ATCGATCG")
        assert 'delta_entropy' in features


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
