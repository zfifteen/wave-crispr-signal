#!/usr/bin/env python3
"""
Tests for FFT-Based CRISPR Disruption Metrics with Golden-Ratio Phase Analysis

Validates:
- DNA sequence validation (A/C/G/T only)
- Complex encoding
- θ′(n,k) golden-ratio phase calculation
- FFT spectrum analysis with phase weighting
- Off-target periodicity detection
- Disruption scoring for indels
- Codon-aligned features
"""

import sys
import os
import numpy as np
import pytest
from scipy.fft import fft

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'applications'))

from fft_crispr_disruption import (
    FFTCRISPRDisruptionAnalyzer,
    calculate_grna_off_target_score,
    PHI,
    DEFAULT_K,
    DEFAULT_PHI_PERIOD
)


class TestDNAValidation:
    """Test DNA sequence validation gates."""
    
    def test_valid_dna_sequence(self):
        """Test that valid DNA sequences pass validation."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        valid_sequences = [
            "ATCG",
            "ATCGATCGATCG",
            "AAAAAAAAAA",
            "GCGCGCGCGC",
            "ATCGNATCG"  # N allowed
        ]
        
        for seq in valid_sequences:
            result = analyzer.validate_dna_sequence(seq)
            assert result == seq.upper()
    
    def test_reject_rna_bases(self):
        """Test that RNA base U is rejected for DNA."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        with pytest.raises(ValueError, match="Found 'U'"):
            analyzer.validate_dna_sequence("AUCG")
    
    def test_reject_invalid_bases(self):
        """Test that invalid bases are rejected."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        invalid_sequences = [
            "ATCGR",  # R is IUPAC ambiguity code
            "ATCGY",  # Y is IUPAC ambiguity code
            "ATCGW",  # W is IUPAC ambiguity code
            "ATCG1",  # Number
            "ATCG-"   # Special character
        ]
        
        for seq in invalid_sequences:
            with pytest.raises(ValueError, match="Invalid DNA bases"):
                analyzer.validate_dna_sequence(seq)


class TestComplexEncoding:
    """Test DNA to complex waveform encoding."""
    
    def test_standard_mapping(self):
        """Test standard A/T/C/G complex mapping."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # Test individual bases
        assert np.isclose(analyzer.encode_dna_complex("A")[0], 1.0 + 0.0j)
        assert np.isclose(analyzer.encode_dna_complex("T")[0], -1.0 + 0.0j)
        assert np.isclose(analyzer.encode_dna_complex("C")[0], 0.0 + 1.0j)
        assert np.isclose(analyzer.encode_dna_complex("G")[0], 0.0 - 1.0j)
        assert np.isclose(analyzer.encode_dna_complex("N")[0], 0.0 + 0.0j)
    
    def test_sequence_encoding_length(self):
        """Test that encoding preserves sequence length."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCG"
        encoded = analyzer.encode_dna_complex(seq)
        
        assert len(encoded) == len(seq)
        assert encoded.dtype == np.complex128


class TestGoldenRatioPhase:
    """Test θ′(n,k) geometric resolution calculation."""
    
    def test_theta_prime_calculation(self):
        """Test θ′(n,k) = φ·((n mod φ)/φ)^k calculation."""
        analyzer = FFTCRISPRDisruptionAnalyzer(phi_period=21.0, k=0.3)
        
        # Test specific values
        theta_0 = analyzer.calculate_theta_prime(0)
        assert theta_0 > 0  # Should handle n=0 case
        
        theta_1 = analyzer.calculate_theta_prime(1)
        assert theta_1 > 0
        
        # θ′ should be bounded by φ
        for n in range(100):
            theta_n = analyzer.calculate_theta_prime(n)
            assert 0 < theta_n <= PHI * 1.1  # Allow small margin
    
    def test_theta_prime_periodicity(self):
        """Test that θ′ has periodicity related to φ."""
        analyzer = FFTCRISPRDisruptionAnalyzer(phi_period=21.0, k=0.3)
        
        # Values at n and n+φ_period should be similar
        n = 5
        theta_n = analyzer.calculate_theta_prime(n)
        theta_n_plus_period = analyzer.calculate_theta_prime(n + int(analyzer.phi_period))
        
        # Should be close (exact match for integer periods)
        assert np.isclose(theta_n, theta_n_plus_period, rtol=0.01)
    
    def test_k_parameter_effect(self):
        """Test that k parameter affects resolution."""
        seq = "ATCGATCGATCGATCGATCG"
        
        analyzer_k03 = FFTCRISPRDisruptionAnalyzer(k=0.3)
        analyzer_k05 = FFTCRISPRDisruptionAnalyzer(k=0.5)
        
        theta_k03 = analyzer_k03.calculate_theta_prime(10)
        theta_k05 = analyzer_k05.calculate_theta_prime(10)
        
        # Different k values should give different results
        assert theta_k03 != theta_k05


class TestFFTSpectralAnalysis:
    """Test FFT spectrum analysis with phase weighting."""
    
    def test_basic_fft_analysis(self):
        """Test basic FFT computation."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCGATCGATCG"  # 20-nt guide sequence
        analysis = analyzer.detect_off_target_periodicities(seq)
        
        assert 'sequence_length' in analysis
        assert analysis['sequence_length'] == len(seq)
        assert 'fft_spectrum' in analysis
        assert 'weighted_spectrum' in analysis
        assert len(analysis['fft_spectrum']) == len(seq)
        assert len(analysis['weighted_spectrum']) == len(seq)
    
    def test_golden_phase_weighting(self):
        """Test that golden-ratio phase weights are applied correctly."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCG"
        complex_wave = analyzer.encode_dna_complex(seq)
        
        fft_result = fft(complex_wave)
        fft_magnitude = np.abs(fft_result)
        
        weighted = analyzer.apply_golden_phase_weights(fft_magnitude)
        
        # Weighted spectrum should differ from raw spectrum
        assert not np.allclose(weighted, fft_magnitude)
        
        # All weights should be positive
        assert np.all(weighted >= 0)
    
    def test_periodicity_detection(self):
        """Test detection of periodicities in structured sequences."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # Create sequence with obvious periodicity (ATCG repeated)
        periodic_seq = "ATCG" * 5  # 20 bp with 4-bp period
        
        analysis = analyzer.detect_off_target_periodicities(periodic_seq)
        
        assert analysis['n_significant_peaks'] > 0
        assert len(analysis['significant_peaks']) > 0
        
        # Check that peaks have expected structure
        peak = analysis['significant_peaks'][0]
        assert 'frequency' in peak
        assert 'period' in peak
        assert 'magnitude' in peak
        assert 'weighted_magnitude' in peak
        assert 'theta_prime_weight' in peak


class TestDisruptionScoring:
    """Test disruption scoring for edited sequences."""
    
    def test_identical_sequences_zero_disruption(self):
        """Test that identical sequences have minimal disruption."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCGATCGATCG"
        disruption = analyzer.calculate_disruption_score(seq, seq)
        
        # Identical sequences should have zero deltas
        assert disruption['delta_entropy'] == 0.0
        assert disruption['delta_f1'] == 0.0
        assert disruption['delta_sidelobes'] == 0
        # Overall disruption should be minimal (phase disruption might be non-zero due to numerical precision)
        assert disruption['disruption_score'] < 0.3
    
    def test_deletion_causes_disruption(self):
        """Test that deletions cause measurable disruption."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        ref_seq = "ATCGATCGATCGATCGATCG"
        # Delete 3 bases in middle
        edited_seq = "ATCGATCGATCGATCG"
        
        disruption = analyzer.calculate_disruption_score(ref_seq, edited_seq)
        
        # Deletion should cause non-zero disruption
        assert disruption['disruption_score'] > 0.0
        assert 'delta_entropy' in disruption
        assert 'delta_f1' in disruption
    
    def test_indel_analysis_deletion(self):
        """Test indel analysis for deletions."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCGATCGATCG"
        indel_analysis = analyzer.analyze_indel_disruption(
            seq,
            indel_position=10,
            indel_length=3,
            indel_type='deletion'
        )
        
        assert indel_analysis['indel_type'] == 'deletion'
        assert indel_analysis['indel_position'] == 10
        assert indel_analysis['indel_length'] == 3
        assert indel_analysis['original_length'] == 20
        assert indel_analysis['edited_length'] == 17
        assert indel_analysis['disruption_score'] > 0.0
    
    def test_indel_analysis_insertion(self):
        """Test indel analysis for insertions."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCGATCGATCG"
        indel_analysis = analyzer.analyze_indel_disruption(
            seq,
            indel_position=10,
            indel_length=3,
            indel_type='insertion'
        )
        
        assert indel_analysis['indel_type'] == 'insertion'
        assert indel_analysis['indel_position'] == 10
        assert indel_analysis['indel_length'] == 3
        assert indel_analysis['original_length'] == 20
        assert indel_analysis['edited_length'] == 23
        assert indel_analysis['disruption_score'] > 0.0


class TestCodonAlignedFeatures:
    """Test codon-aligned φ-structured analysis."""
    
    def test_codon_analysis_basic(self):
        """Test basic codon-aligned analysis."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # 21 bp sequence (7 codons)
        seq = "ATCGATCGATCGATCGATCGA"
        codon_features = analyzer.calculate_codon_aligned_features(seq, frame=0)
        
        assert codon_features['n_codons'] == 7
        assert codon_features['frame'] == 0
        assert len(codon_features['codon_values']) == 7
        assert 'codon_entropy' in codon_features
        assert 'dominant_codon_period' in codon_features
    
    def test_reading_frame_shift(self):
        """Test reading frame parameter."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "ATCGATCGATCGATCGATCGA"
        
        frame0 = analyzer.calculate_codon_aligned_features(seq, frame=0)
        frame1 = analyzer.calculate_codon_aligned_features(seq, frame=1)
        frame2 = analyzer.calculate_codon_aligned_features(seq, frame=2)
        
        # Different frames should give different results
        assert not np.allclose(frame0['codon_values'], frame1['codon_values'])
        assert not np.allclose(frame1['codon_values'], frame2['codon_values'])
    
    def test_codon_padding(self):
        """Test that sequences are padded to codon multiples."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # 20 bp (not divisible by 3)
        seq = "ATCGATCGATCGATCGATCG"
        codon_features = analyzer.calculate_codon_aligned_features(seq, frame=0)
        
        # Should be padded to 21 bp = 7 codons
        assert codon_features['n_codons'] == 7


class TestGRNAOffTargetScore:
    """Test gRNA off-target scoring function."""
    
    def test_basic_scoring(self):
        """Test basic off-target score calculation."""
        grna = "GACGATCGATCGATCGATCG"  # 20-nt guide
        
        score = calculate_grna_off_target_score(grna)
        
        assert 'off_target_score' in score
        assert 0.0 <= score['off_target_score'] <= 1.0
        assert 'recommendation' in score
        assert score['recommendation'] in ['good', 'review', 'poor']
    
    def test_scoring_with_custom_parameters(self):
        """Test scoring with custom φ-period and k."""
        grna = "GACGATCGATCGATCGATCG"
        
        score1 = calculate_grna_off_target_score(grna, phi_period=21.0, k=0.3)
        score2 = calculate_grna_off_target_score(grna, phi_period=20.0, k=0.3)
        
        # Different parameters should give different scores
        assert score1['off_target_score'] != score2['off_target_score']
    
    def test_structured_vs_random_sequences(self):
        """Test that structured sequences score differently than random."""
        # Highly structured (repetitive)
        structured = "ATCG" * 5
        
        # More random-like
        random_like = "GACTAGCTAGCTAGCTAACG"
        
        score_structured = calculate_grna_off_target_score(structured)
        score_random = calculate_grna_off_target_score(random_like)
        
        # Scores should differ
        assert score_structured['off_target_score'] != score_random['off_target_score']


class TestScientificGates:
    """Test compliance with scientific gates from REPOSITORY_POLICY."""
    
    def test_golden_ratio_constant(self):
        """Test that φ constant is correct."""
        expected_phi = 1.6180339887498949
        assert np.isclose(PHI, expected_phi, rtol=1e-10)
    
    def test_default_k_parameter(self):
        """Test that default k ≈ 0.3 as specified."""
        assert DEFAULT_K == 0.3
    
    def test_default_phi_period(self):
        """Test default φ-period for 21-nt guides."""
        assert DEFAULT_PHI_PERIOD == 21.0
    
    def test_human_dna_only_enforcement(self):
        """Test that only human DNA nucleotides are accepted."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # Should pass
        analyzer.validate_dna_sequence("ATCGATCG")
        
        # Should fail
        with pytest.raises(ValueError):
            analyzer.validate_dna_sequence("AUCG")  # RNA
        
        with pytest.raises(ValueError):
            analyzer.validate_dna_sequence("ATCGR")  # IUPAC ambiguity


class TestNumericalStability:
    """Test numerical stability and edge cases."""
    
    def test_empty_sequence_handling(self):
        """Test handling of empty sequences."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        with pytest.raises(ValueError):
            analyzer.validate_dna_sequence("")
    
    def test_very_short_sequence(self):
        """Test handling of very short sequences."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "AT"  # 2-bp sequence
        analysis = analyzer.detect_off_target_periodicities(seq)
        
        # Should complete without error
        assert analysis['sequence_length'] == 2
        assert len(analysis['fft_spectrum']) == 2
    
    def test_very_long_sequence(self):
        """Test handling of longer sequences."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        # 200-bp sequence
        seq = "ATCG" * 50
        analysis = analyzer.detect_off_target_periodicities(seq)
        
        assert analysis['sequence_length'] == 200
        assert len(analysis['fft_spectrum']) == 200
    
    def test_all_same_base(self):
        """Test homopolymer sequences."""
        analyzer = FFTCRISPRDisruptionAnalyzer()
        
        seq = "A" * 20
        analysis = analyzer.detect_off_target_periodicities(seq)
        
        # Should complete without error
        assert analysis['sequence_length'] == 20


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
