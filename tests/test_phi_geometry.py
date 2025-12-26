"""
Unit tests for the φ-geometry module.

Tests cover:
1. Constants and their physical validity
2. DNA sequence validation
3. Phase computation functions
4. Curvature computation functions
5. Feature score functions
6. Baseline comparison functions
"""

import pytest
import numpy as np
import math
import sys
from pathlib import Path

# Add repo root to path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from wave_crispr_signal.features.phi_geometry import (
    # Constants
    PHI,
    DNA_LENGTH_DIAMETER_RATIO,
    DNA_MAJOR_MINOR_SEP_RATIO,
    DNA_HELICAL_PERIOD,
    PHI_MODES,
    VALID_DNA_BASES,
    # Validation
    validate_dna_sequence,
    # Phase functions
    dna_phi_phase,
    dna_phi_phase_vectorized,
    # Curvature functions
    dna_phi_curvature,
    dna_phi_curvature_windowed,
    # Score functions
    phi_phase_score,
    phi_curvature_score,
    compute_phi_features,
    # Baseline functions
    uniform_phase_score,
    random_phase_score,
    simple_gc_content,
)


# ============================================================================
# CONSTANTS TESTS
# ============================================================================

class TestConstants:
    """Test physical constants and their validity."""
    
    def test_phi_is_golden_ratio(self):
        """φ should equal (1 + √5) / 2."""
        expected = (1.0 + math.sqrt(5.0)) / 2.0
        assert abs(PHI - expected) < 1e-10
    
    def test_phi_approximately_1618(self):
        """φ ≈ 1.618."""
        assert 1.617 < PHI < 1.619
    
    def test_dna_length_diameter_ratio(self):
        """DNA L/D ratio should be ≈ 1.6088 (Larsen, 2021)."""
        assert abs(DNA_LENGTH_DIAMETER_RATIO - 1.6088) < 0.001
    
    def test_dna_major_minor_sep_ratio(self):
        """DNA major/minor separation should be ≈ 1.6407."""
        assert abs(DNA_MAJOR_MINOR_SEP_RATIO - 1.6407) < 0.001
    
    def test_dna_helical_period(self):
        """Helical period should be 10 bp."""
        assert DNA_HELICAL_PERIOD == 10
    
    def test_phi_modes_include_fundamental(self):
        """PHI_MODES should include 10 bp fundamental."""
        assert 10.0 in PHI_MODES
    
    def test_phi_modes_include_phi_scaled(self):
        """PHI_MODES should include φ-scaled harmonics."""
        # 10/φ ≈ 6.18, 10*φ ≈ 16.18
        phi_scaled = [m for m in PHI_MODES if abs(m - 10.0/PHI) < 0.1 or abs(m - 10.0*PHI) < 0.1]
        assert len(phi_scaled) >= 2


# ============================================================================
# VALIDATION TESTS
# ============================================================================

class TestValidation:
    """Test DNA sequence validation."""
    
    def test_valid_dna_sequence(self):
        """Valid DNA should pass validation."""
        seq = validate_dna_sequence("ATCGATCGATCG")
        assert seq == "ATCGATCGATCG"
    
    def test_lowercase_converted(self):
        """Lowercase should be converted to uppercase."""
        seq = validate_dna_sequence("atcg")
        assert seq == "ATCG"
    
    def test_n_is_valid(self):
        """N (ambiguous) should be valid."""
        seq = validate_dna_sequence("ATCGN")
        assert seq == "ATCGN"
    
    def test_empty_sequence_raises(self):
        """Empty sequence should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            validate_dna_sequence("")
    
    def test_rna_u_raises(self):
        """RNA U should raise with helpful message."""
        with pytest.raises(ValueError, match="DNA only"):
            validate_dna_sequence("AUCG")
    
    def test_invalid_chars_raise(self):
        """Invalid characters should raise."""
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            validate_dna_sequence("ATCGX")


# ============================================================================
# PHASE FUNCTION TESTS
# ============================================================================

class TestPhaseFunctions:
    """Test φ-phase computation functions."""
    
    def test_phase_at_zero(self):
        """Phase at index 0 should be 0."""
        phase = dna_phi_phase(0)
        assert phase == 0.0
    
    def test_phase_range(self):
        """Phase should always be in [0, 1)."""
        for i in range(100):
            phase = dna_phi_phase(i)
            assert 0.0 <= phase < 1.0
    
    def test_phase_negative_raises(self):
        """Negative index should raise."""
        with pytest.raises(ValueError, match="non-negative"):
            dna_phi_phase(-1)
    
    def test_phase_vectorized_matches_scalar(self):
        """Vectorized phase should match scalar version."""
        indices = np.arange(20)
        phases_vec = dna_phi_phase_vectorized(indices)
        phases_scalar = [dna_phi_phase(i) for i in indices]
        
        np.testing.assert_array_almost_equal(phases_vec, phases_scalar)
    
    def test_phase_vectorized_negative_raises(self):
        """Vectorized phase with negative should raise."""
        with pytest.raises(ValueError, match="non-negative"):
            dna_phi_phase_vectorized(np.array([-1, 0, 1]))


# ============================================================================
# CURVATURE FUNCTION TESTS
# ============================================================================

class TestCurvatureFunctions:
    """Test φ-curvature computation functions."""
    
    def test_curvature_periodic_signal(self):
        """Periodic signal at 10 bp should have high curvature."""
        x = np.sin(2 * np.pi * np.arange(100) / 10)
        curv = dna_phi_curvature(x)
        # Should be in valid range
        assert 0.0 <= curv <= 1.0
        # Should be reasonably high for periodic
        assert curv > 0.3
    
    def test_curvature_random_signal(self):
        """Random signal should have lower curvature than periodic."""
        np.random.seed(42)
        x_periodic = np.sin(2 * np.pi * np.arange(100) / 10)
        x_random = np.random.randn(100)
        
        curv_periodic = dna_phi_curvature(x_periodic)
        curv_random = dna_phi_curvature(x_random)
        
        assert curv_periodic > curv_random
    
    def test_curvature_range(self):
        """Curvature should always be in [0, 1]."""
        for _ in range(10):
            x = np.random.randn(50)
            curv = dna_phi_curvature(x)
            assert 0.0 <= curv <= 1.0
    
    def test_curvature_short_signal(self):
        """Very short signal should return 0."""
        x = np.array([1.0, 2.0, 3.0])
        curv = dna_phi_curvature(x)
        assert curv == 0.0
    
    def test_curvature_constant_signal(self):
        """Constant signal should return 0."""
        x = np.ones(50)
        curv = dna_phi_curvature(x)
        assert curv == 0.0
    
    def test_curvature_windowed_shape(self):
        """Windowed curvature should have correct output shape."""
        x = np.sin(2 * np.pi * np.arange(50) / 10)
        window_size = 21
        curvatures = dna_phi_curvature_windowed(x, window_size)
        
        expected_len = len(x) - window_size + 1
        assert len(curvatures) == expected_len


# ============================================================================
# SCORE FUNCTION TESTS
# ============================================================================

class TestScoreFunctions:
    """Test φ-based score functions."""
    
    def test_phi_phase_score_range(self):
        """Phase score should be in [0, 1]."""
        seqs = ["ATCGATCGATCGATCGATCG", "AAAAAAAAAAAAAAAAAAAA", "GCGCGCGCGCGCGCGCGCGC"]
        for seq in seqs:
            score = phi_phase_score(seq)
            assert 0.0 <= score <= 1.0
    
    def test_phi_curvature_score_range(self):
        """Curvature score should be in [0, 1]."""
        seqs = ["ATCGATCGATCGATCGATCG", "AAAAAAAAAAAAAAAAAAAA", "GCGCGCGCGCGCGCGCGCGC"]
        for seq in seqs:
            score = phi_curvature_score(seq)
            assert 0.0 <= score <= 1.0
    
    def test_phi_curvature_homopolymer_low(self):
        """Homopolymer should have low curvature (no variation)."""
        score = phi_curvature_score("AAAAAAAAAAAAAAAAAAAA")
        assert score == 0.0
    
    def test_compute_phi_features_keys(self):
        """compute_phi_features should return expected keys."""
        features = compute_phi_features("ATCGATCGATCGATCGATCG")
        assert 'phi_phase_score' in features
        assert 'phi_curvature_score' in features
        assert 'phi_combined_score' in features
    
    def test_compute_phi_features_combined(self):
        """Combined score should be geometric mean of phase and curvature."""
        features = compute_phi_features("ATCGATCGATCGATCGATCG")
        # Geometric mean formula (with epsilon for safety)
        eps = 1e-10
        expected = math.sqrt((features['phi_phase_score'] + eps) * 
                            (features['phi_curvature_score'] + eps))
        assert abs(features['phi_combined_score'] - expected) < 1e-6


# ============================================================================
# BASELINE FUNCTION TESTS
# ============================================================================

class TestBaselineFunctions:
    """Test baseline comparison functions."""
    
    def test_uniform_phase_is_constant(self):
        """Uniform phase should return constant value."""
        seqs = ["ATCGATCGATCG", "AAAAAAAAAAAAAAAAAAAA", "GC" * 10]
        scores = [uniform_phase_score(seq) for seq in seqs]
        
        # All should be equal
        assert all(s == scores[0] for s in scores)
        assert scores[0] == 0.5
    
    def test_random_phase_reproducible(self):
        """Random phase should be reproducible with same seed."""
        seq = "ATCGATCGATCGATCGATCG"
        score1 = random_phase_score(seq, seed=42)
        score2 = random_phase_score(seq, seed=42)
        
        assert score1 == score2
    
    def test_random_phase_different_seeds(self):
        """Random phase should differ with different seeds."""
        seq = "ATCGATCGATCGATCGATCG"
        score1 = random_phase_score(seq, seed=42)
        score2 = random_phase_score(seq, seed=123)
        
        # Very unlikely to be exactly equal
        assert score1 != score2
    
    def test_gc_content_pure_gc(self):
        """Pure GC sequence should have GC content = 1."""
        score = simple_gc_content("GCGCGCGCGCGCGCGCGCGC")
        assert score == 1.0
    
    def test_gc_content_pure_at(self):
        """Pure AT sequence should have GC content = 0."""
        score = simple_gc_content("ATATATATATATATATATATAT")
        assert score == 0.0
    
    def test_gc_content_mixed(self):
        """Mixed sequence should have correct GC content."""
        # 10 GC, 12 AT = 22 total, GC fraction = 10/22 ≈ 0.4545
        score = simple_gc_content("GCGCGCGCGCATATATATATAT")
        expected = 10 / 22
        assert abs(score - expected) < 0.01


# ============================================================================
# INTEGRATION TESTS
# ============================================================================

class TestIntegration:
    """Integration tests for complete workflows."""
    
    def test_features_from_package_import(self):
        """Features should be importable from package __init__."""
        from wave_crispr_signal.features import (
            PHI, dna_phi_phase, compute_phi_features
        )
        assert PHI > 1.6
        assert callable(dna_phi_phase)
        assert callable(compute_phi_features)
    
    def test_real_grna_length_sequence(self):
        """Test with realistic 20-nt gRNA sequence."""
        grna = "ATCGATCGATCGATCGATCG"
        features = compute_phi_features(grna)
        
        # All scores should be valid
        assert 0.0 <= features['phi_phase_score'] <= 1.0
        assert 0.0 <= features['phi_curvature_score'] <= 1.0
        assert 0.0 <= features['phi_combined_score'] <= 1.0


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
