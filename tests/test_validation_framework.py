"""
Tests for Wave-CRISPR validation framework

This module tests the synthetic and real DNA validation scripts
to ensure proper functioning for FDA approval validation requirements.
"""

import pytest
import sys
import os
import numpy as np
from pathlib import Path

# Add parent directories for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'applications'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'proof_pack'))

from crispr_guide_designer import CRISPRGuideDesigner


class TestSyntheticValidation:
    """Test synthetic sequence validation functionality."""
    
    def test_sequence_generation_with_gc_content(self):
        """Test that synthetic sequences have correct GC content."""
        from validate_synthetic import generate_sequence_with_gc_content
        
        # Test different GC contents
        for target_gc in [0.25, 0.50, 0.75]:
            seq = generate_sequence_with_gc_content(100, target_gc)
            
            # Check length
            assert len(seq) == 100
            
            # Check GC content (allow 10% tolerance)
            actual_gc = (seq.count('G') + seq.count('C')) / len(seq)
            assert abs(actual_gc - target_gc) < 0.10
            
            # Check only valid bases
            assert all(base in 'ACGT' for base in seq)
    
    def test_wave_metrics_calculation(self):
        """Test wave metrics calculation for DNA sequences."""
        from validate_synthetic import calculate_wave_metrics
        
        designer = CRISPRGuideDesigner()
        
        # Test sequences with different characteristics
        test_sequences = [
            'AAAAAAAAAAAAAAAAAAAA',  # AT-rich
            'GGGGGGGGGGGGGGGGGGGG',  # GC-rich
            'ACGTACGTACGTACGTACGT',  # Balanced
        ]
        
        for seq in test_sequences:
            metrics = calculate_wave_metrics(seq, designer)
            
            # Check all required metrics are present
            assert 'wave_wobble' in metrics
            assert 'spectral_entropy' in metrics
            assert 'sidelobe_count' in metrics
            assert 'gc_content' in metrics
            
            # Check metric values are reasonable
            assert metrics['wave_wobble'] >= 0
            assert metrics['spectral_entropy'] >= 0
            assert metrics['sidelobe_count'] >= 0
            assert 0 <= metrics['gc_content'] <= 1
    
    def test_synthetic_sequence_generation(self):
        """Test generation of synthetic sequence categories."""
        from validate_synthetic import generate_synthetic_sequences
        
        num_sequences = 10
        sequences = generate_synthetic_sequences(num_sequences, seed=42)
        
        # Check all categories are present
        assert 'at_rich' in sequences
        assert 'gc_rich' in sequences
        assert 'balanced' in sequences
        
        # Check correct number of sequences
        assert len(sequences['at_rich']) == num_sequences
        assert len(sequences['gc_rich']) == num_sequences
        assert len(sequences['balanced']) == num_sequences
        
        # Check GC content ranges (use inclusive bounds for consistency)
        for seq in sequences['at_rich']:
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            assert 0.15 <= gc <= 0.40  # AT-rich (20-35% target)
        
        for seq in sequences['gc_rich']:
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            assert 0.60 <= gc <= 0.85  # GC-rich (65-80% target)
        
        for seq in sequences['balanced']:
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            assert 0.40 <= gc <= 0.60  # Balanced (45-55% target)
    
    def test_confidence_interval_calculation(self):
        """Test bootstrap confidence interval calculation."""
        from validate_synthetic import calculate_confidence_intervals
        
        # Test with known distribution
        data = np.random.normal(10.0, 2.0, 100)
        lower, upper = calculate_confidence_intervals(data, confidence=0.95, n_bootstrap=100)
        
        # Check CI bounds are reasonable
        assert lower < np.mean(data) < upper
        assert upper - lower > 0  # CI has positive width


class TestRealDNAValidation:
    """Test real DNA validation functionality."""
    
    def test_mutation_introduction(self):
        """Test point mutation introduction."""
        from validate_real_dna import introduce_mutation
        
        seq = 'ACGTACGT'
        
        # Test valid mutation
        mutated = introduce_mutation(seq, 2, 'T')
        assert mutated == 'ACTTACGT'
        assert len(mutated) == len(seq)
        
        # Test invalid position
        with pytest.raises(ValueError):
            introduce_mutation(seq, 10, 'A')
        
        # Test invalid base
        with pytest.raises(ValueError):
            introduce_mutation(seq, 2, 'X')
    
    def test_wave_disruption_calculation(self):
        """Test wave disruption metric calculation."""
        from validate_real_dna import calculate_wave_disruption
        
        designer = CRISPRGuideDesigner()
        
        original = 'ACGTACGTACGTACGTACGT'
        mutated = 'ACGTACGTACGTACGTACGG'  # T->G at last position
        
        disruption = calculate_wave_disruption(original, mutated, designer)
        
        # Check all metrics are present
        assert 'spectral_distance' in disruption
        assert 'phase_disruption' in disruption
        assert 'entropy_change' in disruption
        assert 'amplitude_disruption' in disruption
        assert 'sidelobe_change' in disruption
        assert 'composite_disruption' in disruption
        
        # Check metrics are non-negative
        assert disruption['spectral_distance'] >= 0
        assert disruption['phase_disruption'] >= 0
        assert disruption['entropy_change'] >= 0
        assert disruption['composite_disruption'] >= 0
    
    def test_mutation_effects_analysis(self):
        """Test mutation effects analysis on a sequence."""
        from validate_real_dna import analyze_mutation_effects
        
        designer = CRISPRGuideDesigner()
        
        # Use a longer test sequence
        sequence = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT' * 3
        seq_name = 'test_sequence'
        
        # Analyze with small number of mutations for speed
        results = analyze_mutation_effects(
            sequence, seq_name, designer,
            num_mutations=10, window_size=30
        )
        
        # Check result structure
        assert 'sequence_name' in results
        assert 'sequence_length' in results
        assert 'num_mutations' in results
        assert 'mutations' in results
        assert 'summary' in results
        
        # Check correct number of mutations
        assert len(results['mutations']) == 10
        
        # Check summary statistics
        summary = results['summary']
        assert 'mean_disruption' in summary
        assert 'std_disruption' in summary
        assert 'median_disruption' in summary
        assert 'max_disruption' in summary


class TestCRISPRCLIEnhancements:
    """Test enhancements to CRISPR CLI."""
    
    def test_biopython_availability(self):
        """Test that Biopython is available."""
        try:
            from Bio import SeqIO
            assert True
        except ImportError:
            pytest.skip("Biopython not installed")
    
    def test_fasta_reading_with_biopython(self):
        """Test FASTA reading with Biopython integration."""
        # This test requires a real FASTA file
        # Skip if not available
        test_fasta = Path('data/test_human_cdna.fasta')
        if not test_fasta.exists():
            pytest.skip("Test FASTA file not found")
        
        from crispr_cli import read_fasta
        
        sequences = read_fasta(str(test_fasta))
        
        # Check that sequences were loaded
        assert len(sequences) > 0
        
        # Check that sequences are uppercase
        for seq in sequences.values():
            assert seq.isupper()
            # Check only valid DNA bases
            assert all(base in 'ACGTN' for base in seq)


class TestValidationPipeline:
    """Test validation pipeline integration."""
    
    def test_designer_initialization(self):
        """Test that CRISPRGuideDesigner initializes correctly."""
        designer = CRISPRGuideDesigner()
        
        assert designer.pam_pattern is not None
        assert designer.guide_length == 20
    
    def test_waveform_generation(self):
        """Test waveform generation for DNA sequences."""
        designer = CRISPRGuideDesigner()
        
        sequence = 'ACGTACGTACGTACGTACGT'
        waveform = designer.build_waveform(sequence)
        
        # Check waveform shape
        assert len(waveform) == len(sequence)
        
        # Check waveform is complex
        assert np.iscomplexobj(waveform)
    
    def test_spectrum_computation(self):
        """Test FFT spectrum computation."""
        designer = CRISPRGuideDesigner()
        
        sequence = 'ACGTACGTACGTACGTACGT'
        waveform = designer.build_waveform(sequence)
        spectrum = designer.compute_spectrum(waveform)
        
        # Check spectrum is real-valued and non-negative
        assert np.all(spectrum >= 0)
        assert len(spectrum) == len(sequence)


class TestStatisticalMethods:
    """Test statistical methods used in validation."""
    
    def test_bootstrap_reproducibility(self):
        """Test that bootstrap results are reproducible with seed."""
        from validate_synthetic import calculate_confidence_intervals
        
        np.random.seed(42)
        data = np.random.normal(10.0, 2.0, 100)
        
        # Calculate CI twice with same seed
        np.random.seed(42)
        ci1 = calculate_confidence_intervals(data, n_bootstrap=100)
        
        np.random.seed(42)
        ci2 = calculate_confidence_intervals(data, n_bootstrap=100)
        
        # Should be identical
        assert ci1[0] == ci2[0]
        assert ci1[1] == ci2[1]
    
    def test_gc_content_calculation(self):
        """Test GC content calculation accuracy."""
        test_cases = [
            ('AAAA', 0.0),
            ('GGGG', 1.0),
            ('ACGT', 0.5),
            ('AACCGGTT', 0.5),
        ]
        
        for seq, expected_gc in test_cases:
            gc = (seq.count('G') + seq.count('C')) / len(seq)
            assert abs(gc - expected_gc) < 1e-10


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
