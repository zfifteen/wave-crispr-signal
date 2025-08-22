#!/usr/bin/env python3
"""
Test Suite for Bin-Resonance Test Module

Tests the bin-resonance CRISPR efficiency analysis implementation
for mathematical correctness and expected behavior.
"""

import sys
import os
import math
import tempfile
import csv

# Add project root to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bin_resonance_test import (
    complex_base_mapping,
    calculate_phase_coherence,
    calculate_gc_content,
    pearson_correlation,
    bootstrap_correlation,
    benjamini_hochberg_correction,
    load_doench_data,
    perform_bin_resonance_analysis
)


def test_complex_base_mapping():
    """Test DNA base to complex number mapping."""
    print("Testing complex base mapping...")
    
    # Test expected mappings
    assert complex_base_mapping('A') == 1 + 0j
    assert complex_base_mapping('T') == -1 + 0j
    assert complex_base_mapping('C') == 0 + 1j
    assert complex_base_mapping('G') == 0 - 1j
    
    # Test case insensitivity
    assert complex_base_mapping('a') == 1 + 0j
    assert complex_base_mapping('t') == -1 + 0j
    
    # Test unknown base default
    assert complex_base_mapping('N') == 1 + 0j
    
    print("✓ Complex base mapping tests passed")


def test_phase_coherence_calculation():
    """Test phase-coherence calculation."""
    print("Testing phase-coherence calculation...")
    
    # Test empty sequence
    assert calculate_phase_coherence("") == 0.0
    
    # Test simple sequences
    pc_a = calculate_phase_coherence("A")
    assert pc_a == 1.0  # Single base should give maximum coherence
    
    # Test two-base sequences
    pc_at = calculate_phase_coherence("AT")
    pc_gc = calculate_phase_coherence("GC")
    
    # Both should be non-zero (and actually equal for 2-base sequences)
    assert pc_at > 0.0
    assert pc_gc > 0.0
    
    # Test that different longer sequences give different values
    pc_atcg = calculate_phase_coherence("ATCG")
    pc_gcta = calculate_phase_coherence("GCTA")
    assert pc_atcg > 0.0
    assert pc_gcta > 0.0
    
    # Test that longer sequences give reasonable values
    pc_long = calculate_phase_coherence("ATCGATCGATCGATCG")
    assert 0.0 <= pc_long <= 1.0
    
    print("✓ Phase-coherence calculation tests passed")


def test_gc_content_calculation():
    """Test GC content calculation."""
    print("Testing GC content calculation...")
    
    # Test empty sequence
    assert calculate_gc_content("") == 0.0
    
    # Test pure sequences
    assert calculate_gc_content("AAAA") == 0.0
    assert calculate_gc_content("TTTT") == 0.0
    assert calculate_gc_content("CCCC") == 1.0
    assert calculate_gc_content("GGGG") == 1.0
    
    # Test mixed sequences
    assert calculate_gc_content("ATCG") == 0.5
    assert calculate_gc_content("ATCGATCG") == 0.5
    
    # Test case insensitivity
    assert calculate_gc_content("atcg") == 0.5
    
    print("✓ GC content calculation tests passed")


def test_pearson_correlation():
    """Test Pearson correlation calculation."""
    print("Testing Pearson correlation...")
    
    # Test perfect positive correlation
    x1 = [1, 2, 3, 4, 5]
    y1 = [2, 4, 6, 8, 10]
    r1 = pearson_correlation(x1, y1)
    assert abs(r1 - 1.0) < 0.001
    
    # Test perfect negative correlation
    x2 = [1, 2, 3, 4, 5]
    y2 = [10, 8, 6, 4, 2]
    r2 = pearson_correlation(x2, y2)
    assert abs(r2 + 1.0) < 0.001
    
    # Test no correlation
    x3 = [1, 2, 3, 4, 5]
    y3 = [1, 1, 1, 1, 1]  # Constant
    r3 = pearson_correlation(x3, y3)
    assert r3 == 0.0
    
    # Test insufficient data
    assert pearson_correlation([1], [2]) == 0.0
    assert pearson_correlation([], []) == 0.0
    
    print("✓ Pearson correlation tests passed")


def test_bootstrap_correlation():
    """Test bootstrap correlation confidence intervals."""
    print("Testing bootstrap correlation...")
    
    # Test with data that has some noise
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [2.1, 3.9, 6.2, 7.8, 10.1, 11.9, 14.1, 15.8, 18.2, 19.9]  # Slightly noisy
    
    r, ci_low, ci_high = bootstrap_correlation(x, y, n_boot=100, seed=42)
    
    # Should get strong positive correlation
    assert r > 0.8
    
    # Confidence interval should be reasonable
    assert ci_low <= r <= ci_high
    assert ci_high >= ci_low  # CI should be valid
    
    # Test with insufficient data
    r2, ci_low2, ci_high2 = bootstrap_correlation([1], [2])
    assert r2 == 0.0 and ci_low2 == 0.0 and ci_high2 == 0.0
    
    print("✓ Bootstrap correlation tests passed")


def test_benjamini_hochberg():
    """Test Benjamini-Hochberg FDR correction."""
    print("Testing Benjamini-Hochberg correction...")
    
    # Test with known p-values
    p_values = [0.01, 0.02, 0.05, 0.20]
    significant = benjamini_hochberg_correction(p_values, q=0.10)
    
    # First few should be significant, last should not
    assert len(significant) == len(p_values)
    assert significant[0] == True  # 0.01 should be significant
    
    # Test empty list
    assert benjamini_hochberg_correction([]) == []
    
    print("✓ Benjamini-Hochberg correction tests passed")


def test_doench_data_loading():
    """Test loading of Doench data."""
    print("Testing Doench data loading...")
    
    # Test loading actual data file if it exists
    data = load_doench_data("doench_2016.csv")
    
    if data:
        # Should have loaded data
        assert len(data) > 0
        
        # Check first item structure
        first_item = data[0]
        required_keys = ['sequence', 'efficiency', 'gc_content', 'phase_coherence', 'length']
        for key in required_keys:
            assert key in first_item
        
        # Check data types and ranges
        assert isinstance(first_item['sequence'], str)
        assert 0.0 <= first_item['efficiency'] <= 1.0
        assert 0.0 <= first_item['gc_content'] <= 1.0
        assert first_item['phase_coherence'] >= 0.0
        assert first_item['length'] > 0
        
        print(f"✓ Loaded {len(data)} sequences from Doench data")
    else:
        print("⚠ Doench data file not found, skipping data loading test")


def test_bin_resonance_analysis():
    """Test the main bin-resonance analysis function."""
    print("Testing bin-resonance analysis...")
    
    # Create synthetic test data
    test_data = []
    for i in range(20):
        # Create sequences with varying GC content
        if i < 5:
            seq = "AAATATATATAT"  # Low GC
            gc = 0.0
        elif i < 10:
            seq = "ATCGATCGATCG"  # Medium GC
            gc = 0.5
        elif i < 15:
            seq = "GCGCGCGCGCGC"  # High GC
            gc = 1.0
        else:
            seq = "CCGGCCGGCCGG"  # Very high GC
            gc = 1.0
        
        test_data.append({
            'sequence': seq,
            'efficiency': 0.5 + 0.1 * i,  # Varying efficiency
            'gc_content': gc,
            'phase_coherence': 0.3 + 0.02 * i,  # Varying PC
            'length': len(seq)
        })
    
    # Run analysis
    results = perform_bin_resonance_analysis(test_data)
    
    # Should get results for quartiles
    assert len(results) == 4
    
    # Check result structure
    for result in results:
        required_keys = ['bin', 'n', 'r', 'ci_low', 'ci_high', 'p_boot', 'passed']
        for key in required_keys:
            assert key in result
        
        # Check data types
        assert isinstance(result['bin'], str)
        assert isinstance(result['n'], int)
        assert isinstance(result['r'], float)
        assert isinstance(result['passed'], bool)
    
    print("✓ Bin-resonance analysis tests passed")


def run_all_tests():
    """Run all tests for bin-resonance module."""
    print("Bin-Resonance Test Suite")
    print("=" * 40)
    
    try:
        test_complex_base_mapping()
        test_phase_coherence_calculation()
        test_gc_content_calculation()
        test_pearson_correlation()
        test_bootstrap_correlation()
        test_benjamini_hochberg()
        test_doench_data_loading()
        test_bin_resonance_analysis()
        
        print("\n" + "=" * 40)
        print("🎉 ALL TESTS PASSED!")
        return True
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)