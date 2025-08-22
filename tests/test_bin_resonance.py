#!/usr/bin/env python3
"""
Test Suite for Bin-Resonance Test Module

Tests the bin-resonance CRISPR efficiency analysis implementation
for mathematical correctness and expected behavior.
"""

import sys
import os
import math

# Add project root to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import from the new reference implementation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin"))

from bin_resonance_test import (
    MAP,
    phase_coherence,
    gc_fraction,
    pearson_r,
    bootstrap_ci_r,
    permutation_pvalue,
    bh_fdr,
    read_data,
    quartile_edges,
    calculate_gc_midpoints,
    calculate_lift_at_k,
    bootstrap_lift_ci,
    fisher_z_transform,
    calculate_heterogeneity_stats,
    run
)


def test_complex_base_mapping():
    """Test DNA base to complex number mapping."""
    print("Testing complex base mapping...")
    
    # Test expected mappings
    assert MAP['A'] == 1 + 0j
    assert MAP['T'] == -1 + 0j
    assert MAP['C'] == 0 + 1j
    assert MAP['G'] == 0 - 1j
    
    # Test that unknown bases default to 0j (not in MAP)
    assert MAP.get('N', 0j) == 0j
    
    print("✓ Complex base mapping tests passed")


def test_phase_coherence_calculation():
    """Test phase-coherence calculation."""
    print("Testing phase-coherence calculation...")
    
    # Test empty sequence
    import math
    assert math.isnan(phase_coherence(""))
    
    # Test simple sequences
    pc_a = phase_coherence("A")
    assert pc_a == 1.0  # Single base should give maximum coherence
    
    # Test two-base sequences
    pc_at = phase_coherence("AT")
    pc_gc = phase_coherence("GC")
    
    # Both should be non-zero (and actually equal for 2-base sequences)
    assert pc_at > 0.0
    assert pc_gc > 0.0
    
    # Test that different longer sequences give different values
    pc_atcg = phase_coherence("ATCG")
    pc_gcta = phase_coherence("GCTA")
    assert pc_atcg > 0.0
    assert pc_gcta > 0.0
    
    # Test that longer sequences give reasonable values
    pc_long = phase_coherence("ATCGATCGATCGATCG")
    assert 0.0 <= pc_long <= 1.0
    
    print("✓ Phase-coherence calculation tests passed")


def test_gc_content_calculation():
    """Test GC content calculation."""
    print("Testing GC content calculation...")
    
    import math
    # Test empty sequence
    assert math.isnan(gc_fraction(""))
    
    # Test pure sequences
    assert gc_fraction("AAAA") == 0.0
    assert gc_fraction("TTTT") == 0.0
    assert gc_fraction("CCCC") == 1.0
    assert gc_fraction("GGGG") == 1.0
    
    # Test mixed sequences
    assert gc_fraction("ATCG") == 0.5
    assert gc_fraction("ATCGATCG") == 0.5
    
    # Test case insensitivity
    assert gc_fraction("atcg") == 0.5
    
    print("✓ GC content calculation tests passed")


def test_pearson_correlation():
    """Test Pearson correlation calculation."""
    print("Testing Pearson correlation...")
    
    # Test perfect positive correlation
    x1 = [1, 2, 3, 4, 5]
    y1 = [2, 4, 6, 8, 10]
    r1 = pearson_r(x1, y1)
    assert abs(r1 - 1.0) < 0.001
    
    # Test perfect negative correlation
    x2 = [1, 2, 3, 4, 5]
    y2 = [10, 8, 6, 4, 2]
    r2 = pearson_r(x2, y2)
    assert abs(r2 + 1.0) < 0.001
    
    # Test no correlation - constant values
    x3 = [1, 2, 3, 4, 5]
    y3 = [1, 1, 1, 1, 1]  # Constant
    r3 = pearson_r(x3, y3)
    import math
    assert math.isnan(r3)  # Division by zero in denominator
    
    print("✓ Pearson correlation tests passed")


def test_bootstrap_correlation():
    """Test bootstrap correlation confidence intervals."""
    print("Testing bootstrap correlation...")
    
    # Test with data that has some noise
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [2.1, 3.9, 6.2, 7.8, 10.1, 11.9, 14.1, 15.8, 18.2, 19.9]  # Slightly noisy
    
    ci_low, ci_high = bootstrap_ci_r(x, y, n_boot=100, seed=42)
    r = pearson_r(x, y)
    
    # Should get strong positive correlation
    assert r > 0.8
    
    # Confidence interval should be reasonable
    assert ci_low <= r <= ci_high
    assert ci_high >= ci_low  # CI should be valid
    
    print("✓ Bootstrap correlation tests passed")


def test_permutation_pvalue():
    """Test permutation p-value calculation."""
    print("Testing permutation p-value...")
    
    # Test with data that has some correlation
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [2.1, 3.9, 6.2, 7.8, 10.1, 11.9, 14.1, 15.8, 18.2, 19.9]  # Slightly noisy
    
    r = pearson_r(x, y)
    p_two = permutation_pvalue(x, y, r, n_perm=100, seed=42, tail="two")
    p_greater = permutation_pvalue(x, y, r, n_perm=100, seed=42, tail="greater")
    
    # P-values should be reasonable (between 0 and 1)
    assert 0.0 <= p_two <= 1.0
    assert 0.0 <= p_greater <= 1.0
    
    # For strong positive correlation, greater should be smaller than two-tailed
    if r > 0:
        assert p_greater <= p_two
    
    print("✓ Permutation p-value tests passed")


def test_benjamini_hochberg():
    """Test Benjamini-Hochberg FDR correction."""
    print("Testing Benjamini-Hochberg correction...")
    
    # Test with known p-values
    p_values = [0.01, 0.02, 0.05, 0.20]
    significant, cutoff = bh_fdr(p_values, alpha=0.05)
    
    # Should get significance decisions and cutoff
    assert len(significant) == len(p_values)
    assert isinstance(cutoff, float)
    
    # Test empty list
    empty_sig, empty_cutoff = bh_fdr([])
    assert empty_sig == []
    import math
    assert math.isnan(empty_cutoff)
    
    print("✓ Benjamini-Hochberg correction tests passed")


def test_doench_data_loading():
    """Test loading of Doench data."""
    print("Testing Doench data loading...")
    
    # Test loading actual data file if it exists
    data = read_data("doench_2016.csv")
    
    if data:
        # Should have loaded data
        assert len(data) > 0
        
        # Check first item structure - should be tuple (sequence, efficiency)
        first_item = data[0]
        assert isinstance(first_item, tuple)
        assert len(first_item) == 2
        sequence, efficiency = first_item
        
        # Check data types and ranges
        assert isinstance(sequence, str)
        assert isinstance(efficiency, float)
        assert 0.0 <= efficiency <= 1.0
        
        print(f"✓ Loaded {len(data)} sequences from Doench data")
    else:
        print("⚠ Doench data file not found, skipping data loading test")


def test_deterministic_output():
    """Test that the analysis produces deterministic output with fixed seeds."""
    print("Testing deterministic output...")
    
    # Create a temporary test file
    import tempfile
    import csv
    
    # Create small test dataset
    test_data = [
        ("ATCGATCGATCGATCG", 0.7),
        ("GCGCGCGCGCGCGCGC", 0.3),
        ("AAATTTCCCGGGAAAA", 0.5),
        ("TTTTGGGGCCCCAAAA", 0.6),
        ("ATCGATCGATCGATCG", 0.8)
    ]
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['sequence', 'efficiency'])
        writer.writerows(test_data)
        temp_file = f.name
    
    try:
        # Run twice with same seed
        import io
        from contextlib import redirect_stdout
        
        # Capture output from first run
        f1 = io.StringIO()
        with redirect_stdout(f1):
            run(temp_file, n_boot=100, n_perm=100, seed=42, tail="two", save_control=False)
        output1 = f1.getvalue()
        
        # Capture output from second run
        f2 = io.StringIO()
        with redirect_stdout(f2):
            run(temp_file, n_boot=100, n_perm=100, seed=42, tail="two", save_control=False)
        output2 = f2.getvalue()
        
        # Outputs should be identical
        assert output1 == output2, "Outputs should be deterministic with same seed"
        
    finally:
        # Clean up
        import os
        os.unlink(temp_file)
    
    print("✓ Deterministic output tests passed")


def test_csv_output_with_quartile_edges():
    """Test CSV output includes quartile edges for auditability."""
    print("Testing CSV output with quartile edges...")
    
    import tempfile
    import csv
    
    # Create test dataset with varying GC content
    test_data = [
        ("AAAAAAAAAAAAAAAA", 0.1),  # Low GC
        ("ATATATATATATATAT", 0.2),  # Low GC
        ("ATCGATCGATCGATCG", 0.5),  # Medium GC
        ("GCGCGCGCGCGCGCGC", 0.8),  # High GC
        ("GGGGGGGGGGGGGGGG", 0.9),  # High GC
        ("CCCCCCCCCCCCCCCC", 0.7),  # High GC
        ("TTTTGGGGCCCCAAAA", 0.6),  # Medium GC
        ("AAATTTCCCGGGAAAA", 0.3),  # Medium GC
    ]
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['sequence', 'efficiency'])
        writer.writerows(test_data)
        temp_file = f.name
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        output_file = f.name
    
    try:
        # Run analysis with CSV output
        run(temp_file, n_boot=50, n_perm=50, seed=42, tail="two", out=output_file, save_control=False)
        
        # Read and verify CSV output
        with open(output_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
        
        # Check that new columns are present
        assert "q1_edge" in header
        assert "q2_edge" in header
        assert "q3_edge" in header
        
        # Check that data has the right number of columns
        expected_cols = len(header)
        for row in rows:
            assert len(row) == expected_cols
        
        # Check that quartile edges are consistent across rows
        q1_values = set(row[header.index("q1_edge")] for row in rows)
        q2_values = set(row[header.index("q2_edge")] for row in rows)
        q3_values = set(row[header.index("q3_edge")] for row in rows)
        
        # All rows should have the same quartile edges
        assert len(q1_values) == 1
        assert len(q2_values) == 1  
        assert len(q3_values) == 1
        
    finally:
        # Clean up
        import os
        os.unlink(temp_file)
        os.unlink(output_file)
    
    print("✓ CSV output with quartile edges tests passed")


def test_negative_control_generation():
    """Test negative control generation with --save_control flag."""
    print("Testing negative control generation...")
    
    import tempfile
    import csv
    import os
    
    # Create test dataset
    test_data = [
        ("AAAAAAAAAAAAAAAA", 0.1),
        ("ATATATATATATATAT", 0.2),
        ("ATCGATCGATCGATCG", 0.5),
        ("GCGCGCGCGCGCGCGC", 0.8),
        ("GGGGGGGGGGGGGGGG", 0.9),
        ("CCCCCCCCCCCCCCCC", 0.7),
        ("TTTTGGGGCCCCAAAA", 0.6),
        ("AAATTTCCCGGGAAAA", 0.3),
    ]
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['sequence', 'efficiency'])
        writer.writerows(test_data)
        temp_file = f.name
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        output_file = f.name
    
    try:
        # Run analysis with control generation
        run(temp_file, n_boot=50, n_perm=50, seed=42, tail="two", out=output_file, save_control=True)
        
        # Check that both files were created
        assert os.path.exists(output_file), "Primary output file should exist"
        
        control_file = output_file.replace('.csv', '_control.csv')
        assert os.path.exists(control_file), "Control output file should exist"
        
        # Read both files and verify they have the same structure
        with open(output_file, 'r') as f:
            primary_reader = csv.reader(f)
            primary_header = next(primary_reader)
            primary_rows = list(primary_reader)
        
        with open(control_file, 'r') as f:
            control_reader = csv.reader(f)
            control_header = next(control_reader)
            control_rows = list(control_reader)
        
        # Headers should be identical
        assert primary_header == control_header
        
        # Should have same number of rows (same quartile structure)
        assert len(primary_rows) == len(control_rows)
        
        # Quartile edges should be the same (same input sequences)
        q1_idx = primary_header.index("q1_edge")
        q2_idx = primary_header.index("q2_edge") 
        q3_idx = primary_header.index("q3_edge")
        
        for p_row, c_row in zip(primary_rows, control_rows):
            assert p_row[q1_idx] == c_row[q1_idx]
            assert p_row[q2_idx] == c_row[q2_idx]
            assert p_row[q3_idx] == c_row[q3_idx]
        
        # Clean up control file
        os.unlink(control_file)
        
    finally:
        # Clean up
        os.unlink(temp_file)
        os.unlink(output_file)
    
    print("✓ Negative control generation tests passed")


def test_lift_calculation():
    """Test lift@k% calculation."""
    print("Testing lift@k% calculation...")
    
    import statistics
    
    # Test case: high PC should correlate with high efficiency for positive lift
    x = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]  # PC descending
    y = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]  # Eff matching
    
    lift10 = calculate_lift_at_k(x, y, 10.0)  # Top 10% = 1 guide
    lift20 = calculate_lift_at_k(x, y, 20.0)  # Top 20% = 2 guides
    
    # Top 10% mean = 0.9, overall mean = 0.5, lift = 0.9/0.5 - 1 = 0.8
    expected_lift10 = (0.9 / statistics.mean(y)) - 1
    assert abs(lift10 - expected_lift10) < 0.01
    
    # Top 20% mean = (0.9 + 0.8)/2 = 0.85, lift = 0.85/0.5 - 1 = 0.7  
    expected_lift20 = (statistics.mean([0.9, 0.8]) / statistics.mean(y)) - 1
    assert abs(lift20 - expected_lift20) < 0.01
    
    # Test empty case
    lift_empty = calculate_lift_at_k([], [], 10.0)
    assert math.isnan(lift_empty)
    
    print("✓ Lift@k% calculation tests passed")


def test_fisher_z_transform():
    """Test Fisher z-transformation."""
    print("Testing Fisher z-transformation...")
    
    
    # Test known values
    z_zero = fisher_z_transform(0.0)
    assert abs(z_zero) < 1e-10
    
    # Fisher z of 0.5 should be about 0.549
    z_half = fisher_z_transform(0.5)
    expected = 0.5 * math.log(1.5 / 0.5)  # 0.5 * ln(3) ≈ 0.549
    assert abs(z_half - expected) < 0.01
    
    # Test edge cases
    z_nan = fisher_z_transform(float("nan"))
    assert math.isnan(z_nan)
    
    z_one = fisher_z_transform(1.0)
    assert math.isnan(z_one)  # Should be NaN at |r| = 1
    
    print("✓ Fisher z-transformation tests passed")


def test_heterogeneity_stats():
    """Test heterogeneity statistics calculation."""
    print("Testing heterogeneity statistics...")
    
    # Test with no heterogeneity (similar correlations)
    correlations = [0.5, 0.52, 0.48, 0.51]
    sample_sizes = [50, 60, 55, 45]
    
    stats_result = calculate_heterogeneity_stats(correlations, sample_sizes)
    
    # Should have low Q and I² for homogeneous correlations
    assert not math.isnan(stats_result["Q"])
    assert not math.isnan(stats_result["I2"])
    assert stats_result["df"] == 3
    
    # Test with high heterogeneity
    het_correlations = [0.8, -0.2, 0.1, 0.7]
    het_stats = calculate_heterogeneity_stats(het_correlations, sample_sizes)
    
    # Should have higher Q for heterogeneous correlations
    assert het_stats["Q"] > stats_result["Q"]
    
    # Test edge cases
    empty_stats = calculate_heterogeneity_stats([], [])
    assert math.isnan(empty_stats["Q"])
    
    single_stats = calculate_heterogeneity_stats([0.5], [50])
    assert math.isnan(single_stats["Q"])  # Need at least 2 for heterogeneity test
    
    print("✓ Heterogeneity statistics tests passed")


def test_gc_midpoints():
    """Test GC midpoint calculation."""
    print("Testing GC midpoint calculation...")
    
    # Test with known quartile edges
    gc_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    q1, q2, q3 = 0.25, 0.5, 0.75
    
    midpoints = calculate_gc_midpoints(q1, q2, q3, gc_values)
    
    # Q1: values ≤ 0.25 -> [0.1, 0.2] -> midpoint 0.15
    # Q2: values 0.25 < x ≤ 0.5 -> [0.3, 0.4, 0.5] -> midpoint 0.4
    # Q3: values 0.5 < x ≤ 0.75 -> [0.6, 0.7] -> midpoint 0.65
    # Q4: values > 0.75 -> [0.8, 0.9, 1.0] -> midpoint 0.9
    
    assert len(midpoints) == 4
    assert abs(midpoints[0] - 0.15) < 0.01  # Q1
    assert abs(midpoints[1] - 0.4) < 0.01   # Q2  
    assert abs(midpoints[2] - 0.65) < 0.01  # Q3
    assert abs(midpoints[3] - 0.9) < 0.01   # Q4
    
    print("✓ GC midpoint calculation tests passed")


def test_enhanced_csv_output():
    """Test CSV output includes all new enhancement metrics."""
    print("Testing enhanced CSV output...")
    
    import tempfile
    import csv
    
    # Create test dataset with sufficient data
    test_data = []
    for i in range(20):
        # Create sequences with varying GC content and some correlation structure
        gc = i / 19.0  # GC from 0 to 1
        if gc < 0.3:
            seq_base = "AAATTTAAATTTAAAT"
            eff = 0.2 + 0.3 * (i % 3)  # Some variation
        elif gc < 0.6:
            seq_base = "ATCGATCGATCGATCG"
            eff = 0.4 + 0.2 * (i % 4)
        else:
            seq_base = "GGCCGGCCGGCCGGCC"
            eff = 0.6 + 0.3 * (i % 3)
        
        test_data.append((seq_base, eff))
    
    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.writer(f)
        writer.writerow(['sequence', 'efficiency'])
        writer.writerows(test_data)
        temp_file = f.name
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        output_file = f.name
    
    try:
        # Run analysis with CSV output
        run(temp_file, n_boot=20, n_perm=20, seed=42, tail="two", out=output_file, save_control=False)
        
        # Read and verify CSV output
        with open(output_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
        
        # Check that new columns are present
        required_cols = [
            "lift_top10", "lift10_ci_low", "lift10_ci_high",
            "lift_top20", "lift20_ci_low", "lift20_ci_high", 
            "gc_midpoint", "Fisher_Q", "I2", "p_heterogeneity"
        ]
        
        for col in required_cols:
            assert col in header, f"Missing column: {col}"
        
        # Check that data has the right number of columns
        expected_cols = len(header)
        for row in rows:
            assert len(row) == expected_cols
        
        # Check that heterogeneity stats are consistent across rows
        q_values = set(row[header.index("Fisher_Q")] for row in rows)
        i2_values = set(row[header.index("I2")] for row in rows)
        
        # All rows should have the same heterogeneity statistics
        assert len(q_values) == 1
        assert len(i2_values) == 1
        
    finally:
        # Clean up
        import os
        os.unlink(temp_file)
        os.unlink(output_file)
    
    print("✓ Enhanced CSV output tests passed")


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
        test_permutation_pvalue()
        test_benjamini_hochberg()
        test_doench_data_loading()
        test_deterministic_output()
        test_csv_output_with_quartile_edges()
        test_negative_control_generation()
        
        # New Z-style enhancement tests
        test_lift_calculation()
        test_fisher_z_transform()
        test_heterogeneity_stats()
        test_gc_midpoints()
        test_enhanced_csv_output()
        
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