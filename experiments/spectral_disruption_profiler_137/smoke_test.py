#!/usr/bin/env python3
"""
Smoke Test for Spectral Disruption Profiler Falsification Experiment

Quick validation test designed to complete in <5 seconds for CI.
Tests core functionality without full statistical analysis.

Author: Z Framework Falsification Team
Date: 2025-11-17
"""

import sys
import time
from pathlib import Path

# Add parent directories for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.spectral_disruption_profiler_137.spectral_disruption_profiler import (
    SpectralDisruptionFalsifier
)


def test_dna_validation():
    """Test DNA sequence validation (human DNA only, A/C/G/T/N)."""
    print("Testing DNA validation...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42)
    
    # Valid sequences
    valid_sequences = [
        "ATCG",
        "GATTACA",
        "ACGTACGTACGT",
        "NNNNATCGNNNN"
    ]
    
    for seq in valid_sequences:
        try:
            validated = falsifier.validate_dna_sequence(seq)
            assert validated == seq.upper()
        except ValueError as e:
            print(f"  ❌ FAILED: Valid sequence rejected: {seq}")
            print(f"     Error: {e}")
            return False
    
    # Invalid sequences (should raise ValueError)
    invalid_sequences = [
        ("ATCGU", "Contains U (RNA)"),
        ("ATCGR", "Contains IUPAC code R"),
        ("ATCG123", "Contains numbers"),
        ("", "Empty sequence")
    ]
    
    for seq, reason in invalid_sequences:
        try:
            falsifier.validate_dna_sequence(seq)
            print(f"  ❌ FAILED: Invalid sequence accepted: {seq} ({reason})")
            return False
        except ValueError:
            pass  # Expected
    
    print("  ✓ DNA validation passed")
    return True


def test_complex_encoding():
    """Test DNA to complex waveform encoding."""
    print("Testing complex encoding...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42)
    
    # Test standard mapping
    seq = "ATCG"
    waveform = falsifier.encode_dna_complex(seq)
    
    expected = [1+0j, -1+0j, 0+1j, 0-1j]
    
    if len(waveform) != len(expected):
        print(f"  ❌ FAILED: Wrong length {len(waveform)} != {len(expected)}")
        return False
    
    for i, (got, exp) in enumerate(zip(waveform, expected)):
        if abs(got - exp) > 1e-10:
            print(f"  ❌ FAILED: Position {i}: {got} != {exp}")
            return False
    
    print("  ✓ Complex encoding passed")
    return True


def test_theta_prime():
    """Test geometric resolution function θ′(n,k)."""
    print("Testing θ′(n,k) geometric resolution...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42, k_parameter=0.3)
    
    # Test specific values
    theta_0 = falsifier.theta_prime(0, 0.3)
    theta_21 = falsifier.theta_prime(21, 0.3)
    
    # θ′(0,k) should be 0 since (0 mod φ)/φ = 0
    if abs(theta_0) > 1e-10:
        print(f"  ❌ FAILED: θ′(0,0.3) = {theta_0}, expected ≈0")
        return False
    
    # θ′(n,k) should be bounded
    for n in [1, 5, 10, 21, 42]:
        theta = falsifier.theta_prime(n, 0.3)
        if not (0 <= theta <= 10):  # Reasonable upper bound
            print(f"  ❌ FAILED: θ′({n},0.3) = {theta} out of bounds")
            return False
    
    print("  ✓ θ′(n,k) passed")
    return True


def test_fft_features():
    """Test FFT feature computation."""
    print("Testing FFT feature computation...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42)
    
    # Test with simple sequence
    seq = "ATCGATCGATCGATCGATCG"
    
    # Compute features with and without phase weighting
    features_weighted = falsifier.compute_fft_features(seq, use_phase_weighting=True)
    features_unweighted = falsifier.compute_fft_features(seq, use_phase_weighting=False)
    
    # Check required keys
    required_keys = [
        'dominant_freq',
        'spectral_entropy',
        'sidelobe_count',
        'gc_content',
        'sequence_length',
        'phase_weighted'
    ]
    
    for key in required_keys:
        if key not in features_weighted:
            print(f"  ❌ FAILED: Missing key '{key}' in weighted features")
            return False
        if key not in features_unweighted:
            print(f"  ❌ FAILED: Missing key '{key}' in unweighted features")
            return False
    
    # Check phase_weighted flag
    if not features_weighted['phase_weighted']:
        print("  ❌ FAILED: phase_weighted flag wrong for weighted features")
        return False
    if features_unweighted['phase_weighted']:
        print("  ❌ FAILED: phase_weighted flag wrong for unweighted features")
        return False
    
    # Check GC content calculation
    expected_gc = 0.5  # 10 G/C out of 20 bases
    if abs(features_weighted['gc_content'] - expected_gc) > 0.01:
        print(f"  ❌ FAILED: GC content {features_weighted['gc_content']} != {expected_gc}")
        return False
    
    print("  ✓ FFT features passed")
    return True


def test_synthetic_generation():
    """Test synthetic sequence generation."""
    print("Testing synthetic sequence generation...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42)
    
    # Generate sequences
    n_seq = 10
    length = 21
    gc_range = (0.4, 0.6)
    
    sequences = falsifier.generate_synthetic_sequences(
        n_sequences=n_seq,
        length=length,
        gc_range=gc_range,
        add_noise=False
    )
    
    # Check count
    if len(sequences) != n_seq:
        print(f"  ❌ FAILED: Generated {len(sequences)} sequences, expected {n_seq}")
        return False
    
    # Check each sequence
    for i, seq in enumerate(sequences):
        # Check length
        if len(seq) != length:
            print(f"  ❌ FAILED: Sequence {i} length {len(seq)} != {length}")
            return False
        
        # Check valid bases
        if not all(base in 'ATCG' for base in seq):
            print(f"  ❌ FAILED: Sequence {i} contains invalid bases: {seq}")
            return False
        
        # Check GC content roughly in range
        gc = (seq.count('G') + seq.count('C')) / len(seq)
        if not (gc_range[0] - 0.1 <= gc <= gc_range[1] + 0.1):
            print(f"  ❌ FAILED: Sequence {i} GC {gc:.2f} out of range {gc_range}")
            return False
    
    print("  ✓ Synthetic generation passed")
    return True


def test_z_framework_integration():
    """Test Z Framework calculator integration."""
    print("Testing Z Framework integration...")
    
    falsifier = SpectralDisruptionFalsifier(seed=42)
    
    # Check Z calculator is initialized
    if falsifier.z_calc is None:
        print("  ❌ FAILED: Z Framework calculator not initialized")
        return False
    
    # Test a simple Z calculation (if available)
    try:
        # This is just checking the calculator exists and has expected methods
        if not hasattr(falsifier.z_calc, 'calculate_z_value'):
            print("  ⚠ WARNING: Z calculator missing calculate_z_value method")
        print("  ✓ Z Framework integration passed")
        return True
    except Exception as e:
        print(f"  ❌ FAILED: Z Framework error: {e}")
        return False


def test_performance():
    """Test performance meets <5s target for smoke test."""
    print("Testing performance...")
    
    start_time = time.time()
    
    falsifier = SpectralDisruptionFalsifier(
        seed=42,
        n_bootstrap=10,  # Minimal for smoke test
        n_permutation=10
    )
    
    # Generate small dataset
    sequences = falsifier.generate_synthetic_sequences(n_sequences=10, length=21)
    
    # Process sequences
    for seq in sequences:
        _ = falsifier.compute_fft_features(seq, use_phase_weighting=True)
        _ = falsifier.compute_fft_features(seq, use_phase_weighting=False)
    
    elapsed = time.time() - start_time
    
    if elapsed > 5.0:
        print(f"  ❌ FAILED: Smoke test took {elapsed:.2f}s > 5s target")
        return False
    
    print(f"  ✓ Performance passed ({elapsed:.2f}s < 5s)")
    return True


def run_smoke_tests():
    """Run all smoke tests."""
    print("="*60)
    print("SPECTRAL DISRUPTION PROFILER - SMOKE TEST")
    print("="*60)
    print()
    
    tests = [
        ("DNA Validation", test_dna_validation),
        ("Complex Encoding", test_complex_encoding),
        ("θ′(n,k) Function", test_theta_prime),
        ("FFT Features", test_fft_features),
        ("Synthetic Generation", test_synthetic_generation),
        ("Z Framework Integration", test_z_framework_integration),
        ("Performance", test_performance)
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ❌ FAILED: Unexpected error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
        print()
    
    print("="*60)
    print(f"SMOKE TEST RESULTS: {passed}/{len(tests)} passed")
    print("="*60)
    
    if failed > 0:
        print(f"❌ {failed} test(s) failed")
        return 1
    else:
        print("✓ All tests passed")
        return 0


if __name__ == '__main__':
    sys.exit(run_smoke_tests())
