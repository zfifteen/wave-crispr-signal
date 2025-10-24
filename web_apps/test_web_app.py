#!/usr/bin/env python3
"""
Test script for DNA Breathing Dynamics Web Application

This script tests the core functionality of the web app without starting the Flask server.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# Import the app components
from web_apps.breathing_dynamics_demo import (
    DNAValidator,
    BreathingDynamicsEncoder,
    ArbitraryEncoder,
    BreathingDynamicsAnalyzer,
    BREATHING_FREQ
)

def test_validator():
    """Test DNA sequence validation"""
    print("Testing DNA Validator...")
    
    # Valid sequence
    is_valid, msg = DNAValidator.validate_sequence("ATGCGATCGATCG")
    assert is_valid, f"Valid sequence failed: {msg}"
    print("✓ Valid DNA sequence accepted")
    
    # Invalid bases
    is_valid, msg = DNAValidator.validate_sequence("ATGCXYZ")
    assert not is_valid, "Invalid bases should be rejected"
    print("✓ Invalid bases rejected")
    
    # Too short
    is_valid, msg = DNAValidator.validate_sequence("ATGC")
    assert not is_valid, "Too short sequence should be rejected"
    print("✓ Too short sequence rejected")
    
    # With U (RNA)
    is_valid, msg = DNAValidator.validate_sequence("AUGC")
    assert not is_valid, "RNA (U) should be rejected"
    print("✓ RNA (U) rejected")
    
    print()

def test_breathing_encoder():
    """Test breathing dynamics encoder"""
    print("Testing Breathing Dynamics Encoder...")
    
    encoder = BreathingDynamicsEncoder()
    
    # Check weights exist
    assert 'A' in encoder.weights
    assert 'T' in encoder.weights
    assert 'C' in encoder.weights
    assert 'G' in encoder.weights
    print("✓ Weights initialized for all bases")
    
    # Check frequency separation
    a_weight = encoder.weights['A']
    c_weight = encoder.weights['C']
    
    # AT should have lower real part (lower frequency)
    # GC should have higher real part (higher frequency)
    assert a_weight.real < c_weight.real, "AT should have lower frequency encoding than GC"
    print(f"✓ Frequency separation: AT={a_weight.real:.2f}, GC={c_weight.real:.2f}")
    
    # Test encoding
    sequence = "ATGC"
    encoded = encoder.encode_sequence(sequence)
    assert len(encoded) == len(sequence)
    print(f"✓ Encoded sequence length correct: {len(encoded)}")
    
    # Get weights info
    info = encoder.get_weights_info()
    assert 'A' in info
    assert info['A']['frequency_hz'] == BREATHING_FREQ['A']
    print("✓ Weights info correct")
    
    print()

def test_arbitrary_encoder():
    """Test arbitrary encoder"""
    print("Testing Arbitrary Encoder...")
    
    encoder = ArbitraryEncoder(seed=42)
    
    # Check weights exist
    assert 'A' in encoder.weights
    print("✓ Arbitrary weights initialized")
    
    # Test encoding
    sequence = "ATGC"
    encoded = encoder.encode_sequence(sequence)
    assert len(encoded) == len(sequence)
    print(f"✓ Encoded sequence length correct: {len(encoded)}")
    
    # Test reproducibility with seed
    encoder2 = ArbitraryEncoder(seed=42)
    assert encoder.weights['A'] == encoder2.weights['A']
    print("✓ Seed reproducibility works")
    
    print()

def test_analyzer():
    """Test sequence analyzer"""
    print("Testing Breathing Dynamics Analyzer...")
    
    analyzer = BreathingDynamicsAnalyzer()
    
    # Test with a real human DNA sequence
    sequence = "ATGCGATCGATCGATCGATCGATCGATCGATCG"
    
    results = analyzer.analyze_sequence(sequence)
    
    # Check results structure
    assert 'sequence_length' in results
    assert 'gc_content' in results
    assert 'base_composition' in results
    assert 'breathing_spectrum' in results
    assert 'arbitrary_spectrum' in results
    assert 'mutation_analysis' in results
    print("✓ Results structure correct")
    
    # Check sequence info
    assert results['sequence_length'] == len(sequence)
    print(f"✓ Sequence length: {results['sequence_length']}")
    
    # Check GC content calculation
    gc_count = sequence.count('G') + sequence.count('C')
    expected_gc = (gc_count / len(sequence)) * 100
    assert abs(results['gc_content'] - expected_gc) < 0.1
    print(f"✓ GC content: {results['gc_content']:.1f}%")
    
    # Check mutation analysis
    mut = results['mutation_analysis']
    assert 'gc_affecting' in mut
    assert 'at_affecting' in mut
    assert 'random_mutation' in mut
    print("✓ Mutation analysis present")
    
    # Check that GC-affecting mutation has data
    if mut['gc_affecting']:
        gc_mut = mut['gc_affecting']
        assert 'z_breathing' in gc_mut
        assert 'z_arbitrary' in gc_mut
        assert 'winner' in gc_mut
        print(f"✓ GC-affecting mutation: {gc_mut['description']}")
        print(f"  Breathing Z: {gc_mut['z_breathing']:.4f}")
        print(f"  Arbitrary Z: {gc_mut['z_arbitrary']:.4f}")
        print(f"  Winner: {gc_mut['winner']}")
    
    print()

def main():
    """Run all tests"""
    print("="*70)
    print("DNA Breathing Dynamics Web App - Test Suite")
    print("="*70)
    print()
    
    try:
        test_validator()
        test_breathing_encoder()
        test_arbitrary_encoder()
        test_analyzer()
        
        print("="*70)
        print("✅ ALL TESTS PASSED")
        print("="*70)
        return 0
    
    except Exception as e:
        print("\n" + "="*70)
        print("❌ TEST FAILED")
        print("="*70)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
