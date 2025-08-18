#!/usr/bin/env python3
"""
Test script for Pain Management Z Framework Application

This script tests the pain management application functionality to ensure
proper integration with the Z Framework and validation of all key features.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from pain_management_application import (
    PainManagementAnalyzer, 
    PainManagementTarget,
    demo_pain_management_analysis,
    format_pain_analysis_results
)
from z_framework import format_mpmath_for_display
import mpmath as mp

def test_pain_management_analyzer():
    """Test the PainManagementAnalyzer initialization and basic functionality"""
    print("\n--- Testing PainManagementAnalyzer Initialization ---")
    
    # Test initialization
    analyzer = PainManagementAnalyzer(precision_dps=30)  # Lower precision for faster testing
    
    print(f"Analyzer precision: {analyzer.precision}")
    print(f"Number of Casgevy targets: {len(analyzer.casgevy_targets)}")
    print(f"Number of JOURNAVX targets: {len(analyzer.journavx_targets)}")
    
    # Verify target initialization
    assert len(analyzer.casgevy_targets) > 0, "Should have Casgevy targets"
    assert len(analyzer.journavx_targets) > 0, "Should have JOURNAVX targets"
    
    print("✓ PainManagementAnalyzer initialization successful")
    return analyzer

def test_prime_curvature_analysis(analyzer):
    """Test prime curvature analysis functionality"""
    print("\n--- Testing Prime Curvature Analysis ---")
    
    # Test with a Casgevy target
    casgevy_target = analyzer.casgevy_targets[0]
    print(f"Testing with Casgevy target: {casgevy_target.name}")
    print(f"Sequence length: {len(casgevy_target.sequence)} bp")
    
    results = analyzer.analyze_prime_curvature(casgevy_target)
    
    # Verify results structure
    required_keys = [
        'target_name', 'target_type', 'sequence_length', 'z_mean', 'z_variance',
        'prime_curvature', 'pain_efficacy_score', 'therapeutic_index',
        'phi_convergence', 'variance_convergence', 'curvature_disruption',
        'binding_prediction'
    ]
    
    for key in required_keys:
        assert key in results, f"Missing key: {key}"
    
    print(f"Target: {results['target_name']}")
    print(f"Pain efficacy score: {format_mpmath_for_display(results['pain_efficacy_score'])}")
    print(f"Therapeutic index: {format_mpmath_for_display(results['therapeutic_index'])}")
    print(f"Prime curvature: {format_mpmath_for_display(results['prime_curvature'])}")
    print(f"Binding prediction: {results['binding_prediction']}")
    
    # Test with a JOURNAVX target
    journavx_target = analyzer.journavx_targets[0]
    print(f"\nTesting with JOURNAVX target: {journavx_target.name}")
    
    journavx_results = analyzer.analyze_prime_curvature(journavx_target)
    print(f"Pain efficacy score: {format_mpmath_for_display(journavx_results['pain_efficacy_score'])}")
    print(f"Therapeutic index: {format_mpmath_for_display(journavx_results['therapeutic_index'])}")
    
    print("✓ Prime curvature analysis successful")
    return results

def test_z5d_predictor(analyzer):
    """Test Z5D predictor functionality"""
    print("\n--- Testing Z5D Predictor ---")
    
    # Test with a shorter sequence first (for speed)
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    target_n = 1000  # Smaller target for testing
    
    print(f"Testing Z5D predictor with sequence length: {len(test_sequence)}")
    print(f"Target N: {target_n}")
    
    results = analyzer.implement_z5d_predictor(test_sequence, target_n=target_n)
    
    # Verify results structure
    required_keys = [
        'sequence_length', 'target_n', 'z5d_dimension_1', 'z5d_dimension_2',
        'z5d_dimension_3', 'z5d_dimension_4', 'z5d_dimension_5',
        'density_boost_achieved', 'density_boost_target', 'density_enhancement_success',
        'confidence_interval_lower', 'confidence_interval_upper', 'statistical_significance'
    ]
    
    for key in required_keys:
        assert key in results, f"Missing key: {key}"
    
    print(f"Extended sequence length: {results['sequence_length']}")
    print(f"Density boost achieved: {format_mpmath_for_display(results['density_boost_achieved'])}x")
    print(f"Density boost target: {results['density_boost_target']}x")
    print(f"Enhancement success: {results['density_enhancement_success']}")
    print(f"Statistical significance: {results['statistical_significance']}")
    
    # Print Z5D dimensions
    print("\nZ5D Dimensions:")
    for i in range(1, 6):
        dim_key = f'z5d_dimension_{i}'
        print(f"  Dimension {i}: {format_mpmath_for_display(results[dim_key])}")
    
    print("✓ Z5D predictor test successful")
    return results

def test_comprehensive_analysis(analyzer):
    """Test comprehensive pain management analysis"""
    print("\n--- Testing Comprehensive Analysis ---")
    
    # Run analysis on a subset of targets for speed
    test_targets = analyzer.casgevy_targets[:2] + analyzer.journavx_targets[:1]
    
    print(f"Running analysis on {len(test_targets)} targets")
    
    results = analyzer.run_comprehensive_pain_analysis(targets=test_targets)
    
    # Verify results structure
    assert 'analysis_summary' in results
    assert 'target_analyses' in results
    assert 'z5d_analyses' in results
    assert 'comparative_metrics' in results
    
    summary = results['analysis_summary']
    print(f"Total targets analyzed: {summary['total_targets']}")
    print(f"Casgevy targets: {summary['casgevy_targets']}")
    print(f"JOURNAVX targets: {summary['journavx_targets']}")
    
    # Verify we have results for each target
    assert len(results['target_analyses']) == len(test_targets)
    assert len(results['z5d_analyses']) == len(test_targets)
    
    # Print comparative metrics
    if results['comparative_metrics']:
        metrics = results['comparative_metrics']
        print(f"Mean efficacy score: {format_mpmath_for_display(metrics['mean_efficacy_score'])}")
        print(f"Mean therapeutic index: {format_mpmath_for_display(metrics['mean_therapeutic_index'])}")
    
    print("✓ Comprehensive analysis test successful")
    return results

def test_target_creation():
    """Test PainManagementTarget creation and validation"""
    print("\n--- Testing Target Creation ---")
    
    # Create a custom target
    custom_target = PainManagementTarget(
        name="test_target",
        sequence="ATCGATCGATCGATCGATCG",
        target_type="casgevy",
        clinical_stage="preclinical",
        binding_affinity=75.0
    )
    
    print(f"Created target: {custom_target.name}")
    print(f"Sequence length: {len(custom_target.sequence)}")
    print(f"Target type: {custom_target.target_type}")
    print(f"Clinical stage: {custom_target.clinical_stage}")
    print(f"Binding affinity: {custom_target.binding_affinity}")
    
    # Test with analyzer
    analyzer = PainManagementAnalyzer(precision_dps=20)
    results = analyzer.analyze_prime_curvature(custom_target)
    
    print(f"Analysis successful for custom target")
    print(f"Pain efficacy score: {format_mpmath_for_display(results['pain_efficacy_score'])}")
    
    print("✓ Target creation test successful")

def test_precision_and_determinism():
    """Test high-precision calculations and deterministic behavior"""
    print("\n--- Testing Precision and Determinism ---")
    
    # Test deterministic behavior
    analyzer1 = PainManagementAnalyzer(precision_dps=30)
    analyzer2 = PainManagementAnalyzer(precision_dps=30)
    
    target = analyzer1.casgevy_targets[0]
    
    results1 = analyzer1.analyze_prime_curvature(target)
    results2 = analyzer2.analyze_prime_curvature(target)
    
    # Compare key metrics
    efficacy_diff = abs(results1['pain_efficacy_score'] - results2['pain_efficacy_score'])
    curvature_diff = abs(results1['prime_curvature'] - results2['prime_curvature'])
    
    print(f"Efficacy score difference: {efficacy_diff}")
    print(f"Prime curvature difference: {curvature_diff}")
    
    # Should be identical for deterministic calculations
    assert efficacy_diff < mp.mpf('1e-10'), "Non-deterministic efficacy calculation"
    assert curvature_diff < mp.mpf('1e-10'), "Non-deterministic curvature calculation"
    
    print("✓ Precision and determinism test successful")

def test_format_results():
    """Test results formatting functionality"""
    print("\n--- Testing Results Formatting ---")
    
    analyzer = PainManagementAnalyzer(precision_dps=20)
    results = analyzer.run_comprehensive_pain_analysis(targets=analyzer.casgevy_targets[:1])
    
    formatted_output = format_pain_analysis_results(results)
    
    # Verify formatted output contains expected sections
    assert "PAIN MANAGEMENT Z FRAMEWORK ANALYSIS RESULTS" in formatted_output
    assert "TARGET ANALYSES:" in formatted_output
    assert "Z5D PREDICTOR ANALYSES:" in formatted_output
    assert "COMPARATIVE METRICS:" in formatted_output
    
    print("Formatted output preview:")
    print("=" * 50)
    print(formatted_output[:500] + "..." if len(formatted_output) > 500 else formatted_output)
    print("=" * 50)
    
    print("✓ Results formatting test successful")

def test_statistical_validation():
    """Test statistical validation features"""
    print("\n--- Testing Statistical Validation ---")
    
    analyzer = PainManagementAnalyzer(precision_dps=30)
    test_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    
    # Test Z5D predictor statistical features
    results = analyzer.implement_z5d_predictor(test_sequence, target_n=500)
    
    # Check confidence intervals
    ci_lower = results['confidence_interval_lower']
    ci_upper = results['confidence_interval_upper']
    boost_ratio = results['density_boost_achieved']
    
    print(f"Density boost: {format_mpmath_for_display(boost_ratio)}")
    print(f"95% CI: [{format_mpmath_for_display(ci_lower)}, {format_mpmath_for_display(ci_upper)}]")
    print(f"Statistical significance: {results['statistical_significance']}")
    
    # Verify CI bounds
    assert ci_lower < boost_ratio < ci_upper, "Boost ratio should be within confidence interval"
    
    print("✓ Statistical validation test successful")

def run_all_tests():
    """Run all test functions"""
    print("=" * 60)
    print("PAIN MANAGEMENT Z FRAMEWORK APPLICATION TESTS")
    print("=" * 60)
    
    try:
        # Basic functionality tests
        analyzer = test_pain_management_analyzer()
        test_prime_curvature_analysis(analyzer)
        test_z5d_predictor(analyzer)
        test_comprehensive_analysis(analyzer)
        
        # Additional tests
        test_target_creation()
        test_precision_and_determinism()
        test_format_results()
        test_statistical_validation()
        
        print("\n" + "=" * 60)
        print("ALL TESTS PASSED SUCCESSFULLY!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

def run_demo():
    """Run the demonstration analysis"""
    print("\n" + "=" * 60)
    print("PAIN MANAGEMENT DEMONSTRATION")
    print("=" * 60)
    
    try:
        results = demo_pain_management_analysis()
        print("\n✓ Demonstration completed successfully")
        return results
    except Exception as e:
        print(f"\n❌ DEMONSTRATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Run tests first
    test_success = run_all_tests()
    
    if test_success:
        # Run demonstration
        demo_results = run_demo()
    else:
        print("Skipping demonstration due to test failures")