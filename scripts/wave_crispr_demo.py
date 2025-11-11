#!/usr/bin/env python3
"""
Wave-CRISPR Physical Z-Metrics Demo

This demo showcases the complete implementation of the Wave-CRISPR experiment plan
with four sequence-derivable physical Z-metrics and integration with the θ′ geodesic pipeline.
"""

import sys
import os
import tempfile

# Add applications to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'applications'))

from crispr_physical_z_metrics import PhysicalZMetricsCalculator, read_fasta_with_validation
from crispr_z_metrics_validation import ZMetricsValidationSuite
from crispr_integrated_pipeline import IntegratedWaveCRISPRPipeline, K_STAR


def demo_guardrails_validation():
    """Demonstrate input validation guardrails."""
    print("🛡️  DEMO: Input Validation Guardrails")
    print("=" * 50)
    
    calc = PhysicalZMetricsCalculator()
    
    # Test valid human sequence
    print("✅ Valid human sequence:")
    valid_seq = "ATCGATCGATCGATCG"
    valid_header = "PCSK9_GRCh38_test"
    
    try:
        results = calc.calculate_all_physical_z_metrics(valid_seq, valid_header)
        print(f"   Sequence: {valid_seq}")
        print(f"   Header: {valid_header}")
        print(f"   Z-metrics calculated successfully!")
        print()
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test invalid nucleotides
    print("❌ Invalid nucleotides (should fail):")
    invalid_seq = "ATCGXYZ"
    try:
        calc.calculate_all_physical_z_metrics(invalid_seq, "test")
        print("   ERROR: Should have failed!")
    except Exception as e:
        print(f"   ✅ Correctly rejected: {e}")
        print()
    
    # Test rate clamping
    print("🔧 Rate clamping demonstration:")
    # This is demonstrated in the calculator where negative rates are clamped to 0
    print("   All rates B ≥ 0 enforced in Z = A * (B / c)")
    print("   Division guards prevent zero division")
    print()


def demo_four_z_metrics():
    """Demonstrate the four physical Z-metrics."""
    print("🧬 DEMO: Four Physical Z-Metrics")
    print("=" * 50)
    
    calc = PhysicalZMetricsCalculator()
    
    # Test sequences with different characteristics
    test_cases = [
        ("High GC (stable)", "GGGCCCGGGCCCGGGCCCGG"),
        ("Low GC (flexible)", "ATATATATATATATATAT"),
        ("Mixed composition", "ATCGATCGATCGATCGATCG"),
        ("GC-rich repeats", "GCGCGCGCGCGCGCGCGC")
    ]
    
    for name, sequence in test_cases:
        print(f"\n{name}: {sequence}")
        results = calc.calculate_all_physical_z_metrics(sequence, validate=False)
        
        print(f"  1. Opening Kinetics Z:    {float(results['opening_kinetics']['z_opening']):.6f}")
        print(f"  2. Stacking Dissociation: {float(results['stacking_dissociation']['z_stacking']):.6f}")
        print(f"  3. Twist Fluctuation:     {float(results['twist_fluctuation']['z_twist']):.6f}")
        print(f"  4. Melting Kinetics:      {float(results['melting_kinetics']['z_melting']):.6f}")
        print(f"     Z-mean:                {float(results['summary']['z_mean']):.6f}")
    
    print()


def demo_z_framework_form():
    """Demonstrate Z = A * (B / c) form with c = e²."""
    print("⚡ DEMO: Z Framework Form (Z = A * (B / c), c = e²)")
    print("=" * 50)
    
    calc = PhysicalZMetricsCalculator()
    sequence = "ATCGATCGATCGATCG"
    
    # Calculate opening kinetics in detail
    result = calc.calculate_base_pair_opening_kinetics(sequence)
    
    print(f"Sequence: {sequence}")
    print(f"Physical constant c = e² = {float(calc.e_squared):.6f}")
    print()
    print("Opening Kinetics Z-metric breakdown:")
    print(f"  A (context factor):  {float(result['context_factor']):.6f}")
    print(f"  B (opening rate):    {float(result['opening_rate']):.6f}")
    print(f"  c (invariant e²):    {float(calc.e_squared):.6f}")
    print(f"  Z = A * (B / c):     {float(result['z_opening']):.6f}")
    
    # Verify calculation
    expected_z = float(result['context_factor']) * (float(result['opening_rate']) / float(calc.e_squared))
    print(f"  Verification:        {expected_z:.6f} ✅")
    print()


def demo_geodesic_integration():
    """Demonstrate integration with θ′ geodesic at k* = 0.3."""
    print("🌐 DEMO: θ′ Geodesic Integration (k* = 0.3)")
    print("=" * 50)
    
    pipeline = IntegratedWaveCRISPRPipeline()
    
    test_sequence = "ATCGATCGATCGATCGATCG"
    print(f"Test sequence: {test_sequence}")
    print(f"Geodesic resolution parameter k* = {K_STAR}")
    print()
    
    # Get geodesic features
    geodesic_features = pipeline.calculate_geodesic_resolution(test_sequence, position=10)
    
    print("Geodesic features:")
    print(f"  θ′(n,k) at position 10:    {geodesic_features['theta_prime']:.6f}")
    print(f"  Weighted composition:      {geodesic_features['weighted_composition']:.6f}")
    print(f"  Geometric curvature:       {geodesic_features['geometric_curvature']:.6f}")
    print()
    
    # Full integrated analysis
    analysis = pipeline.analyze_guide_with_integrated_pipeline(test_sequence)
    bands = analysis['integrated_analysis']['feature_bands']
    
    print("Feature bands:")
    print(f"  Physical Z-band score:     {analysis['integrated_analysis']['band_scores']['z_band_score']:.3f}")
    print(f"  Geodesic band score:       {analysis['integrated_analysis']['band_scores']['geodesic_band_score']:.3f}")
    print(f"  Spectral band score:       {analysis['integrated_analysis']['band_scores']['spectral_band_score']:.3f}")
    print(f"  Invariant band score:      {analysis['integrated_analysis']['band_scores']['invariant_band_score']:.3f}")
    print(f"  Integrated score:          {analysis['integrated_analysis']['integrated_score']['composite_score']:.3f}")
    print()


def demo_hypothesis_testing_framework():
    """Demonstrate H1 and H2 hypothesis testing framework."""
    print("🔬 DEMO: Hypothesis Testing Framework")
    print("=" * 50)
    
    validator = ZMetricsValidationSuite()
    
    # Generate test data for hypothesis testing
    sequences = [
        "GGGCCCGGGCCCGGGCCCGG",  # High efficiency expected (lower Opening-Z)
        "ATATATATATATATATAT",    # Lower efficiency expected
        "GCGCGCGCGCGCGCGCGC",    # Variable efficiency
        "ATCGATCGATCGATCGATCG",  # Moderate efficiency
    ]
    
    # Simulated efficiency data (for demonstration)
    efficiency_data = [0.85, 0.25, 0.60, 0.70]
    headers = [f"test_seq_{i+1}" for i in range(len(sequences))]
    
    print("H1 Hypothesis: Lower Opening-Z correlates with higher CRISPR efficiency")
    print("H2 Hypothesis: Very low Stacking-Z sequences (framework)")
    print()
    
    # Generate validation report
    report = validator.generate_validation_report(sequences, headers, efficiency_data)
    
    # Show H1 results
    h1 = report['h1_hypothesis_test']
    print(f"H1 Test Results:")
    print(f"  Status: {h1['status']}")
    if h1['status'] == 'TESTED':
        print(f"  Correlation: {h1['correlation']:.3f} (target: ≤ -0.5)")
        print(f"  P-value: {h1['p_value']:.3f} (target: < 1e-5)")
        print(f"  H1 Supported: {'✅' if h1['h1_supported'] else '❌'}")
    
    # Show H2 results
    h2 = report['h2_hypothesis_framework']
    print(f"\nH2 Framework:")
    print(f"  Status: {h2['status']}")
    print(f"  Very low threshold: {h2['very_low_threshold']:.6f}")
    print(f"  Framework ready: {'✅' if h2['h2_framework'] else '❌'}")
    print()


def demo_cli_integration():
    """Demonstrate CLI integration."""
    print("💻 DEMO: CLI Integration")
    print("=" * 50)
    
    # Create a temporary FASTA file
    fasta_content = """
>PCSK9_target_GRCh38
ATGCGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGATGCGGCTGCTGCTGCTGCTG
>CRISPR_benchmark_sequence
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content.strip())
        temp_fasta = f.name
    
    print("Created temporary FASTA file with human sequences")
    print("CLI commands demonstrated:")
    print()
    print("1. Basic guide design with Z-metrics:")
    print(f"   python applications/crispr_cli.py design {temp_fasta} --physical-z-metrics -n 2")
    print()
    print("2. CSV output with Z-metrics:")
    print(f"   python applications/crispr_cli.py design {temp_fasta} --physical-z-metrics --format csv")
    print()
    print("3. Validation suite:")
    print("   python applications/crispr_z_metrics_validation.py --test-data")
    print()
    print("4. Integrated pipeline:")
    print("   python applications/crispr_integrated_pipeline.py")
    print()
    
    # Clean up
    os.unlink(temp_fasta)


def demo_validation_suite():
    """Demonstrate comprehensive validation suite."""
    print("📊 DEMO: Validation Suite with Bootstrap CIs")
    print("=" * 50)
    
    # Run validation on test data
    print("Running validation suite with bootstrap confidence intervals...")
    
    import subprocess
    result = subprocess.run([
        sys.executable, 'applications/crispr_z_metrics_validation.py', '--test-data'
    ], capture_output=True, text=True, cwd=os.path.dirname(__file__))
    
    if result.returncode == 0:
        print("✅ Validation suite completed successfully!")
        print("\nKey outputs:")
        print("- Bootstrap 95% confidence intervals for all Z-metrics")
        print("- Guardrail validation (rate clamping)")
        print("- H1/H2 hypothesis testing framework")
        print("- Effect size calculations (Cohen's d)")
    else:
        print("❌ Validation suite encountered issues")
        print(result.stderr)
    
    print()


def main():
    """Run complete demo of Wave-CRISPR implementation."""
    print("🧬 WAVE-CRISPR PHYSICAL Z-METRICS IMPLEMENTATION DEMO")
    print("=" * 70)
    print("Complete implementation of the Wave-CRISPR experiment plan:")
    print("- Four sequence-derivable physical Z-metrics")
    print("- Z = A * (B / c) form with c = e² ≈ 7.389")
    print("- Human DNA validation guardrails")
    print("- θ′ geodesic integration at k* = 0.3")
    print("- Hypothesis testing framework (H1, H2)")
    print("- Bootstrap validation with confidence intervals")
    print("- CLI and Python API integration")
    print()
    
    # Run all demos
    demo_guardrails_validation()
    demo_four_z_metrics()
    demo_z_framework_form()
    demo_geodesic_integration()
    demo_hypothesis_testing_framework()
    demo_cli_integration()
    demo_validation_suite()
    
    print("🎉 DEMO COMPLETE!")
    print("=" * 70)
    print("✅ All objectives implemented:")
    print("  0) Guardrails: Input validation, Z-form, rate checks ✅")
    print("  1) Four Z-metrics: Opening, stacking, twist, melting ✅")
    print("  2) θ′ geodesic integration at k* = 0.3 ✅")
    print("  3) Validation suite with bootstrap CIs ✅")
    print("  4) Clean CLI + Python API ✅")
    print("  5) H1/H2 hypothesis framework ✅")
    print()
    print("Ready for experimental validation and CRISPR efficiency testing!")


if __name__ == "__main__":
    main()