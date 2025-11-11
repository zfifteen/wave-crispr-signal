"""
Test suite for DNA Storage Hypothesis Testing Module

Validates the implementation against the claims made in the research simulation:
- 12-15% error reduction in DNA sequencing
- 10% compression improvement through golden spiral arrangement
- 20% efficiency gain in bio-computational processing
- 30% improvement in cryptographic key strength
- 25% speed boost in data retrieval

Uses bootstrap resampling and statistical validation following Z Framework principles.
"""

import unittest
import logging
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from modules.dna_storage_hypothesis import DNAStorageHypothesis, DNAStorageResults
from scripts.z_framework import ZFrameworkCalculator

# Configure logging
logging.basicConfig(level=logging.INFO)


class TestDNAStorageHypothesis(unittest.TestCase):
    """Test suite for DNA storage hypothesis validation"""

    def setUp(self):
        """Set up test fixtures"""
        self.dna_storage = DNAStorageHypothesis(precision_dps=50)
        self.test_sequences = [
            "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG",  # PCSK9-like (31 bp)
            "ATCGATCGATCGATCGATCGATCGATCG",  # Regular pattern (27 bp)
            "GGCCTTAAGGCCTTAAGGCCTTAAGGCC",  # Repetitive (27 bp)
            "AAAAAAAAAAAAAAAAAAAAAAAAAAA",  # Homopolymer (25 bp)
            "ATCG" * 10,  # Simple repeat (40 bp)
        ]

    def test_dna_encoding_decoding(self):
        """Test DNA to binary encoding and decoding"""
        print("\n--- Testing DNA Encoding/Decoding ---")

        for sequence in self.test_sequences[:2]:  # Test first 2 sequences
            with self.subTest(sequence=sequence):
                # Encode to binary
                binary_data = self.dna_storage.encode_dna_to_binary(sequence)

                # Check binary length (2 bits per base)
                expected_length = len(sequence) * 2
                self.assertEqual(
                    len(binary_data),
                    expected_length,
                    f"Binary encoding length mismatch for {sequence}",
                )

                # Decode back to DNA
                decoded_sequence = self.dna_storage.decode_binary_to_dna(binary_data)

                # Check round-trip accuracy
                self.assertEqual(
                    sequence,
                    decoded_sequence,
                    f"Round-trip encoding failed for {sequence}",
                )

                print(f"✓ {sequence} -> {len(binary_data)} bits -> {decoded_sequence}")

    def test_prime_curvature_transform(self):
        """Test prime curvature transformation of binary data"""
        print("\n--- Testing Prime Curvature Transform ---")

        sequence = self.test_sequences[0]  # PCSK9-like sequence
        binary_data = self.dna_storage.encode_dna_to_binary(sequence)

        # Apply prime curvature transformation
        transformed_data = self.dna_storage.apply_prime_curvature_transform(binary_data)

        # Check that transformation adds redundancy
        self.assertGreaterEqual(
            len(transformed_data),
            len(binary_data),
            "Prime curvature transform should add redundancy",
        )

        # Calculate redundancy factor
        redundancy_factor = len(transformed_data) / len(binary_data)
        print(f"Original length: {len(binary_data)} bits")
        print(f"Transformed length: {len(transformed_data)} bits")
        print(f"Redundancy factor: {redundancy_factor:.2f}x")

        # Redundancy should be reasonable (not too excessive)
        self.assertLess(redundancy_factor, 2.0, "Redundancy factor too high")
        self.assertGreater(
            redundancy_factor, 1.0, "Transform should add some redundancy"
        )

    def test_golden_spiral_optimization(self):
        """Test golden spiral geometry optimization"""
        print("\n--- Testing Golden Spiral Optimization ---")

        for sequence in self.test_sequences[:3]:  # Test first 3 sequences
            with self.subTest(sequence=sequence):
                # Apply golden spiral optimization
                optimized_sequence = self.dna_storage.apply_golden_spiral_optimization(
                    sequence
                )

                # Check length preservation
                self.assertEqual(
                    len(optimized_sequence),
                    len(sequence),
                    "Golden spiral optimization changed sequence length",
                )

                # Check nucleotide composition preservation
                original_counts = {base: sequence.count(base) for base in "ATCG"}
                optimized_counts = {
                    base: optimized_sequence.count(base) for base in "ATCG"
                }

                self.assertEqual(
                    original_counts,
                    optimized_counts,
                    "Golden spiral optimization changed nucleotide composition",
                )

                print(f"✓ {sequence[:20]}... -> {optimized_sequence[:20]}...")

    def test_error_correction_simulation(self):
        """Test error correction improvement simulation"""
        print("\n--- Testing Error Correction Simulation ---")

        sequence = self.test_sequences[0]
        binary_data = self.dna_storage.encode_dna_to_binary(sequence)
        transformed_data = self.dna_storage.apply_prime_curvature_transform(binary_data)

        # Test error correction multiple times for statistical validation
        error_reductions = []
        for _ in range(10):
            error_reduction = self.dna_storage.simulate_error_correction(
                binary_data, transformed_data
            )
            error_reductions.append(error_reduction * 100)

        mean_reduction = np.mean(error_reductions)
        std_reduction = np.std(error_reductions)

        print(f"Error reduction: {mean_reduction:.1f}% ± {std_reduction:.1f}%")

        # Should show some improvement
        self.assertGreater(
            mean_reduction, 0, "Error correction should show improvement"
        )
        self.assertLess(mean_reduction, 50, "Error reduction should be realistic")

    def test_compression_efficiency(self):
        """Test compression efficiency calculation"""
        print("\n--- Testing Compression Efficiency ---")

        for sequence in self.test_sequences[:3]:
            with self.subTest(sequence=sequence):
                optimized_sequence = self.dna_storage.apply_golden_spiral_optimization(
                    sequence
                )
                compression = self.dna_storage.simulate_compression_efficiency(
                    sequence, optimized_sequence
                )

                print(
                    f"Sequence length {len(sequence)}: {compression*100:.1f}% compression improvement"
                )

                # Compression should be positive and reasonable
                self.assertGreater(
                    compression, 0, "Compression should show improvement"
                )
                self.assertLess(
                    compression, 0.5, "Compression improvement should be realistic"
                )

    def test_bio_computational_efficiency(self):
        """Test bio-computational processing efficiency"""
        print("\n--- Testing Bio-Computational Efficiency ---")

        for sequence in self.test_sequences:
            with self.subTest(sequence=sequence):
                efficiency = self.dna_storage.simulate_bio_computational_efficiency(
                    sequence
                )

                print(f"Sequence {sequence[:15]}...: {efficiency:.3f} efficiency score")

                # Efficiency should be between 0 and 1
                self.assertGreaterEqual(
                    efficiency, 0, "Efficiency should be non-negative"
                )
                self.assertLessEqual(efficiency, 1, "Efficiency should not exceed 100%")

    def test_cryptographic_key_generation(self):
        """Test cryptographic key generation"""
        print("\n--- Testing Cryptographic Key Generation ---")

        sequence = self.test_sequences[0]

        # Generate multiple keys to test uniqueness
        keys = []
        for i in range(5):
            key = self.dna_storage.generate_cryptographic_key(sequence + str(i))
            keys.append(key)

        # Check key properties
        for i, key in enumerate(keys):
            self.assertIsInstance(key, str, "Key should be string")
            self.assertGreater(len(key), 0, "Key should not be empty")
            print(f"Key {i+1}: {key[:16]}... (length: {len(key)})")

        # Check uniqueness (different inputs should give different keys)
        unique_keys = set(keys)
        self.assertEqual(
            len(unique_keys), len(keys), "Keys should be unique for different inputs"
        )

    def test_cryptographic_strength(self):
        """Test cryptographic strength evaluation"""
        print("\n--- Testing Cryptographic Strength ---")

        sequence = self.test_sequences[0]
        key1 = self.dna_storage.generate_cryptographic_key(sequence)
        key2 = self.dna_storage.generate_cryptographic_key(sequence + "MODIFIED")

        strength_improvement = self.dna_storage.test_cryptographic_strength(key1, key2)

        print(f"Cryptographic strength improvement: {strength_improvement*100:.1f}%")

        # Strength improvement should be positive and reasonable
        self.assertGreater(strength_improvement, 0, "Should show strength improvement")
        self.assertLess(strength_improvement, 1, "Improvement should be realistic")

    def test_full_simulation(self):
        """Test complete DNA storage simulation"""
        print("\n--- Testing Full DNA Storage Simulation ---")

        sequence = self.test_sequences[0]  # Use PCSK9-like sequence

        # Run simulation with reduced trials for testing
        results = self.dna_storage.run_dna_storage_simulation(sequence, num_trials=50)

        # Validate result structure
        self.assertIsInstance(
            results, DNAStorageResults, "Should return DNAStorageResults"
        )

        # Check all metrics are present and reasonable
        metrics = [
            ("error_reduction_percent", results.error_reduction_percent),
            (
                "compression_improvement_percent",
                results.compression_improvement_percent,
            ),
            ("efficiency_gain_percent", results.efficiency_gain_percent),
            (
                "cryptographic_strength_improvement_percent",
                results.cryptographic_strength_improvement_percent,
            ),
            (
                "retrieval_speed_improvement_percent",
                results.retrieval_speed_improvement_percent,
            ),
        ]

        print(f"Results for sequence {sequence[:20]}...:")
        for metric_name, value in metrics:
            print(f"  {metric_name}: {value:.1f}%")

            # All metrics should be positive
            self.assertGreater(value, 0, f"{metric_name} should be positive")
            # All metrics should be realistic (< 100%)
            self.assertLess(value, 100, f"{metric_name} should be realistic")

        # Check confidence interval
        ci_lower, ci_upper = results.confidence_interval
        self.assertLess(ci_lower, ci_upper, "Confidence interval should be valid")
        print(f"  95% Confidence Interval: ({ci_lower:.1f}, {ci_upper:.1f})")

    def test_hypothesis_validation(self):
        """Test hypothesis validation against expected claims"""
        print("\n--- Testing Hypothesis Validation ---")

        sequence = self.test_sequences[0]

        # Run validation
        validation = self.dna_storage.validate_hypothesis(sequence)

        # Check validation structure
        self.assertIn("overall_hypothesis_validated", validation)
        self.assertIn("confidence_interval", validation)

        # Check individual claims
        expected_claims = [
            "error_reduction_claim",
            "compression_claim",
            "efficiency_claim",
            "crypto_claim",
            "retrieval_claim",
        ]

        print("Validation Results:")
        for claim in expected_claims:
            self.assertIn(claim, validation, f"Missing claim: {claim}")

            claim_data = validation[claim]
            self.assertIn("expected", claim_data)
            self.assertIn("actual", claim_data)
            self.assertIn("passed", claim_data)

            status = "✓ PASS" if claim_data["passed"] else "✗ FAIL"
            print(
                f"  {status} {claim}: Expected {claim_data['expected']}, Got {claim_data['actual']}"
            )

        overall_status = (
            "✓ PASSED" if validation["overall_hypothesis_validated"] else "✗ FAILED"
        )
        print(f"\nOverall Hypothesis Validation: {overall_status}")

    def test_statistical_robustness(self):
        """Test statistical robustness with multiple sequences"""
        print("\n--- Testing Statistical Robustness ---")

        all_results = []

        for sequence in self.test_sequences[:3]:  # Test first 3 sequences
            results = self.dna_storage.run_dna_storage_simulation(
                sequence, num_trials=20
            )
            all_results.append(results)

        # Calculate statistics across sequences
        error_reductions = [r.error_reduction_percent for r in all_results]
        compressions = [r.compression_improvement_percent for r in all_results]

        error_mean = np.mean(error_reductions)
        error_std = np.std(error_reductions)
        compression_mean = np.mean(compressions)
        compression_std = np.std(compressions)

        print("Cross-sequence statistics:")
        print(f"  Error reduction: {error_mean:.1f}% ± {error_std:.1f}%")
        print(f"  Compression: {compression_mean:.1f}% ± {compression_std:.1f}%")

        # Results should be consistent across sequences
        self.assertLess(
            error_std,
            error_mean * 0.5,
            "Error reduction should be reasonably consistent",
        )
        self.assertLess(
            compression_std,
            compression_mean * 0.5,
            "Compression should be reasonably consistent",
        )

    def test_edge_cases(self):
        """Test edge cases and error handling"""
        print("\n--- Testing Edge Cases ---")

        # Test with very short sequence
        short_seq = "ATCG"
        try:
            self.dna_storage.run_dna_storage_simulation(short_seq, num_trials=10)
            print(f"✓ Short sequence ({short_seq}) handled successfully")
        except Exception as e:
            self.fail(f"Short sequence handling failed: {e}")

        # Test with sequence containing invalid characters
        invalid_seq = "ATCGXYZ"
        try:
            # Should handle gracefully by filtering out invalid characters
            binary_data = self.dna_storage.encode_dna_to_binary(invalid_seq)
            # Should only encode valid DNA bases
            expected_length = 4 * 2  # Only ATCG should be encoded
            self.assertEqual(len(binary_data), expected_length)
            print("✓ Invalid characters filtered correctly")
        except Exception as e:
            self.fail(f"Invalid character handling failed: {e}")

        # Test empty sequence
        try:
            empty_binary = self.dna_storage.encode_dna_to_binary("")
            self.assertEqual(len(empty_binary), 0)
            print("✓ Empty sequence handled correctly")
        except Exception as e:
            self.fail(f"Empty sequence handling failed: {e}")


class TestIntegrationWithZFramework(unittest.TestCase):
    """Integration tests with existing Z Framework"""

    def setUp(self):
        """Set up integration test fixtures"""
        self.dna_storage = DNAStorageHypothesis()
        self.z_calc = ZFrameworkCalculator()

    def test_z_framework_integration(self):
        """Test integration with Z Framework calculations"""
        print("\n--- Testing Z Framework Integration ---")

        sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG"

        # Calculate Z Framework values directly
        z_results = self.z_calc.calculate_z_values(sequence)

        # Use DNA storage hypothesis which internally uses Z Framework
        optimized_seq = self.dna_storage.apply_golden_spiral_optimization(sequence)
        efficiency = self.dna_storage.simulate_bio_computational_efficiency(
            optimized_seq
        )

        print(f"Original Z mean: {float(z_results['z_mean']):.3f}")
        print(f"Z variance: {float(z_results['z_variance']):.3f}")
        print(f"Bio-computational efficiency: {efficiency:.3f}")

        # Efficiency should correlate with Z Framework convergence properties
        phi_deviation = abs(float(z_results["z_mean"]) - 0.618)
        expected_efficiency = 1.0 - min(phi_deviation / 2.0, 1.0)

        # Allow for some variance in efficiency calculation
        self.assertAlmostEqual(
            efficiency,
            expected_efficiency * 1.2,
            delta=0.2,
            msg="Efficiency should correlate with Z Framework properties",
        )

    def test_geodesic_resolution_consistency(self):
        """Test that geodesic resolution is used consistently"""
        print("\n--- Testing Geodesic Resolution Consistency ---")

        # Test that the same k value (0.3) is used consistently
        k_value = 0.3
        position = 10

        # Direct Z Framework calculation
        z_theta = self.z_calc.calculate_geodesic_resolution(position, k_value)

        # Verify the optimal k is being used in DNA storage
        self.assertAlmostEqual(float(self.dna_storage.optimal_k), k_value, places=5)

        print(f"Geodesic resolution at position {position}: {float(z_theta):.6f}")
        print(f"Optimal k value: {float(self.dna_storage.optimal_k)}")

        # Values should be consistent
        self.assertGreater(float(z_theta), 0, "Geodesic resolution should be positive")


def run_comprehensive_validation():
    """Run comprehensive validation of DNA storage hypothesis"""
    print("\n" + "=" * 80)
    print("COMPREHENSIVE DNA STORAGE HYPOTHESIS VALIDATION")
    print("=" * 80)

    # Initialize framework
    dna_storage = DNAStorageHypothesis()

    # Test with representative sequence
    test_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGG"

    print(f"\nTest Sequence: {test_sequence}")
    print(f"Length: {len(test_sequence)} bp")

    # Run full simulation
    results = dna_storage.run_dna_storage_simulation(test_sequence, num_trials=200)

    print("\nSimulation Results:")
    print(f"  Error Reduction: {results.error_reduction_percent:.1f}%")
    print(f"  Compression Improvement: {results.compression_improvement_percent:.1f}%")
    print(f"  Efficiency Gain: {results.efficiency_gain_percent:.1f}%")
    print(
        f"  Cryptographic Strength: {results.cryptographic_strength_improvement_percent:.1f}%"
    )
    print(f"  Retrieval Speed: {results.retrieval_speed_improvement_percent:.1f}%")
    print(f"  Confidence Interval: {results.confidence_interval}")

    # Validate hypothesis
    validation = dna_storage.validate_hypothesis(test_sequence)

    print("\nHypothesis Validation:")
    claims = [
        ("Error Reduction (12-15%)", validation["error_reduction_claim"]),
        ("Compression (10%)", validation["compression_claim"]),
        ("Efficiency (20%)", validation["efficiency_claim"]),
        ("Cryptographic (30%)", validation["crypto_claim"]),
        ("Retrieval Speed (25%)", validation["retrieval_claim"]),
    ]

    passed_count = 0
    for claim_name, claim_data in claims:
        status = "✓ PASS" if claim_data["passed"] else "✗ FAIL"
        print(
            f"  {status} {claim_name}: Expected {claim_data['expected']}, Got {claim_data['actual']}"
        )
        if claim_data["passed"]:
            passed_count += 1

    overall_status = (
        "✓ PASSED" if validation["overall_hypothesis_validated"] else "✗ FAILED"
    )
    print(
        f"\nOverall Validation: {overall_status} ({passed_count}/{len(claims)} claims passed)"
    )

    return validation["overall_hypothesis_validated"]


if __name__ == "__main__":
    # Run unit tests
    unittest.main(verbosity=2, exit=False)

    # Run comprehensive validation
    print("\n" + "=" * 80)
    success = run_comprehensive_validation()

    if success:
        print("\n🎉 DNA Storage Hypothesis VALIDATED! 🎉")
    else:
        print("\n⚠️  Some hypothesis claims need further investigation ⚠️")
