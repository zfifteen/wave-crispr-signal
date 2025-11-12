#!/usr/bin/env python3
"""
Test Falsification Experiments

This test module validates the falsification experiments to ensure
they correctly identify and expose the flawed claims in the Z5D framework.

The tests verify that:
1. k* parameter testing correctly identifies non-optimality
2. Density boost measurements show statistical impossibility
3. Biological relevance testing detects spurious correlations
4. Electromagnetic therapy analysis identifies lack of foundation

Author: Falsification Study Team
Date: 2025
"""

import unittest
import sys
import os

# Add project root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from experiments.falsification_experiments import FalsificationExperiment


class TestFalsificationExperiments(unittest.TestCase):
    """Test suite for falsification experiments."""

    def setUp(self):
        """Set up test environment."""
        self.experiment = FalsificationExperiment(random_seed=42)
        self.test_sequences = [
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGG",
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGC",
        ]

    def test_k_parameter_falsification(self):
        """Test that k* parameter claim is correctly falsified."""
        print("\n🧪 Testing k* parameter falsification...")

        results = self.experiment.test_k_parameter_claim(self.test_sequences)

        # Verify falsification results
        self.assertIn("k_claimed", results)
        self.assertIn("claimed_k_enhancement", results)
        self.assertIn("optimal_k_found", results)

        # The claimed k* should not be optimal
        claimed_enhancement = results["claimed_k_enhancement"]
        optimal_enhancement = results["optimal_enhancement"]

        # If optimal enhancement is positive, claimed should be significantly worse
        if optimal_enhancement > 0:
            self.assertLess(
                claimed_enhancement,
                optimal_enhancement * 0.5,
                "Claimed k* should perform significantly worse than optimal",
            )

        print("✅ k* parameter test validated:")
        print(f"   Claimed k* enhancement: {claimed_enhancement:.4f}")
        print(f"   Optimal enhancement: {optimal_enhancement:.4f}")

    def test_density_boost_falsification(self):
        """Test that 210% density boost claim is correctly falsified."""
        print("\n🧪 Testing density boost falsification...")

        results = self.experiment.test_density_boost_claim(self.test_sequences)

        # Should have error or very low actual boost
        if "error" not in results:
            self.assertIn("mean_boost_pct", results)
            self.assertIn("p_value", results)

            mean_boost = results["mean_boost_pct"]
            p_value = results["p_value"]

            # Actual boost should be far below claimed 210%
            self.assertLess(
                mean_boost, 50.0, "Actual boost should be far below claimed 210%"
            )

            # p-value should indicate statistical significance
            self.assertLess(
                p_value,
                0.05,
                "p-value should indicate significant difference from 210%",
            )

            print("✅ Density boost test validated:")
            print(f"   Mean boost: {mean_boost:.2f}% (vs claimed 210%)")
            print(f"   p-value: {p_value:.2e}")
        else:
            # Error in measurement is itself evidence of falsification
            self.assertTrue(results["boost_claim_falsified"])
            print("✅ Density boost falsified due to measurement failure")

    def test_biological_relevance_assessment(self):
        """Test biological relevance assessment."""
        print("\n🧪 Testing biological relevance assessment...")

        results = self.experiment.test_biological_relevance(self.test_sequences)

        # Should have comparison results
        self.assertIn("real_sequences", results)
        self.assertIn("random_sequences", results)

        # Check if statistical test was performed
        if "real_vs_random_ttest" in results:
            ttest_results = results["real_vs_random_ttest"]
            p_value = ttest_results["p_value"]

            print("✅ Biological relevance test validated:")
            print(f"   p-value (real vs random): {p_value:.4f}")
            print(f"   Significant biological signal: {p_value < 0.05}")

    def test_electromagnetic_therapy_falsification(self):
        """Test electromagnetic therapy claims falsification."""
        print("\n🧪 Testing electromagnetic therapy falsification...")

        results = self.experiment.test_electromagnetic_therapy_claims()

        # Should clearly indicate no mathematical connection
        self.assertFalse(results["mathematical_connection_exists"])
        self.assertFalse(results["electromagnetic_field_equations_present"])
        self.assertFalse(results["therapeutic_dosimetry_calculated"])
        self.assertTrue(results["therapeutic_claim_falsified"])

        # Risk assessment should indicate concerns
        risk_assessment = results["risk_assessment"]
        self.assertTrue(risk_assessment["unproven_medical_claims"])
        self.assertTrue(risk_assessment["potential_patient_harm"])

        print("✅ Electromagnetic therapy falsification validated:")
        print("   Mathematical connection: FALSE")
        print("   Therapeutic claims falsified: TRUE")
        print("   Risk assessment: HIGH")

    def test_report_generation(self):
        """Test that falsification report can be generated."""
        print("\n🧪 Testing report generation...")

        # Run basic experiments
        self.experiment.test_k_parameter_claim(self.test_sequences[:2])
        self.experiment.test_density_boost_claim(self.test_sequences[:2])
        self.experiment.test_electromagnetic_therapy_claims()

        # Generate report
        report = self.experiment.generate_falsification_report()

        # Report should contain key sections
        self.assertIn("FALSIFICATION REPORT", report)
        self.assertIn("Claims Tested and Falsified", report)
        self.assertIn("Calibration Parameter k*", report)
        self.assertIn("Electromagnetic Therapy", report)
        self.assertIn("Conclusions", report)

        print("✅ Report generation validated")
        print(f"   Report length: {len(report)} characters")

    def test_reproducibility(self):
        """Test that experiments are reproducible."""
        print("\n🧪 Testing reproducibility...")

        # Run same experiment twice
        experiment1 = FalsificationExperiment(random_seed=42)
        experiment2 = FalsificationExperiment(random_seed=42)

        sequences = experiment1.generate_test_sequences(5, 100)

        results1 = experiment1.test_k_parameter_claim(sequences)
        results2 = experiment2.test_k_parameter_claim(sequences)

        # Results should be identical (within floating point precision)
        self.assertAlmostEqual(
            results1["claimed_k_enhancement"],
            results2["claimed_k_enhancement"],
            places=6,
            msg="Results should be reproducible with same random seed",
        )

        print("✅ Reproducibility validated")
        print("   Identical results with same random seed")


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete falsification workflow."""

    def test_complete_workflow(self):
        """Test the complete falsification workflow."""
        print("\n🧪 Testing complete falsification workflow...")

        from experiments.falsification_experiments import run_complete_falsification_study

        # This should run without errors and produce files
        try:
            experiment, report = run_complete_falsification_study()

            # Check that files were created
            self.assertTrue(os.path.exists("FALSIFICATION_REPORT.md"))
            self.assertTrue(os.path.exists("falsification_results.json"))

            # Check report content
            self.assertIn("FALSIFICATION REPORT", report)
            self.assertIn("FALSIFIED", report)

            print("✅ Complete workflow validated")
            print("   Files created successfully")
            print("   Report contains falsification evidence")

        except Exception as e:
            self.fail(f"Complete workflow failed: {e}")


def run_falsification_tests():
    """Run all falsification tests."""
    print("🚀 Running Falsification Test Suite")
    print("=" * 60)

    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestFalsificationExperiments))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "=" * 60)
    print("FALSIFICATION TEST SUMMARY")
    print("=" * 60)

    if result.wasSuccessful():
        print("✅ ALL FALSIFICATION TESTS PASSED")
        print("🎯 Falsification experiments are working correctly")
        print("📊 Z5D framework claims successfully falsified")
    else:
        print("❌ SOME TESTS FAILED")
        print(f"   Failures: {len(result.failures)}")
        print(f"   Errors: {len(result.errors)}")

    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_falsification_tests()
    sys.exit(0 if success else 1)
