"""
Test suite for CRISPR Guide Designer

This module contains comprehensive tests for the CRISPR guide designer
functionality using signal-theoretic DNA analysis.
"""

import unittest
import numpy as np
import os
import json
from unittest.mock import patch
import sys

# Add applications directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "applications"))

from applications.crispr_guide_designer import CRISPRGuideDesigner
from applications.wave_crispr_metrics import WaveCRISPRMetrics
from applications.crispr_visualization import CRISPRVisualizer


class TestCRISPRGuideDesigner(unittest.TestCase):
    """Test cases for CRISPRGuideDesigner class."""

    def setUp(self):
        """Set up test fixtures."""
        self.designer = CRISPRGuideDesigner()
        self.test_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
        self.test_guide = "GCTGCGGAGACCTGGAGAGA"

    def test_initialization(self):
        """Test designer initialization with different parameters."""
        # Default initialization
        designer = CRISPRGuideDesigner()
        self.assertEqual(designer.guide_length, 20)
        self.assertEqual(designer.pam_pattern, "[ATCG]GG")

        # Custom initialization
        designer_custom = CRISPRGuideDesigner(pam_pattern="TTT[ACG]", guide_length=23)
        self.assertEqual(designer_custom.guide_length, 23)
        self.assertEqual(designer_custom.pam_pattern, "TTT[ACG]")

    def test_build_waveform(self):
        """Test waveform construction from DNA sequence."""
        test_seq = "ATCG"
        wave = self.designer.build_waveform(test_seq)

        # Check output type and length
        self.assertIsInstance(wave, np.ndarray)
        self.assertEqual(len(wave), len(test_seq))
        self.assertTrue(np.iscomplexobj(wave))

        # Test with position scaling
        zn_map = {1: 0.5}
        wave_scaled = self.designer.build_waveform(test_seq, zn_map=zn_map)
        self.assertEqual(len(wave_scaled), len(test_seq))

        # Test empty sequence
        with self.assertRaises(ValueError):
            self.designer.build_waveform("")

    def test_compute_spectrum(self):
        """Test frequency spectrum computation."""
        wave = self.designer.build_waveform(self.test_guide)
        spectrum = self.designer.compute_spectrum(wave)

        # Check output properties
        self.assertIsInstance(spectrum, np.ndarray)
        self.assertEqual(len(spectrum), len(wave))
        self.assertTrue(np.all(spectrum >= 0))  # Magnitude should be non-negative

    def test_normalized_entropy(self):
        """Test normalized entropy calculation."""
        # Test with uniform spectrum
        uniform_spectrum = np.ones(10)
        entropy = self.designer.normalized_entropy(uniform_spectrum)
        expected_entropy = np.log2(10)  # Maximum entropy for uniform distribution
        self.assertAlmostEqual(entropy, expected_entropy, places=5)

        # Test with single peak
        peaked_spectrum = np.zeros(10)
        peaked_spectrum[0] = 1.0
        entropy_peaked = self.designer.normalized_entropy(peaked_spectrum)
        self.assertEqual(entropy_peaked, 0.0)

        # Test with empty spectrum
        self.assertTrue(np.isnan(self.designer.normalized_entropy(np.array([]))))

    def test_count_sidelobes(self):
        """Test sidelobe counting functionality."""
        # Test with known spectrum
        spectrum = np.array([1.0, 0.5, 0.3, 0.1, 0.05])
        sidelobes = self.designer.count_sidelobes(spectrum, threshold_ratio=0.25)

        # Should count peaks above 25% of maximum (1.0 * 0.25 = 0.25)
        expected_count = 3  # 1.0, 0.5, 0.3 are above 0.25
        self.assertEqual(sidelobes, expected_count)

        # Test with different threshold
        sidelobes_strict = self.designer.count_sidelobes(spectrum, threshold_ratio=0.5)
        expected_strict = 1  # Only 1.0 is above 0.5
        self.assertEqual(sidelobes_strict, expected_strict)

    def test_calculate_on_target_score(self):
        """Test on-target efficiency score calculation."""
        score = self.designer.calculate_on_target_score(self.test_guide)

        # Check score range
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

        # Test with different GC contents
        low_gc_guide = "AAATTTAAATTTAAATTTAA"  # 0% GC
        high_gc_guide = "GGCCGGCCGGCCGGCCGGCC"  # 100% GC
        optimal_gc_guide = "ATGCATGCATGCATGCATGC"  # 50% GC

        score_low = self.designer.calculate_on_target_score(low_gc_guide)
        score_high = self.designer.calculate_on_target_score(high_gc_guide)
        score_optimal = self.designer.calculate_on_target_score(optimal_gc_guide)

        # Optimal GC should score higher than extreme GC contents
        self.assertGreater(score_optimal, score_low)
        self.assertGreater(score_optimal, score_high)

    def test_calculate_off_target_risk(self):
        """Test off-target risk calculation."""
        target_seq = self.test_guide

        # Identical sequences should have high risk
        identical_risk = self.designer.calculate_off_target_risk(
            self.test_guide, target_seq, target_seq
        )
        self.assertGreater(identical_risk, 0.8)

        # Completely different sequences should have lower risk
        different_seq = "AAATTTAAATTTAAATTTAA"
        different_risk = self.designer.calculate_off_target_risk(
            self.test_guide, target_seq, different_seq
        )
        self.assertLess(different_risk, identical_risk)

        # Check score range
        self.assertGreaterEqual(different_risk, 0.0)
        self.assertLessEqual(different_risk, 1.0)

    def test_find_pam_sites(self):
        """Test PAM site detection."""
        # Sequence with known PAM sites
        test_seq = "ATCGAGGCTGCGGAGACCTGGAGAGAAAGCAGTGG"  # Contains AGG and TGG
        pam_sites = self.designer.find_pam_sites(test_seq)

        # Check return format
        self.assertIsInstance(pam_sites, list)
        for pos, pam in pam_sites:
            self.assertIsInstance(pos, int)
            self.assertIsInstance(pam, str)
            self.assertEqual(len(pam), 3)  # NGG pattern

        # Check that PAM sequences match pattern
        for pos, pam in pam_sites:
            self.assertTrue(pam.endswith("GG"))
            self.assertIn(pam[0], "ATCG")

    def test_design_guides(self):
        """Test complete guide design workflow."""
        guides = self.designer.design_guides(self.test_sequence, num_guides=3)

        # Check return format
        self.assertIsInstance(guides, list)
        self.assertLessEqual(len(guides), 3)

        # Check guide properties
        for guide in guides:
            self.assertIn("sequence", guide)
            self.assertIn("position", guide)
            self.assertIn("pam_position", guide)
            self.assertIn("pam_sequence", guide)
            self.assertIn("on_target_score", guide)
            self.assertIn("gc_content", guide)
            self.assertIn("length", guide)

            # Validate sequence properties
            self.assertEqual(len(guide["sequence"]), self.designer.guide_length)
            self.assertGreaterEqual(guide["on_target_score"], 0.0)
            self.assertLessEqual(guide["on_target_score"], 1.0)
            self.assertGreaterEqual(guide["gc_content"], 0.0)
            self.assertLessEqual(guide["gc_content"], 1.0)

            # Check sequence validity
            valid_bases = set("ATCG")
            self.assertTrue(all(base in valid_bases for base in guide["sequence"]))

    def test_design_guides_with_phase_weighted_scorecard(self):
        """Test guide design with the phase-weighted scorecard."""
        guides = self.designer.design_guides(
            self.test_sequence, num_guides=3, use_phase_weighted_scorecard=True
        )

        # Check return format
        self.assertIsInstance(guides, list)
        self.assertLessEqual(len(guides), 3)

        # Check guide properties
        for guide in guides:
            self.assertIn("sequence", guide)
            self.assertIn("primary_score", guide)
            self.assertIn("phase_weighted_score", guide)
            self.assertIn("phase_weighted_entropy", guide)
            self.assertIn("phase_weighted_sidelobes", guide)

            # Verify that the primary score is the phase-weighted score
            self.assertEqual(guide["primary_score"], guide["phase_weighted_score"])

        # Check that guides are sorted by primary_score
        scores = [g["primary_score"] for g in guides]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_calculate_phase_weighted_score(self):
        """Test phase-weighted score calculation."""
        score_data = self.designer.calculate_phase_weighted_score(self.test_guide)

        # Check return format
        self.assertIsInstance(score_data, dict)
        required_keys = [
            "phase_weighted_score",
            "phase_weighted_entropy",
            "phase_weighted_sidelobes",
            "phase_weighted_diversity",
        ]
        for key in required_keys:
            self.assertIn(key, score_data)

        # Check score ranges
        self.assertGreaterEqual(score_data["phase_weighted_score"], 0.0)
        self.assertLessEqual(score_data["phase_weighted_score"], 1.0)
        self.assertGreaterEqual(score_data["phase_weighted_entropy"], 0.0)
        self.assertGreaterEqual(score_data["phase_weighted_sidelobes"], 0)

    def test_calculate_phase_weighted_score_edge_cases(self):
        """Test phase-weighted score calculation with edge cases."""
        # Test empty sequence
        with self.assertRaises(ValueError):
            self.designer.calculate_phase_weighted_score("")

        # Test invalid characters
        with self.assertRaises(ValueError):
            self.designer.calculate_phase_weighted_score("ATCGXYZ")

        # Test very short sequence (below minimum of 5 bases)
        # The phase-weighted scorecard requires at least 5 bases
        with self.assertRaises(RuntimeError):
            self.designer.calculate_phase_weighted_score("ATCG")

        # Test minimum valid length (5 bases)
        try:
            score_data = self.designer.calculate_phase_weighted_score("ATCGG")
            self.assertIsInstance(score_data, dict)
        except Exception as e:
            self.fail(f"5-base sequence should not raise exception: {e}")

    def test_csv_export_with_phase_weighted_scores(self):
        """Test CSV export includes phase-weighted scores correctly."""
        from applications.crispr_cli import format_csv

        # Design guides with phase-weighted scoring
        guides = self.designer.design_guides(
            self.test_sequence, num_guides=2, use_phase_weighted_scorecard=True
        )

        results = [
            {
                "sequence_name": "test_sequence",
                "sequence_length": len(self.test_sequence),
                "guides": guides,
            }
        ]

        # Format as CSV
        try:
            csv_output = format_csv(results)
            self.assertIsInstance(csv_output, str)
            self.assertIn("phase_weighted_score", csv_output)
            self.assertIn("phase_weighted_entropy", csv_output)
            self.assertIn("phase_weighted_sidelobes", csv_output)

            # Verify no empty float formatting errors
            lines = csv_output.split("\n")
            self.assertGreater(len(lines), 1)  # At least header + 1 guide
        except Exception as e:
            self.fail(f"CSV formatting should not fail: {e}")

    def test_csv_export_without_phase_weighted_scores(self):
        """Test CSV export handles missing phase-weighted scores gracefully."""
        from applications.crispr_cli import format_csv

        # Design guides WITHOUT phase-weighted scoring
        guides = self.designer.design_guides(
            self.test_sequence, num_guides=2, use_phase_weighted_scorecard=False
        )

        results = [
            {
                "sequence_name": "test_sequence",
                "sequence_length": len(self.test_sequence),
                "guides": guides,
            }
        ]

        # Format as CSV (should handle missing phase-weighted fields)
        try:
            csv_output = format_csv(results)
            self.assertIsInstance(csv_output, str)

            # Verify no errors with missing fields
            lines = csv_output.split("\n")
            self.assertGreater(len(lines), 1)
        except Exception as e:
            self.fail(f"CSV formatting should handle missing fields gracefully: {e}")

    def test_performance_optimization_filters(self):
        """Test that quality filters are applied before expensive calculations."""
        # Create a sequence with guides that will fail quality filters
        # (high GC content and poly-T stretches)
        bad_sequence = "GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG"

        # Design guides - should return empty or very few guides due to filters
        guides = self.designer.design_guides(
            bad_sequence, num_guides=5, use_phase_weighted_scorecard=True
        )

        # Should have filtered out most/all guides before expensive calculations
        # This test mainly ensures no errors occur during filtering
        self.assertIsInstance(guides, list)

    def test_predict_repair_outcomes(self):
        """Test repair outcome predictions."""
        repair = self.designer.predict_repair_outcomes(
            self.test_guide, self.test_sequence[:60]
        )

        # Check return format
        required_keys = [
            "nhej_probability",
            "mmej_probability",
            "hdr_efficiency",
            "spectral_entropy",
            "stability_score",
        ]
        for key in required_keys:
            self.assertIn(key, repair)

        # Check probability ranges
        for prob_key in ["nhej_probability", "mmej_probability", "hdr_efficiency"]:
            self.assertGreaterEqual(repair[prob_key], 0.0)
            self.assertLessEqual(repair[prob_key], 1.0)

        # Check that probabilities are reasonable (not all zero)
        total_prob = (
            repair["nhej_probability"]
            + repair["mmej_probability"]
            + repair["hdr_efficiency"]
        )
        self.assertGreater(total_prob, 0.0)


class TestWaveCRISPRMetrics(unittest.TestCase):
    """Test cases for WaveCRISPRMetrics class."""

    def setUp(self):
        """Set up test fixtures."""
        self.metrics = WaveCRISPRMetrics()
        self.test_sequence = "GCTGCGGAGACCTGGAGAGA"

    def test_calculate_spectral_entropy(self):
        """Test spectral entropy calculation."""
        entropy = self.metrics.calculate_spectral_entropy(self.test_sequence)

        self.assertIsInstance(entropy, float)
        self.assertGreaterEqual(entropy, 0.0)

        # Test with different sequences
        uniform_seq = "AAAAAAAAAAAAAAAAAAA"  # Uniform sequence
        random_seq = "ATCCGTAGACGTTAGCATG"   # More random sequence

        entropy_uniform = self.metrics.calculate_spectral_entropy(uniform_seq)
        entropy_random = self.metrics.calculate_spectral_entropy(random_seq)

        # Both should be valid entropy values
        self.assertIsInstance(entropy_uniform, float)
        self.assertIsInstance(entropy_random, float)
        self.assertGreaterEqual(entropy_uniform, 0.0)
        self.assertGreaterEqual(entropy_random, 0.0)

    def test_calculate_spectral_complexity(self):
        """Test spectral complexity analysis."""
        complexity = self.metrics.calculate_spectral_complexity(self.test_sequence)

        # Check return format
        required_keys = [
            "spectral_entropy",
            "spectral_flatness",
            "spectral_centroid",
            "spectral_rolloff",
            "zero_crossing_rate",
            "peak_count",
        ]
        for key in required_keys:
            self.assertIn(key, complexity)
            self.assertIsInstance(complexity[key], (int, float))

        # Check reasonable ranges
        self.assertGreaterEqual(complexity["spectral_flatness"], 0.0)
        self.assertLessEqual(complexity["spectral_flatness"], 1.0)
        self.assertGreaterEqual(complexity["spectral_rolloff"], 0.0)
        self.assertLessEqual(complexity["spectral_rolloff"], 1.0)

    def test_calculate_harmonic_content(self):
        """Test harmonic content analysis."""
        harmonics = self.metrics.calculate_harmonic_content(self.test_sequence)

        # Check return format
        required_keys = [
            "fundamental_frequency",
            "fundamental_power",
            "harmonic_power",
            "total_harmonic_distortion",
            "harmonic_to_noise_ratio",
            "num_peaks",
        ]
        for key in required_keys:
            self.assertIn(key, harmonics)

        # Check reasonable ranges
        self.assertGreaterEqual(harmonics["fundamental_frequency"], 0.0)
        self.assertLessEqual(harmonics["fundamental_frequency"], 1.0)
        self.assertGreaterEqual(harmonics["fundamental_power"], 0.0)
        self.assertGreaterEqual(harmonics["harmonic_power"], 0.0)
        self.assertGreaterEqual(harmonics["num_peaks"], 0)

    def test_calculate_mutational_sensitivity(self):
        """Test mutational sensitivity analysis."""
        # Test with subset of positions for speed
        positions = list(range(0, len(self.test_sequence), 5))
        sensitivity = self.metrics.calculate_mutational_sensitivity(
            self.test_sequence, positions=positions
        )

        # Check return format
        required_keys = [
            "average_sensitivity",
            "max_sensitivity",
            "sensitivity_variance",
            "hotspot_positions",
            "position_effects",
        ]
        for key in required_keys:
            self.assertIn(key, sensitivity)

        # Check position effects format
        for pos in positions:
            if pos in sensitivity["position_effects"]:
                effect = sensitivity["position_effects"][pos]
                self.assertIn("average_disruption", effect)
                self.assertIn("original_base", effect)
                self.assertIn(effect["original_base"], "ATCG")

    def test_calculate_comprehensive_score(self):
        """Test comprehensive scoring system."""
        target_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        comprehensive = self.metrics.calculate_comprehensive_score(
            self.test_sequence, target_seq
        )

        # Check return format
        required_keys = [
            "composite_score",
            "component_scores",
            "detailed_metrics",
            "weights",
        ]
        for key in required_keys:
            self.assertIn(key, comprehensive)

        # Check score range
        self.assertGreaterEqual(comprehensive["composite_score"], 0.0)
        self.assertLessEqual(comprehensive["composite_score"], 1.0)

        # Check component scores
        components = comprehensive["component_scores"]
        for component, score in components.items():
            self.assertIsInstance(score, (int, float))
            if component != "spectral":  # Spectral can exceed 1 before clipping
                self.assertGreaterEqual(score, 0.0)
                self.assertLessEqual(score, 1.0)

        # Check weights sum to 1
        weights = comprehensive["weights"]
        weight_sum = sum(weights.values())
        self.assertAlmostEqual(weight_sum, 1.0, places=5)


class TestCRISPRVisualization(unittest.TestCase):
    """Test cases for CRISPRVisualization class."""

    def setUp(self):
        """Set up test fixtures."""
        self.designer = CRISPRGuideDesigner()
        self.visualizer = CRISPRVisualizer(self.designer)
        self.test_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
        self.guides = self.designer.design_guides(self.test_sequence, num_guides=3)

    def test_plot_sequence_spectrum(self):
        """Test sequence spectrum plotting."""
        # Mock matplotlib to avoid actual plotting
        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = unittest.mock.Mock()
            mock_axes = [unittest.mock.Mock(), unittest.mock.Mock()]
            mock_subplots.return_value = (mock_fig, mock_axes)

            fig = self.visualizer.plot_sequence_spectrum(self.test_sequence)
            self.assertIsNotNone(fig)
            mock_subplots.assert_called_once()

    def test_plot_guide_scores(self):
        """Test guide score distribution plotting."""
        if not self.guides:
            self.skipTest("No guides generated for testing")

        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = unittest.mock.Mock()
            mock_axes = [
                [unittest.mock.Mock(), unittest.mock.Mock()],
                [unittest.mock.Mock(), unittest.mock.Mock()],
            ]
            mock_subplots.return_value = (mock_fig, mock_axes)

            fig = self.visualizer.plot_guide_scores(self.guides)
            self.assertIsNotNone(fig)

    def test_plot_guide_scores_empty(self):
        """Test guide score plotting with empty guide list."""
        with self.assertRaises(ValueError):
            self.visualizer.plot_guide_scores([])

    @unittest.skip("Mocking issue with matplotlib colorbar")
    def test_spectral_heatmap(self):
        """Test spectral heatmap creation."""
        sequences = {
            "seq1": "ATCGATCGATCGATCGATCG",
            "seq2": "GGCCGGCCGGCCGGCCGGCC",
            "seq3": "AAATTTAAATTTAAATTTAA",
        }

        with patch("matplotlib.pyplot.subplots") as mock_subplots:
            mock_fig = unittest.mock.Mock()
            mock_ax = unittest.mock.Mock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            mock_ax.get_position.return_value.frozen.return_value = unittest.mock.Mock(xmin=0, ymin=0, xmax=10, ymax=10)
            fig = self.visualizer.plot_spectral_heatmap(sequences)
            self.assertIsNotNone(fig)


class TestIntegration(unittest.TestCase):
    """Integration tests for complete workflows."""

    def setUp(self):
        """Set up test fixtures."""
        self.designer = CRISPRGuideDesigner()
        self.metrics = WaveCRISPRMetrics()
        self.test_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

    def test_complete_guide_design_workflow(self):
        """Test complete guide design and analysis workflow."""
        # Design guides
        guides = self.designer.design_guides(self.test_sequence, num_guides=3)
        self.assertGreater(len(guides), 0)

        # Add repair predictions
        for guide in guides:
            context_start = max(0, guide["position"] - 20)
            context_end = min(
                len(self.test_sequence), guide["position"] + guide["length"] + 20
            )
            context_seq = self.test_sequence[context_start:context_end]

            repair = self.designer.predict_repair_outcomes(
                guide["sequence"], context_seq
            )
            guide["repair_outcomes"] = repair

            # Validate repair outcomes
            self.assertIn("nhej_probability", repair)
            self.assertIn("mmej_probability", repair)
            self.assertIn("hdr_efficiency", repair)

        # Benchmark guides
        benchmark = self.metrics.benchmark_guide_set(guides, self.test_sequence)

        # Validate benchmark results
        self.assertIn("num_guides", benchmark)
        self.assertIn("score_statistics", benchmark)
        self.assertIn("top_guides", benchmark)
        self.assertEqual(benchmark["num_guides"], len(guides))

        # Check that top guides are properly ranked
        top_guides = benchmark["top_guides"]
        for i in range(len(top_guides) - 1):
            score1 = top_guides[i]["comprehensive_score"]["composite_score"]
            score2 = top_guides[i + 1]["comprehensive_score"]["composite_score"]
            self.assertGreaterEqual(score1, score2)

    def test_custom_pam_workflow(self):
        """Test workflow with custom PAM pattern."""
        # Test with Cas12a PAM (TTTV)
        cas12a_designer = CRISPRGuideDesigner(pam_pattern="TTT[ACG]", guide_length=23)

        # Create sequence with Cas12a PAM sites
        cas12a_sequence = "ATCGTTTAATCGATCGATCGATCGTTTCATCGATCGATCGATCGTTTGATCGATCGATCG"

        guides = cas12a_designer.design_guides(cas12a_sequence, num_guides=2)

        # Validate guide properties for Cas12a
        for guide in guides:
            self.assertEqual(guide["length"], 23)
            self.assertTrue(guide["pam_sequence"].startswith("TTT"))
            self.assertIn(guide["pam_sequence"][3], "ACG")

    def test_error_handling(self):
        """Test error handling for edge cases."""
        # Empty sequence
        guides_empty = self.designer.design_guides("", num_guides=1)
        self.assertEqual(len(guides_empty), 0)

        # Sequence without PAM sites
        no_pam_seq = "AAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTT"
        guides_no_pam = self.designer.design_guides(no_pam_seq, num_guides=1)
        self.assertEqual(len(guides_no_pam), 0)

        # Invalid characters in sequence
        with self.assertRaises(KeyError):
            invalid_seq = "ATCGXYZ"
            self.designer.build_waveform(invalid_seq)


class TestCLIIntegration(unittest.TestCase):
    """Test command-line interface functionality."""

    def test_cli_imports(self):
        """Test that CLI modules can be imported correctly."""
        try:
            from applications import crispr_cli

            self.assertTrue(hasattr(crispr_cli, "main"))
        except ImportError as e:
            self.fail(f"Failed to import CLI module: {e}")

    def test_json_output_format(self):
        """Test JSON output format structure."""
        designer = CRISPRGuideDesigner()
        guides = designer.design_guides(
            "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGG", num_guides=2
        )

        # Create expected output structure
        result = {
            "sequence_name": "test_sequence",
            "sequence_length": 42,
            "guides": guides,
        }

        # Test JSON serialization
        try:
            json_output = json.dumps(result, indent=2)
            parsed = json.loads(json_output)
            self.assertEqual(parsed["sequence_name"], "test_sequence")
            self.assertEqual(len(parsed["guides"]), len(guides))
        except (TypeError, ValueError) as e:
            self.fail(f"JSON serialization failed: {e}")


def run_performance_benchmark():
    """Run performance benchmark for guide design."""
    import time

    print("\n" + "=" * 60)
    print("PERFORMANCE BENCHMARK")
    print("=" * 60)

    designer = CRISPRGuideDesigner()
    test_sequence = (
        "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"
        * 2
    )

    # Benchmark guide design
    start_time = time.time()
    guides = designer.design_guides(test_sequence, num_guides=10)
    design_time = time.time() - start_time

    print(f"Sequence length: {len(test_sequence)} bp")
    print(f"Guide design time: {design_time:.3f} seconds")
    print(f"Guides found: {len(guides)}")
    print(f"Time per guide: {design_time/max(len(guides), 1):.3f} seconds")

    # Benchmark comprehensive scoring
    if guides:
        metrics = WaveCRISPRMetrics()
        start_time = time.time()
        comprehensive = metrics.calculate_comprehensive_score(
            guides[0]["sequence"], test_sequence
        )
        scoring_time = time.time() - start_time

        print(f"Comprehensive scoring time: {scoring_time:.3f} seconds")
        print(f"Composite score: {comprehensive['composite_score']:.3f}")


if __name__ == "__main__":
    # Run unit tests
    print("Running CRISPR Guide Designer Test Suite")
    print("=" * 50)

    # Create test suite
    test_suite = unittest.TestSuite()

    # Add test classes
    test_classes = [
        TestCRISPRGuideDesigner,
        TestWaveCRISPRMetrics,
        TestCRISPRVisualization,
        TestIntegration,
        TestCLIIntegration,
    ]

    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        test_suite.addTests(tests)

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(test_suite)

    # Print summary
    print("\nTest Summary:")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(
        f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%"
    )

    # Run performance benchmark if tests pass
    if len(result.failures) == 0 and len(result.errors) == 0:
        run_performance_benchmark()

    # Exit with appropriate code
    exit_code = 0 if len(result.failures) == 0 and len(result.errors) == 0 else 1
    sys.exit(exit_code)
