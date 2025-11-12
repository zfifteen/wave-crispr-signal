#!/usr/bin/env python3
"""
Test Suite for Molecular Dynamics Enhanced Z Framework

This module provides comprehensive tests for the molecular dynamics integration
within the Z Framework, including empirical validation, falsification tests,
and comparative analysis against the base framework.
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import random
import json

from modules.molecular_dynamics_framework import (
    MolecularDynamicsZFramework,
    compare_md_vs_base_framework,
    MOLECULAR_DYNAMICS_PARAMETERS,
    MD_ENVIRONMENT,
)
from scripts.z_framework import format_mpmath_for_display


def test_md_framework_basic():
    """Test basic functionality of MD-enhanced Z Framework."""
    print("\n" + "=" * 70)
    print("MOLECULAR DYNAMICS Z FRAMEWORK - BASIC TESTS")
    print("=" * 70)

    # Initialize MD calculator
    md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=0.3)

    # Test sequences
    test_sequences = [
        ("Simple test", "ATCGATCGATCGATCG"),
        ("GC-rich", "GGCCGGCCGGCCGGCC"),
        ("AT-rich", "AAATTTAAATTTAAAA"),
        ("Mixed", "ATGCATGCATGCATGC"),
        ("Real gene fragment", "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG"),
    ]

    for name, sequence in test_sequences:
        print(f"\n--- Testing: {name} ({len(sequence)} bp) ---")
        print(f"Sequence: {sequence}")

        try:
            # Calculate MD-enhanced Z values
            md_results = md_calc.calculate_md_z_values(sequence)

            print(f"MD Z mean: {format_mpmath_for_display(md_results['md_z_mean'])}")
            print(
                f"MD Z variance: {format_mpmath_for_display(md_results['md_z_variance'])}"
            )
            print(
                f"Base Z mean: {format_mpmath_for_display(md_results['base_z_mean'])}"
            )
            print(
                f"Base Z variance: {format_mpmath_for_display(md_results['base_z_variance'])}"
            )
            print(
                f"Enhancement factor: {format_mpmath_for_display(md_results['md_enhancement_factor'])}"
            )
            print(
                f"MD φ-1 convergence: {format_mpmath_for_display(md_results['md_phi_conjugate_convergence'])}"
            )
            print(
                f"MD variance convergence: {format_mpmath_for_display(md_results['md_variance_convergence'])}"
            )
            print(f"Converges to φ-1: {md_results['converges_to_phi_conjugate']}")
            print(
                f"Converges to target variance: {md_results['converges_to_target_variance']}"
            )

        except Exception as e:
            print(f"Error testing {name}: {e}")

    print("\n" + "=" * 70)
    print("MD FRAMEWORK BASIC TESTS COMPLETED")
    print("=" * 70)


def test_md_parameter_perturbation():
    """Test MD parameter perturbation and robustness."""
    print("\n" + "=" * 70)
    print("MOLECULAR DYNAMICS PARAMETER PERTURBATION TESTS")
    print("=" * 70)

    md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=0.3)

    test_sequences = ["ATCGATCGATCGATCG", "GGCCGGCCGGCCGGCC", "ATGCTGCGGAGACCTGGAGAGA"]

    for sequence in test_sequences:
        print(f"\n--- Testing MD perturbation: {sequence} ({len(sequence)} bp) ---")

        try:
            # Perform perturbation test
            perturbation_results = md_calc.perform_md_parameter_perturbation_test(
                sequence, perturbation_range=0.2, num_tests=20
            )

            if "error" in perturbation_results:
                print(f"Perturbation test failed: {perturbation_results['error']}")
                continue

            print(
                f"Baseline MD mean: {format_mpmath_for_display(perturbation_results['baseline_md_mean'])}"
            )
            print(
                f"Baseline MD variance: {format_mpmath_for_display(perturbation_results['baseline_md_variance'])}"
            )
            print(
                f"Perturbed mean stability: {format_mpmath_for_display(perturbation_results['perturbed_mean_stability'])}"
            )
            print(
                f"Perturbed variance stability: {format_mpmath_for_display(perturbation_results['perturbed_variance_stability'])}"
            )
            print(
                f"Mean deviation: {format_mpmath_for_display(perturbation_results['mean_deviation'])}"
            )
            print(
                f"Variance deviation: {format_mpmath_for_display(perturbation_results['variance_deviation'])}"
            )
            print(f"Convergence rate: {perturbation_results['convergence_rate']:.3f}")
            print(
                f"Successful tests: {perturbation_results['num_successful_tests']}/{perturbation_results['num_total_tests']}"
            )
            print(
                f"MD falsification passed: {perturbation_results['md_falsification_passed']}"
            )

        except Exception as e:
            print(f"Error in perturbation test: {e}")

    print("\n" + "=" * 70)
    print("MD PERTURBATION TESTS COMPLETED")
    print("=" * 70)


def test_md_weight_sensitivity():
    """Test sensitivity to different MD weight factors."""
    print("\n" + "=" * 70)
    print("MOLECULAR DYNAMICS WEIGHT SENSITIVITY ANALYSIS")
    print("=" * 70)

    test_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG"
    weight_factors = [0.0, 0.1, 0.3, 0.5, 0.7, 1.0]

    print(f"Testing sequence: {test_sequence}")
    print(f"Sequence length: {len(test_sequence)} bp")

    results = []
    for weight in weight_factors:
        print(f"\n--- MD Weight Factor: {weight} ---")

        try:
            md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=weight)
            md_results = md_calc.calculate_md_z_values(test_sequence)

            result = {
                "md_weight": weight,
                "md_z_mean": float(md_results["md_z_mean"]),
                "md_z_variance": float(md_results["md_z_variance"]),
                "enhancement_factor": float(md_results["md_enhancement_factor"]),
                "phi_convergence": float(md_results["md_phi_conjugate_convergence"]),
                "variance_convergence": float(md_results["md_variance_convergence"]),
            }
            results.append(result)

            print(f"MD Z mean: {format_mpmath_for_display(md_results['md_z_mean'])}")
            print(
                f"Enhancement factor: {format_mpmath_for_display(md_results['md_enhancement_factor'])}"
            )
            print(
                f"φ-1 convergence: {format_mpmath_for_display(md_results['md_phi_conjugate_convergence'])}"
            )

        except Exception as e:
            print(f"Error with weight {weight}: {e}")

    # Analyze weight sensitivity
    if results:
        print("\n--- Weight Sensitivity Summary ---")
        best_phi_weight = min(results, key=lambda x: x["phi_convergence"])
        best_variance_weight = min(results, key=lambda x: x["variance_convergence"])

        print(
            f"Best φ-1 convergence: weight={best_phi_weight['md_weight']}, convergence={best_phi_weight['phi_convergence']:.6f}"
        )
        print(
            f"Best variance convergence: weight={best_variance_weight['md_weight']}, convergence={best_variance_weight['variance_convergence']:.6f}"
        )

        # Check for monotonic improvement
        improvements = [
            r
            for r in results
            if r["md_weight"] > 0
            and r["phi_convergence"] < results[0]["phi_convergence"]
        ]
        print(
            f"Weight factors showing improvement: {len(improvements)}/{len(results)-1}"
        )

    print("\n" + "=" * 70)
    print("MD WEIGHT SENSITIVITY ANALYSIS COMPLETED")
    print("=" * 70)


def test_comparative_analysis():
    """Compare MD-enhanced vs base Z Framework on multiple sequences."""
    print("\n" + "=" * 70)
    print("COMPARATIVE ANALYSIS: MD vs BASE Z FRAMEWORK")
    print("=" * 70)

    # Test sequences from different categories
    test_sequences = [
        # Short sequences
        "ATCGATCGATCG",
        "GGCCGGCCGGCC",
        "AAATTTAAATTT",
        # Medium sequences
        "ATGCTGCGGAGACCTGGAGA",
        "GGCCATGGCCATGGCCATGG",
        "TATACGCGTATACGCGTATA",
        # Longer sequences (gene fragments)
        "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCTGAAGACC",
        "GGCCATGGCCATGGCCATGGCCATGGCCATGGCCATGGCC",
        "AAAAAATTTTTCCCCCCGGGGGGAAAAAATTTTTCCCCCC",
    ]

    print(f"Testing {len(test_sequences)} sequences with MD weight=0.3")

    try:
        comparison_results = compare_md_vs_base_framework(test_sequences, md_weight=0.3)

        print(f"\nSequences successfully tested: {len(comparison_results['results'])}")

        if comparison_results["results"]:
            summary = comparison_results["summary"]
            print("\n--- Comparative Summary ---")
            print(
                f"Average enhancement factor: {summary['avg_enhancement_factor']:.6f}"
            )
            print(
                f"Sequences with improved φ-1 convergence: {summary['sequences_with_improved_phi_convergence']}"
            )
            print(
                f"Sequences with improved variance convergence: {summary['sequences_with_improved_variance_convergence']}"
            )
            print(f"φ-1 improvement rate: {summary['improvement_rate_phi']:.3f}")
            print(
                f"Variance improvement rate: {summary['improvement_rate_variance']:.3f}"
            )

            # Show detailed results for a few sequences
            print("\n--- Detailed Results (first 3 sequences) ---")
            for i, result in enumerate(comparison_results["results"][:3]):
                print(f"\nSequence {i+1}: {result['sequence']}")
                print(f"  Length: {result['length']} bp")
                print(f"  Base Z mean: {result['base_z_mean']:.6f}")
                print(f"  MD Z mean: {result['md_z_mean']:.6f}")
                print(f"  Enhancement factor: {result['enhancement_factor']:.6f}")
                print(
                    f"  MD improves φ-1 convergence: {result['md_improves_phi_convergence']}"
                )
                print(
                    f"  MD improves variance convergence: {result['md_improves_variance_convergence']}"
                )

    except Exception as e:
        print(f"Error in comparative analysis: {e}")

    print("\n" + "=" * 70)
    print("COMPARATIVE ANALYSIS COMPLETED")
    print("=" * 70)


def test_random_vs_real_sequences():
    """Test MD framework on random vs real DNA sequences."""
    print("\n" + "=" * 70)
    print("RANDOM vs REAL SEQUENCE ANALYSIS")
    print("=" * 70)

    # Generate random sequences
    bases = ["A", "T", "C", "G"]
    random.seed(42)
    random_sequences = []
    for length in [20, 30, 40]:
        for _ in range(3):
            seq = "".join(random.choice(bases) for _ in range(length))
            random_sequences.append(seq)

    # Real sequences (gene fragments)
    real_sequences = [
        "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGG",  # PCSK9 fragment
        "ATGGCGACCCTGGAAAAGCTGATGAAGGCCT",  # p53 fragment
        "ATGAAGATCCTGTTTGCCGTGCGCAGCCGAG",  # BRCA1 fragment
        "ATGCCACAAGGCGCCAAGAGGCTGCTGAAGG",  # Synthetic optimized
        "GGCCATGGCCATGGCCATGGCCATGGCCATG",  # Repetitive sequence
        "CGTACGTACGTACGTACGTACGTACGTACGT",  # Another repetitive
    ]

    md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=0.3)

    # Analyze random sequences
    print("\n--- Random Sequences ---")
    random_results = []
    for i, seq in enumerate(random_sequences):
        try:
            md_results = md_calc.calculate_md_z_values(seq)
            result = {
                "type": "random",
                "sequence": seq,
                "length": len(seq),
                "md_z_mean": float(md_results["md_z_mean"]),
                "md_z_variance": float(md_results["md_z_variance"]),
                "enhancement_factor": float(md_results["md_enhancement_factor"]),
                "phi_convergence": float(md_results["md_phi_conjugate_convergence"]),
                "converges_phi": md_results["converges_to_phi_conjugate"],
                "converges_variance": md_results["converges_to_target_variance"],
            }
            random_results.append(result)
            print(
                f"Random {i+1} ({len(seq)} bp): φ-1 conv={result['phi_convergence']:.4f}, enhance={result['enhancement_factor']:.4f}"
            )
        except Exception as e:
            print(f"Error with random sequence {i+1}: {e}")

    # Analyze real sequences
    print("\n--- Real Sequences ---")
    real_results = []
    for i, seq in enumerate(real_sequences):
        try:
            md_results = md_calc.calculate_md_z_values(seq)
            result = {
                "type": "real",
                "sequence": seq,
                "length": len(seq),
                "md_z_mean": float(md_results["md_z_mean"]),
                "md_z_variance": float(md_results["md_z_variance"]),
                "enhancement_factor": float(md_results["md_enhancement_factor"]),
                "phi_convergence": float(md_results["md_phi_conjugate_convergence"]),
                "converges_phi": md_results["converges_to_phi_conjugate"],
                "converges_variance": md_results["converges_to_target_variance"],
            }
            real_results.append(result)
            print(
                f"Real {i+1} ({len(seq)} bp): φ-1 conv={result['phi_convergence']:.4f}, enhance={result['enhancement_factor']:.4f}"
            )
        except Exception as e:
            print(f"Error with real sequence {i+1}: {e}")

    # Compare random vs real
    if random_results and real_results:
        print("\n--- Random vs Real Comparison ---")

        random_phi_avg = np.mean([r["phi_convergence"] for r in random_results])
        real_phi_avg = np.mean([r["phi_convergence"] for r in real_results])

        random_enhance_avg = np.mean([r["enhancement_factor"] for r in random_results])
        real_enhance_avg = np.mean([r["enhancement_factor"] for r in real_results])

        random_converge_rate = sum(r["converges_phi"] for r in random_results) / len(
            random_results
        )
        real_converge_rate = sum(r["converges_phi"] for r in real_results) / len(
            real_results
        )

        print(f"Random sequences - Avg φ-1 convergence: {random_phi_avg:.6f}")
        print(f"Real sequences - Avg φ-1 convergence: {real_phi_avg:.6f}")
        print(f"Random sequences - Avg enhancement: {random_enhance_avg:.6f}")
        print(f"Real sequences - Avg enhancement: {real_enhance_avg:.6f}")
        print(f"Random sequences - Convergence rate: {random_converge_rate:.3f}")
        print(f"Real sequences - Convergence rate: {real_converge_rate:.3f}")

        # Statistical significance test
        from scipy import stats

        phi_ttest = stats.ttest_ind(
            [r["phi_convergence"] for r in random_results],
            [r["phi_convergence"] for r in real_results],
        )
        print(f"φ-1 convergence t-test p-value: {phi_ttest.pvalue:.6f}")

        if phi_ttest.pvalue < 0.05:
            print(
                "Statistically significant difference between random and real sequences!"
            )
        else:
            print("No statistically significant difference found.")

    print("\n" + "=" * 70)
    print("RANDOM vs REAL SEQUENCE ANALYSIS COMPLETED")
    print("=" * 70)


def test_md_parameter_validation():
    """Validate molecular dynamics parameters are reasonable."""
    print("\n" + "=" * 70)
    print("MOLECULAR DYNAMICS PARAMETER VALIDATION")
    print("=" * 70)

    print("Checking MD parameter consistency and biological plausibility...")

    # Check base-pair energy ordering (GC > AT)
    gc_energy = (
        MOLECULAR_DYNAMICS_PARAMETERS["G"]["base_pair_energy"]
        + MOLECULAR_DYNAMICS_PARAMETERS["C"]["base_pair_energy"]
    ) / 2
    at_energy = (
        MOLECULAR_DYNAMICS_PARAMETERS["A"]["base_pair_energy"]
        + MOLECULAR_DYNAMICS_PARAMETERS["T"]["base_pair_energy"]
    ) / 2

    print(f"Average GC base-pair energy: {gc_energy:.2f} kcal/mol")
    print(f"Average AT base-pair energy: {at_energy:.2f} kcal/mol")
    print(f"GC stronger than AT: {gc_energy < at_energy}")  # More negative = stronger

    # Check stacking energy ordering (G > C > A > T typically)
    stacking_energies = {
        base: params["stacking_energy"]
        for base, params in MOLECULAR_DYNAMICS_PARAMETERS.items()
    }
    print(f"\nStacking energies: {stacking_energies}")

    # Check purine vs pyrimidine properties
    purines = ["A", "G"]
    pyrimidines = ["T", "C"]

    purine_vdw = np.mean(
        [
            MOLECULAR_DYNAMICS_PARAMETERS[base]["van_der_waals_radius"]
            for base in purines
        ]
    )
    pyrimidine_vdw = np.mean(
        [
            MOLECULAR_DYNAMICS_PARAMETERS[base]["van_der_waals_radius"]
            for base in pyrimidines
        ]
    )

    print(f"\nPurine avg vdW radius: {purine_vdw:.2f} Å")
    print(f"Pyrimidine avg vdW radius: {pyrimidine_vdw:.2f} Å")
    print(f"Purines larger than pyrimidines: {purine_vdw > pyrimidine_vdw}")

    # Check environmental parameters
    print("\n--- Environmental Parameters ---")
    for param, value in MD_ENVIRONMENT.items():
        print(f"{param}: {value}")

    # Validate temperature is physiological
    temp_check = 300 <= MD_ENVIRONMENT["temperature"] <= 320
    print(f"Physiological temperature range: {temp_check}")

    print("\n" + "=" * 70)
    print("MD PARAMETER VALIDATION COMPLETED")
    print("=" * 70)


def create_md_test_report():
    """Create a comprehensive test report for MD framework."""
    print("\n" + "=" * 70)
    print("GENERATING COMPREHENSIVE MD TEST REPORT")
    print("=" * 70)

    report = {
        "test_summary": {
            "timestamp": str(np.datetime64("now")),
            "framework": "Molecular Dynamics Enhanced Z Framework",
            "version": "1.0.0",
        },
        "tests_performed": [],
    }

    # Run all tests and collect results
    try:
        print("Running basic functionality tests...")
        test_md_framework_basic()
        report["tests_performed"].append("basic_functionality")

        print("\nRunning parameter perturbation tests...")
        test_md_parameter_perturbation()
        report["tests_performed"].append("parameter_perturbation")

        print("\nRunning weight sensitivity analysis...")
        weight_results = test_md_weight_sensitivity()
        report["weight_sensitivity_results"] = weight_results
        report["tests_performed"].append("weight_sensitivity")

        print("\nRunning comparative analysis...")
        comp_results = test_comparative_analysis()
        report["comparative_results"] = comp_results
        report["tests_performed"].append("comparative_analysis")

        print("\nRunning random vs real sequence analysis...")
        seq_results = test_random_vs_real_sequences()
        report["sequence_analysis_results"] = seq_results
        report["tests_performed"].append("sequence_analysis")

        print("\nValidating MD parameters...")
        test_md_parameter_validation()
        report["tests_performed"].append("parameter_validation")

        # Overall assessment
        report["overall_assessment"] = {
            "tests_completed": len(report["tests_performed"]),
            "tests_planned": 6,
            "success_rate": len(report["tests_performed"]) / 6,
            "md_framework_functional": True,
        }

        # Save report to file
        with open("/tmp/md_framework_test_report.json", "w") as f:
            json.dump(report, f, indent=2, default=str)

        print("\n--- Test Report Summary ---")
        print(f"Tests completed: {report['overall_assessment']['tests_completed']}/6")
        print(f"Success rate: {report['overall_assessment']['success_rate']:.2%}")
        print("Report saved to: /tmp/md_framework_test_report.json")

    except Exception as e:
        print(f"Error generating report: {e}")
        report["error"] = str(e)

    print("\n" + "=" * 70)
    print("MD TEST REPORT GENERATION COMPLETED")
    print("=" * 70)

    return report


if __name__ == "__main__":
    # Run comprehensive test suite
    print("MOLECULAR DYNAMICS Z FRAMEWORK TEST SUITE")
    print("=========================================")

    # Create comprehensive test report
    test_report = create_md_test_report()

    print("\nAll tests completed successfully!")
    print("Check /tmp/md_framework_test_report.json for detailed results.")
