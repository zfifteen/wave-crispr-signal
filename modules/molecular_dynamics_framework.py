#!/usr/bin/env python3
"""
Molecular Dynamics Framework for Z Framework Enhancement

This module implements molecular dynamics parameters and calculations
to extend the Z Framework with explicit molecular dynamics considerations.
It includes base-pair energetics, stacking interactions, and steric constraints.

The enhanced framework maintains the core Z = n(Δₙ/Δₘₐₓ) formulation while
incorporating molecular dynamics influences through parameter adjustments.
"""

import mpmath as mp
import numpy as np
from typing import List, Dict
import logging
from z_framework import ZFrameworkCalculator, TARGET_VARIANCE

# Configure high precision for molecular dynamics calculations
mp.dps = 50

# Configure logging
logger = logging.getLogger(__name__)

# Molecular dynamics parameters for DNA bases
# These represent simplified coarse-grained properties for MD simulations
MOLECULAR_DYNAMICS_PARAMETERS = {
    "A": {
        # Base-pair energetics (kcal/mol)
        "base_pair_energy": -2.1,  # A-T base pair energy
        "stacking_energy": -1.2,  # Average stacking with adjacent bases
        "hydrogen_bond_energy": -1.9,  # H-bond contribution
        # Steric constraints (Angstroms)
        "van_der_waals_radius": 3.4,
        "steric_clash_threshold": 2.8,
        "backbone_flexibility": 0.85,  # Normalized flexibility
        # Dynamic properties
        "vibrational_frequency": 1580,  # cm⁻¹ (characteristic vibration)
        "rotational_barrier": 2.3,  # kcal/mol for glycosidic rotation
        "thermal_motion_amplitude": 0.12,  # Angstroms at 300K
        # Electrostatic properties
        "partial_charge": -0.15,
        "dipole_moment": 1.2,  # Debye
        "polarizability": 1.65,  # Cubic Angstroms
    },
    "T": {
        "base_pair_energy": -2.1,  # T-A base pair energy
        "stacking_energy": -1.0,
        "hydrogen_bond_energy": -1.9,
        "van_der_waals_radius": 3.2,
        "steric_clash_threshold": 2.7,
        "backbone_flexibility": 0.88,
        "vibrational_frequency": 1650,
        "rotational_barrier": 2.1,
        "thermal_motion_amplitude": 0.14,
        "partial_charge": 0.12,
        "dipole_moment": 1.1,
        "polarizability": 1.52,
    },
    "C": {
        "base_pair_energy": -3.1,  # C-G base pair energy (stronger)
        "stacking_energy": -1.4,
        "hydrogen_bond_energy": -2.8,  # Triple H-bond
        "van_der_waals_radius": 3.1,
        "steric_clash_threshold": 2.6,
        "backbone_flexibility": 0.82,
        "vibrational_frequency": 1620,
        "rotational_barrier": 2.5,
        "thermal_motion_amplitude": 0.10,
        "partial_charge": -0.18,
        "dipole_moment": 1.4,
        "polarizability": 1.44,
    },
    "G": {
        "base_pair_energy": -3.1,  # G-C base pair energy
        "stacking_energy": -1.5,  # Strongest stacking
        "hydrogen_bond_energy": -2.8,
        "van_der_waals_radius": 3.5,
        "steric_clash_threshold": 2.9,
        "backbone_flexibility": 0.80,
        "vibrational_frequency": 1590,
        "rotational_barrier": 2.7,
        "thermal_motion_amplitude": 0.09,
        "partial_charge": -0.12,
        "dipole_moment": 1.3,
        "polarizability": 1.72,
    },
}

# Temperature and environmental parameters
MD_ENVIRONMENT = {
    "temperature": 310.0,  # K (physiological)
    "ionic_strength": 0.15,  # M (physiological salt)
    "pH": 7.4,  # Physiological pH
    "pressure": 1.0,  # atm
    "dielectric_constant": 78.5,  # Water at 25°C
    "viscosity": 0.89e-3,  # Pa·s (water at 25°C)
}


class MolecularDynamicsZFramework(ZFrameworkCalculator):
    """
    Enhanced Z Framework calculator incorporating molecular dynamics parameters.

    Extends the base ZFrameworkCalculator with molecular dynamics considerations
    including base-pair energetics, stacking interactions, and steric constraints.
    """

    def __init__(self, precision_dps: int = 50, md_weight: float = 0.3):
        """
        Initialize the MD-enhanced Z Framework calculator.

        Args:
            precision_dps: Decimal precision for mpmath calculations
            md_weight: Weight factor for molecular dynamics contribution (0-1)
        """
        super().__init__(precision_dps)
        self.md_weight = mp.mpf(md_weight)
        self.md_params = MOLECULAR_DYNAMICS_PARAMETERS
        self.environment = MD_ENVIRONMENT

        logger.info(
            f"Initialized MD-enhanced Z Framework with weight factor {md_weight}"
        )

    def calculate_md_energy_term(self, sequence: str, position: int) -> mp.mpf:
        """
        Calculate molecular dynamics energy contribution for a specific position.

        Args:
            sequence: DNA sequence string
            position: Position in the sequence (0-indexed)

        Returns:
            MD energy term as high-precision value
        """
        if position >= len(sequence):
            raise ValueError(
                f"Position {position} out of range for sequence length {len(sequence)}"
            )

        base = sequence[position].upper()
        if base not in self.md_params:
            raise ValueError(f"Invalid base: {base}")

        params = self.md_params[base]

        # Base-pair energy contribution
        bp_energy = mp.mpf(params["base_pair_energy"])

        # Stacking interaction with adjacent bases
        stacking_energy = mp.mpf(0)
        if position > 0:
            prev_base = sequence[position - 1].upper()
            if prev_base in self.md_params:
                stacking_energy += mp.mpf(self.md_params[prev_base]["stacking_energy"])

        if position < len(sequence) - 1:
            next_base = sequence[position + 1].upper()
            if next_base in self.md_params:
                stacking_energy += mp.mpf(self.md_params[next_base]["stacking_energy"])

        # Steric constraint contribution
        steric_factor = mp.mpf(1) / (
            mp.mpf(params["van_der_waals_radius"])
            + mp.mpf(params["steric_clash_threshold"])
        )

        # Vibrational contribution (temperature-dependent)
        kb = mp.mpf("1.380649e-23")  # Boltzmann constant
        temp = mp.mpf(self.environment["temperature"])
        thermal_energy = kb * temp
        vib_freq = mp.mpf(params["vibrational_frequency"])
        vib_contribution = thermal_energy * vib_freq / mp.mpf("1000")  # Normalize

        # Combine all MD contributions
        total_md_energy = (
            bp_energy
            + stacking_energy * mp.mpf("0.5")
            + steric_factor * mp.mpf("10")
            + vib_contribution
        )

        return total_md_energy

    def calculate_md_enhanced_delta_n(self, sequence: str, position: int) -> mp.mpf:
        """
        Calculate Δₙ with molecular dynamics enhancement.

        Args:
            sequence: DNA sequence string
            position: Position in the sequence

        Returns:
            MD-enhanced delta value
        """
        # Get base Z Framework delta
        sequence_values = self.map_dna_sequence(sequence)
        base_delta = self.calculate_delta_n(sequence_values, position)

        # Calculate MD energy term
        md_energy = self.calculate_md_energy_term(sequence, position)

        # Calculate geodesic resolution at this position
        theta_prime = self.calculate_geodesic_resolution(position, k=0.3)

        # Combine base delta with MD contribution
        md_enhanced_delta = base_delta * (
            mp.mpf("1") + self.md_weight * md_energy * theta_prime / mp.mpf("10")
        )

        return md_enhanced_delta

    def calculate_md_z_values(self, sequence: str) -> Dict[str, mp.mpf]:
        """
        Calculate Z values with molecular dynamics enhancement.

        Args:
            sequence: DNA sequence string

        Returns:
            Dictionary containing MD-enhanced Z Framework results
        """
        if not sequence:
            raise ValueError("DNA sequence cannot be empty")

        n = len(sequence)
        self.map_dna_sequence(sequence)

        # Calculate MD-enhanced delta values
        md_delta_values = []
        for i in range(n):
            md_delta = self.calculate_md_enhanced_delta_n(sequence, i)
            md_delta_values.append(md_delta)

        # Calculate MD-enhanced delta_max
        md_delta_max = max(abs(delta) for delta in md_delta_values)
        if md_delta_max == 0:
            md_delta_max = mp.mpf("1e-10")  # Avoid division by zero

        # Calculate MD-enhanced Z values: Z = n(Δₙ/Δₘₐₓ)
        md_z_values = [mp.mpf(n) * (delta / md_delta_max) for delta in md_delta_values]

        # Statistical analysis
        md_z_mean = sum(md_z_values) / mp.mpf(len(md_z_values))
        md_z_variance = sum((z - md_z_mean) ** 2 for z in md_z_values) / mp.mpf(
            len(md_z_values)
        )
        md_z_std = mp.sqrt(md_z_variance)

        # Calculate convergence metrics
        phi_conjugate_convergence = abs(md_z_mean - self.phi_conjugate)
        variance_convergence = abs(md_z_variance - TARGET_VARIANCE)

        # Compare with base Z Framework
        base_results = self.calculate_z_values(sequence)
        md_enhancement_factor = abs(md_z_mean - base_results["z_mean"])

        results = {
            "sequence_length": mp.mpf(n),
            "md_weight": self.md_weight,
            "md_z_values": md_z_values,
            "md_delta_values": md_delta_values,
            "md_delta_max": md_delta_max,
            "md_z_mean": md_z_mean,
            "md_z_variance": md_z_variance,
            "md_z_std": md_z_std,
            "md_phi_conjugate_convergence": phi_conjugate_convergence,
            "md_variance_convergence": variance_convergence,
            "md_enhancement_factor": md_enhancement_factor,
            "base_z_mean": base_results["z_mean"],
            "base_z_variance": base_results["z_variance"],
            "converges_to_phi_conjugate": phi_conjugate_convergence < mp.mpf("0.01"),
            "converges_to_target_variance": variance_convergence < mp.mpf("0.01"),
        }

        logger.info(f"Calculated MD-enhanced Z values for sequence of length {n}")
        logger.info(f"MD Z mean: {md_z_mean}, Base Z mean: {base_results['z_mean']}")
        logger.info(f"MD enhancement factor: {md_enhancement_factor}")

        return results

    def perform_md_parameter_perturbation_test(
        self, sequence: str, perturbation_range: float = 0.2, num_tests: int = 50
    ) -> Dict[str, any]:
        """
        Perform falsification tests with MD parameter perturbations.

        Args:
            sequence: DNA sequence to test
            perturbation_range: Range of parameter perturbation (fraction)
            num_tests: Number of perturbation tests to run

        Returns:
            Dictionary containing falsification test results
        """
        import random

        random.seed(42)  # For reproducibility

        # Calculate baseline MD results
        baseline_results = self.calculate_md_z_values(sequence)
        baseline_mean = baseline_results["md_z_mean"]
        baseline_variance = baseline_results["md_z_variance"]

        # Storage for perturbation results
        perturbed_means = []
        perturbed_variances = []
        convergence_preserved = []

        # Save original parameters
        original_params = {}
        for base in self.md_params:
            original_params[base] = self.md_params[base].copy()

        for i in range(num_tests):
            try:
                # Perturb MD parameters
                for base in self.md_params:
                    for param_name in [
                        "base_pair_energy",
                        "stacking_energy",
                        "vibrational_frequency",
                    ]:
                        original_value = original_params[base][param_name]
                        perturbation = random.uniform(
                            -perturbation_range, perturbation_range
                        )
                        self.md_params[base][param_name] = original_value * (
                            1 + perturbation
                        )

                # Calculate perturbed results
                perturbed_results = self.calculate_md_z_values(sequence)
                perturbed_means.append(perturbed_results["md_z_mean"])
                perturbed_variances.append(perturbed_results["md_z_variance"])

                # Check convergence preservation
                phi_convergence = abs(
                    perturbed_results["md_z_mean"] - self.phi_conjugate
                ) < mp.mpf("0.05")
                var_convergence = abs(
                    perturbed_results["md_z_variance"] - TARGET_VARIANCE
                ) < mp.mpf("0.05")
                convergence_preserved.append(phi_convergence and var_convergence)

                # Restore original parameters for next iteration
                for base in self.md_params:
                    self.md_params[base] = original_params[base].copy()

            except Exception as e:
                logger.warning(f"MD perturbation {i} failed: {e}")
                continue

        # Calculate statistics
        if perturbed_means:
            mean_stability = sum(perturbed_means) / mp.mpf(len(perturbed_means))
            variance_stability = sum(perturbed_variances) / mp.mpf(
                len(perturbed_variances)
            )
            convergence_rate = sum(convergence_preserved) / len(convergence_preserved)

            # Deviation from baseline
            mean_deviation = abs(mean_stability - baseline_mean)
            variance_deviation = abs(variance_stability - baseline_variance)

            results = {
                "baseline_md_mean": baseline_mean,
                "baseline_md_variance": baseline_variance,
                "perturbed_mean_stability": mean_stability,
                "perturbed_variance_stability": variance_stability,
                "mean_deviation": mean_deviation,
                "variance_deviation": variance_deviation,
                "convergence_rate": convergence_rate,
                "num_successful_tests": len(perturbed_means),
                "num_total_tests": num_tests,
                "md_falsification_passed": convergence_rate
                > 0.6,  # 60% threshold for MD robustness
            }
        else:
            results = {
                "error": "No successful perturbation tests",
                "md_falsification_passed": False,
            }

        logger.info(
            f"MD falsification test completed: {len(perturbed_means)}/{num_tests} successful"
        )
        if perturbed_means:
            logger.info(f"MD convergence rate: {convergence_rate:.3f}")

        return results


def compare_md_vs_base_framework(
    sequences: List[str], md_weight: float = 0.3
) -> Dict[str, any]:
    """
    Compare MD-enhanced Z Framework against base Z Framework.

    Args:
        sequences: List of DNA sequences to test
        md_weight: Weight factor for MD contribution

    Returns:
        Comparative analysis results
    """
    base_calc = ZFrameworkCalculator(precision_dps=30)
    md_calc = MolecularDynamicsZFramework(precision_dps=30, md_weight=md_weight)

    comparison_results = {
        "sequences_tested": len(sequences),
        "md_weight": md_weight,
        "results": [],
    }

    for seq in sequences:
        if len(seq) < 5:  # Skip very short sequences
            continue

        try:
            base_results = base_calc.calculate_z_values(seq)
            md_results = md_calc.calculate_md_z_values(seq)

            seq_comparison = {
                "sequence": seq,
                "length": len(seq),
                "base_z_mean": float(base_results["z_mean"]),
                "md_z_mean": float(md_results["md_z_mean"]),
                "base_z_variance": float(base_results["z_variance"]),
                "md_z_variance": float(md_results["md_z_variance"]),
                "enhancement_factor": float(md_results["md_enhancement_factor"]),
                "base_phi_convergence": float(
                    base_results["phi_conjugate_convergence"]
                ),
                "md_phi_convergence": float(md_results["md_phi_conjugate_convergence"]),
                "md_improves_phi_convergence": md_results[
                    "md_phi_conjugate_convergence"
                ]
                < base_results["phi_conjugate_convergence"],
                "md_improves_variance_convergence": md_results[
                    "md_variance_convergence"
                ]
                < base_results["variance_convergence"],
            }

            comparison_results["results"].append(seq_comparison)

        except Exception as e:
            logger.warning(f"Comparison failed for sequence {seq}: {e}")
            continue

    # Calculate aggregate statistics
    if comparison_results["results"]:
        results = comparison_results["results"]
        comparison_results["summary"] = {
            "avg_enhancement_factor": np.mean(
                [r["enhancement_factor"] for r in results]
            ),
            "sequences_with_improved_phi_convergence": sum(
                r["md_improves_phi_convergence"] for r in results
            ),
            "sequences_with_improved_variance_convergence": sum(
                r["md_improves_variance_convergence"] for r in results
            ),
            "improvement_rate_phi": sum(
                r["md_improves_phi_convergence"] for r in results
            )
            / len(results),
            "improvement_rate_variance": sum(
                r["md_improves_variance_convergence"] for r in results
            )
            / len(results),
        }

    return comparison_results
