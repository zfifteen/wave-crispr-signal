"""
Validation Module for Invariant Features

Implements bootstrap validation and phase-stability metrics for
the invariant CRISPR features as outlined in the evaluation protocol.
"""

import numpy as np
from typing import Dict, List
import random
from scipy import stats
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from invariant_features import InvariantFeatureSet
from applications.crispr_guide_designer import CRISPRGuideDesigner


class InvariantValidator:
    """
    Validator for invariant features with bootstrap and phase-stability analysis.
    """

    def __init__(self, seed: int = 42):
        """
        Initialize validator.

        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        random.seed(seed)
        np.random.seed(seed)

        self.invariant_features = InvariantFeatureSet()
        self.crispr_designer = CRISPRGuideDesigner()

    def calculate_phase_stability_index(
        self, sequences: List[str], num_bootstrap: int = 100
    ) -> Dict[str, float]:
        """
        Calculate phase-stability index using bootstrap resampling.

        Phase-stability index = 1 - sd(Δ_phase) / mean|Δ_phase|
        Higher values indicate more biologically real phase effects.

        Args:
            sequences: List of DNA sequences for analysis
            num_bootstrap: Number of bootstrap resamples

        Returns:
            Dictionary with phase stability metrics
        """
        delta_phase_entropy_samples = []
        delta_phase_flatness_samples = []
        delta_phase_f1_samples = []

        for _ in range(num_bootstrap):
            # Bootstrap resample sequences
            resampled_sequences = np.random.choice(
                sequences, size=len(sequences), replace=True
            )

            # Calculate phase differences for this resample
            entropy_deltas = []
            flatness_deltas = []
            f1_deltas = []

            for seq in resampled_sequences:
                try:
                    features = self.invariant_features.calculate_complete_feature_set(
                        seq
                    )
                    entropy_deltas.append(features.get("delta_phase_entropy", 0))
                    flatness_deltas.append(features.get("delta_phase_flatness", 0))
                    f1_deltas.append(features.get("delta_phase_f1_magnitude", 0))
                except Exception:
                    continue  # Skip problematic sequences

            if entropy_deltas:
                delta_phase_entropy_samples.append(np.mean(entropy_deltas))
                delta_phase_flatness_samples.append(np.mean(flatness_deltas))
                delta_phase_f1_samples.append(np.mean(f1_deltas))

        # Calculate stability indices
        def calc_stability_index(values):
            values = np.array(values)
            if len(values) == 0 or np.all(values == 0):
                return 0.0

            mean_abs = np.mean(np.abs(values))
            std_val = np.std(values)

            if mean_abs == 0:
                return 0.0

            stability = 1.0 - (std_val / mean_abs)
            return max(0.0, stability)  # Ensure non-negative

        results = {
            "phase_stability_entropy": calc_stability_index(
                delta_phase_entropy_samples
            ),
            "phase_stability_flatness": calc_stability_index(
                delta_phase_flatness_samples
            ),
            "phase_stability_f1": calc_stability_index(delta_phase_f1_samples),
            "num_bootstrap_samples": len(delta_phase_entropy_samples),
            "num_sequences": len(sequences),
        }

        # Overall phase stability (mean of individual metrics)
        individual_stabilities = [
            results["phase_stability_entropy"],
            results["phase_stability_flatness"],
            results["phase_stability_f1"],
        ]
        results["overall_phase_stability"] = np.mean(individual_stabilities)

        return results

    def stratify_mutation_analysis(
        self, sequences: List[str], confidence_level: float = 0.95
    ) -> Dict[str, Dict]:
        """
        Stratify analysis by mutation class (G→C vs others) with confidence intervals.

        Args:
            sequences: List of DNA sequences
            confidence_level: Confidence level for intervals

        Returns:
            Dictionary with stratified analysis results
        """
        gc_transitions = []
        other_transitions = []

        for seq in sequences:
            # Find G positions for G→C analysis
            g_positions = [i for i, base in enumerate(seq) if base == "G"]

            if g_positions:
                # Analyze G→C transition
                features = self.invariant_features.calculate_complete_feature_set(
                    seq, g_positions[0], "C"
                )

                gc_data = {
                    "delta_phase_entropy": features.get(
                        "delta_phase_entropy_change", 0
                    ),
                    "delta_phase_flatness": features.get(
                        "delta_phase_flatness_change", 0
                    ),
                    "delta_phi": features.get("delta_phi", 0),
                    "phase_bit": features.get("phase_bit", 0),
                }
                gc_transitions.append(gc_data)

            # Analyze other transitions (A→T as example)
            a_positions = [i for i, base in enumerate(seq) if base == "A"]
            if a_positions:
                features = self.invariant_features.calculate_complete_feature_set(
                    seq, a_positions[0], "T"
                )

                other_data = {
                    "delta_phase_entropy": features.get(
                        "delta_phase_entropy_change", 0
                    ),
                    "delta_phase_flatness": features.get(
                        "delta_phase_flatness_change", 0
                    ),
                    "delta_phi": features.get("delta_phi", 0),
                    "phase_bit": features.get("phase_bit", 0),
                }
                other_transitions.append(other_data)

        # Calculate confidence intervals
        1.0 - confidence_level

        def calc_ci(data_list, metric):
            if not data_list:
                return {"mean": 0, "ci_lower": 0, "ci_upper": 0, "n": 0}

            values = [d[metric] for d in data_list if metric in d]
            if not values:
                return {"mean": 0, "ci_lower": 0, "ci_upper": 0, "n": 0}

            mean_val = np.mean(values)
            sem = stats.sem(values)
            h = sem * stats.t.ppf((1 + confidence_level) / 2.0, len(values) - 1)

            return {
                "mean": mean_val,
                "ci_lower": mean_val - h,
                "ci_upper": mean_val + h,
                "n": len(values),
                "std": np.std(values),
            }

        results = {
            "gc_transitions": {
                "delta_phase_entropy": calc_ci(gc_transitions, "delta_phase_entropy"),
                "delta_phase_flatness": calc_ci(gc_transitions, "delta_phase_flatness"),
                "delta_phi": calc_ci(gc_transitions, "delta_phi"),
                "phase_bit_distribution": [d["phase_bit"] for d in gc_transitions],
            },
            "other_transitions": {
                "delta_phase_entropy": calc_ci(
                    other_transitions, "delta_phase_entropy"
                ),
                "delta_phase_flatness": calc_ci(
                    other_transitions, "delta_phase_flatness"
                ),
                "delta_phi": calc_ci(other_transitions, "delta_phi"),
                "phase_bit_distribution": [d["phase_bit"] for d in other_transitions],
            },
            "confidence_level": confidence_level,
        }

        # Statistical tests comparing G→C vs others
        if gc_transitions and other_transitions:
            gc_entropy = [d["delta_phase_entropy"] for d in gc_transitions]
            other_entropy = [d["delta_phase_entropy"] for d in other_transitions]

            if len(gc_entropy) > 1 and len(other_entropy) > 1:
                # Welch's t-test (unequal variances)
                t_stat, p_value = stats.ttest_ind(
                    gc_entropy, other_entropy, equal_var=False
                )
                results["statistical_comparison"] = {
                    "t_statistic": t_stat,
                    "p_value": p_value,
                    "significant": p_value < 0.05,
                }

        return results

    def evaluate_guide_performance_lift(
        self,
        sequences: List[str],
        use_invariants: bool = True,
        baseline_scoring: bool = False,
    ) -> Dict[str, float]:
        """
        Evaluate performance lift from invariant features in guide scoring.

        Args:
            sequences: List of target sequences for guide design
            use_invariants: Whether to use invariant features
            baseline_scoring: Whether to compare against baseline

        Returns:
            Dictionary with performance metrics
        """
        invariant_scores = []
        baseline_scores = []

        for seq in sequences:
            try:
                # Get guides with invariant features
                guides_inv = self.crispr_designer.design_guides(
                    seq, num_guides=3, use_invariants=True
                )

                # Get guides with baseline scoring
                guides_base = self.crispr_designer.design_guides(
                    seq, num_guides=3, use_invariants=False
                )

                if guides_inv and guides_base:
                    # Compare top guide scores
                    best_inv_score = guides_inv[0]["comprehensive_score"]
                    best_base_score = guides_base[0]["on_target_score"]

                    invariant_scores.append(best_inv_score)
                    baseline_scores.append(best_base_score)

            except Exception:
                continue  # Skip problematic sequences

        if not invariant_scores or not baseline_scores:
            return {"error": "No valid guides found", "n_sequences": 0}

        # Calculate performance metrics
        mean_inv = np.mean(invariant_scores)
        mean_base = np.mean(baseline_scores)

        # Performance lift calculation
        relative_improvement = (
            (mean_inv - mean_base) / mean_base if mean_base > 0 else 0
        )

        # Statistical significance test
        t_stat, p_value = stats.ttest_rel(invariant_scores, baseline_scores)

        results = {
            "mean_invariant_score": mean_inv,
            "mean_baseline_score": mean_base,
            "relative_improvement": relative_improvement,
            "improvement_percentage": relative_improvement * 100,
            "t_statistic": t_stat,
            "p_value": p_value,
            "significant_improvement": p_value < 0.05 and mean_inv > mean_base,
            "n_sequences": len(invariant_scores),
            "effect_size": (mean_inv - mean_base)
            / np.sqrt((np.var(invariant_scores) + np.var(baseline_scores)) / 2),
        }

        return results

    def comprehensive_validation_report(
        self, sequences: List[str], num_bootstrap: int = 50
    ) -> Dict[str, any]:
        """
        Generate comprehensive validation report for invariant features.

        Args:
            sequences: List of DNA sequences for validation
            num_bootstrap: Number of bootstrap samples

        Returns:
            Comprehensive validation report
        """
        print("🔬 Running Comprehensive Invariant Feature Validation...")
        print("=" * 60)

        report = {
            "metadata": {
                "num_sequences": len(sequences),
                "num_bootstrap_samples": num_bootstrap,
                "seed": self.seed,
            }
        }

        # 1. Phase stability analysis
        print("1. Calculating phase stability indices...")
        phase_stability = self.calculate_phase_stability_index(sequences, num_bootstrap)
        report["phase_stability"] = phase_stability

        # 2. Mutation stratification analysis
        print("2. Performing mutation stratification analysis...")
        mutation_analysis = self.stratify_mutation_analysis(sequences)
        report["mutation_stratification"] = mutation_analysis

        # 3. Guide performance evaluation
        print("3. Evaluating guide performance lift...")
        performance_lift = self.evaluate_guide_performance_lift(sequences)
        report["performance_evaluation"] = performance_lift

        # 4. Summary and recommendations
        print("4. Generating summary and recommendations...")
        summary = self._generate_summary(report)
        report["summary"] = summary

        return report

    def _generate_summary(self, report: Dict) -> Dict[str, any]:
        """Generate summary and recommendations from validation results."""
        summary = {}

        # Phase stability assessment
        phase_stab = report["phase_stability"]["overall_phase_stability"]
        if phase_stab > 0.7:
            phase_quality = "High"
        elif phase_stab > 0.4:
            phase_quality = "Moderate"
        else:
            phase_quality = "Low"

        summary["phase_stability_quality"] = phase_quality
        summary["phase_stability_score"] = phase_stab

        # Performance improvement assessment
        perf_data = report["performance_evaluation"]
        if "improvement_percentage" in perf_data:
            improvement_pct = perf_data["improvement_percentage"]
            is_significant = perf_data.get("significant_improvement", False)

            if improvement_pct > 10 and is_significant:
                performance_assessment = "Significant improvement"
            elif improvement_pct > 5:
                performance_assessment = "Moderate improvement"
            elif improvement_pct > 0:
                performance_assessment = "Slight improvement"
            else:
                performance_assessment = "No improvement"

            summary["performance_assessment"] = performance_assessment
            summary["improvement_percentage"] = improvement_pct
            summary["statistically_significant"] = is_significant

        # Overall recommendation
        if phase_quality == "High" and summary.get("statistically_significant", False):
            recommendation = "Strongly recommend using invariant features"
        elif (
            phase_quality in ["High", "Moderate"]
            or summary.get("improvement_percentage", 0) > 5
        ):
            recommendation = "Recommend using invariant features with validation"
        else:
            recommendation = "Further validation needed before deployment"

        summary["recommendation"] = recommendation

        return summary


def generate_test_sequences(num_sequences: int = 20, seq_length: int = 50) -> List[str]:
    """
    Generate test sequences for validation.

    Args:
        num_sequences: Number of sequences to generate
        seq_length: Length of each sequence

    Returns:
        List of test DNA sequences
    """
    bases = ["A", "T", "C", "G"]
    sequences = []

    for _ in range(num_sequences):
        seq = "".join(random.choices(bases, k=seq_length))
        # Ensure sequence has PAM sites
        if "GG" not in seq:
            # Insert a PAM site
            pos = random.randint(3, len(seq) - 3)
            seq = seq[:pos] + "AGG" + seq[pos + 3 :]
        sequences.append(seq)

    return sequences


if __name__ == "__main__":
    # Example validation run
    validator = InvariantValidator()
    test_sequences = generate_test_sequences(15, 60)

    # Run comprehensive validation
    report = validator.comprehensive_validation_report(test_sequences, num_bootstrap=30)

    # Print summary
    print("\n📊 VALIDATION SUMMARY")
    print("=" * 40)
    for key, value in report["summary"].items():
        print(f"{key.replace('_', ' ').title()}: {value}")
