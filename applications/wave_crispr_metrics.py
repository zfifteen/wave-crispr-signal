"""
Wave CRISPR Metrics Module

This module provides comprehensive metrics and analysis tools for CRISPR guides
using signal-theoretic approaches combined with traditional bioinformatics metrics.
"""

import numpy as np
from scipy.signal import find_peaks
from typing import List, Dict, Optional

try:
    from .crispr_guide_designer import CRISPRGuideDesigner
except ImportError:
    # Handle relative import for direct execution
    import sys

    sys.path.append(".")
    from crispr_guide_designer import CRISPRGuideDesigner


class WaveCRISPRMetrics:
    """Advanced metrics calculator for CRISPR guides using signal theory."""

    def __init__(self, designer: Optional[CRISPRGuideDesigner] = None):
        """Initialize metrics calculator."""
        self.designer = designer or CRISPRGuideDesigner()

    def calculate_spectral_entropy(self, sequence: str, base: int = 2) -> float:
        """
        Calculate spectral entropy of DNA sequence.

        Args:
            sequence: DNA sequence
            base: Logarithm base for entropy calculation

        Returns:
            Spectral entropy value
        """
        wave = self.designer.build_waveform(sequence)
        spectrum = self.designer.compute_spectrum(wave)

        # Normalize spectrum
        ps = spectrum / np.sum(spectrum)
        ps = ps[ps > 0]  # Remove zeros to avoid log(0)

        if len(ps) == 0:
            return 0.0

        return -np.sum(ps * np.log(ps) / np.log(base))

    def calculate_spectral_complexity(self, sequence: str) -> Dict[str, float]:
        """
        Calculate multiple spectral complexity measures.

        Args:
            sequence: DNA sequence

        Returns:
            Dictionary with complexity metrics
        """
        wave = self.designer.build_waveform(sequence)
        spectrum = self.designer.compute_spectrum(wave)

        # Spectral entropy
        entropy = self.calculate_spectral_entropy(sequence)

        # Spectral flatness (Wiener entropy)
        geometric_mean = np.exp(np.mean(np.log(spectrum + 1e-10)))
        arithmetic_mean = np.mean(spectrum)
        spectral_flatness = (
            geometric_mean / arithmetic_mean if arithmetic_mean > 0 else 0
        )

        # Spectral centroid
        freqs = np.arange(len(spectrum))
        spectral_centroid = (
            np.sum(freqs * spectrum) / np.sum(spectrum) if np.sum(spectrum) > 0 else 0
        )

        # Spectral rolloff (95% of energy)
        cumsum = np.cumsum(spectrum)
        total_energy = cumsum[-1]
        rolloff_idx = np.where(cumsum >= 0.95 * total_energy)[0]
        spectral_rolloff = (
            rolloff_idx[0] / len(spectrum) if len(rolloff_idx) > 0 else 1.0
        )

        # Zero crossing rate of spectrum
        zcr = np.sum(np.diff(np.sign(spectrum - np.mean(spectrum))) != 0) / len(
            spectrum
        )

        return {
            "spectral_entropy": entropy,
            "spectral_flatness": spectral_flatness,
            "spectral_centroid": spectral_centroid,
            "spectral_rolloff": spectral_rolloff,
            "zero_crossing_rate": zcr,
            "peak_count": int(self.designer.count_sidelobes(spectrum)),
        }

    def calculate_harmonic_content(self, sequence: str) -> Dict[str, float]:
        """
        Analyze harmonic content of DNA sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Dictionary with harmonic analysis
        """
        wave = self.designer.build_waveform(sequence)
        spectrum = self.designer.compute_spectrum(wave)

        # Find peaks in spectrum
        peaks, properties = find_peaks(spectrum, height=np.max(spectrum) * 0.1)

        # Fundamental frequency (strongest peak)
        fundamental_idx = peaks[np.argmax(spectrum[peaks])] if len(peaks) > 0 else 0
        fundamental_power = (
            spectrum[fundamental_idx] if fundamental_idx < len(spectrum) else 0
        )

        # Harmonic series strength
        harmonics = []
        for h in range(2, 6):  # Check harmonics 2-5
            harmonic_idx = min(fundamental_idx * h, len(spectrum) - 1)
            harmonics.append(spectrum[harmonic_idx])

        harmonic_power = np.sum(harmonics)
        total_power = np.sum(spectrum)

        # Harmonic-to-noise ratio
        noise_power = total_power - fundamental_power - harmonic_power
        hnr = fundamental_power / noise_power if noise_power > 0 else np.inf

        return {
            "fundamental_frequency": fundamental_idx / len(spectrum),
            "fundamental_power": fundamental_power,
            "harmonic_power": harmonic_power,
            "total_harmonic_distortion": (
                harmonic_power / fundamental_power if fundamental_power > 0 else 0
            ),
            "harmonic_to_noise_ratio": hnr,
            "num_peaks": len(peaks),
        }

    def calculate_mutational_sensitivity(
        self, sequence: str, positions: Optional[List[int]] = None
    ) -> Dict:
        """
        Calculate sensitivity to mutations using spectral disruption.

        Args:
            sequence: DNA sequence
            positions: Specific positions to analyze (default: all)

        Returns:
            Dictionary with sensitivity metrics
        """
        if positions is None:
            positions = list(range(len(sequence)))

        base_wave = self.designer.build_waveform(sequence)
        base_spec = self.designer.compute_spectrum(base_wave)
        base_entropy = self.calculate_spectral_entropy(sequence)

        sensitivity_scores = []
        position_effects = {}

        for pos in positions:
            if pos >= len(sequence):
                continue

            pos_scores = []
            original_base = sequence[pos]

            for new_base in "ATCG":
                if new_base != original_base:
                    # Create mutated sequence
                    mutated = list(sequence)
                    mutated[pos] = new_base
                    mutated_seq = "".join(mutated)

                    # Calculate spectral disruption
                    mut_wave = self.designer.build_waveform(mutated_seq)
                    mut_spec = self.designer.compute_spectrum(mut_wave)
                    mut_entropy = self.calculate_spectral_entropy(mutated_seq)

                    # Spectral distance
                    min_len = min(len(base_spec), len(mut_spec))
                    spec_distance = np.linalg.norm(
                        base_spec[:min_len] - mut_spec[:min_len]
                    )

                    # Entropy change
                    entropy_change = abs(mut_entropy - base_entropy)

                    # Combined disruption score
                    disruption = spec_distance + entropy_change * 10
                    pos_scores.append(disruption)

            avg_sensitivity = np.mean(pos_scores)
            sensitivity_scores.append(avg_sensitivity)
            position_effects[pos] = {
                "average_disruption": avg_sensitivity,
                "max_disruption": np.max(pos_scores),
                "min_disruption": np.min(pos_scores),
                "original_base": original_base,
            }

        return {
            "average_sensitivity": np.mean(sensitivity_scores),
            "max_sensitivity": np.max(sensitivity_scores),
            "sensitivity_variance": np.var(sensitivity_scores),
            "hotspot_positions": sorted(
                position_effects.keys(),
                key=lambda x: position_effects[x]["average_disruption"],
                reverse=True,
            )[:5],
            "position_effects": position_effects,
        }

    def calculate_off_target_metrics(
        self, guide_seq: str, target_seq: str, off_targets: List[str]
    ) -> Dict:
        """
        Calculate comprehensive off-target risk metrics.

        Args:
            guide_seq: Guide RNA sequence
            target_seq: On-target sequence
            off_targets: List of potential off-target sequences

        Returns:
            Dictionary with off-target metrics
        """
        target_wave = self.designer.build_waveform(target_seq)
        target_spec = self.designer.compute_spectrum(target_wave)
        # target_entropy = self.calculate_spectral_entropy(target_seq)  # Not used currently

        off_target_risks = []
        spectral_similarities = []
        sequence_similarities = []

        for off_seq in off_targets:
            # Spectral similarity
            off_wave = self.designer.build_waveform(off_seq)
            off_spec = self.designer.compute_spectrum(off_wave)

            min_len = min(len(target_spec), len(off_spec))
            if min_len > 0:
                correlation = np.corrcoef(target_spec[:min_len], off_spec[:min_len])[
                    0, 1
                ]
                correlation = correlation if not np.isnan(correlation) else 0
            else:
                correlation = 0

            spectral_similarities.append(correlation)

            # Sequence similarity (normalized Hamming distance)
            min_seq_len = min(len(target_seq), len(off_seq))
            if min_seq_len > 0:
                seq_sim = (
                    sum(
                        a == b
                        for a, b in zip(target_seq[:min_seq_len], off_seq[:min_seq_len])
                    )
                    / min_seq_len
                )
            else:
                seq_sim = 0

            sequence_similarities.append(seq_sim)

            # Combined risk score
            risk = self.designer.calculate_off_target_risk(
                guide_seq, target_seq, off_seq
            )
            off_target_risks.append(risk)

        return {
            "num_off_targets": len(off_targets),
            "max_risk": np.max(off_target_risks) if off_target_risks else 0,
            "mean_risk": np.mean(off_target_risks) if off_target_risks else 0,
            "high_risk_count": sum(1 for r in off_target_risks if r > 0.7),
            "spectral_similarity_mean": (
                np.mean(spectral_similarities) if spectral_similarities else 0
            ),
            "sequence_similarity_mean": (
                np.mean(sequence_similarities) if sequence_similarities else 0
            ),
            "risk_variance": np.var(off_target_risks) if off_target_risks else 0,
        }

    def calculate_repair_prediction_confidence(
        self, guide_seq: str, target_context: str
    ) -> Dict:
        """
        Calculate confidence metrics for repair outcome predictions.

        Args:
            guide_seq: Guide RNA sequence
            target_context: Target sequence context around cut site

        Returns:
            Dictionary with confidence metrics
        """
        # Get repair predictions
        repair_outcomes = self.designer.predict_repair_outcomes(
            guide_seq, target_context
        )

        # Calculate sequence features that affect repair
        cut_site = len(target_context) // 2

        # Microhomology analysis (affects MMEJ)
        microhomology_count = 0
        for length in range(2, 8):  # 2-7 bp microhomologies
            for i in range(max(0, cut_site - 20), cut_site):
                for j in range(cut_site + 3, min(len(target_context), cut_site + 20)):
                    if i + length <= len(target_context) and j + length <= len(
                        target_context
                    ):
                        if (
                            target_context[i : i + length]
                            == target_context[j : j + length]
                        ):
                            microhomology_count += 1

        # GC content around cut site
        context_window = target_context[
            max(0, cut_site - 10) : min(len(target_context), cut_site + 10)
        ]
        gc_content = (context_window.count("G") + context_window.count("C")) / len(
            context_window
        )

        # Repetitive elements
        repeat_score = 0
        for length in range(3, 8):
            for i in range(len(target_context) - length):
                motif = target_context[i : i + length]
                count = target_context.count(motif)
                if count > 1:
                    repeat_score += count - 1

        # Prediction confidence based on feature reliability
        confidence_factors = {
            "microhomology_reliability": min(
                1.0, microhomology_count / 5.0
            ),  # More microhomology = higher MMEJ confidence
            "gc_reliability": 1.0
            - abs(gc_content - 0.5) * 2,  # Moderate GC = higher reliability
            "repeat_reliability": max(
                0.2, 1.0 - repeat_score / 10.0
            ),  # Fewer repeats = higher reliability
            "spectral_reliability": repair_outcomes[
                "stability_score"
            ],  # Higher stability = higher reliability
        }

        overall_confidence = np.mean(list(confidence_factors.values()))

        return {
            "overall_confidence": overall_confidence,
            "confidence_factors": confidence_factors,
            "microhomology_count": microhomology_count,
            "gc_content": gc_content,
            "repeat_score": repeat_score,
            "predicted_dominant_pathway": max(
                repair_outcomes,
                key=lambda k: (
                    repair_outcomes[k]
                    if k.endswith("_probability") or k.endswith("_efficiency")
                    else 0
                ),
            ),
        }

    def calculate_comprehensive_score(
        self, guide_seq: str, target_seq: str, off_targets: Optional[List[str]] = None
    ) -> Dict:
        """
        Calculate comprehensive guide quality score combining all metrics.

        Args:
            guide_seq: Guide RNA sequence
            target_seq: Target sequence
            off_targets: Optional list of off-target sequences

        Returns:
            Dictionary with comprehensive scoring
        """
        # Basic scores
        on_target_score = self.designer.calculate_on_target_score(guide_seq)

        # Spectral complexity
        complexity = self.calculate_spectral_complexity(guide_seq)

        # Harmonic content
        harmonics = self.calculate_harmonic_content(guide_seq)

        # Mutational sensitivity
        sensitivity = self.calculate_mutational_sensitivity(guide_seq)

        # Off-target analysis
        if off_targets:
            off_target_metrics = self.calculate_off_target_metrics(
                guide_seq, target_seq, off_targets
            )
        else:
            off_target_metrics = {"max_risk": 0, "mean_risk": 0}

        # Repair prediction confidence
        repair_confidence = self.calculate_repair_prediction_confidence(
            guide_seq, target_seq
        )

        # Calculate composite scores
        spectral_score = (
            complexity["spectral_entropy"] / 10.0 * 0.3
            + complexity["spectral_flatness"] * 0.2
            + (1.0 - complexity["spectral_rolloff"]) * 0.2
            + harmonics["harmonic_to_noise_ratio"] / 10.0 * 0.3
        )

        stability_score = 1.0 / (1.0 + sensitivity["average_sensitivity"])

        off_target_score = 1.0 - off_target_metrics["max_risk"]

        # Weighted composite score
        weights = {
            "on_target": 0.35,
            "spectral": 0.25,
            "stability": 0.20,
            "off_target": 0.15,
            "repair_confidence": 0.05,
        }

        composite_score = (
            on_target_score * weights["on_target"]
            + np.clip(spectral_score, 0, 1) * weights["spectral"]
            + stability_score * weights["stability"]
            + off_target_score * weights["off_target"]
            + repair_confidence["overall_confidence"] * weights["repair_confidence"]
        )

        return {
            "composite_score": np.clip(composite_score, 0, 1),
            "component_scores": {
                "on_target": on_target_score,
                "spectral": spectral_score,
                "stability": stability_score,
                "off_target": off_target_score,
                "repair_confidence": repair_confidence["overall_confidence"],
            },
            "detailed_metrics": {
                "complexity": complexity,
                "harmonics": harmonics,
                "sensitivity": sensitivity,
                "off_target": off_target_metrics,
                "repair_confidence": repair_confidence,
            },
            "weights": weights,
        }

    def benchmark_guide_set(self, guides: List[Dict], target_seq: str) -> Dict:
        """
        Benchmark a set of guides and provide comparative analysis.

        Args:
            guides: List of guide dictionaries
            target_seq: Target sequence

        Returns:
            Dictionary with benchmark results
        """
        scores = []
        detailed_results = []

        for guide in guides:
            guide_seq = guide["sequence"]

            # Calculate comprehensive score
            result = self.calculate_comprehensive_score(guide_seq, target_seq)
            scores.append(result["composite_score"])

            detailed_results.append({"guide": guide, "comprehensive_score": result})

        # Statistical analysis
        scores_array = np.array(scores)

        return {
            "num_guides": len(guides),
            "score_statistics": {
                "mean": np.mean(scores_array),
                "std": np.std(scores_array),
                "min": np.min(scores_array),
                "max": np.max(scores_array),
                "median": np.median(scores_array),
                "q25": np.percentile(scores_array, 25),
                "q75": np.percentile(scores_array, 75),
            },
            "top_guides": sorted(
                detailed_results,
                key=lambda x: x["comprehensive_score"]["composite_score"],
                reverse=True,
            )[:5],
            "quality_distribution": {
                "excellent": sum(1 for s in scores if s > 0.8),
                "good": sum(1 for s in scores if 0.6 < s <= 0.8),
                "fair": sum(1 for s in scores if 0.4 < s <= 0.6),
                "poor": sum(1 for s in scores if s <= 0.4),
            },
            "detailed_results": detailed_results,
        }


def main():
    """Example usage of wave CRISPR metrics."""
    # Example sequence
    target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

    # Design guides
    designer = CRISPRGuideDesigner()
    guides = designer.design_guides(target_sequence, num_guides=5)

    # Initialize metrics calculator
    metrics = WaveCRISPRMetrics(designer)

    print("📊 Wave CRISPR Metrics Analysis")
    print("=" * 50)

    # Analyze first guide in detail
    if guides:
        guide = guides[0]
        guide_seq = guide["sequence"]

        print(f"\n🧬 Detailed Analysis: {guide_seq}")
        print("-" * 50)

        # Spectral complexity
        complexity = metrics.calculate_spectral_complexity(guide_seq)
        print(f"Spectral Entropy: {complexity['spectral_entropy']:.3f}")
        print(f"Spectral Flatness: {complexity['spectral_flatness']:.3f}")
        print(f"Spectral Centroid: {complexity['spectral_centroid']:.3f}")

        # Harmonic content
        harmonics = metrics.calculate_harmonic_content(guide_seq)
        print(f"Fundamental Power: {harmonics['fundamental_power']:.3f}")
        print(f"Harmonic Power: {harmonics['harmonic_power']:.3f}")
        print(f"THD: {harmonics['total_harmonic_distortion']:.3f}")

        # Mutational sensitivity
        sensitivity = metrics.calculate_mutational_sensitivity(guide_seq)
        print(f"Average Sensitivity: {sensitivity['average_sensitivity']:.3f}")
        print(f"Hotspot Positions: {sensitivity['hotspot_positions'][:3]}")

        # Comprehensive score
        comprehensive = metrics.calculate_comprehensive_score(
            guide_seq, target_sequence
        )
        print(f"\n🎯 Comprehensive Score: {comprehensive['composite_score']:.3f}")

        print("\nComponent Scores:")
        for component, score in comprehensive["component_scores"].items():
            print(f"  {component}: {score:.3f}")

    # Benchmark all guides
    print(f"\n📈 Benchmarking {len(guides)} guides...")
    benchmark = metrics.benchmark_guide_set(guides, target_sequence)

    print("Score Statistics:")
    stats = benchmark["score_statistics"]
    print(f"  Mean: {stats['mean']:.3f} ± {stats['std']:.3f}")
    print(f"  Range: {stats['min']:.3f} - {stats['max']:.3f}")
    print(f"  Median: {stats['median']:.3f}")

    print("\nQuality Distribution:")
    dist = benchmark["quality_distribution"]
    for quality, count in dist.items():
        print(f"  {quality.capitalize()}: {count} guides")

    print("\n🏆 Top 3 Guides:")
    for i, result in enumerate(benchmark["top_guides"][:3], 1):
        guide = result["guide"]
        score = result["comprehensive_score"]["composite_score"]
        print(f"  {i}. {guide['sequence']} (Score: {score:.3f})")


if __name__ == "__main__":
    main()
