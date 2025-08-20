"""
CRISPR Guide Designer using Signal-Theoretic DNA Analysis

This module implements CRISPR guide RNA design using spectral disruption profiling
based on the wave-crispr-signal methodology with enhanced invariant features.
"""

import numpy as np
from scipy.fft import fft
from scipy.stats import entropy
import re
import sys
import os
from typing import List, Dict, Tuple, Optional

# Add parent directory for invariant features import
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.invariant_features import InvariantFeatureSet


class CRISPRGuideDesigner:
    """CRISPR guide RNA designer using signal-theoretic analysis with invariant features."""

    # Base weights for wave function mapping
    WEIGHTS = {"A": 1 + 0j, "T": -1 + 0j, "C": 0 + 1j, "G": 0 - 1j}

    def __init__(self, pam_pattern: str = "NGG", guide_length: int = 20):
        """
        Initialize CRISPR guide designer.

        Args:
            pam_pattern: PAM sequence pattern (default: NGG for Cas9)
            guide_length: Length of guide RNA (default: 20)
        """
        self.pam_pattern = pam_pattern.replace("N", "[ATCG]")
        self.guide_length = guide_length

        # Initialize invariant feature calculator
        self.invariant_features = InvariantFeatureSet(pam_pattern)

    def build_waveform(
        self, seq: str, d: float = 0.34, zn_map: Optional[Dict[int, float]] = None
    ) -> np.ndarray:
        """
        Build complex waveform from DNA sequence.

        Args:
            seq: DNA sequence
            d: Base spacing parameter (default: 0.34 nm)
            zn_map: Position-dependent scaling factors

        Returns:
            Complex waveform array
        """
        if zn_map is None:
            spacings = [d] * len(seq)
        else:
            spacings = [d * (1 + zn_map.get(i, 0)) for i in range(len(seq))]

        s = np.cumsum(spacings)
        wave = [
            self.WEIGHTS[base] * np.exp(2j * np.pi * s[i]) for i, base in enumerate(seq)
        ]
        return np.array(wave)

    def compute_spectrum(self, waveform: np.ndarray) -> np.ndarray:
        """Compute frequency spectrum of waveform."""
        return np.abs(fft(waveform))

    def normalized_entropy(self, spectrum: np.ndarray) -> float:
        """Calculate normalized entropy of spectrum."""
        ps = spectrum / np.sum(spectrum)
        return entropy(ps, base=2)

    def count_sidelobes(
        self, spectrum: np.ndarray, threshold_ratio: float = 0.25
    ) -> int:
        """Count spectral sidelobes above threshold."""
        peak = np.max(spectrum)
        return np.sum(spectrum > (threshold_ratio * peak))

    def calculate_on_target_score(self, guide_seq: str) -> float:
        """
        Calculate on-target efficiency score using spectral features.

        Args:
            guide_seq: Guide RNA sequence

        Returns:
            On-target efficiency score (0-1)
        """
        wave = self.build_waveform(guide_seq)
        spec = self.compute_spectrum(wave)

        # Spectral features for on-target prediction
        entropy_score = 1.0 - (
            self.normalized_entropy(spec) / 10.0
        )  # Lower entropy = higher efficiency
        sidelobe_score = 1.0 - (
            self.count_sidelobes(spec) / len(spec)
        )  # Fewer sidelobes = higher efficiency

        # GC content factor (optimal ~40-60%)
        gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq)
        gc_score = 1.0 - abs(gc_content - 0.5) * 2  # Penalize extreme GC content

        # Combined score
        return np.clip(
            (entropy_score * 0.4 + sidelobe_score * 0.4 + gc_score * 0.2), 0, 1
        )

    def calculate_comprehensive_score(
        self, sequence: str, include_invariants: bool = True
    ) -> Dict[str, float]:
        """
        Calculate comprehensive guide score including invariant features.

        Args:
            sequence: Guide sequence
            include_invariants: Whether to include invariant features

        Returns:
            Dictionary with comprehensive scoring metrics
        """
        # Traditional spectral features
        wave = self.build_waveform(sequence)
        spec = self.compute_spectrum(wave)

        entropy_val = self.normalized_entropy(spec)
        sidelobes = self.count_sidelobes(spec)
        gc_content = (sequence.count("G") + sequence.count("C")) / len(sequence)

        # Base scoring
        entropy_score = 1.0 - (entropy_val / 10.0)  # Lower entropy → higher efficiency
        sidelobe_score = 1.0 - (
            sidelobes / len(spec)
        )  # Fewer sidelobes → higher efficiency
        gc_score = 1.0 - abs(gc_content - 0.5) * 2  # Optimal GC ~50%

        base_score = entropy_score * 0.4 + sidelobe_score * 0.4 + gc_score * 0.2

        results = {
            "base_score": base_score,
            "entropy_score": entropy_score,
            "sidelobe_score": sidelobe_score,
            "gc_score": gc_score,
            "spectral_entropy": entropy_val,
            "sidelobe_count": sidelobes,
            "gc_content": gc_content,
        }

        # Add invariant features if requested
        if include_invariants:
            invariant_features = self.invariant_features.calculate_complete_feature_set(
                sequence
            )

            # Extract key invariant metrics for scoring
            phase_bit = invariant_features.get("phase_bit", 0)
            delta_phi = invariant_features.get("delta_phi", 1.0)

            # Phase coherence boost (phase-specific optimization)
            phase_boost = (
                1.1 if phase_bit == 0 else 1.0
            )  # Slight preference for phase 0

            # Golden proximity bonus (lower δφ = more stable)
            golden_bonus = max(
                0, 1.0 - delta_phi * 2
            )  # Bonus for proximity to golden ratio

            # Phase-difference stability (consistent phase differences indicate real effects)
            phase_entropy_diff = abs(invariant_features.get("delta_phase_entropy", 0))
            phase_stability = 1.0 - min(
                phase_entropy_diff / 5.0, 1.0
            )  # Normalize to [0,1]

            # Comprehensive score with invariant features
            invariant_score = (
                base_score * phase_boost + golden_bonus * 0.2 + phase_stability * 0.3
            )
            invariant_score = np.clip(
                invariant_score, 0, 2.0
            )  # Allow scores > 1 for exceptional guides

            results.update(
                {
                    "comprehensive_score": invariant_score,
                    "phase_bit": phase_bit,
                    "delta_phi": delta_phi,
                    "golden_bonus": golden_bonus,
                    "phase_stability": phase_stability,
                    "phase_boost": phase_boost,
                    **invariant_features,  # Include all invariant features
                }
            )
        else:
            results["comprehensive_score"] = base_score

        return results

    def analyze_gc_transition_effects(
        self, sequence: str, window_size: int = 25
    ) -> Dict[str, float]:
        """
        Analyze G→C transition effects using phase-coherent validation.

        Args:
            sequence: DNA sequence to analyze
            window_size: Window size for local analysis

        Returns:
            Dictionary with G→C transition analysis results
        """
        results = {}
        g_positions = [i for i, base in enumerate(sequence) if base == "G"]

        if not g_positions:
            return {"gc_transition_score": 0.0, "num_g_sites": 0}

        transition_scores = []
        phase_coherence_scores = []

        for pos in g_positions:
            # Calculate features for G→C mutation
            features = self.invariant_features.calculate_complete_feature_set(
                sequence, pos, "C"
            )

            # Phase coherence check: real G→C effects should be phase-coherent
            delta_phase_entropy = features.get("delta_phase_entropy_change", 0)
            delta_phase_flatness = features.get("delta_phase_flatness_change", 0)

            # G→C transitions should show consistent signed phase differences
            phase_coherence = (
                1.0 - abs(delta_phase_entropy + delta_phase_flatness) / 10.0
            )
            phase_coherence = max(0, phase_coherence)

            # Curvature disruption analysis
            curv_disruption = features.get("delta_curv_structural_complexity", 0)
            transition_impact = abs(curv_disruption) / 10.0  # Normalize

            transition_scores.append(transition_impact)
            phase_coherence_scores.append(phase_coherence)

        results = {
            "gc_transition_score": (
                np.mean(transition_scores) if transition_scores else 0.0
            ),
            "phase_coherence_score": (
                np.mean(phase_coherence_scores) if phase_coherence_scores else 0.0
            ),
            "num_g_sites": len(g_positions),
            "avg_transition_impact": (
                np.mean(transition_scores) if transition_scores else 0.0
            ),
            "gc_transition_consistency": (
                1.0 - np.std(transition_scores) if len(transition_scores) > 1 else 1.0
            ),
        }

        return results

    def calculate_off_target_risk(
        self, guide_seq: str, target_seq: str, off_target_seq: str
    ) -> float:
        """
        Calculate off-target risk using spectral signature comparison.

        Args:
            guide_seq: Guide RNA sequence
            target_seq: On-target genomic sequence
            off_target_seq: Potential off-target sequence

        Returns:
            Off-target risk score (0-1)
        """
        target_wave = self.build_waveform(target_seq)
        off_target_wave = self.build_waveform(off_target_seq)

        target_spec = self.compute_spectrum(target_wave)
        off_target_spec = self.compute_spectrum(off_target_wave)

        # Spectral similarity using cross-correlation
        min_len = min(len(target_spec), len(off_target_spec))
        corr = np.corrcoef(target_spec[:min_len], off_target_spec[:min_len])[0, 1]

        # Sequence similarity (Hamming distance)
        min_seq_len = min(len(target_seq), len(off_target_seq))
        seq_similarity = (
            sum(
                a == b
                for a, b in zip(target_seq[:min_seq_len], off_target_seq[:min_seq_len])
            )
            / min_seq_len
        )

        # Combined risk score (high correlation + high sequence similarity = high risk)
        spectral_risk = max(0, corr) if not np.isnan(corr) else 0
        return np.clip((spectral_risk * 0.6 + seq_similarity * 0.4), 0, 1)

    def find_pam_sites(self, sequence: str) -> List[Tuple[int, str]]:
        """
        Find PAM sites in sequence.

        Args:
            sequence: DNA sequence to search

        Returns:
            List of (position, PAM_sequence) tuples
        """
        pam_sites = []
        pattern = re.compile(self.pam_pattern)

        for match in pattern.finditer(sequence.upper()):
            pam_sites.append((match.start(), match.group()))

        return pam_sites

    def design_guides(
        self, sequence: str, num_guides: int = 5, use_invariants: bool = True
    ) -> List[Dict]:
        """
        Design CRISPR guides for a target sequence using comprehensive scoring.

        Args:
            sequence: Target DNA sequence
            num_guides: Number of top guides to return
            use_invariants: Whether to use invariant features for scoring

        Returns:
            List of guide dictionaries with scores and positions
        """
        guides = []
        pam_sites = self.find_pam_sites(sequence)

        for pam_pos, pam_seq in pam_sites:
            # Extract guide sequence (upstream of PAM)
            guide_start = max(0, pam_pos - self.guide_length)
            guide_end = pam_pos

            if guide_end - guide_start >= self.guide_length:
                guide_seq = sequence[guide_start:guide_end].upper()

                # Calculate comprehensive scores
                if use_invariants:
                    score_data = self.calculate_comprehensive_score(
                        guide_seq, include_invariants=True
                    )
                    primary_score = score_data["comprehensive_score"]
                else:
                    score_data = self.calculate_comprehensive_score(
                        guide_seq, include_invariants=False
                    )
                    primary_score = score_data["base_score"]

                # G→C transition analysis
                gc_analysis = self.analyze_gc_transition_effects(guide_seq)

                # Basic quality filters
                gc_content = score_data["gc_content"]
                has_poly_t = "TTTT" in guide_seq  # Avoid poly-T stretches

                if 0.2 <= gc_content <= 0.8 and not has_poly_t:
                    guide_data = {
                        "sequence": guide_seq,
                        "position": guide_start,
                        "pam_position": pam_pos,
                        "pam_sequence": pam_seq,
                        "comprehensive_score": primary_score,
                        "on_target_score": score_data[
                            "base_score"
                        ],  # Keep traditional score for comparison
                        "gc_content": gc_content,
                        "length": len(guide_seq),
                        **score_data,  # Include all scoring metrics
                        **gc_analysis,  # Include G→C analysis
                    }
                    guides.append(guide_data)

        # Sort by comprehensive score and return top guides
        sort_key = "comprehensive_score" if use_invariants else "on_target_score"
        guides.sort(key=lambda x: x[sort_key], reverse=True)
        return guides[:num_guides]

    def predict_repair_outcomes(self, guide_seq: str, target_seq: str) -> Dict:
        """
        Predict CRISPR repair outcomes using spectral stability analysis.

        Args:
            guide_seq: Guide RNA sequence
            target_seq: Target genomic sequence around cut site

        Returns:
            Dictionary with repair outcome predictions
        """
        # Find cut site (typically 3 bp upstream of PAM)
        cut_site = len(target_seq) // 2  # Assume cut site is in middle

        # Analyze sequence context around cut site
        context_seq = target_seq[max(0, cut_site - 10) : cut_site + 10]
        wave = self.build_waveform(context_seq)
        spec = self.compute_spectrum(wave)

        entropy_val = self.normalized_entropy(spec)
        sidelobes = self.count_sidelobes(spec)

        # Predict repair pathway bias
        # Low entropy → NHEJ bias
        # High sidelobes → MMEJ bias
        nhej_bias = 1.0 - (entropy_val / 10.0)
        mmej_bias = sidelobes / len(spec)
        hdr_efficiency = max(0, 1.0 - nhej_bias - mmej_bias)

        return {
            "nhej_probability": np.clip(nhej_bias, 0, 1),
            "mmej_probability": np.clip(mmej_bias, 0, 1),
            "hdr_efficiency": np.clip(hdr_efficiency, 0, 1),
            "spectral_entropy": entropy_val,
            "stability_score": 1.0 - entropy_val / 10.0,
        }


def main():
    """Example usage of CRISPR guide designer."""
    # Example target sequence (PCSK9 exon 1)
    target_sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

    designer = CRISPRGuideDesigner()

    print("🧬 CRISPR Guide Designer - Signal-Theoretic Analysis")
    print("=" * 60)

    # Design guides
    guides = designer.design_guides(target_sequence, num_guides=5)

    print(f"Target sequence length: {len(target_sequence)} bp")
    print(f"Found {len(guides)} high-quality guide candidates:\n")

    for i, guide in enumerate(guides, 1):
        print(f"Guide #{i}:")
        print(f"  Sequence: {guide['sequence']}")
        print(f"  Position: {guide['position']}-{guide['position'] + guide['length']}")
        print(f"  PAM: {guide['pam_sequence']} (pos {guide['pam_position']})")
        print(f"  On-target score: {guide['on_target_score']:.3f}")
        print(f"  GC content: {guide['gc_content']:.1%}")

        # Predict repair outcomes
        context_start = max(0, guide["position"] - 20)
        context_end = min(
            len(target_sequence), guide["position"] + guide["length"] + 20
        )
        context_seq = target_sequence[context_start:context_end]

        repair = designer.predict_repair_outcomes(guide["sequence"], context_seq)
        print(
            f"  Repair outcomes - NHEJ: {repair['nhej_probability']:.2f}, MMEJ: {repair['mmej_probability']:.2f}, HDR: {repair['hdr_efficiency']:.2f}"
        )
        print()


if __name__ == "__main__":
    main()
