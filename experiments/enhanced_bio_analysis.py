"""
Enhanced Bio vs Arbitrary Analysis

This script extends the original bio_v_arbitrary.py with additional biological encoding strategies
and real CRISPR sequence patterns to provide a more comprehensive test of biological anchoring.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import random
from typing import List
import sys
import os

# Add applications directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "applications"))

# Import our main analysis
from modules.bio_v_arbitrary import EmpiricalValidator, ArbitraryEncoder


class BiologicalEncoder3:
    """Third bio-anchored encoding using thermodynamic and kinetic properties"""

    def __init__(self):
        # Based on nearest-neighbor thermodynamic parameters and CRISPR cutting efficiency
        # This encoding focuses on DNA melting temperatures and binding kinetics

        # Melting temperature contribution (normalized)
        melting_contrib = {"A": -1.0, "T": -1.0, "C": 1.0, "G": 1.0}  # GC stabilizes

        # CRISPR cutting preference (empirical data suggests certain patterns)
        crispr_pref = {
            "A": 0.8,
            "T": 0.9,
            "C": 1.2,
            "G": 1.1,
        }  # T slightly favored, C strongly

        # Combine into complex weights
        self.bio_weights = {}
        for base in ["A", "T", "C", "G"]:
            real_part = melting_contrib[base] * 8.0  # Scale to similar magnitude
            imag_part = crispr_pref[base] * 3.0
            self.bio_weights[base] = real_part + imag_part * 1j

    def encode_sequence(self, sequence: str) -> np.ndarray:
        """Encode sequence using thermodynamic properties"""
        encoded = []
        for i, base in enumerate(sequence):
            if base in self.bio_weights:
                # Apply position-based phase modulation with thermodynamic context
                phase = 2 * np.pi * i / len(sequence)

                # Add thermodynamic context from neighboring bases
                thermo_weight = 1.0
                if i > 0 and i < len(sequence) - 1:
                    # Simple nearest-neighbor effect
                    gc_neighbors = sum(
                        1 for b in [sequence[i - 1], sequence[i + 1]] if b in "GC"
                    )
                    thermo_weight = 1.0 + 0.15 * gc_neighbors  # GC neighbors stabilize

                encoded_base = (
                    self.bio_weights[base] * thermo_weight * np.exp(1j * phase)
                )
                encoded.append(encoded_base)
            else:
                encoded.append(0 + 0j)
        return np.array(encoded)


class RealCRISPRSequencesEnhanced:
    """Generator for real CRISPR guide sequences and targets"""

    def __init__(self):
        # Known high-efficiency CRISPR guides (from literature)
        self.high_efficiency_guides = [
            "GCTGCGGAGACCTGGAGAGA",  # PCSK9 guide
            "GTCGAGTGGTTCTGGTAGAG",  # EMX1 guide
            "GCTGTCCTGCAGGAACATCT",  # FANCF guide
            "GTGGCCGAGATTCTGGATGG",  # CD55 guide
            "GCGGGGCCACGGGTCCTCGT",  # AAVS1 guide
            "GGTGGAGGAGGAGGAGCTGG",  # High GC content
            "GATATTCATGAGATGGAATG",  # Balanced composition
            "GTTTCTCGGGAGCTGAATGG",  # Medium efficiency
        ]

        # Target sequences (extended context)
        self.target_contexts = [
            "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAG",
            "CTAGTCGAGTGGTTCTGGTAGAGCGTGAAACCTGTGATCCTGGAAGCTGGGATGCTGGAGTGCTGAAG",
            "GGTGCTGTCCTGCAGGAACATCTGGCGGTGATGGTATCTGGAACGATGGAGACGCTGAAGGTTGAGGT",
            "CTGGTGGCCGAGATTCTGGATGGCTGAACGATGGAGACGCTGAAGGTTGAGGTTCTGAACGCTGTGAA",
        ]

    def get_real_sequences(self, n_sequences: int = 50) -> List[str]:
        """Get a mix of real and realistic CRISPR sequences"""

        # Add real guides (repeated as needed)
        real_guides = []
        for i in range(n_sequences // 2):
            guide = self.high_efficiency_guides[i % len(self.high_efficiency_guides)]
            real_guides.append(guide)

        # Add derived sequences (mutations of real guides)
        derived_guides = []
        for i in range(n_sequences - len(real_guides)):
            base_guide = random.choice(self.high_efficiency_guides)
            # Make 1-2 mutations
            mutated = list(base_guide)
            n_muts = random.randint(1, 2)
            for _ in range(n_muts):
                pos = random.randint(0, len(mutated) - 1)
                new_base = random.choice([b for b in "ATCG" if b != mutated[pos]])
                mutated[pos] = new_base
            derived_guides.append("".join(mutated))

        return real_guides + derived_guides


def run_enhanced_analysis():
    """Run enhanced analysis with real CRISPR sequences and additional bio encoder"""
    print("=" * 80)
    print("ENHANCED BIO vs ARBITRARY ANALYSIS")
    print("=" * 80)
    print("Including real CRISPR sequences and thermodynamic encoding...")
    print()

    # Initialize components
    validator = EmpiricalValidator()
    validator.bio_encoder3 = BiologicalEncoder3()
    real_crispr = RealCRISPRSequencesEnhanced()

    # Test with real CRISPR sequences
    print("Generating real CRISPR sequences...")
    real_sequences = real_crispr.get_real_sequences(80)

    print("Analyzing with 3 biological encoding methods...")
    bio_results1 = validator.analyze_sequence_set(
        real_sequences, validator.bio_encoder, "Bio_v1_MultiProp"
    )
    bio_results2 = validator.analyze_sequence_set(
        real_sequences, validator.bio_encoder2, "Bio_v2_Pairing"
    )
    bio_results3 = validator.analyze_sequence_set(
        real_sequences, validator.bio_encoder3, "Bio_v3_Thermo"
    )

    print("Analyzing arbitrary encodings...")
    arbitrary_results = []
    for trial in range(10):
        arbitrary_encoder = ArbitraryEncoder(seed=trial + 100)  # Different seeds
        arb_result = validator.analyze_sequence_set(
            real_sequences, arbitrary_encoder, f"Arbitrary_Trial_{trial}"
        )
        arbitrary_results.append(arb_result)

    # Statistical analysis
    bio1_z = bio_results1["z_scores"]
    bio2_z = bio_results2["z_scores"]
    bio3_z = bio_results3["z_scores"]
    arb_z_means = [result["mean_z"] for result in arbitrary_results]

    # Comprehensive statistical tests
    t_stat1, p_val1 = stats.ttest_ind(bio1_z, arb_z_means)
    t_stat2, p_val2 = stats.ttest_ind(bio2_z, arb_z_means)
    t_stat3, p_val3 = stats.ttest_ind(bio3_z, arb_z_means)

    # Effect sizes
    def cohens_d(group1, group2):
        pooled_std = np.sqrt(
            (
                (len(group1) - 1) * np.var(group1, ddof=1)
                + (len(group2) - 1) * np.var(group2, ddof=1)
            )
            / (len(group1) + len(group2) - 2)
        )
        return (np.mean(group1) - np.mean(group2)) / pooled_std

    d1 = cohens_d(bio1_z, arb_z_means)
    d2 = cohens_d(bio2_z, arb_z_means)
    d3 = cohens_d(bio3_z, arb_z_means)

    # ANOVA to test if any biological method differs from arbitrary
    all_bio_scores = np.concatenate([bio1_z, bio2_z, bio3_z])
    all_arb_scores = np.array(arb_z_means)

    # Prepare data for ANOVA
    bio_group = np.ones(len(all_bio_scores))
    arb_group = np.zeros(len(all_arb_scores))

    np.concatenate([all_bio_scores, all_arb_scores])
    np.concatenate([bio_group, arb_group])

    f_stat, f_pval = stats.f_oneway(all_bio_scores, all_arb_scores)

    # Print results
    print("=" * 80)
    print("ENHANCED RESULTS SUMMARY")
    print("=" * 80)

    print(f"Bio Method 1 (Multi-prop): {np.mean(bio1_z):.4f} ± {np.std(bio1_z):.4f}")
    print(f"Bio Method 2 (Pairing):    {np.mean(bio2_z):.4f} ± {np.std(bio2_z):.4f}")
    print(f"Bio Method 3 (Thermo):     {np.mean(bio3_z):.4f} ± {np.std(bio3_z):.4f}")
    print(
        f"Arbitrary Mean:            {np.mean(arb_z_means):.4f} ± {np.std(arb_z_means):.4f}"
    )
    print()

    print("Individual t-tests vs Arbitrary:")
    for i, (method, t_stat, p_val, cohen_d) in enumerate(
        [
            ("Bio v1", t_stat1, p_val1, d1),
            ("Bio v2", t_stat2, p_val2, d2),
            ("Bio v3", t_stat3, p_val3, d3),
        ],
        1,
    ):
        sig = "YES" if p_val < 0.05 else "NO"
        print(f"  {method}: t={t_stat:.3f}, p={p_val:.4f}, d={cohen_d:.3f}, Sig={sig}")

    print("\nOverall ANOVA:")
    print(f"  F-statistic: {f_stat:.4f}")
    print(f"  p-value:     {f_pval:.6f}")
    print(f"  Significant: {'YES' if f_pval < 0.05 else 'NO'}")

    # Generate enhanced plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Distribution comparison
    ax1.hist(bio1_z, bins=15, alpha=0.5, label="Bio v1", color="blue", density=True)
    ax1.hist(bio2_z, bins=15, alpha=0.5, label="Bio v2", color="green", density=True)
    ax1.hist(bio3_z, bins=15, alpha=0.5, label="Bio v3", color="purple", density=True)
    ax1.hist(
        arb_z_means, bins=10, alpha=0.7, label="Arbitrary", color="red", density=True
    )
    ax1.set_xlabel("Z-Score")
    ax1.set_ylabel("Density")
    ax1.set_title("Enhanced Z-Score Distributions")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Box plot
    data = [bio1_z, bio2_z, bio3_z, arb_z_means]
    ax2.boxplot(data, tick_labels=["Bio v1", "Bio v2", "Bio v3", "Arbitrary"])
    ax2.set_ylabel("Z-Score")
    ax2.set_title("Z-Score Comparison (Real CRISPR Sequences)")
    ax2.grid(True, alpha=0.3)

    # Plot 3: Mean comparison
    means = [np.mean(bio1_z), np.mean(bio2_z), np.mean(bio3_z), np.mean(arb_z_means)]
    stds = [np.std(bio1_z), np.std(bio2_z), np.std(bio3_z), np.std(arb_z_means)]
    colors = ["blue", "green", "purple", "red"]

    bars = ax3.bar(range(4), means, yerr=stds, capsize=5, color=colors, alpha=0.7)
    ax3.set_xticks(range(4))
    ax3.set_xticklabels(["Bio v1", "Bio v2", "Bio v3", "Arbitrary"])
    ax3.set_ylabel("Mean Z-Score")
    ax3.set_title("Mean Z-Scores with Error Bars")
    ax3.grid(True, alpha=0.3)

    # Plot 4: Effect sizes
    effect_sizes = [d1, d2, d3]
    p_values = [p_val1, p_val2, p_val3]
    colors_bio = ["blue", "green", "purple"]

    bars = ax4.bar(range(3), effect_sizes, color=colors_bio, alpha=0.7)
    ax4.set_xticks(range(3))
    ax4.set_xticklabels(["Bio v1", "Bio v2", "Bio v3"])
    ax4.set_ylabel("Cohen's d (vs Arbitrary)")
    ax4.set_title("Effect Sizes: Biological vs Arbitrary")
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0, color="black", linestyle="--", alpha=0.5)

    # Add significance indicators
    for i, (bar, p_val) in enumerate(zip(bars, p_values)):
        if p_val < 0.001:
            sig_text = "***"
        elif p_val < 0.01:
            sig_text = "**"
        elif p_val < 0.05:
            sig_text = "*"
        else:
            sig_text = "ns"

        height = bar.get_height()
        y_pos = height + 0.05 if height > 0 else height - 0.1
        ax4.text(
            bar.get_x() + bar.get_width() / 2.0,
            y_pos,
            sig_text,
            ha="center",
            va="bottom" if height > 0 else "top",
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig("enhanced_bio_analysis.png", dpi=300, bbox_inches="tight")
    print("\nEnhanced analysis plot saved to: enhanced_bio_analysis.png")

    # Final interpretation
    print("\n" + "=" * 80)
    print("ENHANCED CONCLUSION")
    print("=" * 80)

    bio_better_count = sum(
        1
        for mean_bio, mean_arb in [
            (np.mean(bio1_z), np.mean(arb_z_means)),
            (np.mean(bio2_z), np.mean(arb_z_means)),
            (np.mean(bio3_z), np.mean(arb_z_means)),
        ]
        if mean_bio > mean_arb
    )

    sig_count = sum(1 for p in [p_val1, p_val2, p_val3] if p < 0.05)

    print("Summary:")
    print(f"- {bio_better_count}/3 biological methods outperformed arbitrary encodings")
    print(f"- {sig_count}/3 comparisons showed statistically significant differences")
    print(f"- Overall ANOVA p-value: {f_pval:.6f}")

    if f_pval < 0.05:
        print(
            "\n✓ STATISTICAL EVIDENCE for differences between biological and arbitrary encodings"
        )
    else:
        print("\n✗ NO STRONG STATISTICAL EVIDENCE for biological anchoring superiority")

    print("\nThe analysis using real CRISPR sequences provides a more realistic test")
    print("of biological anchoring in spectral DNA encoding. Results suggest that")
    if bio_better_count >= 2:
        print("biological properties may have some advantages in this context.")
    else:
        print("arbitrary encodings remain competitive with biological approaches.")

    return {
        "bio_results": [bio_results1, bio_results2, bio_results3],
        "arbitrary_results": arbitrary_results,
        "statistics": {
            "f_statistic": f_stat,
            "f_pvalue": f_pval,
            "individual_tests": [
                {
                    "method": "Bio_v1",
                    "t_stat": t_stat1,
                    "p_val": p_val1,
                    "cohens_d": d1,
                },
                {
                    "method": "Bio_v2",
                    "t_stat": t_stat2,
                    "p_val": p_val2,
                    "cohens_d": d2,
                },
                {
                    "method": "Bio_v3",
                    "t_stat": t_stat3,
                    "p_val": p_val3,
                    "cohens_d": d3,
                },
            ],
        },
    }


if __name__ == "__main__":
    results = run_enhanced_analysis()
