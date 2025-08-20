#!/usr/bin/env python3
"""
Generate Synthetic Datasets for Z Framework Validation

Creates realistic synthetic datasets for pain management analysis validation.
All data is synthetic and for research validation purposes only.

Usage:
    python generate_synthetic_data.py
"""

import numpy as np
import pandas as pd
import mpmath as mp
from pathlib import Path
import logging

# Set reproducible seed
np.random.seed(42)
mp.dps = 50

# Constants
PHI = mp.mpf("1.618033988749894848204586834365638117720309179805762862135")


def generate_neural_spikes_data(n_samples=1000):
    """Generate synthetic neural spike data for DRG neuron analysis."""

    # Simulate different pain conditions
    conditions = ["baseline", "neuropathic", "inflammatory", "treated"]

    data = []
    for i in range(n_samples):
        condition = np.random.choice(conditions)

        # Generate spike train parameters based on condition
        if condition == "baseline":
            spike_rate = np.random.normal(5.0, 1.0)  # Hz
            burst_ratio = np.random.normal(0.2, 0.05)
        elif condition == "neuropathic":
            spike_rate = np.random.normal(15.0, 3.0)  # Elevated
            burst_ratio = np.random.normal(0.7, 0.1)  # High bursting
        elif condition == "inflammatory":
            spike_rate = np.random.normal(12.0, 2.5)
            burst_ratio = np.random.normal(0.5, 0.08)
        else:  # treated
            spike_rate = np.random.normal(6.0, 1.5)  # Reduced
            burst_ratio = np.random.normal(0.25, 0.06)

        # Generate sequence representation
        sequence_length = np.random.randint(50, 200)
        # Simple DNA-like encoding for mathematical analysis
        bases = ["A", "T", "G", "C"]
        sequence = "".join(np.random.choice(bases, sequence_length))

        # Encode sequence numerically
        numeric_sequence = [ord(base) - ord("A") + 1 for base in sequence]

        # Generate pain score (0-10 scale)
        if condition == "baseline":
            pain_score = np.random.normal(1.0, 0.5)
        elif condition == "neuropathic":
            pain_score = np.random.normal(8.5, 1.0)
        elif condition == "inflammatory":
            pain_score = np.random.normal(7.0, 1.2)
        else:  # treated
            pain_score = np.random.normal(2.5, 0.8)

        pain_score = np.clip(pain_score, 0, 10)

        data.append(
            {
                "sample_id": f"DRG_{i:04d}",
                "condition": condition,
                "spike_rate_hz": max(0, spike_rate),
                "burst_ratio": np.clip(burst_ratio, 0, 1),
                "sequence": sequence,
                "sequence_numeric": str(numeric_sequence),  # Store as string for CSV
                "pain_score": pain_score,
                "therapeutic_response": 1 if condition == "treated" else 0,
            }
        )

    return pd.DataFrame(data)


def generate_nav18_panel_data(n_samples=500):
    """Generate synthetic Nav1.8 binding affinity panel data."""

    data = []
    for i in range(n_samples):
        # Generate compound properties
        molecular_weight = np.random.normal(350, 50)  # Typical small molecule
        logp = np.random.normal(3.2, 0.8)  # Lipophilicity

        # Generate Nav1.8 selectivity data
        nav18_affinity = np.random.lognormal(3.0, 0.8)  # nM range
        nav17_affinity = nav18_affinity * np.random.lognormal(
            2.0, 0.5
        )  # Less selective
        nav16_affinity = nav18_affinity * np.random.lognormal(
            2.5, 0.6
        )  # Even less selective

        # Calculate selectivity ratio
        selectivity_ratio = nav17_affinity / nav18_affinity

        # Generate sequence representation for the compound
        # Use a pseudo-sequence based on chemical features
        sequence_length = np.random.randint(20, 50)
        bases = ["A", "T", "G", "C"]

        # Weight sequence based on compound properties
        if nav18_affinity < 50:  # High affinity compounds
            weights = [0.3, 0.2, 0.3, 0.2]  # Favor A and G
        else:  # Lower affinity
            weights = [0.25, 0.25, 0.25, 0.25]  # Random

        sequence = "".join(np.random.choice(bases, sequence_length, p=weights))
        numeric_sequence = [ord(base) - ord("A") + 1 for base in sequence]

        # Classification: high affinity if < 100 nM
        high_affinity = 1 if nav18_affinity < 100 else 0

        data.append(
            {
                "compound_id": f"COMP_{i:04d}",
                "molecular_weight": molecular_weight,
                "logp": logp,
                "nav18_affinity_nm": nav18_affinity,
                "nav17_affinity_nm": nav17_affinity,
                "nav16_affinity_nm": nav16_affinity,
                "selectivity_ratio": selectivity_ratio,
                "sequence": sequence,
                "sequence_numeric": str(numeric_sequence),
                "high_affinity": high_affinity,
            }
        )

    return pd.DataFrame(data)


def generate_bcl11a_edits_data(n_samples=300):
    """Generate synthetic BCL11A CRISPR edit data with HbF induction scores."""

    data = []
    for i in range(n_samples):
        # Generate gRNA properties
        grna_length = np.random.choice([20, 21, 22])  # Standard gRNA lengths

        # Generate target sequence
        bases = ["A", "T", "G", "C"]
        grna_sequence = "".join(np.random.choice(bases, grna_length))

        # Calculate GC content
        gc_content = (grna_sequence.count("G") + grna_sequence.count("C")) / len(
            grna_sequence
        )

        # Generate off-target score (lower is better)
        off_target_score = np.random.exponential(0.3)  # Most should be low

        # Generate on-target efficiency (CFD score-like)
        cfd_score = np.random.beta(2, 1) * 100  # Skewed toward higher values

        # Generate HbF induction based on realistic parameters
        # Better gRNAs (high CFD, low off-target, good GC content) induce more HbF
        base_hbf = 5.0  # Baseline HbF %

        # Efficiency factors
        efficiency_factor = (cfd_score / 100) * 0.8  # Max 80% boost from efficiency
        gc_factor = min(1.0, max(0.2, 2 * abs(gc_content - 0.5)))  # Optimal ~50% GC
        off_target_penalty = np.exp(-off_target_score * 2)  # Exponential penalty

        hbf_induction = base_hbf + (
            30 * efficiency_factor * gc_factor * off_target_penalty
        )
        hbf_induction += np.random.normal(0, 2)  # Add noise
        hbf_induction = np.clip(hbf_induction, 0, 95)  # Realistic range

        # Encode sequence numerically
        numeric_sequence = [ord(base) - ord("A") + 1 for base in grna_sequence]

        # Success classification (>20% HbF considered successful)
        successful_edit = 1 if hbf_induction > 20 else 0

        data.append(
            {
                "grna_id": f"gRNA_{i:04d}",
                "grna_sequence": grna_sequence,
                "grna_length": grna_length,
                "gc_content": gc_content,
                "cfd_score": cfd_score,
                "off_target_score": off_target_score,
                "hbf_induction_percent": hbf_induction,
                "sequence_numeric": str(numeric_sequence),
                "successful_edit": successful_edit,
            }
        )

    return pd.DataFrame(data)


def main():
    """Generate all synthetic datasets."""

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Create output directory
    output_dir = Path(__file__).parent
    output_dir.mkdir(exist_ok=True)

    logger.info("Generating synthetic neural spikes data...")
    neural_data = generate_neural_spikes_data(1000)
    neural_data.to_csv(output_dir / "neural_spikes.csv", index=False)
    logger.info(f"Saved {len(neural_data)} neural spike samples")

    logger.info("Generating synthetic Nav1.8 panel data...")
    nav18_data = generate_nav18_panel_data(500)
    nav18_data.to_csv(output_dir / "nav1.8_panel.csv", index=False)
    logger.info(f"Saved {len(nav18_data)} Nav1.8 compound samples")

    logger.info("Generating synthetic BCL11A edit data...")
    bcl11a_data = generate_bcl11a_edits_data(300)
    bcl11a_data.to_csv(output_dir / "bcl11a_edits.csv", index=False)
    logger.info(f"Saved {len(bcl11a_data)} BCL11A edit samples")

    logger.info("All synthetic datasets generated successfully!")
    logger.info("Remember: This is synthetic data for research validation only")


if __name__ == "__main__":
    main()
