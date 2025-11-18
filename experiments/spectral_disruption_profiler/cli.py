"""
Command-Line Interface for Spectral Disruption Profiler

Provides batch processing and analysis commands for CRISPR guide sequences.

Usage:
    python cli.py score --reference REF.fa --mutant MUT.fa --output results.json
    python cli.py batch --input guides.csv --output scores.csv
    python cli.py analyze --sequences data.csv --output analysis/
"""

import argparse
import json
import csv
import sys
import os
from pathlib import Path
from typing import List, Dict
import logging

# Add parent directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from experiments.spectral_disruption_profiler.encoding import validate_dna_sequence
from experiments.spectral_disruption_profiler.analysis import (
    analyze_disruption,
    batch_analyze_disruption,
)
from experiments.spectral_disruption_profiler.scoring import (
    compute_composite_score,
    score_with_confidence,
    auto_optimize_k,
,
)
from experiments.spectral_disruption_profiler.detection import (
    detect_off_targets,
    compute_gc_resonance,
,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def load_sequences_from_csv(
    filepath: str, ref_col: str = "reference", mut_col: str = "mutant"
) -> List[Dict[str, str]]:
    """Load sequence pairs from CSV file."""
    sequences = []

    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if ref_col in row and mut_col in row:
                sequences.append(
                    {
                        "reference": row[ref_col],
                        "mutant": row[mut_col],
                        "id": row.get("id", f"seq_{len(sequences)}"),
                    }
                )

    logger.info(f"Loaded {len(sequences)} sequence pairs from {filepath}")
    return sequences


def load_sequence_from_fasta(filepath: str) -> str:
    """Load single sequence from FASTA file."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Skip header line(s) starting with '>'
    sequence_lines = [line.strip() for line in lines if not line.startswith(">")]
    sequence = "".join(sequence_lines)

    return sequence


def cmd_score(args):
    """Score a single sequence pair."""
    logger.info("=== Spectral Disruption Profiler: Score Command ===")

    # Load sequences
    if args.reference.endswith(".fa") or args.reference.endswith(".fasta"):
        reference_seq = load_sequence_from_fasta(args.reference)
    else:
        reference_seq = args.reference

    if args.mutant.endswith(".fa") or args.mutant.endswith(".fasta"):
        mutant_seq = load_sequence_from_fasta(args.mutant)
    else:
        mutant_seq = args.mutant

    logger.info(
        f"Reference: {reference_seq[:50]}..."
        if len(reference_seq) > 50
        else f"Reference: {reference_seq}"
    )
    logger.info(
        f"Mutant: {mutant_seq[:50]}..."
        if len(mutant_seq) > 50
        else f"Mutant: {mutant_seq}"
    )

    # Validate sequences
    try:
        validate_dna_sequence(reference_seq, is_rna=args.rna)
        validate_dna_sequence(mutant_seq, is_rna=args.rna)
    except ValueError as e:
        logger.error(f"Sequence validation failed: {e}")
        return 1

    # Auto-optimize k if requested
    if args.auto_k:
        logger.info("Auto-optimizing k parameter...")
        k = auto_optimize_k(reference_seq, is_rna=args.rna)
        logger.info(f"Optimized k: {k:.3f}")
    else:
        k = args.k

    # Analyze disruption
    features = analyze_disruption(
        mutant_seq, reference_seq, is_rna=args.rna, phi=args.phi, k=k
    )

    # Compute composite score
    score = compute_composite_score(features)

    # Detect off-targets
    off_target_flags = detect_off_targets([features])

    # Prepare results
    results = {
        "reference_sequence": reference_seq,
        "mutant_sequence": mutant_seq,
        "features": features,
        "composite_score": score,
        "off_target_flag": off_target_flags[0],
        "parameters": {
            "phi": args.phi,
            "k": k,
            "is_rna": args.rna,
        },
    }

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to {args.output}")
    logger.info(f"Composite score: {score:.4f}")
    logger.info(f"Off-target flag: {off_target_flags[0]}")

    return 0


def cmd_batch(args):
    """Batch process multiple sequences."""
    logger.info("=== Spectral Disruption Profiler: Batch Command ===")

    # Load sequences
    sequence_pairs = load_sequences_from_csv(args.input)

    if len(sequence_pairs) == 0:
        logger.error("No sequences loaded from input file")
        return 1

    # Extract sequences
    references = [pair["reference"] for pair in sequence_pairs]
    mutants = [pair["mutant"] for pair in sequence_pairs]
    ids = [pair["id"] for pair in sequence_pairs]

    # Analyze disruption
    logger.info(f"Analyzing {len(sequence_pairs)} sequence pairs...")
    features_list = batch_analyze_disruption(
        mutants, references, is_rna=args.rna, phi=args.phi, k=args.k
    )

    # Compute scores with confidence intervals
    if args.bootstrap > 0:
        logger.info(f"Computing scores with {args.bootstrap} bootstrap resamples...")
        score_results = score_with_confidence(
            features_list, n_bootstrap=args.bootstrap, seed=args.seed
        )

    # Compute individual scores
    scores = [compute_composite_score(f) for f in features_list]

    # Detect off-targets
    off_target_flags = detect_off_targets(features_list)

    # Compute GC resonance
    gc_resonance = compute_gc_resonance(
        features_list, n_permutations=args.permutations, seed=args.seed
    )

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "id",
                "reference",
                "mutant",
                "delta_f1",
                "delta_entropy",
                "delta_sidelobes",
                "gc_content",
                "composite_score",
                "off_target_flag",
            ],
        )
        writer.writeheader()

        for i, (seq_id, ref, mut, features, score, flag) in enumerate(
            zip(ids, references, mutants, features_list, scores, off_target_flags)
        ):
            writer.writerow(
                {
                    "id": seq_id,
                    "reference": ref,
                    "mutant": mut,
                    "delta_f1": features.get("delta_f1", 0.0),
                    "delta_entropy": features.get("delta_entropy", 0.0),
                    "delta_sidelobes": features.get("delta_sidelobes", 0),
                    "gc_content": features.get("gc_content", 0.0),
                    "composite_score": score,
                    "off_target_flag": flag,
                }
            )

    logger.info(f"Results saved to {args.output}")

    # Save summary statistics
    summary_path = output_path.parent / f"{output_path.stem}_summary.json"
    summary = {
        "n_sequences": len(sequence_pairs),
        "n_off_targets": sum(off_target_flags),
        "gc_resonance": gc_resonance,
        "parameters": {
            "phi": args.phi,
            "k": args.k,
            "is_rna": args.rna,
            "bootstrap": args.bootstrap,
            "permutations": args.permutations,
            "seed": args.seed,
        },
    }

    if args.bootstrap > 0:
        summary["score_statistics"] = score_results

    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Summary saved to {summary_path}")

    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Spectral Disruption Profiler for CRISPR Guide Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to execute")

    # Score command
    score_parser = subparsers.add_parser("score", help="Score a single sequence pair")
    score_parser.add_argument(
        "--reference", required=True, help="Reference sequence or FASTA file"
    )
    score_parser.add_argument(
        "--mutant", required=True, help="Mutant sequence or FASTA file"
    )
    score_parser.add_argument(
        "--output", default="results.json", help="Output JSON file"
    )
    score_parser.add_argument(
        "--rna", action="store_true", help="Use RNA encoding (U instead of T)"
    )
    score_parser.add_argument(
        "--phi", type=float, default=21.0, help="Geometric period (default: 21)"
    )
    score_parser.add_argument(
        "--k", type=float, default=0.3, help="Curvature parameter (default: 0.3)"
    )
    score_parser.add_argument(
        "--auto-k", action="store_true", help="Auto-optimize k parameter"
    )

    # Batch command
    batch_parser = subparsers.add_parser(
        "batch", help="Batch process multiple sequences"
    )
    batch_parser.add_argument(
        "--input",
        required=True,
        help="Input CSV file with reference and mutant columns",
    )
    batch_parser.add_argument("--output", default="scores.csv", help="Output CSV file")
    batch_parser.add_argument("--rna", action="store_true", help="Use RNA encoding")
    batch_parser.add_argument(
        "--phi", type=float, default=21.0, help="Geometric period"
    )
    batch_parser.add_argument(
        "--k", type=float, default=0.3, help="Curvature parameter"
    )
    batch_parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Bootstrap resamples (default: 1000)",
    )
    batch_parser.add_argument(
        "--permutations",
        type=int,
        default=1000,
        help="Permutations for GC resonance (default: 1000)",
    )
    batch_parser.add_argument(
        "--seed", type=int, default=42, help="Random seed (default: 42)"
    )

    args = parser.parse_args()

    if args.command == "score":
        return cmd_score(args)
    elif args.command == "batch":
        return cmd_batch(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
