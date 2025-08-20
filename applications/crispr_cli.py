"""
CRISPR CLI - Command Line Interface for CRISPR Guide Design

This module provides a command-line interface for the CRISPR guide designer
using signal-theoretic DNA analysis.
"""

import argparse
import sys
import json
from typing import List, Dict, Optional

try:
    from .crispr_guide_designer import CRISPRGuideDesigner
except ImportError:
    # Handle relative import for direct execution
    import sys

    sys.path.append(".")
    from crispr_guide_designer import CRISPRGuideDesigner


def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Read sequences from FASTA file.

    Args:
        filepath: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Add last sequence
        if current_name:
            sequences[current_name] = "".join(current_seq)

    return sequences


def write_output(
    results: List[Dict], output_file: Optional[str] = None, format: str = "json"
):
    """
    Write results to file or stdout.

    Args:
        results: Guide design results
        output_file: Output file path (None for stdout)
        format: Output format ('json', 'csv', 'tsv')
    """
    if format == "json":
        output = json.dumps(results, indent=2)
    elif format == "csv":
        output = format_csv(results)
    elif format == "tsv":
        output = format_tsv(results)
    else:
        raise ValueError(f"Unsupported format: {format}")

    if output_file:
        with open(output_file, "w") as f:
            f.write(output)
        print(f"Results written to {output_file}")
    else:
        print(output)


def format_csv(results: List[Dict]) -> str:
    """Format results as CSV."""
    if not results:
        return ""

    headers = [
        "sequence_name",
        "guide_sequence",
        "position",
        "pam_sequence",
        "on_target_score",
        "gc_content",
        "nhej_prob",
        "mmej_prob",
        "hdr_eff",
    ]

    lines = [",".join(headers)]

    for result in results:
        for guide in result.get("guides", []):
            repair = guide.get("repair_outcomes", {})
            row = [
                result["sequence_name"],
                guide["sequence"],
                str(guide["position"]),
                guide["pam_sequence"],
                f"{guide['on_target_score']:.3f}",
                f"{guide['gc_content']:.3f}",
                f"{repair.get('nhej_probability', 0):.3f}",
                f"{repair.get('mmej_probability', 0):.3f}",
                f"{repair.get('hdr_efficiency', 0):.3f}",
            ]
            lines.append(",".join(row))

    return "\n".join(lines)


def format_tsv(results: List[Dict]) -> str:
    """Format results as TSV."""
    return format_csv(results).replace(",", "\t")


def design_guides_command(args):
    """Handle guide design command."""
    designer = CRISPRGuideDesigner(pam_pattern=args.pam, guide_length=args.length)

    # Read input sequences
    if args.input.endswith(".fa") or args.input.endswith(".fasta"):
        sequences = read_fasta(args.input)
    else:
        # Treat as raw sequence
        sequences = {"input_sequence": args.input.upper()}

    results = []

    for seq_name, sequence in sequences.items():
        print(
            f"Designing guides for {seq_name} ({len(sequence)} bp)...", file=sys.stderr
        )

        # Design guides
        guides = designer.design_guides(sequence, num_guides=args.num_guides)

        # Add repair outcome predictions if requested
        if args.predict_repair:
            for guide in guides:
                context_start = max(0, guide["position"] - 20)
                context_end = min(
                    len(sequence), guide["position"] + guide["length"] + 20
                )
                context_seq = sequence[context_start:context_end]
                repair = designer.predict_repair_outcomes(
                    guide["sequence"], context_seq
                )
                guide["repair_outcomes"] = repair

        results.append(
            {
                "sequence_name": seq_name,
                "sequence_length": len(sequence),
                "guides": guides,
            }
        )

    # Write output
    write_output(results, args.output, args.format)


def score_guide_command(args):
    """Handle guide scoring command."""
    designer = CRISPRGuideDesigner(pam_pattern=args.pam, guide_length=args.length)

    guide_seq = args.guide.upper()

    # Calculate scores
    on_target_score = designer.calculate_on_target_score(guide_seq)

    # Calculate GC content
    gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq)

    result = {
        "guide_sequence": guide_seq,
        "length": len(guide_seq),
        "on_target_score": on_target_score,
        "gc_content": gc_content,
    }

    # Add off-target analysis if provided
    if args.off_target:
        off_target_risk = designer.calculate_off_target_risk(
            guide_seq, args.guide, args.off_target.upper()
        )
        result["off_target_risk"] = off_target_risk

    # Add repair predictions if target context provided
    if args.target_context:
        repair = designer.predict_repair_outcomes(
            guide_seq, args.target_context.upper()
        )
        # Convert numpy types to Python types for JSON serialization
        repair = {k: float(v) if hasattr(v, "item") else v for k, v in repair.items()}
        result["repair_outcomes"] = repair

    # Output result
    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Results written to {args.output}")
    else:
        print(json.dumps(result, indent=2))


def batch_score_command(args):
    """Handle batch scoring command."""
    designer = CRISPRGuideDesigner(pam_pattern=args.pam, guide_length=args.length)

    # Read guide sequences
    guides = []
    with open(args.guides_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                guides.append(line.upper())

    results = []

    for i, guide_seq in enumerate(guides):
        print(f"Scoring guide {i+1}/{len(guides)}...", file=sys.stderr)

        on_target_score = designer.calculate_on_target_score(guide_seq)
        gc_content = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq)

        result = {
            "guide_sequence": guide_seq,
            "length": len(guide_seq),
            "on_target_score": on_target_score,
            "gc_content": gc_content,
        }

        results.append(result)

    # Write output
    write_output(results, args.output, args.format)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="CRISPR Guide Designer CLI - Signal-Theoretic DNA Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Design guides for a sequence
  python crispr_cli.py design ATGCGATCGATCGATCG

  # Design guides from FASTA file
  python crispr_cli.py design target.fasta -o results.json

  # Score a specific guide
  python crispr_cli.py score GACGATCGATCGATCGATCG --target-context AATGACGATCGATCGATCGATCGTGG

  # Batch score guides from file
  python crispr_cli.py batch-score guides.txt -o scores.csv -f csv
        """,
    )

    # Global arguments
    parser.add_argument("--pam", default="NGG", help="PAM pattern (default: NGG)")
    parser.add_argument(
        "--length", type=int, default=20, help="Guide length (default: 20)"
    )
    parser.add_argument("--output", "-o", help="Output file (default: stdout)")
    parser.add_argument(
        "--format",
        "-f",
        choices=["json", "csv", "tsv"],
        default="json",
        help="Output format (default: json)",
    )

    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Design command
    design_parser = subparsers.add_parser("design", help="Design CRISPR guides")
    design_parser.add_argument("input", help="Input sequence or FASTA file")
    design_parser.add_argument(
        "--num-guides",
        "-n",
        type=int,
        default=5,
        help="Number of guides to return (default: 5)",
    )
    design_parser.add_argument(
        "--predict-repair", action="store_true", help="Predict repair outcomes"
    )

    # Score command
    score_parser = subparsers.add_parser("score", help="Score a specific guide")
    score_parser.add_argument("guide", help="Guide RNA sequence to score")
    score_parser.add_argument(
        "--off-target", help="Off-target sequence for risk analysis"
    )
    score_parser.add_argument(
        "--target-context", help="Target context for repair prediction"
    )

    # Batch score command
    batch_parser = subparsers.add_parser(
        "batch-score", help="Score multiple guides from file"
    )
    batch_parser.add_argument(
        "guides_file", help="File with guide sequences (one per line)"
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        if args.command == "design":
            design_guides_command(args)
        elif args.command == "score":
            score_guide_command(args)
        elif args.command == "batch-score":
            batch_score_command(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
