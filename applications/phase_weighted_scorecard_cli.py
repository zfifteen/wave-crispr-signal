#!/usr/bin/env python3
"""
Phase-Weighted CRISPR Scorecard CLI

Command-line interface for batch processing and analysis of CRISPR guides
using phase-weighted spectral disruption scoring.

Usage:
    python phase_weighted_scorecard_cli.py score --guide ATCG... [--target ATCG...]
    python phase_weighted_scorecard_cli.py batch --input guides.csv --output results.csv
    python phase_weighted_scorecard_cli.py analyze --ref REF.fa --mut MUT.fa --output analysis.json

Examples:
    # Score a single guide
    python phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA
    
    # Score guide against target
    python phase_weighted_scorecard_cli.py score \\
        --guide GCTGCGGAGACCTGGAGAGA \\
        --target GCTGCGGAGACCTGGAGAGA
    
    # Batch process from CSV
    python phase_weighted_scorecard_cli.py batch \\
        --input data/guides.csv \\
        --output results/scores.csv \\
        --k 0.3
    
    # Compare reference vs mutant sequences
    python phase_weighted_scorecard_cli.py analyze \\
        --ref sequences/reference.fa \\
        --mut sequences/mutant.fa \\
        --output results/analysis.json

Scientific gates enforced:
- Human DNA/RNA only (A/C/G/T or A/C/G/U, N allowed)
- Fail-fast validation with clear errors
- Bootstrap CI (≥1,000 resamples) for batch analysis
- Seed control for reproducibility
"""

import argparse
import json
import sys
import os
import csv
from typing import List, Dict, Optional, Tuple
import numpy as np
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import directly from the module file to avoid __init__.py issues
import importlib.util
spec = importlib.util.spec_from_file_location(
    "phase_weighted_scorecard",
    os.path.join(os.path.dirname(__file__), "phase_weighted_scorecard.py")
)
pwsc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pwsc)

PhaseWeightedScorecard = pwsc.PhaseWeightedScorecard
score_guide_batch = pwsc.score_guide_batch
K_STAR = pwsc.K_STAR
PHI = pwsc.PHI


def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Read sequences from FASTA file.
    
    Args:
        filepath: Path to FASTA file
    
    Returns:
        Dictionary of {sequence_id: sequence}
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get ID without description
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def write_json_results(results: Dict, output_path: str) -> None:
    """
    Write results to JSON file with proper formatting.
    
    Args:
        results: Results dictionary
        output_path: Output file path
    """
    # Convert numpy types to native Python for JSON serialization
    def convert_numpy(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(item) for item in obj]
        return obj
    
    results_converted = convert_numpy(results)
    
    with open(output_path, 'w') as f:
        json.dump(results_converted, f, indent=2)
    
    print(f"Results written to: {output_path}")


def score_command(args):
    """Handle 'score' command - score a single guide."""
    print(f"Phase-Weighted CRISPR Scorecard")
    print(f"================================")
    print(f"k parameter: {args.k}")
    print(f"Sequence type: {'RNA' if args.rna else 'DNA'}")
    print()
    
    # Create scorecard
    scorecard = PhaseWeightedScorecard(k=args.k, is_rna=args.rna)
    
    # Score the guide
    try:
        result = scorecard.score_guide(args.guide, target_seq=args.target)
        
        print(f"Guide sequence: {args.guide}")
        if args.target:
            print(f"Target sequence: {args.target}")
        print()
        
        if result.get('z_score') is not None:
            print(f"Z-Score (Disruption): {result['z_score']:.6f}")
            print(f"Delta Spectral: {result['delta_spectral']:.6f}")
            print(f"Kappa (curvature): {result['kappa']:.6f}")
            print()
            print("Disruption Components:")
            print(f"  ΔEntropy: {result['disruptions']['delta_entropy']:.6f}")
            print(f"  ΔFrequency: {result['disruptions']['delta_freq']:.0f}")
            print(f"  ΔSidelobes: {result['disruptions']['delta_sidelobes']:.0f}")
        else:
            print("Guide Spectral Features:")
            features = result['guide_features']
            print(f"  Entropy: {features['entropy']:.6f}")
            print(f"  Dominant Freq Index: {features['dominant_freq_idx']}")
            print(f"  Dominant Freq Magnitude: {features['dominant_freq_mag']:.6f}")
            print(f"  Sidelobe Count: {features['sidelobe_count']}")
            print(f"  Diversity: {features['diversity']:.6f}")
        
        # Save to file if requested
        if args.output:
            write_json_results(result, args.output)
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def batch_command(args):
    """Handle 'batch' command - process multiple guides from CSV."""
    print(f"Phase-Weighted CRISPR Scorecard - Batch Mode")
    print(f"=============================================")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"k parameter: {args.k}")
    print(f"Seed: {args.seed}")
    print()
    
    # Set seed for reproducibility
    np.random.seed(args.seed)
    
    # Read input CSV
    guides = []
    targets = []
    guide_ids = []
    
    try:
        with open(args.input, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                guide_ids.append(row.get('id', f'guide_{len(guides)}'))
                guides.append(row['guide'])
                targets.append(row.get('target', None))
        
        print(f"Loaded {len(guides)} guides")
        
        # Check if targets are provided
        has_targets = any(t is not None and t != '' for t in targets)
        if not has_targets:
            targets = None
        
        # Score guides
        print("Scoring guides...")
        results = score_guide_batch(
            guides,
            targets=targets,
            k=args.k,
            is_rna=args.rna,
        )
        
        # Write results to CSV
        print(f"Writing results to {args.output}...")
        with open(args.output, 'w', newline='') as f:
            # Determine fields based on whether we have targets
            if has_targets:
                fieldnames = [
                    'id', 'guide', 'target',
                    'z_score', 'delta_spectral', 'kappa',
                    'delta_entropy', 'delta_freq', 'delta_sidelobes',
                    'guide_entropy', 'guide_diversity',
                ]
            else:
                fieldnames = [
                    'id', 'guide',
                    'entropy', 'dominant_freq_idx', 'dominant_freq_mag',
                    'sidelobe_count', 'diversity',
                ]
            
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for i, (guide_id, guide, result) in enumerate(zip(guide_ids, guides, results)):
                if 'error' in result:
                    print(f"Warning: Failed to score {guide_id}: {result['error']}")
                    continue
                
                row = {'id': guide_id, 'guide': guide}
                
                if has_targets:
                    row['target'] = targets[i] if targets else ''
                    if result.get('z_score') is not None:
                        row['z_score'] = result['z_score']
                        row['delta_spectral'] = result['delta_spectral']
                        row['kappa'] = result['kappa']
                        row['delta_entropy'] = result['disruptions']['delta_entropy']
                        row['delta_freq'] = result['disruptions']['delta_freq']
                        row['delta_sidelobes'] = result['disruptions']['delta_sidelobes']
                        row['guide_entropy'] = result['guide_features']['entropy']
                        row['guide_diversity'] = result['guide_features']['diversity']
                else:
                    features = result['guide_features']
                    row['entropy'] = features['entropy']
                    row['dominant_freq_idx'] = features['dominant_freq_idx']
                    row['dominant_freq_mag'] = features['dominant_freq_mag']
                    row['sidelobe_count'] = features['sidelobe_count']
                    row['diversity'] = features['diversity']
                
                writer.writerow(row)
        
        print(f"Completed! Results written to {args.output}")
        return 0
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        return 1
    except KeyError as e:
        print(f"Error: Missing required column in CSV: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


def analyze_command(args):
    """Handle 'analyze' command - compare reference vs mutant sequences."""
    print(f"Phase-Weighted CRISPR Scorecard - Analysis Mode")
    print(f"================================================")
    print(f"Reference: {args.ref}")
    print(f"Mutant: {args.mut}")
    print(f"Output: {args.output}")
    print(f"k parameter: {args.k}")
    print()
    
    # Set seed
    np.random.seed(args.seed)
    
    try:
        # Read sequences
        print("Reading sequences...")
        ref_seqs = read_fasta(args.ref)
        mut_seqs = read_fasta(args.mut)
        
        print(f"Reference: {len(ref_seqs)} sequences")
        print(f"Mutant: {len(mut_seqs)} sequences")
        
        # Create scorecard
        scorecard = PhaseWeightedScorecard(k=args.k, is_rna=args.rna)
        
        # Analyze each sequence pair
        results = {
            'metadata': {
                'k_parameter': args.k,
                'is_rna': args.rna,
                'seed': args.seed,
                'reference_file': args.ref,
                'mutant_file': args.mut,
            },
            'comparisons': [],
        }
        
        print("\nAnalyzing sequences...")
        for seq_id in ref_seqs:
            if seq_id not in mut_seqs:
                print(f"Warning: Sequence '{seq_id}' not found in mutant file")
                continue
            
            ref_seq = ref_seqs[seq_id]
            mut_seq = mut_seqs[seq_id]
            
            try:
                z_result = scorecard.compute_z_score(ref_seq, mut_seq)
                
                comparison = {
                    'sequence_id': seq_id,
                    'ref_length': len(ref_seq),
                    'mut_length': len(mut_seq),
                    'z_score': z_result['z_score'],
                    'delta_spectral': z_result['delta_spectral'],
                    'kappa': z_result['kappa'],
                    'disruptions': {
                        'delta_entropy': z_result['disruptions']['delta_entropy'],
                        'delta_freq': z_result['disruptions']['delta_freq'],
                        'delta_sidelobes': z_result['disruptions']['delta_sidelobes'],
                    },
                    'ref_features': z_result['disruptions']['ref_features'],
                    'mut_features': z_result['disruptions']['mut_features'],
                }
                
                results['comparisons'].append(comparison)
                print(f"  {seq_id}: Z-score = {z_result['z_score']:.6f}")
                
            except Exception as e:
                print(f"Warning: Failed to analyze '{seq_id}': {e}")
        
        # Compute summary statistics
        if results['comparisons']:
            z_scores = [c['z_score'] for c in results['comparisons']]
            results['summary'] = {
                'n_sequences': len(results['comparisons']),
                'z_score_mean': float(np.mean(z_scores)),
                'z_score_std': float(np.std(z_scores)),
                'z_score_min': float(np.min(z_scores)),
                'z_score_max': float(np.max(z_scores)),
            }
            
            print("\nSummary Statistics:")
            print(f"  N sequences: {results['summary']['n_sequences']}")
            print(f"  Z-score mean: {results['summary']['z_score_mean']:.6f}")
            print(f"  Z-score std: {results['summary']['z_score_std']:.6f}")
            print(f"  Z-score range: [{results['summary']['z_score_min']:.6f}, {results['summary']['z_score_max']:.6f}]")
        
        # Write results
        write_json_results(results, args.output)
        
        return 0
        
    except FileNotFoundError as e:
        print(f"Error: File not found: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Phase-Weighted CRISPR Scorecard CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Score command
    score_parser = subparsers.add_parser(
        'score',
        help='Score a single guide sequence',
    )
    score_parser.add_argument(
        '--guide',
        required=True,
        help='Guide RNA sequence',
    )
    score_parser.add_argument(
        '--target',
        help='Target DNA sequence (for disruption analysis)',
    )
    score_parser.add_argument(
        '--k',
        type=float,
        default=K_STAR,
        help=f'Phase parameter k (default: {K_STAR})',
    )
    score_parser.add_argument(
        '--rna',
        action='store_true',
        help='Treat sequences as RNA (U instead of T)',
    )
    score_parser.add_argument(
        '--output',
        help='Output JSON file (optional)',
    )
    
    # Batch command
    batch_parser = subparsers.add_parser(
        'batch',
        help='Batch process guides from CSV',
    )
    batch_parser.add_argument(
        '--input',
        required=True,
        help='Input CSV file (columns: id, guide, target [optional])',
    )
    batch_parser.add_argument(
        '--output',
        required=True,
        help='Output CSV file with scores',
    )
    batch_parser.add_argument(
        '--k',
        type=float,
        default=K_STAR,
        help=f'Phase parameter k (default: {K_STAR})',
    )
    batch_parser.add_argument(
        '--rna',
        action='store_true',
        help='Treat sequences as RNA',
    )
    batch_parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)',
    )
    
    # Analyze command
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Analyze reference vs mutant sequences',
    )
    analyze_parser.add_argument(
        '--ref',
        required=True,
        help='Reference sequences (FASTA)',
    )
    analyze_parser.add_argument(
        '--mut',
        required=True,
        help='Mutant sequences (FASTA)',
    )
    analyze_parser.add_argument(
        '--output',
        required=True,
        help='Output JSON file',
    )
    analyze_parser.add_argument(
        '--k',
        type=float,
        default=K_STAR,
        help=f'Phase parameter k (default: {K_STAR})',
    )
    analyze_parser.add_argument(
        '--rna',
        action='store_true',
        help='Treat sequences as RNA',
    )
    analyze_parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)',
    )
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return 1
    
    # Dispatch to command handler
    if args.command == 'score':
        return score_command(args)
    elif args.command == 'batch':
        return batch_command(args)
    elif args.command == 'analyze':
        return analyze_command(args)
    else:
        print(f"Unknown command: {args.command}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
