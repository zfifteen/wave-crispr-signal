#!/usr/bin/env python3
"""
Example Usage: FFT-Based CRISPR Disruption Metrics with Golden-Ratio Phase

This script demonstrates practical applications of the FFT-based analyzer for:
1. Scoring gRNA off-target risk
2. Detecting periodicities in edited sequences
3. Quantifying disruption from indels
4. Analyzing codon-aligned patterns

Author: Z Framework CRISPR Signal Team
"""

from fft_crispr_disruption import (
    FFTCRISPRDisruptionAnalyzer,
    calculate_grna_off_target_score
)


def example_1_score_grna_guides():
    """Example 1: Score multiple gRNA candidates for off-target risk."""
    print("="*70)
    print("EXAMPLE 1: Score gRNA Candidates for Off-Target Risk")
    print("="*70)
    
    # Candidate gRNA sequences (20-nt guides)
    candidates = {
        "Guide-A": "GACGATCGATCGATCGATCG",
        "Guide-B": "ATCGTAGCTACGATCGTAGC",
        "Guide-C": "GGGGAAAACCCCTTTTAAAA",  # Highly structured
        "Guide-D": "GATCAGTCGACTGACTGACT",
    }
    
    print("\nScoring gRNA candidates...")
    print("-"*70)
    
    results = []
    for name, sequence in candidates.items():
        score = calculate_grna_off_target_score(sequence)
        results.append((name, score))
        
        print(f"\n{name}: {sequence}")
        print(f"  Off-Target Score: {score['off_target_score']:.4f}")
        print(f"  Recommendation: {score['recommendation'].upper()}")
        print(f"  Significant Peaks: {score['n_significant_peaks']}")
    
    # Sort by score (higher is better)
    results.sort(key=lambda x: x[1]['off_target_score'], reverse=True)
    
    print("\n" + "-"*70)
    print("RANKING (Best to Worst):")
    print("-"*70)
    for i, (name, score) in enumerate(results, 1):
        print(f"{i}. {name}: {score['off_target_score']:.4f} ({score['recommendation']})")


def example_2_detect_periodicities():
    """Example 2: Detect off-target periodicities in a sequence."""
    print("\n" + "="*70)
    print("EXAMPLE 2: Detect Off-Target Periodicities")
    print("="*70)
    
    analyzer = FFTCRISPRDisruptionAnalyzer()
    
    sequence = "ATCGATCGATCGATCGATCG"  # 4-bp periodicity
    
    print(f"\nAnalyzing: {sequence}")
    print("-"*70)
    
    analysis = analyzer.detect_off_target_periodicities(sequence)
    
    print(f"\nSequence Length: {analysis['sequence_length']} bp")
    print(f"Significant Peaks Detected: {analysis['n_significant_peaks']}")
    
    if analysis['n_significant_peaks'] > 0:
        print("\nTop Periodicities:")
        for i, peak in enumerate(analysis['significant_peaks'][:3], 1):
            print(f"\n  Peak {i}:")
            print(f"    Period: {peak['period']:.2f} bp")
            print(f"    Frequency: {peak['frequency']:.4f} cycles/bp")
            print(f"    Magnitude: {peak['magnitude']:.3f}")
            print(f"    Weighted Magnitude: {peak['weighted_magnitude']:.3f}")
            print(f"    θ′ Weight: {peak['theta_prime_weight']:.3f}")


def example_3_quantify_indel_disruption():
    """Example 3: Quantify disruption from CRISPR-induced indels."""
    print("\n" + "="*70)
    print("EXAMPLE 3: Quantify Disruption from CRISPR Indels")
    print("="*70)
    
    analyzer = FFTCRISPRDisruptionAnalyzer()
    
    reference = "ATCGATCGATCGATCGATCG"
    
    print(f"\nReference Sequence: {reference}")
    print("-"*70)
    
    # Test different indel scenarios
    scenarios = [
        {
            'name': '1-bp deletion',
            'position': 10,
            'length': 1,
            'type': 'deletion'
        },
        {
            'name': '3-bp deletion',
            'position': 10,
            'length': 3,
            'type': 'deletion'
        },
        {
            'name': '3-bp insertion',
            'position': 10,
            'length': 3,
            'type': 'insertion'
        }
    ]
    
    for scenario in scenarios:
        print(f"\n{scenario['name'].upper()}:")
        print(f"  Position: {scenario['position']}")
        print(f"  Length: {scenario['length']} bp")
        
        result = analyzer.analyze_indel_disruption(
            reference,
            scenario['position'],
            scenario['length'],
            scenario['type']
        )
        
        print(f"\n  Disruption Metrics:")
        print(f"    Overall Score: {result['disruption_score']:.4f}")
        print(f"    ΔEntropy: {result['delta_entropy']:.4f}")
        print(f"    Δf₁ (relative): {result['delta_f1_relative']:.4f}")
        print(f"    Phase Disruption: {result['phase_disruption']:.4f}")
        
        # Interpretation
        if result['disruption_score'] < 0.5:
            severity = "Low"
        elif result['disruption_score'] < 1.0:
            severity = "Moderate"
        else:
            severity = "High"
        print(f"    Severity: {severity}")


def example_4_compare_before_after_editing():
    """Example 4: Compare reference vs edited sequence."""
    print("\n" + "="*70)
    print("EXAMPLE 4: Compare Reference vs Edited Sequence")
    print("="*70)
    
    analyzer = FFTCRISPRDisruptionAnalyzer()
    
    reference = "ATCGATCGATCGATCGATCG"
    edited = "ATCGATCGATCGATCG"  # 3-bp deletion
    
    print(f"\nReference: {reference}")
    print(f"Edited:    {edited}")
    print("-"*70)
    
    disruption = analyzer.calculate_disruption_score(reference, edited)
    
    print("\nDisruption Analysis:")
    print(f"  Composite Score: {disruption['disruption_score']:.4f}")
    
    print("\n  Component Metrics:")
    print(f"    Reference Entropy: {disruption['reference_entropy']:.4f}")
    print(f"    Edited Entropy: {disruption['edited_entropy']:.4f}")
    print(f"    ΔEntropy: {disruption['delta_entropy']:.4f}")
    
    print(f"\n    Δf₁: {disruption['delta_f1']:.4f}")
    print(f"    Δf₁ (relative): {disruption['delta_f1_relative']:.4f}")
    
    print(f"\n    Reference Peaks: {disruption['reference_peaks']}")
    print(f"    Edited Peaks: {disruption['edited_peaks']}")
    print(f"    ΔPeaks: {disruption['delta_sidelobes']}")
    
    print(f"\n    Phase Disruption: {disruption['phase_disruption']:.4f}")


def example_5_codon_aligned_analysis():
    """Example 5: Analyze codon-aligned patterns."""
    print("\n" + "="*70)
    print("EXAMPLE 5: Codon-Aligned φ-Structured Analysis")
    print("="*70)
    
    analyzer = FFTCRISPRDisruptionAnalyzer()
    
    # 21-bp sequence (7 codons)
    sequence = "ATCGATCGATCGATCGATCGA"
    
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bp (expected: 7 codons)")
    print("-"*70)
    
    # Analyze all three reading frames
    for frame in [0, 1, 2]:
        print(f"\nReading Frame {frame}:")
        
        codon_features = analyzer.calculate_codon_aligned_features(
            sequence,
            frame=frame
        )
        
        print(f"  Codons Analyzed: {codon_features['n_codons']}")
        print(f"  Codon Entropy: {codon_features['codon_entropy']:.4f}")
        print(f"  Dominant Period: {codon_features['dominant_codon_period']:.2f} codons")
        
        # Show first 3 codon values
        codon_vals = codon_features['codon_values'][:3]
        print(f"  First 3 Codons: {[f'{v:.3f}' for v in codon_vals]}")


def example_6_batch_analysis():
    """Example 6: Batch analysis of multiple guides."""
    print("\n" + "="*70)
    print("EXAMPLE 6: Batch Analysis of Multiple Guides")
    print("="*70)
    
    # Simulate a set of designed guides
    guides = [
        "GACGATCGATCGATCGATCG",
        "ATCGTAGCTACGATCGTAGC",
        "GATCAGTCGACTGACTGACT",
        "TACGATCGACTGATCGATCG",
        "CGTAGCTACGATCGTAGCTA",
    ]
    
    print(f"\nAnalyzing {len(guides)} gRNA candidates...")
    print("-"*70)
    
    results = []
    for i, guide in enumerate(guides, 1):
        score = calculate_grna_off_target_score(guide)
        results.append({
            'id': f'Guide-{i}',
            'sequence': guide,
            'score': score['off_target_score'],
            'recommendation': score['recommendation'],
            'peaks': score['n_significant_peaks']
        })
    
    # Sort by score
    results.sort(key=lambda x: x['score'], reverse=True)
    
    print(f"\n{'Rank':<6} {'Guide ID':<12} {'Score':<10} {'Peaks':<8} {'Recommendation':<15}")
    print("-"*70)
    
    for rank, result in enumerate(results, 1):
        print(f"{rank:<6} {result['id']:<12} {result['score']:<10.4f} "
              f"{result['peaks']:<8} {result['recommendation']:<15}")
    
    # Statistics
    scores = [r['score'] for r in results]
    print("\n" + "-"*70)
    print(f"Mean Score: {sum(scores)/len(scores):.4f}")
    print(f"Best Score: {max(scores):.4f}")
    print(f"Worst Score: {min(scores):.4f}")


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("FFT-BASED CRISPR DISRUPTION METRICS")
    print("Golden-Ratio Phase Analysis Examples")
    print("="*70)
    
    # Run examples
    example_1_score_grna_guides()
    example_2_detect_periodicities()
    example_3_quantify_indel_disruption()
    example_4_compare_before_after_editing()
    example_5_codon_aligned_analysis()
    example_6_batch_analysis()
    
    print("\n" + "="*70)
    print("EXAMPLES COMPLETE")
    print("="*70)
    print("\nFor more details, see docs/FFT_GOLDEN_RATIO_CRISPR.md")


if __name__ == "__main__":
    main()
