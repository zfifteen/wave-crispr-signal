#!/usr/bin/env python3
"""
Minimal Bin-Resonance Test for Human CRISPR Efficiency

This module implements a simple spectral phase-coherence feature analysis
for sgRNA efficiency prediction using the Z Framework methodology. The test
examines phase-coherence correlation with on-target efficiency within
GC-content quartile bins.

Mathematical Foundation:
- Phase-coherence: PC = (1/L) * |∑(n=1 to L) w_n * e^(2πin/L)|
- Complex base mapping: A=1+0j, T=-1+0j, C=0+1j, G=0-1j
- Statistical validation: Bootstrap confidence intervals with multiple testing correction

Dependencies: Python stdlib only (csv, math, cmath, random, statistics)
"""

import csv
import math
import cmath
import random
import statistics
from typing import List, Tuple, Dict, Any


def complex_base_mapping(base: str) -> complex:
    """Map DNA bases to complex numbers for phase-coherence analysis.
    
    Args:
        base: Single DNA base character (A, T, C, G)
        
    Returns:
        Complex number representation of the base
    """
    mapping = {
        'A': 1 + 0j,
        'T': -1 + 0j, 
        'C': 0 + 1j,
        'G': 0 - 1j
    }
    return mapping.get(base.upper(), 1 + 0j)  # Default to A for unknown bases


def calculate_phase_coherence(sequence: str) -> float:
    """Calculate phase-coherence feature for a DNA sequence.
    
    Phase-coherence quantifies the spectral organization of base composition
    using complex exponential weighting across sequence positions.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Phase-coherence value (0.0 to 1.0)
    """
    if not sequence:
        return 0.0
        
    L = len(sequence)
    if L == 0:
        return 0.0
    
    # Calculate weighted sum with complex exponential phase factors
    coherence_sum = 0.0 + 0.0j
    for n in range(L):
        base = sequence[n]
        w_n = complex_base_mapping(base)
        phase_factor = cmath.exp(2j * math.pi * (n + 1) / L)
        coherence_sum += w_n * phase_factor
    
    # Return magnitude normalized by sequence length
    return abs(coherence_sum) / L


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content fraction for a sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as fraction (0.0 to 1.0)
    """
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return gc_count / len(sequence)


def pearson_correlation(x: List[float], y: List[float]) -> float:
    """Calculate Pearson correlation coefficient.
    
    Args:
        x: First variable values
        y: Second variable values
        
    Returns:
        Pearson correlation coefficient (-1.0 to 1.0)
    """
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    
    n = len(x)
    mean_x = statistics.mean(x)
    mean_y = statistics.mean(y)
    
    # Calculate numerator and denominators
    numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    sum_sq_x = sum((x[i] - mean_x) ** 2 for i in range(n))
    sum_sq_y = sum((y[i] - mean_y) ** 2 for i in range(n))
    
    denominator = math.sqrt(sum_sq_x * sum_sq_y)
    
    if denominator == 0:
        return 0.0
    
    return numerator / denominator


def bootstrap_correlation(x: List[float], y: List[float], n_boot: int = 1000, 
                         seed: int = 42) -> Tuple[float, float, float]:
    """Calculate bootstrap confidence interval for correlation coefficient.
    
    Args:
        x: First variable values
        y: Second variable values  
        n_boot: Number of bootstrap samples
        seed: Random seed for reproducibility
        
    Returns:
        Tuple of (correlation, ci_low, ci_high)
    """
    if len(x) != len(y) or len(x) < 2:
        return 0.0, 0.0, 0.0
    
    random.seed(seed)
    n = len(x)
    
    # Original correlation
    original_r = pearson_correlation(x, y)
    
    # Bootstrap resampling
    boot_correlations = []
    for _ in range(n_boot):
        # Resample with replacement
        indices = [random.randint(0, n - 1) for _ in range(n)]
        boot_x = [x[i] for i in indices]
        boot_y = [y[i] for i in indices]
        
        boot_r = pearson_correlation(boot_x, boot_y)
        boot_correlations.append(boot_r)
    
    # Calculate 95% confidence interval
    boot_correlations.sort()
    ci_low = boot_correlations[int(0.025 * n_boot)]
    ci_high = boot_correlations[int(0.975 * n_boot)]
    
    return original_r, ci_low, ci_high


def benjamini_hochberg_correction(p_values: List[float], q: float = 0.10) -> List[bool]:
    """Apply Benjamini-Hochberg FDR correction for multiple testing.
    
    Args:
        p_values: List of p-values to correct
        q: False discovery rate threshold
        
    Returns:
        List of boolean values indicating significance after correction
    """
    if not p_values:
        return []
    
    n = len(p_values)
    # Create list of (p_value, original_index) and sort by p-value
    indexed_p = [(p_values[i], i) for i in range(n)]
    indexed_p.sort()
    
    # Apply BH correction
    significant = [False] * n
    for rank, (p_val, orig_idx) in enumerate(indexed_p):
        critical_value = (rank + 1) / n * q
        if p_val <= critical_value:
            significant[orig_idx] = True
        else:
            break  # All subsequent p-values will also fail
    
    return significant


def load_doench_data(filepath: str = "doench_2016.csv") -> List[Dict[str, Any]]:
    """Load Doench 2016 sgRNA efficiency data.
    
    Args:
        filepath: Path to CSV file with sequence and efficiency columns
        
    Returns:
        List of dictionaries with sequence data and calculated features
    """
    data = []
    
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sequence = row['sequence'].strip().upper()
                try:
                    efficiency = float(row['efficiency'])
                except (ValueError, KeyError):
                    continue
                
                # Calculate features
                gc_content = calculate_gc_content(sequence)
                phase_coherence = calculate_phase_coherence(sequence)
                
                data.append({
                    'sequence': sequence,
                    'efficiency': efficiency,
                    'gc_content': gc_content,
                    'phase_coherence': phase_coherence,
                    'length': len(sequence)
                })
                
    except FileNotFoundError:
        print(f"Error: Could not find data file '{filepath}'")
        return []
    
    return data


def perform_bin_resonance_analysis(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Perform bin-resonance analysis on sgRNA efficiency data.
    
    Splits data into GC content quartiles and analyzes phase-coherence 
    correlation with efficiency within each bin using bootstrap confidence
    intervals and multiple testing correction.
    
    Args:
        data: List of sequence data with calculated features
        
    Returns:
        List of analysis results for each GC quartile bin
    """
    if not data:
        return []
    
    # Extract GC contents for quartile calculation
    gc_contents = [item['gc_content'] for item in data]
    gc_contents.sort()
    
    n = len(gc_contents)
    q1_threshold = gc_contents[n // 4]
    q2_threshold = gc_contents[n // 2] 
    q3_threshold = gc_contents[3 * n // 4]
    
    # Define quartile bins
    bins = [
        ('Q1', lambda gc: gc <= q1_threshold),
        ('Q2', lambda gc: q1_threshold < gc <= q2_threshold),
        ('Q3', lambda gc: q2_threshold < gc <= q3_threshold),
        ('Q4', lambda gc: gc > q3_threshold)
    ]
    
    results = []
    p_values = []
    
    for bin_name, bin_filter in bins:
        # Filter data for this bin
        bin_data = [item for item in data if bin_filter(item['gc_content'])]
        
        if len(bin_data) < 3:  # Need minimum samples for meaningful analysis
            results.append({
                'bin': bin_name,
                'n': len(bin_data),
                'r': 0.0,
                'ci_low': 0.0,
                'ci_high': 0.0,
                'p_boot': 1.0,
                'passed': False
            })
            p_values.append(1.0)
            continue
        
        # Extract features for correlation analysis
        phase_coherences = [item['phase_coherence'] for item in bin_data]
        efficiencies = [item['efficiency'] for item in bin_data]
        
        # Calculate correlation with bootstrap CI
        r, ci_low, ci_high = bootstrap_correlation(phase_coherences, efficiencies)
        
        # Calculate p-value (proportion of bootstrap samples <= 0)
        boot_rs = []
        n_samples = len(bin_data)
        for _ in range(1000):
            indices = [random.randint(0, n_samples - 1) for _ in range(n_samples)]
            boot_pc = [phase_coherences[i] for i in indices]
            boot_eff = [efficiencies[i] for i in indices]
            boot_rs.append(pearson_correlation(boot_pc, boot_eff))
        
        p_boot = sum(1 for br in boot_rs if br <= r) / len(boot_rs)
        p_values.append(p_boot)
        
        results.append({
            'bin': bin_name,
            'n': len(bin_data),
            'r': r,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'p_boot': p_boot,
            'passed': False  # Will be updated after FDR correction
        })
    
    # Apply multiple testing correction
    fdr_significant = benjamini_hochberg_correction(p_values, q=0.10)
    
    # Update passed status based on success criteria
    for i, result in enumerate(results):
        r = result['r']
        ci_low = result['ci_low']
        
        # Success criteria: r > 0 AND ci_low > 0 AND r >= 0.10 AND FDR significant
        success_criteria = (
            r > 0 and 
            ci_low > 0 and 
            r >= 0.10 and 
            fdr_significant[i]
        )
        
        result['passed'] = success_criteria
    
    return results


def sanity_check_with_shuffled_labels(data: List[Dict[str, Any]]) -> float:
    """Perform sanity check by analyzing correlation with shuffled efficiency labels.
    
    Args:
        data: List of sequence data with calculated features
        
    Returns:
        Global correlation with shuffled labels (should be ~0)
    """
    if not data:
        return 0.0
    
    phase_coherences = [item['phase_coherence'] for item in data]
    efficiencies = [item['efficiency'] for item in data]
    
    # Shuffle efficiency labels
    random.seed(42)
    shuffled_efficiencies = efficiencies.copy()
    random.shuffle(shuffled_efficiencies)
    
    return pearson_correlation(phase_coherences, shuffled_efficiencies)


def save_results_csv(results: List[Dict[str, Any]], filepath: str = "bin_resonance_results.csv"):
    """Save analysis results to CSV file.
    
    Args:
        results: List of analysis results from perform_bin_resonance_analysis
        filepath: Output CSV file path
    """
    if not results:
        print("No results to save.")
        return
    
    fieldnames = ['bin', 'n', 'r', 'ci_low', 'ci_high', 'p_boot', 'passed']
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in results:
            # Format numerical values for cleaner output
            formatted_result = {
                'bin': result['bin'],
                'n': result['n'],
                'r': f"{result['r']:.4f}",
                'ci_low': f"{result['ci_low']:.4f}",
                'ci_high': f"{result['ci_high']:.4f}",
                'p_boot': f"{result['p_boot']:.4f}",
                'passed': result['passed']
            }
            writer.writerow(formatted_result)
    
    print(f"Results saved to {filepath}")


def main():
    """Main execution function for bin-resonance test."""
    print("Minimal Bin-Resonance Test (Human CRISPR)")
    print("=" * 50)
    
    # Load data
    print("Loading Doench 2016 data...")
    data = load_doench_data()
    
    if not data:
        print("Error: No data loaded. Exiting.")
        return
    
    print(f"Loaded {len(data)} sequences")
    
    # Perform analysis
    print("Performing bin-resonance analysis...")
    results = perform_bin_resonance_analysis(data)
    
    # Display results
    print("\nResults by GC Quartile:")
    print("-" * 60)
    for result in results:
        status = "PASS" if result['passed'] else "FAIL"
        print(f"Bin {result['bin']}: n={result['n']}, r={result['r']:.4f} "
              f"[{result['ci_low']:.4f}, {result['ci_high']:.4f}], "
              f"p={result['p_boot']:.4f}, {status}")
    
    # Sanity check
    print("\nSanity check with shuffled labels...")
    shuffled_r = sanity_check_with_shuffled_labels(data)
    print(f"Global correlation with shuffled labels: {shuffled_r:.4f} (should be ~0)")
    
    # Save results
    save_results_csv(results)
    
    # Overall assessment
    any_passed = any(result['passed'] for result in results)
    print(f"\nOverall result: {'PASS' if any_passed else 'FAIL'}")
    if any_passed:
        print("At least one GC quartile shows significant positive correlation.")
    else:
        print("No GC quartile meets the success criteria.")


if __name__ == "__main__":
    main()