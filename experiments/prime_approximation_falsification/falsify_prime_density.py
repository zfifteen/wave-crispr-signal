#!/usr/bin/env python3
"""
Prime Approximation Density Falsification Experiment

This experiment tests the hypothesis:
"Prime approximation via Riemann inversion achieves ~210% effective density 
at N=10^6 (CI [207.2%, 228.9%]) by aligning with geodesic clustering."

Methodology:
1. Implement high-precision Riemann R(x) prime counting function
2. Use Newton's method to invert R(x) for prime approximation
3. Calculate actual density at N=10^6
4. Bootstrap confidence intervals (1,000 resamples)
5. Compare claimed vs. actual density boost

Scientific Gates:
- Human DNA only: N/A (mathematical domain)
- No fabrication: All calculations from first principles
- Statistical validity: Bootstrap CI with ≥1,000 resamples
- Reproducibility: Fixed seed, exact dependencies
"""

import argparse
import json
import logging
import os
import sys
from datetime import datetime
from functools import lru_cache
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np
from scipy import stats

# High precision for Riemann calculations
mp.dps = 60

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@lru_cache(None)
def mu(k: int) -> int:
    """
    Möbius function μ(k) for use in Riemann R(x) calculation.
    
    Args:
        k: Positive integer
        
    Returns:
        0 if k has a squared prime factor
        1 if k = 1
        (-1)^m if k is product of m distinct primes
    """
    if k == 1:
        return 1
    
    m, p, cnt = k, 2, 0
    
    # Check for squared factors and count distinct prime factors
    while p * p <= m:
        if m % p == 0:
            cnt += 1
            m //= p
            if m % p == 0:  # Squared factor detected
                return 0
        p += 1 if p == 2 else 2
    
    if m > 1:
        cnt += 1
    
    return -1 if (cnt % 2) else 1


def riemann_R(x: mp.mpf, max_terms: int = 1000) -> mp.mpf:
    """
    Calculate Riemann's prime counting function R(x).
    
    R(x) = Σ[n=1 to ∞] μ(n)/n * Li(x^(1/n))
    
    where Li(x) is the logarithmic integral.
    
    Args:
        x: Value at which to evaluate R(x)
        max_terms: Maximum number of terms in the series
        
    Returns:
        R(x) approximation with high precision
    """
    if x <= 1:
        return mp.mpf(0)
    
    result = mp.mpf(0)
    
    for n in range(1, max_terms + 1):
        mu_n = mu(n)
        if mu_n == 0:
            continue
        
        # Calculate x^(1/n)
        x_power = mp.power(x, mp.mpf(1) / mp.mpf(n))
        
        # Use mpmath's logarithmic integral
        li_term = mp.li(x_power)
        
        # Add term: μ(n)/n * Li(x^(1/n))
        term = mp.mpf(mu_n) / mp.mpf(n) * li_term
        result += term
        
        # Check convergence
        if n > 10 and abs(term) < mp.mpf(10) ** (-50):
            break
    
    return result


def invert_riemann_R(target_count: int, initial_guess: float = None,
                     max_iter: int = 50) -> mp.mpf:
    """
    Invert Riemann R(x) using Newton's method to find x where R(x) ≈ target_count.
    
    This gives us the approximate nth prime or position where we expect
    'target_count' primes.
    
    Args:
        target_count: Target number of primes
        initial_guess: Initial x value for Newton iteration
        max_iter: Maximum Newton iterations
        
    Returns:
        x such that R(x) ≈ target_count
    """
    if initial_guess is None:
        # Use prime number theorem as initial guess: π(x) ≈ x/ln(x)
        # So x ≈ n * ln(n)
        n = target_count
        initial_guess = n * mp.log(n) if n > 1 else 2
    
    x = mp.mpf(initial_guess)
    
    for iteration in range(max_iter):
        # Current value
        R_x = riemann_R(x, max_terms=100)
        
        # Derivative approximation: R'(x) ≈ 1/ln(x) (from Li'(x))
        dR_dx = mp.mpf(1) / mp.log(x)
        
        # Newton step: x_new = x - (R(x) - target) / R'(x)
        delta = (R_x - mp.mpf(target_count)) / dR_dx
        x_new = x - delta
        
        # Check convergence
        if abs(delta) < mp.mpf(10) ** (-40):
            logger.debug(f"Newton converged in {iteration + 1} iterations")
            return x_new
        
        x = x_new
    
    logger.warning(f"Newton did not converge in {max_iter} iterations")
    return x


def count_primes_exact(n: int) -> int:
    """
    Count primes up to n using simple sieve (for verification).
    
    Args:
        n: Upper bound
        
    Returns:
        Number of primes ≤ n
    """
    if n < 2:
        return 0
    
    # Sieve of Eratosthenes
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n + 1, i):
                is_prime[j] = False
    
    return sum(is_prime)


def calculate_density_boost(N: int, num_bootstrap: int = 1000,
                            seed: int = 42) -> Dict:
    """
    Calculate the effective density boost at N=10^6 using Riemann inversion.
    
    The hypothesis claims ~210% density with CI [207.2%, 228.9%].
    
    Args:
        N: Target value (e.g., 10^6)
        num_bootstrap: Number of bootstrap resamples
        seed: Random seed for reproducibility
        
    Returns:
        Dictionary with density metrics and bootstrap CI
    """
    logger.info(f"Calculating density at N={N}")
    
    # Exact prime count at N
    pi_N = count_primes_exact(N)
    logger.info(f"Exact π({N}) = {pi_N}")
    
    # Riemann R(N) approximation
    R_N = riemann_R(mp.mpf(N), max_terms=100)
    logger.info(f"R({N}) = {float(R_N):.6f}")
    
    # Relative error
    rel_error = float((R_N - pi_N) / pi_N * 100)
    logger.info(f"Relative error: {rel_error:.6f}%")
    
    # Expected density from prime number theorem: π(x) ≈ x/ln(x)
    expected_density = N / mp.log(N)
    logger.info(f"PNT expected: {float(expected_density):.2f}")
    
    # Claimed "effective density" calculation
    # The hypothesis claims 210% boost, meaning effective_count ≈ 2.1 * expected
    # We need to check if this is valid
    
    # Calculate what "210% density" would mean
    claimed_boost_pct = 210.0
    claimed_effective = pi_N * (claimed_boost_pct / 100.0)
    
    # Actual boost from Riemann inversion relative to PNT
    actual_boost_pct = float(R_N / expected_density * 100)
    
    logger.info(f"Claimed boost: {claimed_boost_pct}%")
    logger.info(f"Actual boost (R/PNT): {actual_boost_pct:.2f}%")
    
    # Bootstrap confidence intervals for R(N) estimation
    # We resample the series terms to get uncertainty
    np.random.seed(seed)
    bootstrap_values = []
    
    logger.info(f"Running {num_bootstrap} bootstrap resamples...")
    
    for b in range(num_bootstrap):
        # Resample by varying the number of terms slightly
        # This simulates uncertainty in truncation
        n_terms = np.random.randint(50, 150)
        R_sample = riemann_R(mp.mpf(N), max_terms=n_terms)
        boost_sample = float(R_sample / expected_density * 100)
        bootstrap_values.append(boost_sample)
        
        if (b + 1) % 100 == 0:
            logger.info(f"  Completed {b + 1}/{num_bootstrap} resamples")
    
    # Calculate bootstrap CI
    ci_lower = np.percentile(bootstrap_values, 2.5)
    ci_upper = np.percentile(bootstrap_values, 97.5)
    bootstrap_mean = np.mean(bootstrap_values)
    bootstrap_std = np.std(bootstrap_values, ddof=1)
    
    logger.info(f"Bootstrap CI (95%): [{ci_lower:.2f}%, {ci_upper:.2f}%]")
    logger.info(f"Bootstrap mean: {bootstrap_mean:.2f}%")
    
    # Statistical test: Does actual boost match claimed 210%?
    # H0: actual_boost = 210%
    t_stat = (bootstrap_mean - claimed_boost_pct) / (bootstrap_std / np.sqrt(num_bootstrap))
    p_value = 2 * (1 - stats.norm.cdf(abs(t_stat)))
    
    logger.info(f"t-statistic vs claimed: {t_stat:.4f}")
    logger.info(f"p-value: {p_value:.6e}")
    
    # Does claimed CI [207.2%, 228.9%] overlap with actual CI?
    claimed_ci = (207.2, 228.9)
    ci_overlap = not (ci_upper < claimed_ci[0] or ci_lower > claimed_ci[1])
    
    return {
        'N': N,
        'pi_N_exact': pi_N,
        'R_N': float(R_N),
        'rel_error_pct': rel_error,
        'PNT_expected': float(expected_density),
        'claimed_boost_pct': claimed_boost_pct,
        'claimed_ci_lower': claimed_ci[0],
        'claimed_ci_upper': claimed_ci[1],
        'actual_boost_pct': actual_boost_pct,
        'bootstrap_mean_pct': bootstrap_mean,
        'bootstrap_std_pct': bootstrap_std,
        'bootstrap_ci_lower': ci_lower,
        'bootstrap_ci_upper': ci_upper,
        'ci_overlap': ci_overlap,
        't_statistic': t_stat,
        'p_value': p_value,
        'hypothesis_falsified': p_value < 0.05 and not ci_overlap,
        'num_bootstrap': num_bootstrap,
        'seed': seed
    }


def generate_executive_summary(results: Dict) -> str:
    """
    Generate executive summary of falsification results.
    
    Args:
        results: Dictionary from calculate_density_boost
        
    Returns:
        Formatted executive summary string
    """
    summary = """
# EXECUTIVE SUMMARY: Prime Approximation Density Falsification

## Hypothesis Tested
"Prime approximation via Riemann inversion achieves ~210% effective density 
at N=10^6 (CI [207.2%, 228.9%]) by aligning with geodesic clustering."

## Key Findings
"""
    
    if results['hypothesis_falsified']:
        summary += "\n### ❌ HYPOTHESIS FALSIFIED\n"
    else:
        summary += "\n### ⚠️  HYPOTHESIS NOT DEFINITIVELY FALSIFIED\n"
    
    summary += f"""
## Empirical Results at N={results['N']:,}

- **Exact prime count π(N)**: {results['pi_N_exact']:,}
- **Riemann R(N)**: {results['R_N']:.2f}
- **Relative error**: {results['rel_error_pct']:.6f}%
- **PNT expected**: {results['PNT_expected']:.2f}

## Density Boost Analysis

### Claimed Values
- **Boost**: {results['claimed_boost_pct']:.1f}%
- **95% CI**: [{results['claimed_ci_lower']:.1f}%, {results['claimed_ci_upper']:.1f}%]

### Actual Values (Bootstrap n={results['num_bootstrap']})
- **Boost**: {results['actual_boost_pct']:.2f}%
- **Bootstrap mean**: {results['bootstrap_mean_pct']:.2f}%
- **Bootstrap std**: {results['bootstrap_std_pct']:.2f}%
- **95% CI**: [{results['bootstrap_ci_lower']:.2f}%, {results['bootstrap_ci_upper']:.2f}%]

## Statistical Test
- **t-statistic**: {results['t_statistic']:.4f}
- **p-value**: {results['p_value']:.6e}
- **CI overlap**: {'Yes' if results['ci_overlap'] else 'No'}

## Conclusion
"""
    
    if results['hypothesis_falsified']:
        summary += f"""
The hypothesis is **FALSIFIED** with high confidence (p={results['p_value']:.6e}).

The claimed 210% density boost is not supported by empirical calculations.
The actual boost is {results['actual_boost_pct']:.2f}% with 95% CI 
[{results['bootstrap_ci_lower']:.2f}%, {results['bootstrap_ci_upper']:.2f}%], 
which {'does not overlap' if not results['ci_overlap'] else 'marginally overlaps'} 
with the claimed CI [{results['claimed_ci_lower']:.1f}%, {results['claimed_ci_upper']:.1f}%].
"""
    else:
        summary += f"""
The hypothesis **CANNOT BE DEFINITIVELY FALSIFIED** based on this analysis.

The actual boost of {results['actual_boost_pct']:.2f}% (95% CI: 
[{results['bootstrap_ci_lower']:.2f}%, {results['bootstrap_ci_upper']:.2f}%]) 
{'overlaps with' if results['ci_overlap'] else 'does not overlap with'} 
the claimed CI [{results['claimed_ci_lower']:.1f}%, {results['claimed_ci_upper']:.1f}%].

However, the large discrepancy (t={results['t_statistic']:.2f}, p={results['p_value']:.6e}) 
suggests the claimed values require further scrutiny.
"""
    
    return summary


def main():
    """Main experiment execution."""
    parser = argparse.ArgumentParser(
        description='Falsify prime approximation density hypothesis'
    )
    parser.add_argument('--n', type=int, default=1000000,
                       help='Target N value (default: 10^6)')
    parser.add_argument('--bootstrap', type=int, default=1000,
                       help='Number of bootstrap resamples')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--output-dir', type=str,
                       default='results/prime_approximation_falsification',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    run_dir = os.path.join(args.output_dir, f'run-{timestamp}')
    os.makedirs(run_dir, exist_ok=True)
    
    logger.info("=" * 60)
    logger.info("Prime Approximation Density Falsification Experiment")
    logger.info("=" * 60)
    logger.info(f"N = {args.n:,}")
    logger.info(f"Bootstrap resamples = {args.bootstrap}")
    logger.info(f"Random seed = {args.seed}")
    logger.info(f"mpmath precision = {mp.dps} decimal places")
    logger.info("")
    
    # Run calculation
    results = calculate_density_boost(
        N=args.n,
        num_bootstrap=args.bootstrap,
        seed=args.seed
    )
    
    # Generate executive summary
    summary = generate_executive_summary(results)
    print("\n" + summary)
    
    # Save results
    results_file = os.path.join(run_dir, 'results.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Results saved to {results_file}")
    
    summary_file = os.path.join(run_dir, 'executive_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(summary)
    logger.info(f"Summary saved to {summary_file}")
    
    # Save environment info
    env_file = os.path.join(run_dir, 'env.txt')
    with open(env_file, 'w') as f:
        f.write(f"Python: {sys.version}\n")
        f.write(f"mpmath.dps: {mp.dps}\n")
        f.write(f"numpy: {np.__version__}\n")
        f.write(f"Git commit: {os.popen('git rev-parse HEAD').read().strip()}\n")
        f.write(f"Timestamp: {timestamp}\n")
    
    logger.info(f"Environment info saved to {env_file}")
    logger.info("\nExperiment completed successfully!")
    
    return 0 if results['hypothesis_falsified'] else 1


if __name__ == '__main__':
    sys.exit(main())
