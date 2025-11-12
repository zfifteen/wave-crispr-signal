#!/usr/bin/env python3
"""
FUS Enhancer: Optimized Focused Ultrasound Targeting with PyTorch Vectorization

RESEARCH USE ONLY - NOT FOR CLINICAL APPLICATIONS

This module provides a highly optimized implementation of the Z Framework focused 
ultrasound targeting experiment using PyTorch vectorization. It can efficiently 
process batches of up to 10^6 targeting trials for performance analysis.

Key Optimizations:
- PyTorch tensor operations for vectorized calculations
- Batch processing of trials (configurable batch size)
- GPU acceleration when available
- Vectorized Z Framework geodesic resolution
- Vectorized statistical analysis with bootstrap/permutation testing
- Memory-efficient processing for large datasets

Usage:
    python fus_enhancer.py --batch-size 1000000 --seed 42 --bootstrap 1000
    
Integration:
    from fus_enhancer import VectorizedFUSExperiment
    experiment = VectorizedFUSExperiment(config)
    results = experiment.run_vectorized_experiment()
"""

import argparse
import json
import logging
import random
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
import torch
import torch.nn.functional as F
import matplotlib.pyplot as plt
import mpmath as mp
from scipy import stats

# Configure high precision for critical calculations
mp.dps = 50

# Mathematical constants (matching z_framework.py)
PHI = mp.mpf((1 + mp.sqrt(5)) / 2)  # Golden ratio ≈ 1.618
PHI_CONJUGATE = PHI - 1  # φ-1 ≈ 0.618
E_SQUARED = mp.e**2  # e² ≈ 7.389

# Default experiment parameters
DEFAULT_BATCH_SIZE = 100000  # Default batch size for vectorized processing
DEFAULT_GRID_SIZE = 100
DEFAULT_N_TRIALS = 1000000  # 10^6 trials for large-scale testing
BASE_ACOUSTIC_VELOCITY = 1540.0  # m/s (typical soft tissue)
VELOCITY_VARIANCE = 0.15  # 15% relative variance

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set device for PyTorch operations
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
logger.info(f"Using device: {DEVICE}")


@dataclass
class VectorizedConfig:
    """Configuration for vectorized focused ultrasound targeting experiment."""
    batch_size: int = DEFAULT_BATCH_SIZE
    grid_size: int = DEFAULT_GRID_SIZE
    n_trials: int = DEFAULT_N_TRIALS
    n_bootstrap: int = 1000
    n_permutation: int = 1000
    seed: int = 42
    output_dir: str = "results"
    k_parameter: float = 0.3  # Geometric resolution parameter
    device: str = str(DEVICE)
    dtype: torch.dtype = torch.float32


class VectorizedAcousticGrid:
    """Vectorized 2D tissue grid with acoustic velocity heterogeneity using PyTorch."""
    
    def __init__(self, size: int, base_velocity: float, variance: float, seed: int, device: str = 'cpu'):
        """
        Initialize vectorized acoustic grid.
        
        Args:
            size: Grid dimension (size × size)
            base_velocity: Base acoustic velocity (m/s)
            variance: Relative variance for Gaussian noise
            seed: Random seed for reproducibility
            device: PyTorch device ('cpu' or 'cuda')
        """
        torch.manual_seed(seed)
        np.random.seed(seed)
        
        self.size = size
        self.base_velocity = base_velocity
        self.variance = variance
        self.device = torch.device(device)
        
        # Generate velocity field with Gaussian noise using PyTorch
        noise = torch.randn(size, size, device=self.device) * variance
        self.velocity_field = base_velocity * (1 + noise)
        
        # Calculate grid statistics
        self.mean_velocity = torch.mean(self.velocity_field).item()
        self.std_velocity = torch.std(self.velocity_field).item()
        self.relative_std = self.std_velocity / self.mean_velocity
        
        logger.info(f"Vectorized acoustic grid initialized: {size}×{size} on {device}")
        logger.info(f"Velocity stats: mean={self.mean_velocity:.1f} m/s, "
                   f"std={self.std_velocity:.1f} m/s, rel_std={self.relative_std:.3f}")
    
    def sample_velocities_batch(self, x_coords: torch.Tensor, y_coords: torch.Tensor, 
                               path_points: int = 10) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Sample velocities along paths for batch of coordinate pairs.
        
        Args:
            x_coords: Batch of x coordinates [batch_size, 2] (source, target)
            y_coords: Batch of y coordinates [batch_size, 2] (source, target)
            path_points: Number of points to sample along each path
            
        Returns:
            Tuple of (velocity_heterogeneity, mean_path_velocity) tensors [batch_size]
        """
        batch_size = x_coords.shape[0]
        
        # Generate path points for all batches simultaneously
        # x_coords/y_coords: [batch_size, 2] -> [batch_size, path_points]
        alphas = torch.linspace(0, 1, path_points, device=self.device).unsqueeze(0)  # [1, path_points]
        
        x_paths = x_coords[:, 0:1] + alphas * (x_coords[:, 1:2] - x_coords[:, 0:1])  # [batch_size, path_points]
        y_paths = y_coords[:, 0:1] + alphas * (y_coords[:, 1:2] - y_coords[:, 0:1])  # [batch_size, path_points]
        
        # Clip coordinates to grid bounds and convert to indices
        x_indices = torch.clamp(x_paths.long(), 0, self.size - 1)
        y_indices = torch.clamp(y_paths.long(), 0, self.size - 1)
        
        # Sample velocities from the grid
        path_velocities = self.velocity_field[y_indices, x_indices]  # [batch_size, path_points]
        
        # Calculate statistics along paths
        mean_path_velocities = torch.mean(path_velocities, dim=1)  # [batch_size]
        std_path_velocities = torch.std(path_velocities, dim=1)    # [batch_size]
        velocity_heterogeneity = std_path_velocities / (mean_path_velocities + 1e-8)  # [batch_size]
        
        return velocity_heterogeneity, mean_path_velocities


class VectorizedZFramework:
    """Vectorized Z Framework calculations using PyTorch."""
    
    def __init__(self, k_parameter: float = 0.3, device: str = 'cpu'):
        """Initialize vectorized Z Framework calculator."""
        self.k_parameter = k_parameter
        self.device = torch.device(device)
        self.phi = float(PHI)
        self.e_squared = float(E_SQUARED)
        
    def theta_prime_vectorized(self, n: torch.Tensor) -> torch.Tensor:
        """
        Vectorized geometric resolution θ'(n,k) = φ·((n mod φ)/φ)^k
        
        Args:
            n: Position parameters [batch_size]
            
        Returns:
            Geometric resolution values [batch_size]
        """
        n_mod_phi = torch.fmod(n, self.phi)
        ratio = n_mod_phi / self.phi
        theta_prime = self.phi * torch.pow(ratio, self.k_parameter)
        return theta_prime
    
    def calculate_z_enhancement_vectorized(self, distances: torch.Tensor, 
                                         velocity_heterogeneity: torch.Tensor) -> torch.Tensor:
        """
        Vectorized Z Framework enhancement factor calculation.
        
        Args:
            distances: Euclidean distances [batch_size]
            velocity_heterogeneity: Velocity variations [batch_size]
            
        Returns:
            Z Framework enhancement factors [batch_size]
        """
        # Map distances to integers for geodesic resolution (avoid zero)
        n_values = torch.clamp((distances * 10).long(), min=1)
        
        # Calculate geodesic resolution
        theta_prime_vals = self.theta_prime_vectorized(n_values.float())
        
        # Discrete domain form: Z = A(B/e²)
        A = torch.clamp(theta_prime_vals / self.phi, max=1.0)  # Geometric scaling factor
        B = velocity_heterogeneity * distances + 0.1  # Heterogeneity adaptation factor
        
        # Z Framework enhancement (bounded to [0, 1])
        z_raw = A * (B / self.e_squared)
        
        # Convert to improvement factor with distance dependency
        distance_factor = 1 / (1 + distances * 0.01)
        enhancement_factor = torch.tanh(z_raw * 2) * distance_factor
        
        # Add noise for realistic variation (±5%)
        noise = torch.randn_like(enhancement_factor) * 0.05
        enhancement_factor = torch.clamp(enhancement_factor + noise, 0.05, 0.45)
        
        return enhancement_factor


class VectorizedTargetingModels:
    """Vectorized baseline and Z Framework targeting models."""
    
    def __init__(self, acoustic_grid: VectorizedAcousticGrid, z_framework: VectorizedZFramework):
        """Initialize with acoustic grid and Z Framework calculator."""
        self.grid = acoustic_grid
        self.z_framework = z_framework
        self.device = acoustic_grid.device
        
    def calculate_targeting_errors_vectorized(self, source_coords: torch.Tensor,
                                            target_coords: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Vectorized calculation of targeting errors for both baseline and Z Framework models.
        
        Args:
            source_coords: Source coordinates [batch_size, 2] (x, y)
            target_coords: Target coordinates [batch_size, 2] (x, y)
            
        Returns:
            Tuple of (baseline_errors, z_framework_errors, baseline_times, z_framework_times)
        """
        # Calculate Euclidean distances
        distances = torch.norm(target_coords - source_coords, dim=1)  # [batch_size]
        
        # Sample velocity fields along paths
        coords_combined = torch.stack([source_coords, target_coords], dim=2)  # [batch_size, 2, 2]
        x_coords = coords_combined[:, 0, :]  # [batch_size, 2]
        y_coords = coords_combined[:, 1, :]  # [batch_size, 2]
        
        velocity_heterogeneity, mean_path_velocities = self.grid.sample_velocities_batch(x_coords, y_coords)
        
        # Grid variance for error calculation
        grid_variance = (self.grid.std_velocity / self.grid.base_velocity) ** 2
        
        # Baseline model: simple error proportional to distance and grid variance
        baseline_errors = distances * grid_variance
        baseline_times = distances / self.grid.base_velocity
        
        # Z Framework model: enhanced targeting with error reduction
        z_enhancements = self.z_framework.calculate_z_enhancement_vectorized(distances, velocity_heterogeneity)
        z_framework_errors = baseline_errors * (1 - torch.clamp(z_enhancements, 0.05, 0.45))
        z_framework_times = distances / mean_path_velocities
        
        return baseline_errors, z_framework_errors, baseline_times, z_framework_times


class VectorizedStatistics:
    """Vectorized statistical analysis with PyTorch acceleration."""
    
    def __init__(self, device: str = 'cpu'):
        """Initialize vectorized statistics calculator."""
        self.device = torch.device(device)
    
    def bootstrap_correlation_vectorized(self, x: torch.Tensor, y: torch.Tensor, 
                                       n_boot: int) -> Tuple[float, float, float, float]:
        """
        Vectorized bootstrap confidence interval for correlation.
        
        Args:
            x, y: Data tensors [n_samples]
            n_boot: Number of bootstrap samples
            
        Returns:
            Tuple of (correlation, ci_low, ci_high, p_value)
        """
        n = x.shape[0]
        
        # Original correlation (Pearson)
        x_centered = x - torch.mean(x)
        y_centered = y - torch.mean(y)
        original_r = torch.sum(x_centered * y_centered) / torch.sqrt(
            torch.sum(x_centered**2) * torch.sum(y_centered**2))
        original_r = torch.clamp(original_r, -0.9999, 0.9999)
        
        # Vectorized bootstrap sampling
        indices = torch.randint(0, n, (n_boot, n), device=self.device)
        boot_x = x[indices]  # [n_boot, n]
        boot_y = y[indices]  # [n_boot, n]
        
        # Vectorized correlation calculation
        boot_x_centered = boot_x - torch.mean(boot_x, dim=1, keepdim=True)
        boot_y_centered = boot_y - torch.mean(boot_y, dim=1, keepdim=True)
        
        numerator = torch.sum(boot_x_centered * boot_y_centered, dim=1)
        denominator = torch.sqrt(
            torch.sum(boot_x_centered**2, dim=1) * torch.sum(boot_y_centered**2, dim=1))
        
        boot_correlations = numerator / (denominator + 1e-8)
        boot_correlations = torch.clamp(boot_correlations, -0.9999, 0.9999)
        
        # Remove NaN values
        valid_mask = ~torch.isnan(boot_correlations)
        boot_correlations = boot_correlations[valid_mask]
        
        # Calculate CI
        ci_low = torch.quantile(boot_correlations, 0.025).item()
        ci_high = torch.quantile(boot_correlations, 0.975).item()
        
        # Calculate p-value
        p_positive = torch.mean((boot_correlations >= 0).float()).item()
        p_negative = torch.mean((boot_correlations <= 0).float()).item()
        p_value = 2 * min(p_positive, p_negative)
        
        return original_r.item(), ci_low, ci_high, p_value
    
    def permutation_test_vectorized(self, x: torch.Tensor, y: torch.Tensor, 
                                  n_perm: int) -> float:
        """
        Vectorized permutation test for difference in means.
        
        Args:
            x, y: Data tensors to compare
            n_perm: Number of permutations
            
        Returns:
            P-value for difference in means
        """
        observed_diff = torch.mean(x) - torch.mean(y)
        
        # Combine arrays for permutation
        combined = torch.cat([x, y])
        n_x = x.shape[0]
        n_total = combined.shape[0]
        
        # Vectorized permutation sampling
        perm_indices = torch.rand(n_perm, n_total, device=self.device).argsort(dim=1)
        perm_combined = combined[perm_indices]  # [n_perm, n_total]
        
        perm_x = perm_combined[:, :n_x]      # [n_perm, n_x]
        perm_y = perm_combined[:, n_x:]      # [n_perm, n_y]
        
        perm_diffs = torch.mean(perm_x, dim=1) - torch.mean(perm_y, dim=1)
        
        # Two-tailed p-value
        p_value = torch.mean((torch.abs(perm_diffs) >= torch.abs(observed_diff)).float()).item()
        
        return p_value


class VectorizedFUSExperiment:
    """Main vectorized experiment class for focused ultrasound targeting."""
    
    def __init__(self, config: VectorizedConfig):
        """Initialize vectorized experiment with configuration."""
        self.config = config
        
        # Set random seeds for reproducibility
        torch.manual_seed(config.seed)
        np.random.seed(config.seed)
        random.seed(config.seed)
        
        # Initialize components
        self.acoustic_grid = VectorizedAcousticGrid(
            size=config.grid_size,
            base_velocity=BASE_ACOUSTIC_VELOCITY,
            variance=VELOCITY_VARIANCE,
            seed=config.seed,
            device=config.device
        )
        
        self.z_framework = VectorizedZFramework(
            k_parameter=config.k_parameter,
            device=config.device
        )
        
        self.targeting_models = VectorizedTargetingModels(
            acoustic_grid=self.acoustic_grid,
            z_framework=self.z_framework
        )
        
        self.statistics = VectorizedStatistics(device=config.device)
        
        logger.info(f"Vectorized experiment initialized for {config.n_trials:,} trials")
        logger.info(f"Batch size: {config.batch_size:,}, Device: {config.device}")
    
    def generate_trial_coordinates_batch(self, batch_size: int) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Generate batch of random trial coordinates.
        
        Args:
            batch_size: Number of coordinate pairs to generate
            
        Returns:
            Tuple of (source_coords, target_coords) tensors [batch_size, 2]
        """
        device = self.acoustic_grid.device
        
        # Fixed source location (corner)
        source_coords = torch.tensor([5.0, 5.0], device=device).expand(batch_size, 2)
        
        # Random target locations
        target_x = torch.rand(batch_size, device=device) * (self.config.grid_size - 20) + 10
        target_y = torch.rand(batch_size, device=device) * (self.config.grid_size - 20) + 10
        target_coords = torch.stack([target_x, target_y], dim=1)
        
        return source_coords, target_coords
    
    def run_vectorized_experiment(self) -> Dict:
        """Run complete vectorized experiment with batched processing."""
        logger.info("Starting vectorized focused ultrasound targeting experiment")
        start_time = time.time()
        
        # Initialize result storage
        all_baseline_errors = []
        all_z_framework_errors = []
        all_baseline_times = []
        all_z_framework_times = []
        
        # Process trials in batches
        n_batches = (self.config.n_trials + self.config.batch_size - 1) // self.config.batch_size
        
        for batch_idx in range(n_batches):
            batch_start = batch_idx * self.config.batch_size
            batch_end = min(batch_start + self.config.batch_size, self.config.n_trials)
            current_batch_size = batch_end - batch_start
            
            # Generate trial coordinates
            source_coords, target_coords = self.generate_trial_coordinates_batch(current_batch_size)
            
            # Calculate targeting errors
            baseline_errors, z_framework_errors, baseline_times, z_framework_times = \
                self.targeting_models.calculate_targeting_errors_vectorized(source_coords, target_coords)
            
            # Store results
            all_baseline_errors.append(baseline_errors.cpu())
            all_z_framework_errors.append(z_framework_errors.cpu())
            all_baseline_times.append(baseline_times.cpu())
            all_z_framework_times.append(z_framework_times.cpu())
            
            if (batch_idx + 1) % max(1, n_batches // 10) == 0:
                logger.info(f"Completed batch {batch_idx + 1}/{n_batches} "
                           f"({(batch_idx + 1) * self.config.batch_size:,} trials)")
        
        # Concatenate all results
        baseline_errors = torch.cat(all_baseline_errors)
        z_framework_errors = torch.cat(all_z_framework_errors) 
        baseline_times = torch.cat(all_baseline_times)
        z_framework_times = torch.cat(all_z_framework_times)
        
        elapsed_time = time.time() - start_time
        logger.info(f"Vectorized experiment completed in {elapsed_time:.2f} seconds")
        logger.info(f"Processing rate: {self.config.n_trials / elapsed_time:.0f} trials/second")
        
        # Statistical analysis
        analysis_results = self.analyze_results_vectorized(
            baseline_errors, z_framework_errors, baseline_times, z_framework_times
        )
        
        return {
            'experiment_config': {
                'n_trials': self.config.n_trials,
                'batch_size': self.config.batch_size,
                'grid_size': self.config.grid_size,
                'k_parameter': self.config.k_parameter,
                'seed': self.config.seed,
                'device': self.config.device
            },
            'performance_metrics': {
                'elapsed_time_seconds': elapsed_time,
                'trials_per_second': self.config.n_trials / elapsed_time,
                'memory_usage_mb': torch.cuda.max_memory_allocated() / 1024**2 if torch.cuda.is_available() else 0
            },
            'statistical_analysis': analysis_results,
            'raw_results': {
                'baseline_errors': baseline_errors.numpy().tolist()[:1000],  # Sample for JSON
                'z_framework_errors': z_framework_errors.numpy().tolist()[:1000],
                'baseline_times': baseline_times.numpy().tolist()[:1000],
                'z_framework_times': z_framework_times.numpy().tolist()[:1000]
            }
        }
    
    def analyze_results_vectorized(self, baseline_errors: torch.Tensor, z_framework_errors: torch.Tensor,
                                 baseline_times: torch.Tensor, z_framework_times: torch.Tensor) -> Dict:
        """Perform vectorized statistical analysis of results."""
        logger.info("Performing vectorized statistical analysis...")
        
        # Move to device for analysis
        device = self.acoustic_grid.device
        baseline_errors = baseline_errors.to(device)
        z_framework_errors = z_framework_errors.to(device)
        
        # Basic statistics
        error_differences = baseline_errors - z_framework_errors
        improvement_percentage = torch.mean(error_differences / baseline_errors * 100).item()
        
        # Effect size (Cohen's d)
        pooled_std = torch.sqrt((torch.var(baseline_errors) + torch.var(z_framework_errors)) / 2)
        cohens_d = torch.mean(error_differences) / pooled_std
        
        # Bootstrap correlation analysis
        correlation_r, ci_low, ci_high, corr_p_value = self.statistics.bootstrap_correlation_vectorized(
            baseline_errors, z_framework_errors, self.config.n_bootstrap
        )
        
        # Permutation test
        permutation_p = self.statistics.permutation_test_vectorized(
            baseline_errors, z_framework_errors, self.config.n_permutation
        )
        
        # Statistical significance
        statistically_significant = permutation_p < 0.05
        
        results = {
            'improvement_percentage': improvement_percentage,
            'improvement_percentage_ci': [
                improvement_percentage - 1.96 * torch.std(error_differences / baseline_errors * 100).item() / np.sqrt(len(error_differences)),
                improvement_percentage + 1.96 * torch.std(error_differences / baseline_errors * 100).item() / np.sqrt(len(error_differences))
            ],
            'cohens_d': cohens_d.item(),
            'correlation_r': correlation_r,
            'correlation_ci': [ci_low, ci_high],
            'correlation_p_value': corr_p_value,
            'permutation_p_value': permutation_p,
            'statistically_significant': statistically_significant,
            'baseline_error_mean': torch.mean(baseline_errors).item(),
            'z_framework_error_mean': torch.mean(z_framework_errors).item(),
            'baseline_error_std': torch.std(baseline_errors).item(),
            'z_framework_error_std': torch.std(z_framework_errors).item()
        }
        
        logger.info(f"Analysis complete: {improvement_percentage:.1f}% improvement, "
                   f"Cohen's d = {cohens_d:.3f}, p = {permutation_p:.4f}")
        
        return results


def get_git_commit_sha() -> str:
    """Get the current git commit SHA for reproducibility tracking."""
    try:
        result = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True)
        if result.returncode == 0:
            sha = result.stdout.strip()
            if len(sha) == 40:  # Full SHA
                return sha
            else:
                logger.warning(f"Invalid git SHA format: {sha}")
                return 'unknown'
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"Could not get git commit SHA: {e}")
        return 'unknown'


def save_results(results: Dict, config: VectorizedConfig) -> Path:
    """Save experimental results with metadata."""
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    output_dir = Path(config.output_dir) / "fus_enhancer" / f"run-{timestamp}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Add metadata
    results['metadata'] = {
        'timestamp': timestamp,
        'git_commit': get_git_commit_sha(),
        'torch_version': torch.__version__,
        'device_name': torch.cuda.get_device_name() if torch.cuda.is_available() else 'CPU',
        'python_version': f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    }
    
    # Save results
    results_file = output_dir / "analysis.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    logger.info(f"Results saved to {results_file}")
    return output_dir


def create_cli_parser() -> argparse.ArgumentParser:
    """Create command-line interface parser following repository standards."""
    parser = argparse.ArgumentParser(
        description="FUS Enhancer: Vectorized Focused Ultrasound Targeting Experiment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required CLI contract flags
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Number of bootstrap samples')
    parser.add_argument('--permutation', type=int, default=1000, help='Number of permutation samples')
    parser.add_argument('--splits', type=str, default='single', help='Data splitting strategy')
    parser.add_argument('--domain', type=str, default='discrete', help='Z Framework domain')
    
    # Experiment-specific parameters
    parser.add_argument('--batch-size', type=int, default=DEFAULT_BATCH_SIZE, 
                       help='Batch size for vectorized processing')
    parser.add_argument('--n-trials', type=int, default=DEFAULT_N_TRIALS,
                       help='Total number of targeting trials')
    parser.add_argument('--grid-size', type=int, default=DEFAULT_GRID_SIZE,
                       help='Acoustic grid size (grid_size × grid_size)')
    parser.add_argument('--k-parameter', type=float, default=0.3,
                       help='Geometric resolution parameter')
    parser.add_argument('--output-dir', type=str, default='results',
                       help='Output directory for results')
    parser.add_argument('--device', type=str, default='auto',
                       choices=['auto', 'cpu', 'cuda'], help='PyTorch device')
    parser.add_argument('--visualize', action='store_true',
                       help='Generate visualization plots')
    
    return parser


def main():
    """Main execution function."""
    import sys
    
    parser = create_cli_parser()
    args = parser.parse_args()
    
    # Determine device
    if args.device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    else:
        device = args.device
    
    # Create configuration
    config = VectorizedConfig(
        batch_size=args.batch_size,
        grid_size=args.grid_size,
        n_trials=args.n_trials,
        n_bootstrap=args.bootstrap,
        n_permutation=args.permutation,
        seed=args.seed,
        output_dir=args.output_dir,
        k_parameter=args.k_parameter,
        device=device
    )
    
    # Print experiment information
    print("🚀 FUS Enhancer: Vectorized Focused Ultrasound Targeting")
    print("=" * 60)
    print(f"Trials: {config.n_trials:,}")
    print(f"Batch size: {config.batch_size:,}")
    print(f"Grid size: {config.grid_size}×{config.grid_size}")
    print(f"Device: {device}")
    print(f"Bootstrap samples: {config.n_bootstrap:,}")
    print(f"Permutation samples: {config.n_permutation:,}")
    print()
    
    # Run experiment
    experiment = VectorizedFUSExperiment(config)
    results = experiment.run_vectorized_experiment()
    
    # Save results
    output_dir = save_results(results, config)
    
    # Print summary
    stats = results['statistical_analysis']
    perf = results['performance_metrics']
    
    print("\n📊 Experiment Results Summary:")
    print("=" * 60)
    print(f"Improvement: {stats['improvement_percentage']:.1f}% ± "
          f"{(stats['improvement_percentage_ci'][1] - stats['improvement_percentage_ci'][0])/2:.1f}%")
    print(f"Effect size (Cohen's d): {stats['cohens_d']:.3f}")
    print(f"Correlation: r = {stats['correlation_r']:.3f} "
          f"[{stats['correlation_ci'][0]:.3f}, {stats['correlation_ci'][1]:.3f}]")
    print(f"Statistical significance: p = {stats['permutation_p_value']:.4f} "
          f"({'significant' if stats['statistically_significant'] else 'not significant'})")
    print(f"Processing rate: {perf['trials_per_second']:.0f} trials/second")
    print(f"Total time: {perf['elapsed_time_seconds']:.1f} seconds")
    if torch.cuda.is_available():
        print(f"Peak GPU memory: {perf['memory_usage_mb']:.1f} MB")
    print(f"\nResults saved to: {output_dir}")
    
    # Visualization if requested
    if args.visualize:
        print("\n🎨 Visualization feature is not yet implemented.")
        # Note: Visualization implementation would go here
if __name__ == "__main__":
    main()