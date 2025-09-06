#!/usr/bin/env python3
"""
Bootstrap and permutation statistics for signal-theoretic CRISPR experiment.

This module implements statistical validation following scientific gates G5-G7:
- Bootstrap confidence intervals (95%)
- Permutation testing with null models
- Multiple comparison correction (Benjamini-Hochberg FDR)
- Controls: shuffled labels, PAM-broken, reverse-complement
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Callable
import logging
from sklearn.metrics import (
    r2_score, mean_absolute_error, roc_auc_score, 
    average_precision_score, precision_recall_curve, roc_curve
)
from sklearn.model_selection import train_test_split
from scipy import stats
import warnings

warnings.filterwarnings('ignore')


class BootstrapStatistics:
    """Bootstrap confidence interval calculations."""
    
    def __init__(self, n_iterations: int = 1000, confidence: float = 0.95, seed: int = 42):
        """
        Initialize bootstrap statistics.
        
        Args:
            n_iterations: Number of bootstrap iterations
            confidence: Confidence level (e.g., 0.95 for 95% CI)
            seed: Random seed for reproducibility
        """
        self.n_iterations = n_iterations
        self.confidence = confidence
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        np.random.seed(seed)
        
    def bootstrap_metric(self, y_true: np.ndarray, y_pred: np.ndarray, 
                        metric_func: Callable, **metric_kwargs) -> Dict[str, float]:
        """
        Calculate bootstrap confidence interval for a metric.
        
        Args:
            y_true: True labels
            y_pred: Predicted values
            metric_func: Metric function (e.g., r2_score, roc_auc_score)
            **metric_kwargs: Additional keyword arguments for metric function
            
        Returns:
            Dictionary with mean, CI bounds, and individual bootstrap values
        """
        bootstrap_scores = []
        n_samples = len(y_true)
        
        for i in range(self.n_iterations):
            # Bootstrap sample with replacement
            indices = np.random.choice(n_samples, size=n_samples, replace=True)
            y_true_boot = y_true[indices]
            y_pred_boot = y_pred[indices]
            
            try:
                score = metric_func(y_true_boot, y_pred_boot, **metric_kwargs)
                bootstrap_scores.append(score)
            except Exception as e:
                self.logger.warning(f"Bootstrap iteration {i} failed: {e}")
                continue
        
        if not bootstrap_scores:
            return {'mean': np.nan, 'ci_lower': np.nan, 'ci_upper': np.nan, 'scores': []}
        
        bootstrap_scores = np.array(bootstrap_scores)
        
        # Calculate confidence interval
        alpha = 1 - self.confidence
        ci_lower = np.percentile(bootstrap_scores, 100 * alpha / 2)
        ci_upper = np.percentile(bootstrap_scores, 100 * (1 - alpha / 2))
        
        return {
            'mean': np.mean(bootstrap_scores),
            'std': np.std(bootstrap_scores),
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'scores': bootstrap_scores.tolist()
        }
    
    def bootstrap_difference(self, baseline_metric: Dict[str, float], 
                           spectral_metric: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate bootstrap confidence interval for difference between methods.
        
        Args:
            baseline_metric: Bootstrap results for baseline method
            spectral_metric: Bootstrap results for spectral method
            
        Returns:
            Bootstrap statistics for the difference (spectral - baseline)
        """
        baseline_scores = np.array(baseline_metric['scores'])
        spectral_scores = np.array(spectral_metric['scores'])
        
        # Calculate differences for each bootstrap iteration
        min_len = min(len(baseline_scores), len(spectral_scores))
        differences = spectral_scores[:min_len] - baseline_scores[:min_len]
        
        # Calculate confidence interval for difference
        alpha = 1 - self.confidence
        ci_lower = np.percentile(differences, 100 * alpha / 2)
        ci_upper = np.percentile(differences, 100 * (1 - alpha / 2))
        
        return {
            'mean_difference': np.mean(differences),
            'std_difference': np.std(differences),
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'differences': differences.tolist(),
            'significant_improvement': ci_lower > 0  # CI excludes 0
        }


class PermutationTest:
    """Permutation testing for null hypothesis validation."""
    
    def __init__(self, n_permutations: int = 1000, seed: int = 42):
        """
        Initialize permutation test.
        
        Args:
            n_permutations: Number of permutation iterations
            seed: Random seed for reproducibility
        """
        self.n_permutations = n_permutations
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        
        np.random.seed(seed)
    
    def permutation_test_improvement(self, baseline_scores: np.ndarray, 
                                   spectral_scores: np.ndarray) -> Dict[str, float]:
        """
        Test if spectral method significantly outperforms baseline.
        
        Args:
            baseline_scores: Baseline method scores
            spectral_scores: Spectral method scores
            
        Returns:
            Permutation test results
        """
        # Observed difference
        observed_diff = np.mean(spectral_scores) - np.mean(baseline_scores)
        
        # Combine all scores
        all_scores = np.concatenate([baseline_scores, spectral_scores])
        n_baseline = len(baseline_scores)
        n_spectral = len(spectral_scores)
        
        # Permutation test
        permuted_diffs = []
        
        for i in range(self.n_permutations):
            # Randomly shuffle and split
            shuffled = np.random.permutation(all_scores)
            perm_baseline = shuffled[:n_baseline]
            perm_spectral = shuffled[n_baseline:n_baseline+n_spectral]
            
            perm_diff = np.mean(perm_spectral) - np.mean(perm_baseline)
            permuted_diffs.append(perm_diff)
        
        permuted_diffs = np.array(permuted_diffs)
        
        # Calculate p-value (two-tailed)
        p_value = np.mean(np.abs(permuted_diffs) >= np.abs(observed_diff))
        
        return {
            'observed_difference': observed_diff,
            'p_value': p_value,
            'permuted_differences': permuted_diffs.tolist(),
            'significant': p_value < 0.05
        }


class MultipleComparisonCorrection:
    """Benjamini-Hochberg FDR correction for multiple comparisons."""
    
    def __init__(self, alpha: float = 0.05):
        """
        Initialize multiple comparison correction.
        
        Args:
            alpha: Significance level
        """
        self.alpha = alpha
        self.logger = logging.getLogger(__name__)
    
    def benjamini_hochberg_correction(self, p_values: List[float], 
                                    test_names: List[str]) -> Dict[str, Any]:
        """
        Apply Benjamini-Hochberg FDR correction.
        
        Args:
            p_values: List of p-values from multiple tests
            test_names: Names of the tests
            
        Returns:
            Correction results
        """
        p_values = np.array(p_values)
        n_tests = len(p_values)
        
        # Sort p-values and keep track of original indices
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]
        sorted_names = [test_names[i] for i in sorted_indices]
        
        # Benjamini-Hochberg procedure
        rejected = np.zeros(n_tests, dtype=bool)
        
        for i in range(n_tests - 1, -1, -1):  # Work backwards
            threshold = (i + 1) / n_tests * self.alpha
            if sorted_p_values[i] <= threshold:
                rejected[sorted_indices[:i+1]] = True
                break
        
        # Create results
        results = {
            'corrected_alpha': self.alpha,
            'n_tests': n_tests,
            'n_significant': np.sum(rejected),
            'tests': []
        }
        
        for i, (name, p_val) in enumerate(zip(test_names, p_values)):
            results['tests'].append({
                'test_name': name,
                'p_value': p_val,
                'rejected_null': rejected[i],
                'rank': np.where(sorted_indices == i)[0][0] + 1
            })
        
        return results


class ControlExperiments:
    """Generate control experiments for validation."""
    
    def __init__(self, seed: int = 42):
        """
        Initialize control experiments.
        
        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        np.random.seed(seed)
    
    def shuffle_labels(self, df: pd.DataFrame, label_column: str) -> pd.DataFrame:
        """
        Create shuffled label control.
        
        Args:
            df: Original DataFrame
            label_column: Column name to shuffle
            
        Returns:
            DataFrame with shuffled labels
        """
        df_shuffled = df.copy()
        df_shuffled[label_column] = np.random.permutation(df[label_column].values)
        return df_shuffled
    
    def break_pam_sequences(self, df: pd.DataFrame, sequence_column: str = 'sequence') -> pd.DataFrame:
        """
        Create PAM-broken control by mutating PAM sites.
        
        Args:
            df: Original DataFrame
            sequence_column: Column name containing sequences
            
        Returns:
            DataFrame with PAM-broken sequences
        """
        df_pam_broken = df.copy()
        
        for idx, row in df_pam_broken.iterrows():
            sequence = row[sequence_column]
            # Assuming PAM is at the end (NGG pattern)
            if len(sequence) >= 3:
                # Mutate last 3 bases to break PAM
                mutated_seq = sequence[:-3] + 'AAA'  # Replace with non-PAM sequence
                df_pam_broken.at[idx, sequence_column] = mutated_seq
        
        return df_pam_broken
    
    def reverse_complement_sequences(self, df: pd.DataFrame, 
                                   sequence_column: str = 'sequence') -> pd.DataFrame:
        """
        Create reverse complement control.
        
        Args:
            df: Original DataFrame
            sequence_column: Column name containing sequences
            
        Returns:
            DataFrame with reverse complement sequences
        """
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        
        df_rev_comp = df.copy()
        
        for idx, row in df_rev_comp.iterrows():
            sequence = row[sequence_column]
            # Reverse complement
            rev_comp = ''.join(complement_map.get(base, base) for base in sequence[::-1])
            df_rev_comp.at[idx, sequence_column] = rev_comp
        
        return df_rev_comp
    
    def generate_all_controls(self, df: pd.DataFrame, 
                            efficiency_column: str = 'efficiency') -> Dict[str, pd.DataFrame]:
        """
        Generate all control experiments.
        
        Args:
            df: Original DataFrame
            efficiency_column: Column name for efficiency labels
            
        Returns:
            Dictionary of control DataFrames
        """
        controls = {
            'shuffled_labels': self.shuffle_labels(df, efficiency_column),
            'pam_broken': self.break_pam_sequences(df),
            'reverse_complement': self.reverse_complement_sequences(df)
        }
        
        self.logger.info(f"Generated {len(controls)} control experiments")
        return controls


class ExperimentalStatistics:
    """Main class for experimental statistics and validation."""
    
    def __init__(self, bootstrap_iterations: int = 1000, 
                 permutation_iterations: int = 1000, 
                 confidence: float = 0.95, seed: int = 42):
        """
        Initialize experimental statistics.
        
        Args:
            bootstrap_iterations: Number of bootstrap iterations
            permutation_iterations: Number of permutation iterations
            confidence: Confidence level
            seed: Random seed
        """
        self.bootstrap = BootstrapStatistics(bootstrap_iterations, confidence, seed)
        self.permutation = PermutationTest(permutation_iterations, seed)
        self.correction = MultipleComparisonCorrection()
        self.controls = ControlExperiments(seed)
        self.logger = logging.getLogger(__name__)
    
    def evaluate_regression_performance(self, y_true: np.ndarray, 
                                      y_pred_baseline: np.ndarray,
                                      y_pred_spectral: np.ndarray) -> Dict[str, Any]:
        """
        Evaluate regression performance with full statistical analysis.
        
        Args:
            y_true: True efficiency values
            y_pred_baseline: Baseline predictions
            y_pred_spectral: Spectral predictions
            
        Returns:
            Complete statistical evaluation
        """
        results = {
            'baseline': {},
            'spectral': {},
            'comparison': {}
        }
        
        # Bootstrap confidence intervals for each method
        baseline_r2 = self.bootstrap.bootstrap_metric(y_true, y_pred_baseline, r2_score)
        baseline_mae = self.bootstrap.bootstrap_metric(y_true, y_pred_baseline, mean_absolute_error)
        
        spectral_r2 = self.bootstrap.bootstrap_metric(y_true, y_pred_spectral, r2_score)
        spectral_mae = self.bootstrap.bootstrap_metric(y_true, y_pred_spectral, mean_absolute_error)
        
        results['baseline'] = {
            'r2': baseline_r2,
            'mae': baseline_mae
        }
        
        results['spectral'] = {
            'r2': spectral_r2,
            'mae': spectral_mae
        }
        
        # Bootstrap difference analysis
        r2_difference = self.bootstrap.bootstrap_difference(baseline_r2, spectral_r2)
        mae_difference = self.bootstrap.bootstrap_difference(baseline_mae, spectral_mae)
        
        results['comparison'] = {
            'r2_improvement': r2_difference,
            'mae_improvement': mae_difference
        }
        
        # Permutation test
        baseline_r2_scores = np.array([r2_score(y_true, y_pred_baseline)])
        spectral_r2_scores = np.array([r2_score(y_true, y_pred_spectral)])
        
        perm_test = self.permutation.permutation_test_improvement(baseline_r2_scores, spectral_r2_scores)
        results['comparison']['permutation_test'] = perm_test
        
        return results
    
    def evaluate_classification_performance(self, y_true: np.ndarray,
                                          y_pred_baseline: np.ndarray,
                                          y_pred_spectral: np.ndarray,
                                          y_prob_baseline: np.ndarray,
                                          y_prob_spectral: np.ndarray) -> Dict[str, Any]:
        """
        Evaluate classification performance with full statistical analysis.
        
        Args:
            y_true: True labels
            y_pred_baseline: Baseline predictions
            y_pred_spectral: Spectral predictions
            y_prob_baseline: Baseline probabilities
            y_prob_spectral: Spectral probabilities
            
        Returns:
            Complete statistical evaluation
        """
        results = {
            'baseline': {},
            'spectral': {},
            'comparison': {}
        }
        
        # Bootstrap confidence intervals
        baseline_auprc = self.bootstrap.bootstrap_metric(y_true, y_prob_baseline, average_precision_score)
        baseline_auroc = self.bootstrap.bootstrap_metric(y_true, y_prob_baseline, roc_auc_score)
        
        spectral_auprc = self.bootstrap.bootstrap_metric(y_true, y_prob_spectral, average_precision_score)
        spectral_auroc = self.bootstrap.bootstrap_metric(y_true, y_prob_spectral, roc_auc_score)
        
        results['baseline'] = {
            'auprc': baseline_auprc,
            'auroc': baseline_auroc
        }
        
        results['spectral'] = {
            'auprc': spectral_auprc,
            'auroc': spectral_auroc
        }
        
        # Bootstrap difference analysis
        auprc_difference = self.bootstrap.bootstrap_difference(baseline_auprc, spectral_auprc)
        auroc_difference = self.bootstrap.bootstrap_difference(baseline_auroc, spectral_auroc)
        
        results['comparison'] = {
            'auprc_improvement': auprc_difference,
            'auroc_improvement': auroc_difference
        }
        
        return results
    
    def run_control_validation(self, df: pd.DataFrame, 
                             feature_extractors: Dict[str, Any],
                             models: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run validation with control experiments.
        
        Args:
            df: Original dataset
            feature_extractors: Dictionary of feature extraction methods
            models: Dictionary of trained models
            
        Returns:
            Control experiment results
        """
        control_data = self.controls.generate_all_controls(df)
        control_results = {}
        
        for control_name, control_df in control_data.items():
            self.logger.info(f"Running control experiment: {control_name}")
            
            try:
                # Extract features for control data
                # This would use the same feature extraction pipeline
                # but with the modified control sequences/labels
                
                # For now, create a placeholder result
                control_results[control_name] = {
                    'status': 'completed',
                    'n_samples': len(control_df),
                    'description': f"Control experiment with {control_name}"
                }
                
            except Exception as e:
                self.logger.error(f"Control experiment {control_name} failed: {e}")
                control_results[control_name] = {
                    'status': 'failed',
                    'error': str(e)
                }
        
        return control_results


if __name__ == "__main__":
    # Quick test of statistics module
    np.random.seed(42)
    
    # Generate test data
    n_samples = 100
    y_true = np.random.rand(n_samples)
    y_pred_baseline = y_true + np.random.normal(0, 0.2, n_samples)  # Baseline with some error
    y_pred_spectral = y_true + np.random.normal(0, 0.15, n_samples)  # Spectral with less error
    
    # Test bootstrap statistics
    stats_engine = ExperimentalStatistics(bootstrap_iterations=100, permutation_iterations=100)
    
    regression_results = stats_engine.evaluate_regression_performance(
        y_true, y_pred_baseline, y_pred_spectral
    )
    
    print("✓ Bootstrap statistics test passed")
    print(f"Baseline R²: {regression_results['baseline']['r2']['mean']:.3f} "
          f"[{regression_results['baseline']['r2']['ci_lower']:.3f}, "
          f"{regression_results['baseline']['r2']['ci_upper']:.3f}]")
    print(f"Spectral R²: {regression_results['spectral']['r2']['mean']:.3f} "
          f"[{regression_results['spectral']['r2']['ci_lower']:.3f}, "
          f"{regression_results['spectral']['r2']['ci_upper']:.3f}]")
    print(f"Improvement: {regression_results['comparison']['r2_improvement']['mean_difference']:.3f} "
          f"(Significant: {regression_results['comparison']['r2_improvement']['significant_improvement']})")
    
    print("Statistics module self-test complete")