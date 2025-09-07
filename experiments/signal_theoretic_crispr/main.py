#!/usr/bin/env python3
"""
Main Signal-Theoretic CRISPR Experiment Runner

This is the primary entry point for the two-arm benchmark experiment comparing:
- Arm A: Batch-scan baseline (Spring Batch style)
- Arm B: WAVE-CRISPR spectral-geodesic features

Enforces all scientific gates and generates comprehensive results.
"""

import os
import sys
import argparse
import logging
import json
import yaml
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd

# Import experiment modules
from .validation import HumanFASTAValidator, DeterminismValidator, validate_experiment_setup
from .baseline import BaselinePipeline, BaselineModels
from .spectral import SpectralFeatureExtractor, SpectralModels
from .statistics import ExperimentalStatistics


class SignalTheoreticExperiment:
    """Main experiment controller."""
    
    def __init__(self, config_path: str, output_dir: str, seed: int = 42):
        """
        Initialize experiment.
        
        Args:
            config_path: Path to experiment manifest
            output_dir: Output directory for results
            seed: Random seed for reproducibility
        """
        self.config_path = config_path
        self.output_dir = Path(output_dir)
        self.seed = seed
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        self.logger = logging.getLogger(__name__)
        
        # Load configuration
        self.config = self.load_config()
        
        # Initialize components
        self.validator = HumanFASTAValidator()
        self.det_validator = DeterminismValidator(seed)
        self.det_validator.set_seeds()
        
        # Initialize pipelines
        self.baseline_pipeline = BaselinePipeline(
            pam_pattern=self.config['parameters']['pam_sequence'],
            seed=seed
        )
        self.spectral_extractor = SpectralFeatureExtractor(
            k=self.config['parameters']['geodesic_k'],
            seed=seed
        )
        
        # Initialize models
        self.baseline_models = BaselineModels(seed)
        self.spectral_models = SpectralModels(seed)
        
        # Initialize statistics
        self.statistics = ExperimentalStatistics(
            bootstrap_iterations=self.config['parameters']['bootstrap_iterations'],
            permutation_iterations=self.config['parameters']['permutation_iterations'],
            seed=seed
        )
        
        # Results storage
        self.results = {
            'experiment_info': {
                'name': self.config['name'],
                'version': self.config['version'],
                'timestamp': datetime.now().isoformat(),
                'config': self.config,
                'reproducibility': self.det_validator.get_reproducibility_info()
            },
            'validation_report': {},
            'baseline_results': {},
            'spectral_results': {},
            'comparison_results': {},
            'control_results': {}
        }
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(self.output_dir / 'experiment.log')
            ]
        )
        
    def load_config(self) -> Dict[str, Any]:
        """Load experiment configuration."""
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def validate_setup(self) -> Dict[str, Any]:
        """
        Validate complete experiment setup according to scientific gates.
        
        Returns:
            Validation report
        """
        self.logger.info("Validating experiment setup...")
        
        try:
            validation_report = validate_experiment_setup(self.config_path)
            self.results['validation_report'] = validation_report
            
            # Check if any datasets failed validation
            failed_datasets = [
                name for name, status in validation_report['gates_status'].items()
                if status.startswith('FAILED')
            ]
            
            if failed_datasets:
                raise ValueError(f"Dataset validation failed for: {failed_datasets}")
            
            self.logger.info("✓ All validation gates passed")
            return validation_report
            
        except Exception as e:
            self.logger.error(f"Setup validation failed: {e}")
            raise
    
    def load_and_validate_data(self) -> Dict[str, pd.DataFrame]:
        """
        Load and validate all datasets.
        
        Returns:
            Dictionary of validated datasets
        """
        self.logger.info("Loading and validating datasets...")
        
        datasets = {}
        
        for dataset_config in self.config['datasets']['primary']:
            dataset_name = dataset_config['name']
            dataset_path = dataset_config['path']
            
            self.logger.info(f"Loading dataset: {dataset_name}")
            
            try:
                # Validate dataset
                validation_report = self.validator.validate_dataset(
                    dataset_path,
                    {'organism': dataset_config['organism']}
                )
                
                # Load data
                if dataset_path.endswith('.csv'):
                    df = pd.read_csv(dataset_path)
                else:
                    raise ValueError(f"Unsupported format: {dataset_path}")
                
                # Additional validation for required columns
                if dataset_config['type'] == 'on_target_efficiency':
                    required_cols = ['sequence', 'efficiency']
                elif dataset_config['type'] == 'off_target_classification':
                    required_cols = ['sequence', 'off_target']
                else:
                    required_cols = ['sequence']
                
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    raise ValueError(f"Missing columns in {dataset_name}: {missing_cols}")
                
                # Validate sequences
                for idx, row in df.iterrows():
                    self.validator.validate_guide_sequence(row['sequence'])
                
                datasets[dataset_name] = df
                self.logger.info(f"✓ Dataset {dataset_name} loaded: {len(df)} sequences")
                
            except Exception as e:
                self.logger.error(f"Failed to load dataset {dataset_name}: {e}")
                raise
        
        return datasets
    
    def run_baseline_arm(self, datasets: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """
        Run Arm A: Batch-scan baseline processing.
        
        Args:
            datasets: Validated datasets
            
        Returns:
            Baseline results
        """
        self.logger.info("Running Arm A: Batch-scan baseline...")
        
        baseline_results = {
            'features_extracted': {},
            'model_performance': {},
            'predictions': {}
        }
        
        for dataset_name, df in datasets.items():
            self.logger.info(f"Processing {dataset_name} with baseline pipeline...")
            
            try:
                # Extract baseline features
                processed_df = self.baseline_pipeline.process_data(df)
                baseline_results['features_extracted'][dataset_name] = {
                    'n_sequences': len(processed_df),
                    'n_features': len([col for col in processed_df.columns if pd.api.types.is_numeric_dtype(processed_df[col])]),
                    'feature_names': [col for col in processed_df.columns if pd.api.types.is_numeric_dtype(processed_df[col])]
                }
                
                # Train models if labels are available
                if 'efficiency' in processed_df.columns:
                    # On-target regression
                    train_metrics = self.baseline_models.train_on_target(processed_df)
                    baseline_results['model_performance'][f'{dataset_name}_on_target'] = train_metrics
                    
                    # Make predictions
                    predictions = self.baseline_models.predict_on_target(processed_df)
                    baseline_results['predictions'][f'{dataset_name}_on_target'] = predictions.tolist()
                
                if 'off_target' in processed_df.columns:
                    # Off-target classification
                    train_metrics = self.baseline_models.train_off_target(processed_df)
                    baseline_results['model_performance'][f'{dataset_name}_off_target'] = train_metrics
                    
                    # Make predictions
                    pred_labels, pred_probs = self.baseline_models.predict_off_target(processed_df)
                    baseline_results['predictions'][f'{dataset_name}_off_target'] = {
                        'labels': pred_labels.tolist(),
                        'probabilities': pred_probs.tolist()
                    }
                
                # Store processed data for comparison
                baseline_results[f'{dataset_name}_processed'] = processed_df
                
            except Exception as e:
                self.logger.error(f"Baseline processing failed for {dataset_name}: {e}")
                raise
        
        self.logger.info("✓ Baseline arm completed")
        return baseline_results
    
    def run_spectral_arm(self, datasets: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """
        Run Arm B: WAVE-CRISPR spectral-geodesic processing.
        
        Args:
            datasets: Validated datasets
            
        Returns:
            Spectral results
        """
        self.logger.info("Running Arm B: WAVE-CRISPR spectral-geodesic...")
        
        spectral_results = {
            'features_extracted': {},
            'model_performance': {},
            'predictions': {},
            'encoding_info': {
                'bio_encoder': self.spectral_extractor.bio_encoder.get_encoding_info(),
                'arb_encoder': self.spectral_extractor.arb_encoder.get_encoding_info(),
                'geodesic_k': self.spectral_extractor.k
            }
        }
        
        for dataset_name, df in datasets.items():
            self.logger.info(f"Processing {dataset_name} with spectral pipeline...")
            
            try:
                # Extract spectral features
                spectral_features_list = []
                
                for idx, row in df.iterrows():
                    guide_data = {
                        'guide_sequence': row['sequence'],
                        'original_index': idx
                    }
                    
                    # Add context window if available (in practice, would extract from hg38)
                    # For now, use the sequence itself
                    guide_data['context_sequence'] = row['sequence']
                    
                    # Extract spectral features
                    features = self.spectral_extractor.extract_features(guide_data)
                    features.update(guide_data)
                    
                    # Add labels if available
                    if 'efficiency' in row:
                        features['efficiency'] = row['efficiency']
                    if 'off_target' in row:
                        features['off_target'] = row['off_target']
                    
                    spectral_features_list.append(features)
                
                processed_df = pd.DataFrame(spectral_features_list)
                
                spectral_results['features_extracted'][dataset_name] = {
                    'n_sequences': len(processed_df),
                    'n_features': len([col for col in processed_df.columns if pd.api.types.is_numeric_dtype(processed_df[col])]),
                    'bio_features': len([col for col in processed_df.columns if col.startswith('bio_')]),
                    'arb_features': len([col for col in processed_df.columns if col.startswith('arb_')]),
                    'z_features': len([col for col in processed_df.columns if col.startswith('z_')])
                }
                
                # Train models if labels are available
                if 'efficiency' in processed_df.columns:
                    # On-target regression
                    train_metrics = self.spectral_models.train_on_target(processed_df)
                    spectral_results['model_performance'][f'{dataset_name}_on_target'] = train_metrics
                    
                    # Make predictions
                    predictions = self.spectral_models.predict_on_target(processed_df)
                    spectral_results['predictions'][f'{dataset_name}_on_target'] = predictions.tolist()
                
                # Store processed data for comparison
                spectral_results[f'{dataset_name}_processed'] = processed_df
                
            except Exception as e:
                self.logger.error(f"Spectral processing failed for {dataset_name}: {e}")
                raise
        
        self.logger.info("✓ Spectral arm completed")
        return spectral_results
    
    def compare_arms(self, baseline_results: Dict[str, Any], 
                    spectral_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compare performance between baseline and spectral arms.
        
        Args:
            baseline_results: Results from baseline arm
            spectral_results: Results from spectral arm
            
        Returns:
            Comparison results with statistical analysis
        """
        self.logger.info("Comparing baseline vs spectral performance...")
        
        comparison_results = {
            'statistical_tests': {},
            'improvement_summary': {},
            'hypothesis_testing': {}
        }
        
        # Find common datasets with predictions
        common_predictions = set(baseline_results['predictions'].keys()) & set(spectral_results['predictions'].keys())
        
        for pred_key in common_predictions:
            self.logger.info(f"Comparing predictions for: {pred_key}")
            
            try:
                # Get true labels
                # Parse prediction key: 'dataset_name_task_type'
                parts = pred_key.split('_')
                if len(parts) >= 3:
                    dataset_name = '_'.join(parts[:-2])  # Everything except last two parts
                    task_type = '_'.join(parts[-2:])     # Last two parts: e.g., 'on_target'
                else:
                    # Fallback for simpler naming
                    dataset_name = parts[0]
                    task_type = '_'.join(parts[1:]) if len(parts) > 1 else ''
                
                self.logger.info(f"Processing dataset: {dataset_name}, task: {task_type}")
                
                if task_type == 'on_target':
                    # On-target regression comparison
                    # Get processed data from either baseline or spectral results
                    processed_data_key = f'{dataset_name}_processed'
                    if processed_data_key in baseline_results:
                        processed_df = baseline_results[processed_data_key]
                    elif processed_data_key in spectral_results:
                        processed_df = spectral_results[processed_data_key]
                    else:
                        raise ValueError(f"No processed data found for {dataset_name}")
                    
                    y_true = processed_df['efficiency'].values
                    y_pred_baseline = np.array(baseline_results['predictions'][pred_key])
                    y_pred_spectral = np.array(spectral_results['predictions'][pred_key])
                    
                    # Statistical evaluation
                    regression_stats = self.statistics.evaluate_regression_performance(
                        y_true, y_pred_baseline, y_pred_spectral
                    )
                    comparison_results['statistical_tests'][pred_key] = regression_stats
                    
                    # Check hypothesis: improvement ≥ 0.05
                    r2_improvement = regression_stats['comparison']['r2_improvement']
                    meets_threshold = r2_improvement['mean_difference'] >= 0.05
                    significant = r2_improvement['significant_improvement']
                    
                    comparison_results['hypothesis_testing'][pred_key] = {
                        'improvement_threshold_met': meets_threshold,
                        'statistically_significant': significant,
                        'hypothesis_h1_supported': meets_threshold and significant
                    }
                
            except Exception as e:
                self.logger.error(f"Comparison failed for {pred_key}: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
                continue
        
        # Generate improvement summary
        improvements = []
        for test_key, test_results in comparison_results['statistical_tests'].items():
            if 'r2_improvement' in test_results['comparison']:
                improvement = test_results['comparison']['r2_improvement']['mean_difference']
                improvements.append(improvement)
        
        if improvements:
            comparison_results['improvement_summary'] = {
                'mean_r2_improvement': np.mean(improvements),
                'median_r2_improvement': np.median(improvements),
                'min_r2_improvement': np.min(improvements),
                'max_r2_improvement': np.max(improvements),
                'n_comparisons': len(improvements)
            }
        
        self.logger.info("✓ Arm comparison completed")
        return comparison_results
    
    def run_control_experiments(self, datasets: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        """
        Run control experiments for validation.
        
        Args:
            datasets: Original datasets
            
        Returns:
            Control experiment results
        """
        self.logger.info("Running control experiments...")
        
        control_results = {}
        
        for dataset_name, df in datasets.items():
            if 'efficiency' in df.columns:  # Only run controls on labeled data
                try:
                    dataset_controls = self.statistics.run_control_validation(
                        df, 
                        {
                            'baseline': self.baseline_pipeline,
                            'spectral': self.spectral_extractor
                        },
                        {
                            'baseline': self.baseline_models,
                            'spectral': self.spectral_models
                        }
                    )
                    control_results[dataset_name] = dataset_controls
                    
                except Exception as e:
                    self.logger.error(f"Control experiments failed for {dataset_name}: {e}")
                    control_results[dataset_name] = {'error': str(e)}
        
        self.logger.info("✓ Control experiments completed")
        return control_results
    
    def save_results(self):
        """Save all results to output directory."""
        self.logger.info("Saving results...")
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save main results as JSON
        results_file = self.output_dir / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        
        # Save detailed CSV results if available
        for arm in ['baseline', 'spectral']:
            arm_results = self.results[f'{arm}_results']
            for dataset_name in arm_results:
                if dataset_name.endswith('_processed'):
                    df = arm_results[dataset_name]
                    csv_file = self.output_dir / f'{arm}_{dataset_name}.csv'
                    df.to_csv(csv_file, index=False)
        
        # Save environment information
        env_file = self.output_dir / 'env.txt'
        with open(env_file, 'w') as f:
            f.write("# Experiment Environment\n")
            f.write(f"Timestamp: {datetime.now().isoformat()}\n")
            f.write(f"Git commit: {self.det_validator.git_commit}\n")
            f.write(f"Random seed: {self.seed}\n")
            f.write("\n# Python Environment\n")
            for key, value in self.det_validator.environment_info.items():
                f.write(f"{key}: {value}\n")
        
        # Save execution log
        log_file = self.output_dir / 'log.txt'
        # Log file is already being written by logging handler
        
        self.logger.info(f"✓ Results saved to {self.output_dir}")
    
    def run_complete_experiment(self) -> Dict[str, Any]:
        """
        Run the complete two-arm benchmark experiment.
        
        Returns:
            Complete experiment results
        """
        self.logger.info("Starting signal-theoretic CRISPR experiment...")
        
        try:
            # 1. Validate setup
            self.validate_setup()
            
            # 2. Load and validate data
            datasets = self.load_and_validate_data()
            
            # 3. Run baseline arm
            baseline_results = self.run_baseline_arm(datasets)
            self.results['baseline_results'] = baseline_results
            
            # 4. Run spectral arm
            spectral_results = self.run_spectral_arm(datasets)
            self.results['spectral_results'] = spectral_results
            
            # 5. Compare arms
            comparison_results = self.compare_arms(baseline_results, spectral_results)
            self.results['comparison_results'] = comparison_results
            
            # 6. Run control experiments
            control_results = self.run_control_experiments(datasets)
            self.results['control_results'] = control_results
            
            # 7. Save results
            self.save_results()
            
            self.logger.info("✓ Experiment completed successfully")
            return self.results
            
        except Exception as e:
            self.logger.error(f"Experiment failed: {e}")
            raise


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Signal-Theoretic CRISPR Experiment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--config', 
        default='experiments/signal_theoretic_crispr/manifest.yml',
        help='Path to experiment configuration file'
    )
    parser.add_argument(
        '--output',
        default=f'results/signal_theoretic_crispr/run-{datetime.now().strftime("%Y%m%d-%H%M%S")}',
        help='Output directory for results'
    )
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Bootstrap iterations')
    parser.add_argument('--permutation', type=int, default=1000, help='Permutation iterations')
    parser.add_argument('--splits', default='split-by-gene', help='Data splitting strategy')
    parser.add_argument('--domain', default='discrete', help='Domain type')
    parser.add_argument('--pam', default='NGG', help='PAM sequence pattern')
    
    args = parser.parse_args()
    
    # Run experiment
    experiment = SignalTheoreticExperiment(
        config_path=args.config,
        output_dir=args.output,
        seed=args.seed
    )
    
    results = experiment.run_complete_experiment()
    
    # Print summary
    print("\n" + "="*60)
    print("EXPERIMENT SUMMARY")
    print("="*60)
    
    if 'improvement_summary' in results.get('comparison_results', {}):
        summary = results['comparison_results']['improvement_summary']
        print(f"Mean R² improvement: {summary['mean_r2_improvement']:.4f}")
        print(f"Number of comparisons: {summary['n_comparisons']}")
    
    # Check hypothesis
    hypothesis_results = results.get('comparison_results', {}).get('hypothesis_testing', {})
    h1_supported = any(
        test.get('hypothesis_h1_supported', False) 
        for test in hypothesis_results.values()
    )
    
    print(f"\nHypothesis H₁ (≥0.05 improvement) supported: {h1_supported}")
    print(f"Results saved to: {args.output}")
    print("="*60)


if __name__ == "__main__":
    main()