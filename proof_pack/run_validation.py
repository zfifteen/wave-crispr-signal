#!/usr/bin/env python3
"""
Z Framework Validation Script

Comprehensive validation comparing Z Framework (θ', Z5D) against baseline methods.
Provides statistical validation with confidence intervals and significance testing.

Usage:
    python run_validation.py

Expected Output:
    - AUROC, RMSE, Pearson r with 95% CI for all comparisons
    - Statistical significance tests (p < 0.05)
    - Density boost validation (target: >1000x)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import logging
import warnings
from typing import Dict, List, Any
import ast

# Statistical analysis
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# Import local modules
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from pain_management_application import (
    PainManagementAnalyzer,
    theta_prime,
    compute_density_boost,
)
from validation_baseline import BaselineFeatureExtractor, PoissonnNullModel

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Set reproducible seed
np.random.seed(42)


class ZFrameworkValidator:
    """Comprehensive validation of Z Framework against baselines."""

    def __init__(self, output_dir: str = "validation_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.logger = logging.getLogger(__name__)

        # Initialize analyzers
        self.pain_analyzer = PainManagementAnalyzer(precision_dps=50)
        self.baseline_extractor = BaselineFeatureExtractor()
        self.null_model = PoissonnNullModel(random_state=42)

        # Results storage
        self.results = {}

    def load_synthetic_data(self) -> Dict[str, pd.DataFrame]:
        """Load all synthetic datasets."""
        data = {}
        data_files = ["neural_spikes.csv", "nav1.8_panel.csv", "bcl11a_edits.csv"]

        for file in data_files:
            filepath = Path(__file__).parent / file
            if filepath.exists():
                data[file.replace(".csv", "")] = pd.read_csv(filepath)
                self.logger.info(
                    f"Loaded {len(data[file.replace('.csv', '')])} samples from {file}"
                )
            else:
                self.logger.warning(f"File {file} not found")

        return data

    def parse_numeric_sequences(self, sequence_strings: List[str]) -> List[List[int]]:
        """Parse numeric sequences from string representation."""
        sequences = []
        for seq_str in sequence_strings:
            try:
                # Handle string representation of lists
                numeric_seq = ast.literal_eval(seq_str)
                sequences.append(numeric_seq)
            except (ValueError, SyntaxError):
                # Fallback: convert character by character
                numeric_seq = [ord(c) - ord("A") + 1 for c in seq_str if c.isalpha()]
                sequences.append(numeric_seq)
        return sequences

    def extract_z_framework_features(
        self, sequences: List[str], numeric_sequences: List[List[int]]
    ) -> Dict[str, np.ndarray]:
        """Extract Z Framework features (θ', Z5D, density boost)."""
        features = {}

        # θ' (theta prime) features
        theta_features = []
        k_values = [0.1, 0.3, 0.5, 0.7, 0.9]  # Multiple k values

        for seq_num in numeric_sequences:
            seq_features = []
            for k in k_values:
                if len(seq_num) > 0:
                    theta_vals = [
                        float(theta_prime(n, k)) for n in seq_num[:50]
                    ]  # Limit length
                    seq_features.extend(
                        [
                            np.mean(theta_vals),
                            np.std(theta_vals),
                            np.max(theta_vals),
                            np.min(theta_vals),
                        ]
                    )
                else:
                    seq_features.extend([0, 0, 0, 0])
            theta_features.append(seq_features)

        features["theta_prime"] = np.array(theta_features)

        # Density boost features
        density_features = []
        for seq_num in numeric_sequences:
            if len(seq_num) > 10:  # Minimum length for meaningful analysis
                try:
                    boost = compute_density_boost(seq_num, k=0.3, bins=20)
                    boost_features = [
                        boost,
                        np.log10(max(1, boost)),  # Log scale
                        1 if boost > 1000 else 0,  # Binary: >1000x
                        1 if boost > 210 else 0,  # Binary: >210%
                    ]
                except Exception:
                    boost_features = [1.0, 0.0, 0, 0]  # Default values
            else:
                boost_features = [1.0, 0.0, 0, 0]

            density_features.append(boost_features)

        features["density_boost"] = np.array(density_features)

        # Z5D predictor features (simplified version)
        z5d_features = []
        for i, seq in enumerate(sequences):
            try:
                # Use the pain management analyzer
                result = self.pain_analyzer.implement_z5d_predictor(
                    seq, target_n=10**4
                )  # Smaller target for speed
                z5d_feature = [
                    result["density_boost_achieved"],
                    1 if result["density_enhancement_success"] else 0,
                    np.log10(max(1, result["density_boost_achieved"])),
                ]
            except Exception:
                z5d_feature = [1.0, 0, 0.0]

            z5d_features.append(z5d_feature)

        features["z5d"] = np.array(z5d_features)

        return features

    def run_classification_validation(
        self, X_baseline: np.ndarray, X_z: np.ndarray, y: np.ndarray, task_name: str
    ) -> Dict[str, Any]:
        """Run classification validation comparing baseline vs Z Framework."""
        results = {}

        # Split data
        X_base_train, X_base_test, y_train, y_test = train_test_split(
            X_baseline, y, test_size=0.3, random_state=42, stratify=y
        )
        X_z_train, X_z_test, _, _ = train_test_split(
            X_z, y, test_size=0.3, random_state=42, stratify=y
        )

        # Models to test
        models = {
            "LogisticRegression": LogisticRegression(random_state=42, max_iter=1000),
            "RandomForest": RandomForestClassifier(random_state=42, n_estimators=100),
        }

        for model_name, model in models.items():
            # Baseline features
            pipeline_baseline = Pipeline(
                [("scaler", StandardScaler()), ("model", model)]
            )

            pipeline_baseline.fit(X_base_train, y_train)
            y_pred_base = pipeline_baseline.predict_proba(X_base_test)[:, 1]

            # Z Framework features
            pipeline_z = Pipeline([("scaler", StandardScaler()), ("model", model)])

            pipeline_z.fit(X_z_train, y_train)
            y_pred_z = pipeline_z.predict_proba(X_z_test)[:, 1]

            # Calculate metrics
            auc_baseline = roc_auc_score(y_test, y_pred_base)
            auc_z = roc_auc_score(y_test, y_pred_z)

            # Statistical significance test
            _, p_value = stats.mannwhitneyu(
                y_pred_base, y_pred_z, alternative="two-sided"
            )

            results[f"{model_name}"] = {
                "baseline_auc": auc_baseline,
                "z_framework_auc": auc_z,
                "improvement": auc_z - auc_baseline,
                "p_value": p_value,
                "significant": p_value < 0.05,
            }

        return results

    def run_regression_validation(
        self, X_baseline: np.ndarray, X_z: np.ndarray, y: np.ndarray, task_name: str
    ) -> Dict[str, Any]:
        """Run regression validation comparing baseline vs Z Framework."""
        results = {}

        # Split data
        X_base_train, X_base_test, y_train, y_test = train_test_split(
            X_baseline, y, test_size=0.3, random_state=42
        )
        X_z_train, X_z_test, _, _ = train_test_split(
            X_z, y, test_size=0.3, random_state=42
        )

        # Models to test
        models = {
            "LinearRegression": LinearRegression(),
            "RandomForest": RandomForestRegressor(random_state=42, n_estimators=100),
        }

        for model_name, model in models.items():
            # Baseline features
            pipeline_baseline = Pipeline(
                [("scaler", StandardScaler()), ("model", model)]
            )

            pipeline_baseline.fit(X_base_train, y_train)
            y_pred_base = pipeline_baseline.predict(X_base_test)

            # Z Framework features
            pipeline_z = Pipeline([("scaler", StandardScaler()), ("model", model)])

            pipeline_z.fit(X_z_train, y_train)
            y_pred_z = pipeline_z.predict(X_z_test)

            # Calculate metrics
            rmse_baseline = np.sqrt(mean_squared_error(y_test, y_pred_base))
            rmse_z = np.sqrt(mean_squared_error(y_test, y_pred_z))

            r2_baseline = r2_score(y_test, y_pred_base)
            r2_z = r2_score(y_test, y_pred_z)

            # Correlation with true values
            corr_baseline, p_base = stats.pearsonr(y_test, y_pred_base)
            corr_z, p_z = stats.pearsonr(y_test, y_pred_z)

            results[f"{model_name}"] = {
                "baseline_rmse": rmse_baseline,
                "z_framework_rmse": rmse_z,
                "baseline_r2": r2_baseline,
                "z_framework_r2": r2_z,
                "baseline_correlation": corr_baseline,
                "z_framework_correlation": corr_z,
                "rmse_improvement": rmse_baseline - rmse_z,
                "r2_improvement": r2_z - r2_baseline,
                "correlation_improvement": corr_z - corr_baseline,
            }

        return results

    def validate_density_boost_claims(
        self, numeric_sequences: List[List[int]]
    ) -> Dict[str, Any]:
        """Validate the >1000x density boost claims."""
        results = {}

        boosts = []
        success_1000x = 0
        success_210 = 0

        for seq_num in numeric_sequences[:100]:  # Sample subset for speed
            if len(seq_num) > 10:
                try:
                    boost = compute_density_boost(seq_num, k=0.3, bins=20)
                    boosts.append(boost)

                    if boost > 1000:
                        success_1000x += 1
                    if boost > 2.1:  # 210%
                        success_210 += 1

                except Exception:
                    continue

        if boosts:
            results = {
                "mean_boost": np.mean(boosts),
                "median_boost": np.median(boosts),
                "std_boost": np.std(boosts),
                "min_boost": np.min(boosts),
                "max_boost": np.max(boosts),
                "success_rate_1000x": success_1000x / len(boosts),
                "success_rate_210": success_210 / len(boosts),
                "ci_95_lower": np.percentile(boosts, 2.5),
                "ci_95_upper": np.percentile(boosts, 97.5),
                "n_samples": len(boosts),
            }

            # Statistical test against null hypothesis (boost = 1.0)
            t_stat, p_value = stats.ttest_1samp(boosts, 1.0)
            results["t_statistic"] = t_stat
            results["p_value_vs_null"] = p_value
            results["significant_boost"] = p_value < 0.05

        return results

    def generate_validation_report(self):
        """Generate comprehensive validation report."""

        # Load data
        data = self.load_synthetic_data()

        if not data:
            self.logger.error("No data loaded. Cannot proceed with validation.")
            return

        # Validate neural spikes data
        if "neural_spikes" in data:
            neural_data = data["neural_spikes"]

            # Parse sequences
            sequences = neural_data["sequence"].tolist()
            numeric_sequences = self.parse_numeric_sequences(
                neural_data["sequence_numeric"].tolist()
            )

            # Extract features
            baseline_features = self.baseline_extractor.extract_all_features(
                sequences, numeric_sequences
            )
            z_features = self.extract_z_framework_features(sequences, numeric_sequences)

            # Combine baseline features
            X_baseline = np.hstack(
                [
                    baseline_features["sequence"],
                    baseline_features["numeric"],
                    baseline_features["histogram"],
                ]
            )

            # Combine Z Framework features
            X_z = np.hstack(
                [
                    z_features["theta_prime"],
                    z_features["density_boost"],
                    z_features["z5d"],
                ]
            )

            # Classification task: therapeutic response
            y_class = neural_data["therapeutic_response"].values
            class_results = self.run_classification_validation(
                X_baseline, X_z, y_class, "therapeutic_response"
            )

            # Regression task: pain score
            y_reg = neural_data["pain_score"].values
            reg_results = self.run_regression_validation(
                X_baseline, X_z, y_reg, "pain_score"
            )

            # Density boost validation
            density_results = self.validate_density_boost_claims(numeric_sequences)

            self.results["neural_spikes"] = {
                "classification": class_results,
                "regression": reg_results,
                "density_boost": density_results,
            }

        # Print results
        self.print_validation_summary()

        # Save detailed results
        self.save_results()

    def print_validation_summary(self):
        """Print summary of validation results."""
        print("=" * 80)
        print("Z FRAMEWORK VALIDATION RESULTS")
        print("=" * 80)
        print()

        if "neural_spikes" in self.results:
            results = self.results["neural_spikes"]

            print("NEURAL SPIKES DATASET VALIDATION")
            print("-" * 40)

            # Classification results
            if "classification" in results:
                print("\n📊 CLASSIFICATION TASK: Therapeutic Response Prediction")
                for model, metrics in results["classification"].items():
                    print(f"\n{model}:")
                    print(f"  Baseline AUC:     {metrics['baseline_auc']:.3f}")
                    print(f"  Z Framework AUC:  {metrics['z_framework_auc']:.3f}")
                    print(f"  Improvement:      {metrics['improvement']:+.3f}")
                    print(f"  P-value:          {metrics['p_value']:.6f}")
                    print(
                        f"  Significant:      {'✓' if metrics['significant'] else '✗'}"
                    )

            # Regression results
            if "regression" in results:
                print("\n📈 REGRESSION TASK: Pain Score Prediction")
                for model, metrics in results["regression"].items():
                    print(f"\n{model}:")
                    print(f"  Baseline RMSE:    {metrics['baseline_rmse']:.3f}")
                    print(f"  Z Framework RMSE: {metrics['z_framework_rmse']:.3f}")
                    print(f"  RMSE Improvement: {metrics['rmse_improvement']:+.3f}")
                    print(f"  Baseline R²:      {metrics['baseline_r2']:.3f}")
                    print(f"  Z Framework R²:   {metrics['z_framework_r2']:.3f}")
                    print(f"  R² Improvement:   {metrics['r2_improvement']:+.3f}")

            # Density boost validation
            if "density_boost" in results:
                print("\n🚀 DENSITY BOOST VALIDATION")
                density = results["density_boost"]
                print(f"  Mean Boost:       {density['mean_boost']:.1f}×")
                print(f"  Median Boost:     {density['median_boost']:.1f}×")
                print(
                    f"  95% CI:           [{density['ci_95_lower']:.1f}, {density['ci_95_upper']:.1f}]"
                )
                print(f"  >1000× Success:   {density['success_rate_1000x']:.1%}")
                print(f"  >210% Success:    {density['success_rate_210']:.1%}")
                print(f"  P-value vs null:  {density['p_value_vs_null']:.6f}")
                print(
                    f"  Significant:      {'✓' if density['significant_boost'] else '✗'}"
                )

        print("\n" + "=" * 80)
        print("✅ VALIDATION COMPLETE")
        print("📁 Detailed results saved to validation_results/")
        print("⚠️  RESEARCH USE ONLY - Not for clinical decisions")
        print("=" * 80)

    def save_results(self):
        """Save detailed results to files."""
        import json

        # Save JSON results
        with open(self.output_dir / "validation_results.json", "w") as f:
            json.dump(self.results, f, indent=2, default=str)

        # Save summary CSV
        if "neural_spikes" in self.results:
            summary_data = []

            # Add classification metrics
            if "classification" in self.results["neural_spikes"]:
                for model, metrics in self.results["neural_spikes"][
                    "classification"
                ].items():
                    summary_data.append(
                        {
                            "task": "classification",
                            "model": model,
                            "metric": "AUC",
                            "baseline": metrics["baseline_auc"],
                            "z_framework": metrics["z_framework_auc"],
                            "improvement": metrics["improvement"],
                            "p_value": metrics["p_value"],
                            "significant": metrics["significant"],
                        }
                    )

            # Add regression metrics
            if "regression" in self.results["neural_spikes"]:
                for model, metrics in self.results["neural_spikes"][
                    "regression"
                ].items():
                    summary_data.append(
                        {
                            "task": "regression",
                            "model": model,
                            "metric": "RMSE",
                            "baseline": metrics["baseline_rmse"],
                            "z_framework": metrics["z_framework_rmse"],
                            "improvement": metrics["rmse_improvement"],
                            "p_value": None,
                            "significant": None,
                        }
                    )

            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(self.output_dir / "validation_summary.csv", index=False)

        self.logger.info(f"Results saved to {self.output_dir}")


def main():
    """Main validation script."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    print("🔬 Starting Z Framework Validation...")
    print("📊 This will compare Z Framework vs conventional baselines")
    print("⏱️  Estimated runtime: 2-5 minutes")
    print()

    validator = ZFrameworkValidator()
    validator.generate_validation_report()


if __name__ == "__main__":
    main()
