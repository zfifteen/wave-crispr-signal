#!/usr/bin/env python3
"""
Smoke test for signal-theoretic CRISPR experiment.

Quick validation that all components work correctly with minimal data.
FAIL-FAST: Requires hg38 reference to be available.
"""

import sys
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import yaml

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from experiments.signal_theoretic_crispr import (
    HumanFASTAValidator,
    BaselinePipeline, 
    SpectralFeatureExtractor,
    ExperimentalStatistics
)
from experiments.signal_theoretic_crispr.main import require_hg38


def main():
    """Run all smoke tests."""
    print("Signal-Theoretic CRISPR Experiment - Smoke Test")
    print("=" * 60)
    
    # Fail fast if hg38 not present
    try:
        hg38_path = require_hg38()
        print(f"✓ hg38 reference validated: {hg38_path}")
    except Exception as e:
        print(f"[FAIL-FAST] hg38 not available: {e}", file=sys.stderr)
        print("Set HG38_FA environment variable to path of hg38 FASTA file.", file=sys.stderr)
        sys.exit(2)
    
    tests = [
        ("Validation", test_validation),
        ("Baseline Pipeline", test_baseline_pipeline),
        ("Spectral Pipeline", test_spectral_pipeline),
        ("Statistics", test_statistics),
        ("Integration", test_integration)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        if test_func():
            passed += 1
        else:
            print(f"  ✗ {test_name} FAILED")
    
    print("\n" + "=" * 60)
    print(f"SMOKE TEST SUMMARY: {passed}/{total} tests passed")
    
    if passed == total:
        print("✓ All smoke tests PASSED - Experiment ready for deployment")
        return True
    else:
        print("✗ Some smoke tests FAILED - Review implementation")
        return False


def create_test_data():
    """Create minimal test dataset."""
    # Create test sequences with known patterns and genomic coordinates
    test_data = pd.DataFrame({
        'sequence': [
            'ATCGATCGATCGATCGATGG',  # 20mer + NGG PAM-like ending
            'GGGGGGGGGGGGGGGGGGGG',  # High GC
            'AAAAAAAAAAAAAAAAAAAA',  # Low GC
            'CCCCCCCCCCCCCCCCCCCC',  # High GC
            'TATATATATATATATATATAT',  # Alternating
            'CGCGCGCGCGCGCGCGCGCG',  # Alternating high GC
        ],
        'efficiency': [0.72, 0.85, 0.45, 0.91, 0.59, 0.83],
        'chromosome': ['chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
        'position': [1000, 1500, 2000, 1000, 1500, 2000]
    })
    
    return test_data


def create_test_manifest(temp_dir):
    """Create test experiment manifest."""
    manifest = {
        'name': 'signal_theoretic_crispr_smoke_test',
        'version': '1.0.0',
        'description': 'Smoke test for signal-theoretic CRISPR experiment',
        'metadata': {
            'experiment_type': 'two_arm_benchmark',
            'domain': 'crispr_genomics',
            'organism': 'homo_sapiens'
        },
        'parameters': {
            'default_seed': 42,
            'bootstrap_iterations': 10,  # Minimal for speed
            'permutation_iterations': 10,
            'window_size': 201,
            'geodesic_k': 0.3,
            'golden_ratio': 1.618033988749895,
            'pam_sequence': 'NGG'
        },
        'datasets': {
            'primary': [{
                'name': 'test_data',
                'path': str(temp_dir / 'test_data.csv'),
                'type': 'on_target_efficiency',
                'organism': 'homo_sapiens'
            }]
        }
    }
    
    return manifest


def test_validation():
    """Test validation module."""
    print("Testing validation module...")
    
    validator = HumanFASTAValidator()
    
    # Test valid sequences
    try:
        validator.validate_guide_sequence("ATCGATCGATCGATCGAT", "DNA")
        validator.validate_guide_sequence("AUCGAUCGAUCGAUCGAU", "RNA")
        print("  ✓ Valid sequences accepted")
    except Exception as e:
        print(f"  ✗ Valid sequence validation failed: {e}")
        return False
    
    # Test invalid sequences
    try:
        validator.validate_guide_sequence("ATCGATCGATCGATCGAX", "DNA")
        print("  ✗ Invalid sequence should have been rejected")
        return False
    except:
        print("  ✓ Invalid sequences correctly rejected")
    
    return True


def test_baseline_pipeline():
    """Test baseline pipeline."""
    print("Testing baseline pipeline...")
    
    test_data = create_test_data()
    pipeline = BaselinePipeline(seed=42)
    
    try:
        # Test processing
        processed_df = pipeline.process_data(test_data)
        
        # Check features were extracted
        feature_cols = [col for col in processed_df.columns if pd.api.types.is_numeric_dtype(processed_df[col])]
        if len(feature_cols) < 5:
            print(f"  ✗ Too few features extracted: {len(feature_cols)}")
            return False
        
        print(f"  ✓ Extracted {len(feature_cols)} baseline features")
        
        # Test model training
        from experiments.signal_theoretic_crispr.baseline import BaselineModels
        models = BaselineModels(seed=42)
        metrics = models.train_on_target(processed_df)
        
        if 'train_r2' not in metrics:
            print("  ✗ Model training failed")
            return False
        
        print(f"  ✓ Model training completed (R² = {metrics['train_r2']:.3f})")
        return True
        
    except Exception as e:
        print(f"  ✗ Baseline pipeline failed: {e}")
        return False


def test_spectral_pipeline():
    """Test spectral pipeline."""
    print("Testing spectral pipeline...")
    
    test_data = create_test_data()
    extractor = SpectralFeatureExtractor(seed=42)
    
    try:
        # Test feature extraction with coordinates
        test_guide = {
            'guide_sequence': test_data.iloc[0]['sequence'],
            'chromosome': test_data.iloc[0]['chromosome'],
            'position': test_data.iloc[0]['position']
        }
        
        features = extractor.extract_features(test_guide)
        
        # Check for expected feature types
        bio_features = [k for k in features.keys() if k.startswith('bio_')]
        arb_features = [k for k in features.keys() if k.startswith('arb_')]
        z_features = [k for k in features.keys() if k.startswith('z_')]
        
        if not (bio_features and arb_features and z_features):
            print(f"  ✗ Missing feature types: bio={len(bio_features)}, arb={len(arb_features)}, z={len(z_features)}")
            return False
        
        print(f"  ✓ Extracted spectral features: bio={len(bio_features)}, arb={len(arb_features)}, z={len(z_features)}")
        return True
        
    except Exception as e:
        print(f"  ✗ Spectral pipeline failed: {e}")
        return False


def test_statistics():
    """Test statistics module."""
    print("Testing statistics module...")
    
    try:
        stats = ExperimentalStatistics(bootstrap_iterations=10, permutation_iterations=10, seed=42)
        
        # Generate test data
        np.random.seed(42)
        n_samples = 20
        y_true = np.random.rand(n_samples)
        y_pred_baseline = y_true + np.random.normal(0, 0.2, n_samples)
        y_pred_spectral = y_true + np.random.normal(0, 0.15, n_samples)
        
        # Test regression evaluation
        results = stats.evaluate_regression_performance(y_true, y_pred_baseline, y_pred_spectral)
        
        required_keys = ['baseline', 'spectral', 'comparison']
        if not all(key in results for key in required_keys):
            print(f"  ✗ Missing result keys: {list(results.keys())}")
            return False
        
        # Check bootstrap CI structure
        baseline_r2 = results['baseline']['r2']
        required_ci_keys = ['mean', 'ci_lower', 'ci_upper']
        if not all(key in baseline_r2 for key in required_ci_keys):
            print(f"  ✗ Missing CI keys: {list(baseline_r2.keys())}")
            return False
        
        print(f"  ✓ Statistics calculation completed")
        print(f"    Baseline R²: {baseline_r2['mean']:.3f} [{baseline_r2['ci_lower']:.3f}, {baseline_r2['ci_upper']:.3f}]")
        return True
        
    except Exception as e:
        print(f"  ✗ Statistics module failed: {e}")
        return False


def test_integration():
    """Test integrated workflow."""
    print("Testing integrated workflow...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        
        try:
            # Create test data
            test_data = create_test_data()
            data_file = temp_dir / 'test_data.csv'
            test_data.to_csv(data_file, index=False)
            
            # Create test manifest
            manifest = create_test_manifest(temp_dir)
            manifest_file = temp_dir / 'manifest.yml'
            with open(manifest_file, 'w') as f:
                yaml.dump(manifest, f)
            
            # Test complete workflow (minimal version)
            from experiments.signal_theoretic_crispr.main import SignalTheoreticExperiment
            
            output_dir = temp_dir / 'output'
            experiment = SignalTheoreticExperiment(
                config_path=str(manifest_file),
                output_dir=str(output_dir),
                seed=42
            )
            
            # Test individual components
            datasets = experiment.load_and_validate_data()
            if 'test_data' not in datasets:
                print("  ✗ Data loading failed")
                return False
            
            baseline_results = experiment.run_baseline_arm(datasets)
            if not baseline_results['features_extracted']:
                print("  ✗ Baseline arm failed")
                return False
            
            spectral_results = experiment.run_spectral_arm(datasets)
            if not spectral_results['features_extracted']:
                print("  ✗ Spectral arm failed")
                return False
            
            comparison_results = experiment.compare_arms(baseline_results, spectral_results)
            if not comparison_results['statistical_tests']:
                print("  ✗ Comparison failed")
                return False
            
            print("  ✓ Integrated workflow completed successfully")
            return True
            
        except Exception as e:
            print(f"  ✗ Integration test failed: {e}")
            import traceback
            traceback.print_exc()
            return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)