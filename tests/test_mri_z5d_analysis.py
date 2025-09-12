#!/usr/bin/env python3
"""
Test Suite for Z5D MRI Analysis Module

Tests the Z5D geodesic analysis functionality for signal processing,
ensuring mathematical correctness and statistical validity.

Test Categories:
1. Core Z5D geodesic calculations
2. Signal pattern analysis
3. Statistical validation methods
4. Integration with Z Framework
5. Reproducibility and precision
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pytest

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from experiments.mri_z5d_analysis import (
        Z5DGeodeskAnalyzer,
        Z5DAnalysisResult,
        StatisticalSummary,
        generate_synthetic_mri_signals,
        load_dicom_signals
    )
    IMPORT_SUCCESS = True
except ImportError as e:
    print(f"Import error: {e}")
    IMPORT_SUCCESS = False


class TestZ5DGeodeskAnalyzer:
    """Test the Z5D Geodesic Analyzer."""
    
    def test_import_success(self):
        """Test that the module imports successfully."""
        assert IMPORT_SUCCESS, "Failed to import Z5D MRI analysis module"
    
    def test_analyzer_initialization(self):
        """Test analyzer initialization."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer(seed=42, precision_dps=20)
        assert analyzer.seed == 42
        assert analyzer.precision_dps == 20
    
    def test_theta_prime_geodesic_calculation(self):
        """Test Z5D theta prime geodesic calculation."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer()
        
        # Test basic calculation
        result = analyzer.theta_prime_geodesic(0.5)
        assert isinstance(result, float)
        assert 0.0 <= result <= 2.0  # Reasonable bounds
        
        # Test with different k parameter
        result_k1 = analyzer.theta_prime_geodesic(0.5, k=0.1)
        result_k2 = analyzer.theta_prime_geodesic(0.5, k=0.2)
        assert result_k1 != result_k2  # Different k should give different results
        
        # Test edge cases
        result_zero = analyzer.theta_prime_geodesic(0.0)
        assert result_zero >= 0.0
        
        result_one = analyzer.theta_prime_geodesic(1.0)
        assert result_one >= 0.0
    
    def test_signal_pattern_analysis(self):
        """Test signal pattern analysis functionality."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Create test signal
        signal = np.array([0.1, 0.5, 0.8, 0.6, 0.3, 0.7, 0.9, 0.4, 0.2])
        
        result = analyzer.analyze_signal_pattern(signal, "test_signal")
        
        # Validate result structure
        assert isinstance(result, Z5DAnalysisResult)
        assert result.sample_id == "test_signal"
        assert result.n_points == len(signal)
        assert isinstance(result.theta_prime_mean, float)
        assert isinstance(result.theta_prime_std, float)
        assert isinstance(result.geodesic_correlation, float)
        assert isinstance(result.focal_accuracy, float)
        assert result.processing_time_ms > 0
        assert result.classification in ["high-coherence", "moderate-coherence", "low-coherence"]
        
        # Test with different signal types
        # High coherence signal
        high_coherence_signal = np.ones(10) * 0.8
        result_high = analyzer.analyze_signal_pattern(high_coherence_signal)
        
        # Low coherence signal (random)
        np.random.seed(42)
        low_coherence_signal = np.random.rand(10)
        result_low = analyzer.analyze_signal_pattern(low_coherence_signal)
        
        # High coherence should have lower std
        assert result_high.theta_prime_std <= result_low.theta_prime_std
    
    def test_bootstrap_correlation_ci(self):
        """Test bootstrap confidence interval calculation."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Create correlated data
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])  # Perfect correlation
        
        ci_low, ci_high = analyzer.bootstrap_correlation_ci(x, y, n_bootstrap=100)
        
        assert isinstance(ci_low, float)
        assert isinstance(ci_high, float)
        assert ci_low <= ci_high
        assert ci_high > 0.8  # Should be high for perfectly correlated data
        
        # Test with uncorrelated data
        x_uncorr = np.array([1, 2, 3, 4, 5])
        y_uncorr = np.array([5, 1, 4, 2, 3])  # No clear correlation
        
        ci_low_uncorr, ci_high_uncorr = analyzer.bootstrap_correlation_ci(x_uncorr, y_uncorr, n_bootstrap=100)
        
        # Uncorrelated data should have wider CI around zero
        assert abs(ci_low_uncorr) < abs(ci_low)
        assert abs(ci_high_uncorr) < abs(ci_high)
    
    def test_permutation_test(self):
        """Test permutation test for correlation significance."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Perfectly correlated data
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        
        p_value = analyzer.permutation_test(x, y, n_permutation=100)
        
        assert isinstance(p_value, float)
        assert 0.0 <= p_value <= 1.0
        assert p_value < 0.1  # Should be significant for perfectly correlated data
        
        # Random uncorrelated data
        np.random.seed(42)
        x_rand = np.random.randn(50)
        y_rand = np.random.randn(50)
        
        p_value_rand = analyzer.permutation_test(x_rand, y_rand, n_permutation=100)
        assert p_value_rand > p_value  # Random data should have higher p-value
    
    def test_statistical_validation(self):
        """Test comprehensive statistical validation."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Create mock results
        results = []
        for i in range(10):
            result = Z5DAnalysisResult(
                sample_id=f"test_{i}",
                n_points=100,
                theta_prime_mean=0.5 + 0.1 * i,
                theta_prime_std=0.1,
                geodesic_correlation=0.7 + 0.02 * i,
                focal_accuracy=0.6 + 0.03 * i,
                processing_time_ms=10.0,
                k_parameter=0.04449,
                classification="moderate-coherence"
            )
            results.append(result)
        
        stats_summary = analyzer.statistical_validation(results, n_bootstrap=50, n_permutation=50)
        
        assert isinstance(stats_summary, StatisticalSummary)
        assert isinstance(stats_summary.pearson_r, float)
        assert isinstance(stats_summary.pearson_p, float)
        assert isinstance(stats_summary.bootstrap_ci_low, float)
        assert isinstance(stats_summary.bootstrap_ci_high, float)
        assert isinstance(stats_summary.effect_size_cohens_d, float)
        assert isinstance(stats_summary.permutation_p, float)
        assert stats_summary.n_bootstrap == 50
        assert stats_summary.n_permutation == 50
        
        # With insufficient data
        stats_insufficient = analyzer.statistical_validation(results[:1])
        assert stats_insufficient.pearson_r == 0.0
        assert stats_insufficient.pearson_p == 1.0


class TestSyntheticSignalGeneration:
    """Test synthetic MRI signal generation."""
    
    def test_generate_synthetic_signals(self):
        """Test synthetic signal generation."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        signals = generate_synthetic_mri_signals(n_samples=10, signal_length=50, seed=42)
        
        assert len(signals) == 10
        assert all(len(signal) == 50 for signal in signals)
        assert all(isinstance(signal, np.ndarray) for signal in signals)
        assert all(np.all(signal >= 0) for signal in signals)  # All signals should be non-negative
        
        # Test reproducibility
        signals1 = generate_synthetic_mri_signals(n_samples=5, signal_length=20, seed=123)
        signals2 = generate_synthetic_mri_signals(n_samples=5, signal_length=20, seed=123)
        
        for s1, s2 in zip(signals1, signals2):
            np.testing.assert_array_almost_equal(s1, s2)


class TestDicomSignalLoading:
    """Test DICOM signal loading functionality."""
    
    def test_load_dicom_signals(self):
        """Test DICOM signal loading."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
        
        # Test DICOM loading with the actual data directories
        cervical_dir = "data/MRI__CERVICAL_SPINE_W_O_CONT/DICOM"
        thoracic_dir = "data/MRI__THORACIC_SPINE_W_O_CONT/DICOM"
        
        # Only test if DICOM directories exist
        if not (os.path.exists(cervical_dir) and os.path.exists(thoracic_dir)):
            pytest.skip("DICOM data directories not found")
        
        signals = load_dicom_signals(
            cervical_dir=cervical_dir,
            thoracic_dir=thoracic_dir,
            signal_length=32,  # Small for testing
            max_files_per_series=2  # Limit for testing
        )
        
        assert len(signals) > 0, "Should load at least some signals"
        
        for signal in signals:
            assert isinstance(signal, np.ndarray)
            assert signal.shape == (32,), f"Expected shape (32,), got {signal.shape}"
            assert np.all(signal >= 0), "Signal values should be non-negative"
            assert np.all(signal <= 1), "Signal values should be normalized to [0,1]"
            assert not np.all(signal == signal[0]), "Signal should not be constant"
    
    def test_dicom_loading_missing_directories(self):
        """Test DICOM loading with missing directories."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
        
        # Test with non-existent directories
        signals = load_dicom_signals(
            cervical_dir="/nonexistent/path1",
            thoracic_dir="/nonexistent/path2",
            signal_length=16,
            max_files_per_series=1
        )
        
        # Should return fallback signals
        assert len(signals) >= 2, "Should provide fallback signals when no DICOM files found"
    
    def test_dicom_vs_synthetic_compatibility(self):
        """Test that DICOM and synthetic signals produce compatible analysis results."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
        
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Test with synthetic signals
        synthetic_signals = generate_synthetic_mri_signals(n_samples=3, signal_length=16, seed=42)
        synthetic_results = []
        for i, signal in enumerate(synthetic_signals):
            result = analyzer.analyze_signal_pattern(signal, f"synthetic_{i}")
            synthetic_results.append(result)
        
        # Test with DICOM signals if available, otherwise use fallback
        cervical_dir = "data/MRI__CERVICAL_SPINE_W_O_CONT/DICOM" 
        thoracic_dir = "data/MRI__THORACIC_SPINE_W_O_CONT/DICOM"
        
        dicom_signals = load_dicom_signals(
            cervical_dir=cervical_dir,
            thoracic_dir=thoracic_dir,
            signal_length=16,
            max_files_per_series=1
        )
        
        dicom_results = []
        for i, signal in enumerate(dicom_signals[:3]):  # Limit to 3 for comparison
            result = analyzer.analyze_signal_pattern(signal, f"dicom_{i}")
            dicom_results.append(result)
        
        # Both should produce valid results
        assert len(synthetic_results) > 0
        assert len(dicom_results) > 0
        
        # Results should have the same structure
        for result in synthetic_results + dicom_results:
            assert isinstance(result, Z5DAnalysisResult)
            assert isinstance(result.theta_prime_mean, float)
            assert isinstance(result.focal_accuracy, float)
            assert result.classification in ['high-coherence', 'moderate-coherence', 'low-coherence']


class TestIntegrationAndReproducibility:
    """Test integration and reproducibility aspects."""
    
    def test_reproducibility_with_seed(self):
        """Test that results are reproducible with same seed."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        # Same configuration
        analyzer1 = Z5DGeodeskAnalyzer(seed=42)
        analyzer2 = Z5DGeodeskAnalyzer(seed=42)
        
        # Same test signal
        signal = np.array([0.1, 0.5, 0.8, 0.6, 0.3])
        
        result1 = analyzer1.analyze_signal_pattern(signal, "test")
        result2 = analyzer2.analyze_signal_pattern(signal, "test")
        
        # Results should be identical
        assert abs(result1.theta_prime_mean - result2.theta_prime_mean) < 1e-10
        assert abs(result1.theta_prime_std - result2.theta_prime_std) < 1e-10
        assert abs(result1.geodesic_correlation - result2.geodesic_correlation) < 1e-10
        assert result1.classification == result2.classification
    
    def test_z_framework_integration(self):
        """Test integration with Z Framework."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        analyzer = Z5DGeodeskAnalyzer()
        
        # Verify Z Framework calculator is initialized
        assert analyzer.z_calc is not None
        assert hasattr(analyzer.z_calc, 'calculate_z_value')
    
    def test_precision_parameters(self):
        """Test high-precision parameter handling."""
        if not IMPORT_SUCCESS:
            pytest.skip("Module import failed")
            
        # Test different precision levels
        analyzer_low = Z5DGeodeskAnalyzer(precision_dps=10)
        analyzer_high = Z5DGeodeskAnalyzer(precision_dps=30)
        
        # Both should work without errors
        result_low = analyzer_low.theta_prime_geodesic(0.5)
        result_high = analyzer_high.theta_prime_geodesic(0.5)
        
        assert isinstance(result_low, float)
        assert isinstance(result_high, float)
        
        # Results should be close but may differ slightly due to precision
        assert abs(result_low - result_high) < 0.01


def run_smoke_test():
    """Run a quick smoke test for CI."""
    print("Running Z5D MRI Analysis smoke test...")
    
    if not IMPORT_SUCCESS:
        print("❌ Import failed")
        return False
    
    try:
        # Quick functional test with both synthetic and DICOM signals
        analyzer = Z5DGeodeskAnalyzer(seed=42)
        
        # Test 1: Synthetic signals
        print("Testing synthetic signal generation...")
        synthetic_signals = generate_synthetic_mri_signals(n_samples=3, signal_length=16, seed=42)
        
        synthetic_results = []
        for i, signal in enumerate(synthetic_signals):
            result = analyzer.analyze_signal_pattern(signal, f"smoke_synthetic_{i}")
            synthetic_results.append(result)
        
        # Test 2: DICOM signals (if available)
        print("Testing DICOM signal loading...")
        dicom_signals = load_dicom_signals(
            signal_length=16,
            max_files_per_series=2
        )
        
        dicom_results = []
        for i, signal in enumerate(dicom_signals[:3]):
            result = analyzer.analyze_signal_pattern(signal, f"smoke_dicom_{i}")
            dicom_results.append(result)
        
        # Quick validation with combined results
        all_results = synthetic_results + dicom_results
        stats = analyzer.statistical_validation(all_results, n_bootstrap=10, n_permutation=10)
        
        print(f"✓ Processed {len(synthetic_results)} synthetic signals")
        print(f"✓ Processed {len(dicom_results)} DICOM signals")
        print(f"✓ Statistical validation complete (r={stats.pearson_r:.3f})")
        print("✓ Z5D MRI Analysis smoke test passed")
        return True
        
    except Exception as e:
        print(f"❌ Smoke test failed: {e}")
        return False


if __name__ == "__main__":
    # Run as standalone test
    import sys
    
    # Check for smoke test argument
    if len(sys.argv) > 1 and "--smoke" in sys.argv:
        success = run_smoke_test()
        sys.exit(0 if success else 1)
    
    # Run pytest if available, otherwise run basic tests
    try:
        import pytest
        print("Running tests with pytest...")
        pytest.main([__file__, "-v"])
    except ImportError:
        print("Running basic tests...")
        
        # Manual test execution
        test_class = TestZ5DGeodeskAnalyzer()
        
        tests = [
            test_class.test_import_success,
            test_class.test_analyzer_initialization,
            test_class.test_theta_prime_geodesic_calculation,
            test_class.test_signal_pattern_analysis,
            test_class.test_bootstrap_correlation_ci,
            test_class.test_permutation_test,
            test_class.test_statistical_validation
        ]
        
        passed = 0
        total = len(tests)
        
        for test in tests:
            try:
                test()
                print(f"✓ {test.__name__}")
                passed += 1
            except Exception as e:
                print(f"❌ {test.__name__}: {e}")
        
        print(f"\nResults: {passed}/{total} tests passed")
        
        if passed == total:
            print("🎉 All tests passed!")
        else:
            print("💥 Some tests failed!")
            sys.exit(1)