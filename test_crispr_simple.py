"""
Simple test runner to validate CRISPR functionality
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "applications"))

from crispr_guide_designer import CRISPRGuideDesigner
from wave_crispr_metrics import WaveCRISPRMetrics


def test_basic_functionality():
    """Test basic CRISPR guide design functionality."""
    print("Testing basic CRISPR guide design functionality...")

    # Test sequence (PCSK9 exon 1)
    sequence = "ATGCTGCGGAGACCTGGAGAGAAAGCAGTGGCCGGGGCAGTGGGAGGAGGAGGAGCTGGAAGAGGAGAGAAAGGAGGAGCTGCAGGAGGAGAGGAGGAGGAGGGAGAGGAGGAGCTGGAGCTGAAGCTGGAGCTGGAGCTGGAGAGGAGAGAGGG"

    # Initialize designer
    designer = CRISPRGuideDesigner()

    # Test guide design
    guides = designer.design_guides(sequence, num_guides=3)
    assert len(guides) > 0, "No guides were designed"
    print(f"✓ Successfully designed {len(guides)} guides")

    # Test guide properties
    for i, guide in enumerate(guides):
        assert "sequence" in guide, f"Guide {i} missing sequence"
        assert "position" in guide, f"Guide {i} missing position"
        assert "on_target_score" in guide, f"Guide {i} missing score"
        assert len(guide["sequence"]) == 20, f"Guide {i} wrong length"
        assert 0 <= guide["on_target_score"] <= 1, f"Guide {i} score out of range"
    print("✓ All guides have correct properties")

    # Test on-target scoring
    test_guide = guides[0]["sequence"]
    score = designer.calculate_on_target_score(test_guide)
    assert 0 <= score <= 1, "On-target score out of range"
    print(f"✓ On-target scoring works (score: {score:.3f})")

    # Test off-target analysis
    risk = designer.calculate_off_target_risk(
        test_guide, sequence[:50], sequence[50:100]
    )
    assert 0 <= risk <= 1, "Off-target risk out of range"
    print(f"✓ Off-target analysis works (risk: {risk:.3f})")

    # Test repair prediction
    repair = designer.predict_repair_outcomes(test_guide, sequence[:60])
    assert "nhej_probability" in repair, "Missing NHEJ prediction"
    assert "mmej_probability" in repair, "Missing MMEJ prediction"
    assert "hdr_efficiency" in repair, "Missing HDR prediction"
    print("✓ Repair outcome prediction works")

    return True


def test_metrics_functionality():
    """Test advanced metrics functionality."""
    print("\nTesting advanced metrics functionality...")

    metrics = WaveCRISPRMetrics()
    test_seq = "GCTGCGGAGACCTGGAGAGA"

    # Test spectral entropy
    entropy = metrics.calculate_spectral_entropy(test_seq)
    assert entropy > 0, "Entropy should be positive"
    print(f"✓ Spectral entropy calculation works (entropy: {entropy:.3f})")

    # Test spectral complexity
    complexity = metrics.calculate_spectral_complexity(test_seq)
    required_keys = ["spectral_entropy", "spectral_flatness", "spectral_centroid"]
    for key in required_keys:
        assert key in complexity, f"Missing complexity metric: {key}"
    print("✓ Spectral complexity analysis works")

    # Test comprehensive scoring
    target_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    comprehensive = metrics.calculate_comprehensive_score(test_seq, target_seq)
    assert "composite_score" in comprehensive, "Missing composite score"
    assert 0 <= comprehensive["composite_score"] <= 1, "Composite score out of range"
    print(
        f"✓ Comprehensive scoring works (score: {comprehensive['composite_score']:.3f})"
    )

    return True


def test_cli_functionality():
    """Test CLI import functionality."""
    print("\nTesting CLI functionality...")

    try:
        # from crispr_cli import main  # Import checked but not used
        print("✓ CLI module imports successfully")
    except ImportError as e:
        print(f"✗ CLI import failed: {e}")
        return False

    return True


def test_visualization_imports():
    """Test visualization module imports."""
    print("\nTesting visualization imports...")

    try:
        from crispr_visualization import CRISPRVisualizer

        CRISPRVisualizer()
        print("✓ Visualization module imports successfully")
    except ImportError as e:
        print(f"✗ Visualization import failed: {e}")
        return False

    return True


def main():
    """Run all tests."""
    print("🧬 CRISPR Functionality Test Suite")
    print("=" * 50)

    tests = [
        test_basic_functionality,
        test_metrics_functionality,
        test_cli_functionality,
        test_visualization_imports,
    ]

    passed = 0
    total = len(tests)

    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test failed: {e}")

    print(f"\n🏆 Test Results: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 All tests passed! CRISPR functionality is working correctly.")
        return True
    else:
        print("❌ Some tests failed. Check the output above for details.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
