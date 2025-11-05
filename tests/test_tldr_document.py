"""
Test suite for WAVE_CRISPR_SIGNAL_TLDR.md document validation.

This test ensures the TL;DR document contains all required sections and 
properly references repository artifacts.
"""

import os
import re


def test_tldr_document_exists():
    """Test that the TL;DR document exists."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    assert os.path.exists(tldr_path), "WAVE_CRISPR_SIGNAL_TLDR.md must exist in docs/"


def test_tldr_contains_required_sections():
    """Test that the TL;DR document contains all required sections."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for required sections
    required_sections = [
        "## 0) TL;DR",
        "### What",
        "### Why It Matters",
        "### Key Validated Findings",
        "### Hypotheses to Test Next",
        "### Repository Artifacts",
        "### Reproducibility",
        "### Safety Note"
    ]
    
    for section in required_sections:
        assert section in content, f"TL;DR must contain section: {section}"


def test_tldr_contains_validated_findings():
    """Test that the TL;DR document contains the three key validated findings."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for the three validated findings
    assert "15% density enhancement" in content, "Must include density enhancement finding"
    assert "κ_geo ≈ 0.3" in content, "Must reference κ_geo parameter"
    assert "<0.01% error" in content, "Must include spectral disruption accuracy"
    assert "R² +0.171" in content, "Must include prediction lift finding"
    
    # Check for [Validated] tags
    validated_count = content.count("[Validated]")
    assert validated_count >= 3, f"Must have at least 3 [Validated] tags, found {validated_count}"


def test_tldr_contains_hypothesis_tags():
    """Test that hypotheses are properly tagged."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for [Hypothesis] tags
    hypothesis_count = content.count("[Hypothesis]")
    assert hypothesis_count >= 3, f"Must have at least 3 [Hypothesis] tags, found {hypothesis_count}"


def test_tldr_references_repo_artifacts():
    """Test that the TL;DR document references key repository artifacts."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for key artifact references
    key_artifacts = [
        "EMPIRICAL_FINDINGS_REPORT.md",
        "TOPOLOGICAL_ANALYSIS.md",
        "proof_pack/run_validation.py",
        "proof_pack/quick_validation_demo.py",
        "applications/crispr_cli.py"
    ]
    
    for artifact in key_artifacts:
        assert artifact in content, f"TL;DR must reference artifact: {artifact}"


def test_tldr_has_reproducibility_commands():
    """Test that the TL;DR document includes reproducibility commands."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for command examples
    assert "python proof_pack/run_validation.py" in content, "Must include validation command"
    assert "python proof_pack/quick_validation_demo.py" in content, "Must include quick demo command"
    assert "python bin/bin_resonance_test.py" in content, "Must include resonance test command"


def test_tldr_has_safety_note():
    """Test that the TL;DR document includes proper safety disclaimer."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for safety language
    assert "analytic method" in content.lower(), "Must mention this is an analytic method"
    assert "not a wet-lab protocol" in content.lower() or "not wet-lab" in content.lower(), \
        "Must clarify this is not a wet-lab protocol"
    assert "complement" in content.lower(), "Must indicate it complements other methods"


def test_tldr_has_metadata():
    """Test that the TL;DR document has proper metadata footer."""
    tldr_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'docs',
        'WAVE_CRISPR_SIGNAL_TLDR.md'
    )
    
    with open(tldr_path, 'r') as f:
        content = f.read()
    
    # Check for metadata
    assert "Document Status" in content, "Must include document status"
    assert "Issue Reference" in content, "Must include issue reference"
    assert "#123" in content, "Must reference issue #123"
    assert "Last Updated" in content, "Must include last updated date"


if __name__ == "__main__":
    # Run all tests
    test_tldr_document_exists()
    print("✓ TL;DR document exists")
    
    test_tldr_contains_required_sections()
    print("✓ TL;DR contains all required sections")
    
    test_tldr_contains_validated_findings()
    print("✓ TL;DR contains validated findings")
    
    test_tldr_contains_hypothesis_tags()
    print("✓ TL;DR contains hypothesis tags")
    
    test_tldr_references_repo_artifacts()
    print("✓ TL;DR references repository artifacts")
    
    test_tldr_has_reproducibility_commands()
    print("✓ TL;DR includes reproducibility commands")
    
    test_tldr_has_safety_note()
    print("✓ TL;DR includes safety note")
    
    test_tldr_has_metadata()
    print("✓ TL;DR includes metadata")
    
    print("\n✅ All TL;DR document tests passed!")
