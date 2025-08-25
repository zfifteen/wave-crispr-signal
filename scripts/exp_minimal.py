#!/usr/bin/env python3
"""
Minimal experimental script for smoke testing biological vs arbitrary encodings.

This script demonstrates preference for local FASTA windows with no network calls,
suitable for CI environments and quick validation of the encoding framework.
"""

import os
import sys
from typing import Optional


# Prefer local FASTA windows; no network in smoke test.
CANDIDATE_FASTAS = [
    "data/ANO1_rs2186797_GRCh38.fasta",
    "data/MYH11_rs6498573_GRCh38.fasta",
    "data/ELN_rs11770437_GRCh38.fasta",
]


def load_first_fasta() -> str:
    """Load the first available local FASTA file."""
    import json
    
    for p in CANDIDATE_FASTAS:
        if os.path.exists(p):
            with open(p, "r", encoding="utf-8") as fh:
                content = fh.read().strip()
                
                # Try JSON format first (our actual format)
                try:
                    data = json.loads(content)
                    if "seq" in data:
                        seq = data["seq"].upper()
                        if not seq:
                            raise RuntimeError(f"Empty sequence in JSON file {p}")
                        print(f"[info] Using local FASTA (JSON format): {p}")
                        return seq
                except json.JSONDecodeError:
                    pass
                
                # Fall back to standard FASTA format
                lines = content.split('\n')
                seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
                if not seq:
                    raise RuntimeError(f"Empty FASTA sequence in {p}")
                print(f"[info] Using local FASTA: {p}")
                return seq
    raise SystemExit("[error] No local FASTA found in data/. Expected one of: " + ", ".join(CANDIDATE_FASTAS))


def create_ref_alt_sequences(ref_window: str, snp_position: Optional[int] = None) -> tuple[str, str]:
    """
    Create REF and ALT sequences from a reference window by substituting alleles.
    
    Args:
        ref_window: Reference sequence
        snp_position: Position to substitute (default: middle of sequence)
    
    Returns:
        Tuple of (REF sequence, ALT sequence)
    """
    if snp_position is None:
        snp_position = len(ref_window) // 2
        
    if snp_position >= len(ref_window):
        raise ValueError(f"SNP position {snp_position} is beyond sequence length {len(ref_window)}")
    
    # Get the current nucleotide at the position
    current_allele = ref_window[snp_position]
    
    # Choose an alternative allele (different from current)
    allele_map = {"A": "G", "T": "C", "G": "A", "C": "T"}
    alt_allele = allele_map.get(current_allele, "A")
    
    # Create ALT sequence by substituting at the position
    alt_window = ref_window[:snp_position] + alt_allele + ref_window[snp_position + 1:]
    
    return ref_window, alt_window


def basic_sequence_analysis(sequence: str) -> dict:
    """Perform basic sequence analysis for smoke testing."""
    analysis = {
        "length": len(sequence),
        "gc_content": (sequence.count("G") + sequence.count("C")) / len(sequence) * 100,
        "nucleotide_counts": {
            "A": sequence.count("A"),
            "T": sequence.count("T"),
            "G": sequence.count("G"),
            "C": sequence.count("C"),
        }
    }
    return analysis


def main():
    """Main smoke test function."""
    print("🧬 CRISPR Signal Analysis - Minimal Smoke Test")
    print("=" * 50)
    
    try:
        # Load reference window from local FASTA (no network calls)
        ref_window = load_first_fasta()
        print(f"[info] Loaded reference sequence: {len(ref_window)} bp")
        
        # Create REF/ALT sequences by substituting allele at known offset
        ref_seq, alt_seq = create_ref_alt_sequences(ref_window)
        
        print(f"[info] Created REF sequence: {len(ref_seq)} bp")
        print(f"[info] Created ALT sequence: {len(alt_seq)} bp")
        
        # Basic analysis for both sequences
        ref_analysis = basic_sequence_analysis(ref_seq)
        alt_analysis = basic_sequence_analysis(alt_seq)
        
        print("\n📊 Reference Sequence Analysis:")
        print(f"  Length: {ref_analysis['length']} bp")
        print(f"  GC Content: {ref_analysis['gc_content']:.1f}%")
        print(f"  Nucleotide counts: {ref_analysis['nucleotide_counts']}")
        
        print("\n📊 Alternative Sequence Analysis:")
        print(f"  Length: {alt_analysis['length']} bp")
        print(f"  GC Content: {alt_analysis['gc_content']:.1f}%")
        print(f"  Nucleotide counts: {alt_analysis['nucleotide_counts']}")
        
        # Verify sequences differ by exactly one nucleotide
        differences = sum(1 for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a)
        print(f"\n🔍 Sequence differences: {differences} nucleotides")
        
        if differences == 1:
            print("✓ Successfully created valid REF/ALT pair")
        else:
            print(f"⚠ Warning: Expected 1 difference, found {differences}")
        
        print("\n✓ Smoke test completed successfully - no network calls made")
        return True
        
    except Exception as e:
        print(f"\n❌ Smoke test failed: {e}")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)