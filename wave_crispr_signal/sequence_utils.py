"""Shared nucleotide validation and complex encoding utilities.

This module is the authoritative source for active codepaths that need:
- strict DNA/RNA validation
- deterministic complex encoding for phase-weighted scoring
"""

from __future__ import annotations

from typing import Dict
import numpy as np


_DNA_VALID = frozenset("ACGTN")
_RNA_VALID = frozenset("ACGUN")


def validate_dna_sequence(seq: str, allow_rna: bool = False) -> str:
    """Validate DNA/RNA sequence and return normalized uppercase sequence."""
    seq_upper = seq.upper()

    if allow_rna:
        invalid_chars = set(seq_upper) - _RNA_VALID
        if invalid_chars:
            raise ValueError(
                f"Invalid RNA sequence: contains {invalid_chars}. "
                f"Only A, C, G, U, N allowed for RNA sequences."
            )
        if "T" in seq_upper:
            raise ValueError(
                "Invalid RNA sequence: contains 'T'. RNA sequences should use 'U' instead of 'T'."
            )
    else:
        invalid_chars = set(seq_upper) - _DNA_VALID
        if invalid_chars:
            raise ValueError(
                f"Invalid DNA sequence: contains {invalid_chars}. "
                f"Only A, C, G, T, N allowed for DNA sequences."
            )
        if "U" in seq_upper:
            raise ValueError(
                "Invalid DNA sequence: contains 'U'. DNA sequences should use 'T' instead of 'U'."
            )

    return seq_upper


def encode_complex_phase_weighted(seq: str, is_rna: bool = False) -> np.ndarray:
    """Encode nucleotides into complex values used by phase-weighted scorecard."""
    seq_upper = validate_dna_sequence(seq, allow_rna=is_rna)

    if is_rna:
        mapping: Dict[str, complex] = {
            "A": 1.0 + 0.0j,
            "U": -1.0 + 0.0j,
            "C": 0.0 + 1.0j,
            "G": 0.0 - 1.0j,
            "N": 0.0 + 0.0j,
        }
    else:
        mapping = {
            "A": 1.0 + 0.0j,
            "T": -1.0 + 0.0j,
            "C": 0.0 + 1.0j,
            "G": 0.0 - 1.0j,
            "N": 0.0 + 0.0j,
        }

    return np.array([mapping[base] for base in seq_upper], dtype=np.complex128)
