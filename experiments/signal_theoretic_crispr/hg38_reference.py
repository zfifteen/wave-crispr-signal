#!/usr/bin/env python3
"""
hg38 Reference Genome Interface for CRISPR experiment.

This module provides functionality to extract genomic windows from the hg38 reference
for contextual analysis of CRISPR guide sequences.

FAIL-FAST MODE: No fallback sequences - requires valid hg38 reference.
"""

import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any
from pyfaidx import Fasta

ALLOWED = set("ACGTN")


class Hg38Extractor:
    """Safe hg38 genomic window extractor with fail-fast validation."""
    
    def __init__(self, fasta_path: str):
        """
        Initialize hg38 extractor.
        
        Args:
            fasta_path: Path to hg38 FASTA file
            
        Raises:
            RuntimeError: If FASTA cannot be opened or indexed
        """
        self.fasta_path = fasta_path
        self.logger = logging.getLogger(__name__)
        
        try:
            self.fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
            self.logger.info(f"Loaded hg38 reference: {fasta_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to open hg38 FASTA at {fasta_path}: {e}")
    
    def window_201(self, chrom: str, pos: int, L: int = 201) -> str:
        """
        Extract safe 201-bp genomic window with clipping and padding.
        
        Args:
            chrom: Chromosome name (e.g., 'chr1', '1')
            pos: 1-based genomic position
            L: Window length (must be odd)
            
        Returns:
            Genomic sequence window of exact length L
            
        Raises:
            ValueError: If window length is even
            RuntimeError: If chromosome not found or sequence invalid
        """
        if L % 2 == 0:
            raise ValueError(f"Context length must be odd, got {L}")
        
        half = L // 2
        
        # Normalize chromosome name
        chrom = self._normalize_chromosome(chrom)
        
        if chrom not in self.fa:
            raise RuntimeError(f"Chromosome {chrom} not found in hg38 reference")
        
        chrom_len = len(self.fa[chrom])
        
        # Calculate safe boundaries with clipping
        start = max(1, pos - half)
        end = min(chrom_len, pos + half)
        
        # Extract sequence
        seq = str(self.fa[chrom][start-1:end])
        
        # Calculate padding needed
        pad_left = max(0, (pos - half) - start + 1) if (pos - half) < 1 else 0
        pad_right = max(0, L - (len(seq) + pad_left))
        
        # Apply padding with N
        if pad_left > 0:
            seq = ("N" * pad_left) + seq
        if pad_right > 0:
            seq = seq + ("N" * pad_right)
        
        # Validate final length
        if len(seq) != L:
            raise RuntimeError(f"Expected {L} nt, got {len(seq)} for {chrom}:{pos}")
        
        # Validate nucleotides
        if set(seq) - ALLOWED:
            invalid_chars = set(seq) - ALLOWED
            raise RuntimeError(f"Non-ACGTN characters in hg38 window: {invalid_chars}")
        
        return seq
    
    def _normalize_chromosome(self, chromosome: str) -> str:
        """
        Normalize chromosome name for FASTA access.
        
        Args:
            chromosome: Input chromosome name
            
        Returns:
            Normalized chromosome name found in FASTA
            
        Raises:
            RuntimeError: If chromosome cannot be found in any format
        """
        # Remove 'chr' prefix if present for checking
        chrom_no_prefix = chromosome.replace('chr', '')
        
        # Get available chromosomes
        available_chroms = list(self.fa.keys())
        
        # Try both with and without 'chr' prefix
        for candidate in [chromosome, f'chr{chrom_no_prefix}', chrom_no_prefix]:
            if candidate in available_chroms:
                return candidate
        
        # If not found, raise error
        raise RuntimeError(f"Chromosome {chromosome} not found. Available: {available_chroms[:10]}...")


class HG38Reference:
    """Interface for hg38 reference genome access with fail-fast validation."""
    
    def __init__(self, reference_path: Optional[str] = None):
        """
        Initialize hg38 reference.
        
        Args:
            reference_path: Path to hg38 FASTA file. If None, uses HG38_FA env var.
            
        Raises:
            RuntimeError: If hg38 reference cannot be loaded
        """
        self.reference_path = reference_path
        self.logger = logging.getLogger(__name__)
        
        if not self.reference_path:
            self.reference_path = os.environ.get("HG38_FA")
        
        if not self.reference_path:
            raise RuntimeError("HG38_FA environment variable not set. Provide path to hg38 FASTA.")
        
        # Initialize extractor (will raise if invalid)
        self.extractor = Hg38Extractor(self.reference_path)
    
    def is_available(self) -> bool:
        """Check if hg38 reference is available."""
        return self.extractor is not None
    
    def get_window(self, chromosome: str, position: int, window_size: int = 201) -> str:
        """
        Extract genomic window around a position.
        
        Args:
            chromosome: Chromosome name (e.g., 'chr1', '1')
            position: 1-based genomic position
            window_size: Size of window to extract
            
        Returns:
            Genomic sequence window
            
        Raises:
            RuntimeError: If window extraction fails
        """
        try:
            return self.extractor.window_201(chromosome, position, window_size)
        except Exception as e:
            raise RuntimeError(f"Failed to extract window at {chromosome}:{position}: {e}")
    
    def get_reference_info(self) -> Dict[str, Any]:
        """Get information about the loaded reference."""
        if not self.extractor:
            return {'available': False, 'mode': 'failed'}
        
        try:
            chroms = list(self.extractor.fa.keys())
            return {
                'available': True,
                'path': self.reference_path,
                'mode': 'hg38',
                'chromosomes': chroms[:10],  # First 10 chromosomes
                'total_chromosomes': len(chroms)
            }
        except Exception as e:
            return {
                'available': False,
                'path': self.reference_path,
                'mode': 'error',
                'error': str(e)
            }


# Global instance for easy access
_hg38_reference = None

def get_hg38_reference(reference_path: Optional[str] = None) -> HG38Reference:
    """
    Get global hg38 reference instance.
    
    Args:
        reference_path: Path to hg38 FASTA file
        
    Returns:
        HG38Reference instance
        
    Raises:
        RuntimeError: If hg38 reference cannot be loaded
    """
    global _hg38_reference
    
    if _hg38_reference is None:
        _hg38_reference = HG38Reference(reference_path)
    
    return _hg38_reference


def extract_genomic_context(guide_data: Dict[str, Any], window_size: int = 201) -> str:
    """
    Extract genomic context for a guide sequence.
    
    Args:
        guide_data: Dictionary containing guide information
        window_size: Size of genomic window to extract
        
    Returns:
        Genomic context sequence
        
    Raises:
        RuntimeError: If hg38 reference is not available or extraction fails
    """
    hg38 = get_hg38_reference()
    
    # Try to extract from genomic coordinates if available
    if 'chromosome' in guide_data and 'position' in guide_data:
        return hg38.get_window(
            guide_data['chromosome'], 
            guide_data['position'], 
            window_size
        )
    
    # If no coordinates available, this is a critical failure
    # No fallback modes allowed per review requirements
    raise RuntimeError(
        f"No genomic coordinates available for guide sequence. "
        f"Required: 'chromosome' and 'position' in guide_data. "
        f"Available keys: {list(guide_data.keys())}"
    )