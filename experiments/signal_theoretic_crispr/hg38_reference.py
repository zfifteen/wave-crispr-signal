#!/usr/bin/env python3
"""
hg38 Reference Genome Interface for CRISPR experiment.

This module provides functionality to extract genomic windows from the hg38 reference
for contextual analysis of CRISPR guide sequences.
"""

import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any
import warnings

try:
    from pyfaidx import Fasta
    PYFAIDX_AVAILABLE = True
except ImportError:
    PYFAIDX_AVAILABLE = False
    warnings.warn("pyfaidx not available - hg38 anchoring will use fallback mode")


class HG38Reference:
    """Interface for hg38 reference genome access."""
    
    def __init__(self, reference_path: Optional[str] = None):
        """
        Initialize hg38 reference.
        
        Args:
            reference_path: Path to hg38 FASTA file. If None, will attempt to locate.
        """
        self.reference_path = reference_path
        self.logger = logging.getLogger(__name__)
        self.fasta = None
        
        # Try to initialize reference
        self._initialize_reference()
    
    def _initialize_reference(self):
        """Initialize the reference genome access."""
        if not PYFAIDX_AVAILABLE:
            self.logger.warning("pyfaidx not available - using fallback mode")
            return
        
        # Try to find hg38 reference file
        candidate_paths = []
        
        if self.reference_path:
            candidate_paths.append(self.reference_path)
        
        # Check common locations
        data_dir = Path("data/get_hg38")
        candidate_paths.extend([
            data_dir / "hg38.fa" / "hg38.fa",
            data_dir / "hg38.fa",
            "/home/runner/work/wave-crispr-signal/wave-crispr-signal/data/get_hg38/hg38.fa/hg38.fa",
            "hg38.fa"
        ])
        
        for path in candidate_paths:
            if os.path.exists(path):
                try:
                    self.fasta = Fasta(str(path))
                    self.reference_path = str(path)
                    self.logger.info(f"Loaded hg38 reference: {path}")
                    return
                except Exception as e:
                    self.logger.warning(f"Failed to load {path}: {e}")
                    continue
        
        self.logger.warning("No hg38 reference found - using fallback mode")
    
    def is_available(self) -> bool:
        """Check if hg38 reference is available."""
        return self.fasta is not None
    
    def get_window(self, chromosome: str, position: int, window_size: int = 201,
                   guide_sequence: Optional[str] = None) -> str:
        """
        Extract genomic window around a position.
        
        Args:
            chromosome: Chromosome name (e.g., 'chr1', '1')
            position: 1-based genomic position
            window_size: Size of window to extract
            guide_sequence: Guide sequence for fallback mode
            
        Returns:
            Genomic sequence window
        """
        if not self.fasta:
            return self._fallback_window(guide_sequence, window_size)
        
        try:
            # Normalize chromosome name
            chrom = self._normalize_chromosome(chromosome)
            
            # Calculate window boundaries
            half_window = window_size // 2
            start = max(1, position - half_window)
            end = start + window_size - 1
            
            # Extract sequence
            sequence = str(self.fasta[chrom][start-1:end]).upper()
            
            # Validate sequence
            if not sequence or len(sequence) < window_size // 2:
                self.logger.warning(f"Short sequence extracted at {chrom}:{position}")
                return self._fallback_window(guide_sequence, window_size)
            
            # Check for valid nucleotides only
            valid_chars = set('ACGTN')
            if not set(sequence) <= valid_chars:
                invalid_chars = set(sequence) - valid_chars
                self.logger.warning(f"Invalid characters in sequence: {invalid_chars}")
                # Replace invalid characters with N
                sequence = ''.join(c if c in valid_chars else 'N' for c in sequence)
            
            return sequence
            
        except Exception as e:
            self.logger.error(f"Failed to extract window at {chromosome}:{position}: {e}")
            return self._fallback_window(guide_sequence, window_size)
    
    def _normalize_chromosome(self, chromosome: str) -> str:
        """Normalize chromosome name for FASTA access."""
        # Remove 'chr' prefix if present for checking
        chrom_no_prefix = chromosome.replace('chr', '')
        
        # Try both with and without 'chr' prefix
        available_chroms = list(self.fasta.keys()) if self.fasta else []
        
        for chrom in [chromosome, f'chr{chrom_no_prefix}', chrom_no_prefix]:
            if chrom in available_chroms:
                return chrom
        
        # Default to provided name
        return chromosome
    
    def _fallback_window(self, guide_sequence: Optional[str], window_size: int) -> str:
        """
        Fallback method when hg38 reference is not available.
        
        This creates a synthetic context using the guide sequence as the center.
        """
        if not guide_sequence:
            # Create a synthetic sequence with balanced nucleotide content
            return self._create_synthetic_window(window_size)
        
        guide_len = len(guide_sequence)
        half_window = window_size // 2
        
        if guide_len >= window_size:
            # Guide is longer than window, take center portion
            start = (guide_len - window_size) // 2
            return guide_sequence[start:start + window_size]
        
        # Extend guide sequence symmetrically
        extend_needed = (window_size - guide_len) // 2
        
        # Create flanking sequences with balanced composition
        left_flank = self._create_balanced_sequence(extend_needed)
        right_flank = self._create_balanced_sequence(window_size - guide_len - extend_needed)
        
        return left_flank + guide_sequence + right_flank
    
    def _create_synthetic_window(self, size: int) -> str:
        """Create a synthetic genomic window with balanced nucleotide composition."""
        import random
        
        # Use a fixed seed for reproducibility
        random.seed(42)
        nucleotides = 'ACGT'
        
        # Create balanced composition (25% each nucleotide)
        sequence = []
        for i in range(size):
            sequence.append(nucleotides[i % 4])
        
        # Shuffle to avoid regular patterns
        random.shuffle(sequence)
        return ''.join(sequence)
    
    def _create_balanced_sequence(self, size: int) -> str:
        """Create a balanced nucleotide sequence of given size."""
        if size <= 0:
            return ""
        
        import random
        random.seed(42)
        
        nucleotides = 'ACGT'
        sequence = []
        
        for i in range(size):
            sequence.append(nucleotides[i % 4])
        
        random.shuffle(sequence)
        return ''.join(sequence)
    
    def get_reference_info(self) -> Dict[str, Any]:
        """Get information about the loaded reference."""
        info = {
            'available': self.is_available(),
            'path': self.reference_path,
            'mode': 'hg38' if self.fasta else 'fallback'
        }
        
        if self.fasta:
            info['chromosomes'] = list(self.fasta.keys())[:10]  # First 10 chromosomes
            info['total_chromosomes'] = len(self.fasta.keys())
        
        return info


# Global instance for easy access
_hg38_reference = None

def get_hg38_reference(reference_path: Optional[str] = None) -> HG38Reference:
    """Get global hg38 reference instance."""
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
    """
    hg38 = get_hg38_reference()
    
    # Try to extract from genomic coordinates if available
    if 'chromosome' in guide_data and 'position' in guide_data:
        return hg38.get_window(
            guide_data['chromosome'], 
            guide_data['position'], 
            window_size,
            guide_data.get('guide_sequence')
        )
    
    # Fallback to guide sequence
    guide_sequence = guide_data.get('guide_sequence', guide_data.get('sequence'))
    return hg38.get_window('chr1', 1000000, window_size, guide_sequence)