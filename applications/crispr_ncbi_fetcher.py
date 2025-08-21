"""
CRISPR NCBI Sequence Fetcher for WAVE-CRISPR Experiment

This module provides BioPython-based NCBI sequence fetching functionality
for reproducible sequence analysis in the WAVE-CRISPR framework.
"""

import time
import logging
from typing import List, Dict, Optional
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# NCBI API rate limiting
NCBI_DELAY = 0.34  # seconds between requests (respect NCBI guidelines)


class CRISPRNCBIFetcher:
    """
    BioPython-based NCBI sequence fetcher for WAVE-CRISPR experiments.
    
    Implements reproducible sequence fetching with error handling and
    rate limiting for empirical validation studies.
    """
    
    def __init__(self, email: str = "research@example.com", tool: str = "wave-crispr-signal"):
        """
        Initialize NCBI fetcher.
        
        Args:
            email: Email for NCBI API (required by NCBI)
            tool: Tool name for NCBI API tracking
        """
        Entrez.email = email
        Entrez.tool = tool
        self.fetch_count = 0
        
    def fetch_sequence(self, accession_id: str, rettype: str = "fasta") -> Optional[SeqRecord]:
        """
        Fetch single sequence from NCBI.
        
        Args:
            accession_id: NCBI accession ID (e.g., "NM_007294.4")
            rettype: Return type (fasta, gb, etc.)
            
        Returns:
            SeqRecord object or None if fetch fails
        """
        try:
            # Rate limiting
            if self.fetch_count > 0:
                time.sleep(NCBI_DELAY)
            
            logger.info(f"Fetching sequence: {accession_id}")
            
            with Entrez.efetch(
                db="nucleotide", 
                id=accession_id, 
                rettype=rettype, 
                retmode="text"
            ) as handle:
                record = SeqIO.read(handle, rettype)
                
            self.fetch_count += 1
            logger.info(f"Successfully fetched {accession_id}, length: {len(record.seq)}")
            return record
            
        except Exception as e:
            logger.error(f"Failed to fetch {accession_id}: {e}")
            return None
    
    def batch_fetch_sequences(self, accession_ids: List[str], 
                            max_retries: int = 3) -> Dict[str, Optional[SeqRecord]]:
        """
        Batch fetch multiple sequences with retry logic.
        
        Args:
            accession_ids: List of NCBI accession IDs
            max_retries: Maximum retry attempts per sequence
            
        Returns:
            Dictionary mapping accession_id to SeqRecord (or None if failed)
        """
        results = {}
        
        for accession_id in accession_ids:
            retry_count = 0
            record = None
            
            while retry_count < max_retries and record is None:
                record = self.fetch_sequence(accession_id)
                if record is None:
                    retry_count += 1
                    if retry_count < max_retries:
                        wait_time = NCBI_DELAY * (2 ** retry_count)  # Exponential backoff
                        logger.warning(f"Retrying {accession_id} in {wait_time}s (attempt {retry_count + 1})")
                        time.sleep(wait_time)
            
            results[accession_id] = record
            
        return results
    
    def fetch_test_sequences(self) -> Dict[str, Optional[SeqRecord]]:
        """
        Fetch standard test sequences for WAVE-CRISPR validation.
        
        Returns:
            Dictionary with test sequences
        """
        test_ids = [
            "NM_007294.4",  # BRCA1
            "NM_000546.6",  # TP53
        ]
        
        logger.info("Fetching standard test sequences for WAVE-CRISPR validation")
        return self.batch_fetch_sequences(test_ids)
    
    def validate_sequence_integrity(self, record: SeqRecord) -> Dict[str, any]:
        """
        Validate sequence integrity for WAVE-CRISPR analysis.
        
        Args:
            record: SeqRecord to validate
            
        Returns:
            Validation metrics dictionary
        """
        if record is None:
            return {"valid": False, "error": "No sequence record"}
        
        seq_str = str(record.seq).upper()
        
        # Basic validation metrics
        validation = {
            "valid": True,
            "length": len(seq_str),
            "gc_content": (seq_str.count('G') + seq_str.count('C')) / len(seq_str) if len(seq_str) > 0 else 0,
            "n_content": seq_str.count('N') / len(seq_str) if len(seq_str) > 0 else 0,
            "valid_bases": all(base in 'ATCGN' for base in seq_str),
            "min_length_ok": len(seq_str) >= 20,  # Minimum for CRISPR analysis
        }
        
        # Overall validity check
        validation["valid"] = (
            validation["valid_bases"] and 
            validation["min_length_ok"] and 
            validation["n_content"] < 0.1  # Less than 10% N bases
        )
        
        return validation


def demo_ncbi_fetcher():
    """Demonstration of NCBI fetcher functionality."""
    print("WAVE-CRISPR NCBI Fetcher Demo")
    print("=" * 40)
    
    # Initialize fetcher
    fetcher = CRISPRNCBIFetcher()
    
    # Fetch test sequences
    sequences = fetcher.fetch_test_sequences()
    
    # Validate and display results
    for accession_id, record in sequences.items():
        validation = fetcher.validate_sequence_integrity(record)
        
        print(f"\nSequence: {accession_id}")
        print(f"Valid: {validation['valid']}")
        if validation['valid']:
            print(f"Length: {validation['length']} bp")
            print(f"GC Content: {validation['gc_content']:.3f}")
            print(f"First 50 bases: {str(record.seq)[:50]}")


if __name__ == "__main__":
    demo_ncbi_fetcher()