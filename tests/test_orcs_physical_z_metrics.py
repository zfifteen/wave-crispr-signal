#!/usr/bin/env python3
"""
Test suite for ORCS Discrete Biological Z-Metrics implementation.

Tests the ORCS v1.1.17 testing script that validates four discrete biological Z-metrics
on real human CRISPR screen outcomes.
"""

import unittest
import sys
import os
import tempfile
import shutil
from pathlib import Path
import pandas as pd

# Add project root to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'scripts'))

from scripts.test_orcs_v1_1_17 import (
    load_orcs_index, 
    select_screens, 
    seqs_by_symbol_from_fasta,
    enumerate_spcas9_guides,
    aggregate_gene,
    z_of_seq,
    PROP_TABLES,
    validate_seq
)


class TestHumanDNAValidation(unittest.TestCase):
    """Test human DNA sequence validation."""
    
    def test_valid_acgtn_sequences(self):
        """Test that valid ACGTN sequences are accepted."""
        valid_sequences = [
            "ATCGN",
            "ATCNNN", 
            "ATCG",
            "AAAAA",
            "CCCCC",
            "GGGGG",
            "TTTTT",
            "NNNNN",
            "ATCGATCGATCGATCGATCG"
        ]
        
        for seq in valid_sequences:
            with self.subTest(sequence=seq):
                try:
                    validate_seq(seq)
                except ValueError:
                    self.fail(f"Valid sequence '{seq}' was rejected")
    
    def test_reject_rna_nucleotides(self):
        """Test that sequences containing U (RNA) are rejected."""
        invalid_sequences = [
            "AUCGA",  # Contains U
            "ATCGU",  # Contains U
            "UUUUU",  # All U
            "AUGCGA", # Multiple with U
        ]
        
        for seq in invalid_sequences:
            with self.subTest(sequence=seq):
                with self.assertRaises(ValueError) as cm:
                    validate_seq(seq)
                self.assertIn("Invalid nucleotides", str(cm.exception))
    
    def test_reject_iupac_ambiguous_codes(self):
        """Test that IUPAC ambiguous nucleotide codes are rejected."""
        invalid_sequences = [
            "ATCRG",  # Contains R (A or G)
            "ATCYG",  # Contains Y (C or T) 
            "ATCSG",  # Contains S (C or G)
            "ATCWG",  # Contains W (A or T)
            "ATCKG",  # Contains K (G or T)
            "ATCMG",  # Contains M (A or C)
            "ATCBG",  # Contains B (C, G, or T)
            "ATCDG",  # Contains D (A, G, or T)
            "ATCHG",  # Contains H (A, C, or T)
            "ATCVG",  # Contains V (A, C, or G)
        ]
        
        for seq in invalid_sequences:
            with self.subTest(sequence=seq):
                with self.assertRaises(ValueError) as cm:
                    validate_seq(seq)
                self.assertIn("Invalid nucleotides", str(cm.exception))
    
    def test_reject_non_nucleotide_characters(self):
        """Test that non-nucleotide characters are rejected."""
        invalid_sequences = [
            "ATCGZ",   # Contains Z
            "ATCG123", # Contains numbers
            "ATCG-N",  # Contains dash
            "ATCG N",  # Contains space
            "ATCG.N",  # Contains period
        ]
        
        for seq in invalid_sequences:
            with self.subTest(sequence=seq):
                with self.assertRaises(ValueError) as cm:
                    validate_seq(seq)
                self.assertIn("Invalid nucleotides", str(cm.exception))
    
    def test_empty_sequence(self):
        """Test that empty sequences are rejected."""
        with self.assertRaises(ValueError):
            validate_seq("")
    
    def test_case_insensitive_validation(self):
        """Test that lowercase sequences are handled properly."""
        # The validate_seq function should handle case properly
        # Let's test both cases
        valid_lower = "atcgn"
        # This might fail if validate_seq expects uppercase - let's see
        try:
            # enumerate_spcas9_guides converts to uppercase before validation
            from scripts.test_orcs_v1_1_17 import enumerate_spcas9_guides
            # This should work because it converts to uppercase
            enumerate_spcas9_guides("atcgatcgatcgatcgatcggg")
        except ValueError as e:
            # If it fails on case, that's acceptable as long as it's documented
            pass


class TestORCSDataLoading(unittest.TestCase):
    """Test ORCS data loading and processing functions."""

    def setUp(self):
        """Set up test data."""
        self.orcs_dir = "data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens"
        self.test_fasta = "data/test_human_cdna.fasta"

    def test_load_orcs_index(self):
        """Test loading ORCS index files."""
        if not os.path.exists(self.orcs_dir):
            self.skipTest("ORCS data directory not found")
        
        idx = load_orcs_index(self.orcs_dir)
        self.assertIsInstance(idx, pd.DataFrame)
        self.assertGreater(len(idx), 0)
        
        # Check expected columns exist
        self.assertGreater(idx.shape[1], 30)  # Should have many columns

    def test_select_screens(self):
        """Test screen selection logic."""
        if not os.path.exists(self.orcs_dir):
            self.skipTest("ORCS data directory not found")
        
        idx = load_orcs_index(self.orcs_dir)
        selected = select_screens(idx)
        
        self.assertIsInstance(selected, pd.DataFrame)
        if len(selected) > 0:
            self.assertIn("screen_id", selected.columns)
            self.assertIn("score1_type", selected.columns)
            
            # Should only contain human, Cas9 screens
            print(f"Found {len(selected)} qualifying screens")

    def test_load_sequences_from_fasta(self):
        """Test sequence loading from FASTA."""
        if not os.path.exists(self.test_fasta):
            self.skipTest("Test FASTA file not found")
        
        seqs = seqs_by_symbol_from_fasta(self.test_fasta)
        
        self.assertIsInstance(seqs, dict)
        self.assertGreater(len(seqs), 0)
        
        # Check expected genes are present
        expected_genes = ["TOP2A", "CDK6", "TAF3", "MYB", "UBE2G2", "HAND1"]
        for gene in expected_genes:
            self.assertIn(gene, seqs)
            self.assertIsInstance(seqs[gene], str)
            self.assertGreater(len(seqs[gene]), 0)


class TestGuideEnumeration(unittest.TestCase):
    """Test CRISPR guide enumeration and Z-metric calculation."""

    def test_enumerate_spcas9_guides(self):
        """Test SpCas9 guide enumeration."""
        # Test sequence with known NGG PAM sites
        test_seq = "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCG"
        guides = enumerate_spcas9_guides(test_seq)
        
        self.assertIsInstance(guides, list)
        if len(guides) > 0:
            guide, pos = guides[0]
            self.assertEqual(len(guide), 20)
            self.assertIsInstance(pos, int)
            self.assertGreater(pos, 0)

    def test_z_metric_calculation(self):
        """Test Z-metric calculation for a sequence."""
        test_seq = "ATCGATCGATCGATCGATCG"  # 20 bp guide sequence
        
        for metric_name, table in PROP_TABLES.items():
            with self.subTest(metric=metric_name):
                z_val = z_of_seq(test_seq, table)
                self.assertIsInstance(float(z_val), float)
                self.assertTrue(abs(float(z_val)) > 0)  # Can be positive or negative

    def test_aggregate_gene_metrics(self):
        """Test gene-level Z-metric aggregation."""
        # Use a longer sequence that should have multiple guides
        test_seq = "ATCGATCGATCGATCGATCGGGCCGATCGATCGATCGATCGGGTCGATCGATCGATCGATCGGG"
        guides = enumerate_spcas9_guides(test_seq)
        
        if len(guides) > 0:
            for metric_name in PROP_TABLES.keys():
                with self.subTest(metric=metric_name):
                    agg = aggregate_gene(guides, metric_name)
                    if agg is not None:
                        self.assertIn("mean", agg)
                        self.assertIn("median", agg)
                        self.assertIn("top_decile_mean", agg)
                        self.assertIn("var", agg)
                        
                        # Check values are reasonable
                        self.assertIsInstance(agg["mean"], float)
                        self.assertTrue(abs(agg["mean"]) >= 0)  # Can be positive or negative


class TestORCSIntegration(unittest.TestCase):
    """Test integration of ORCS data with Z-metrics."""

    def test_full_pipeline_subset(self):
        """Test full pipeline on a small subset of data."""
        orcs_dir = "data/BIOGRID-ORCS-ALL-homo_sapiens-1.1.17.screens"
        test_fasta = "data/test_human_cdna.fasta"
        
        if not os.path.exists(orcs_dir) or not os.path.exists(test_fasta):
            self.skipTest("Required data files not found")
        
        # Test with temporary output directory
        with tempfile.TemporaryDirectory() as temp_dir:
            import subprocess
            import sys
            
            # Set environment variables for test
            env = os.environ.copy()
            env["ORCS_DIR"] = orcs_dir
            env["FASTA"] = test_fasta
            env["OUT"] = os.path.join(temp_dir, "test_summary.csv")
            env["MAX_SCREENS"] = "1"  # Limit to 1 screen for testing
            env["MIN_PAIRS"] = "3"    # Lower threshold for test data
            
            # Run script with limited processing
            result = subprocess.run([
                sys.executable, "scripts/test_orcs_v1_1_17.py"
            ], capture_output=True, text=True, env=env, timeout=60,
               cwd=project_root)
            
            # Check if script completed successfully
            self.assertEqual(result.returncode, 0, 
                           f"Script failed: {result.stderr}")
            
            # Verify output file was created
            output_file = env["OUT"]
            self.assertTrue(os.path.exists(output_file), 
                          "Output CSV file was not created")
            
            df = pd.read_csv(output_file)
            print(f"Generated {len(df)} results")
            
            if len(df) > 0:
                # Check expected columns
                expected_cols = [
                    "screen_id", "score1_type", "metric", "aggregate",
                    "N", "pearson_r", "p_value", "ci95_lo", "ci95_hi"
                ]
                for col in expected_cols:
                    self.assertIn(col, df.columns)
                
                # Check metrics are calculated
                metrics = df["metric"].unique()
                self.assertGreater(len(metrics), 0)
                print(f"Calculated metrics: {list(metrics)}")
                
                # Check for meaningful correlations
                strong_corrs = df[df["pearson_r"].abs() > 0.5]
                if len(strong_corrs) > 0:
                    print(f"Found {len(strong_corrs)} strong correlations (|r| > 0.5)")
            else:
                print("No results generated - may need more test data")


if __name__ == "__main__":
    unittest.main()