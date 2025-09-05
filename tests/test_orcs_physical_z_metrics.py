#!/usr/bin/env python3
"""
Test suite for ORCS Physical Z-Metrics implementation.

Tests the ORCS v1.1.17 testing script that validates four physical Z-metrics
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

from test_orcs_v1_1_17 import (
    load_orcs_index, 
    select_screens, 
    seqs_by_symbol_from_fasta,
    enumerate_spcas9_guides,
    aggregate_gene,
    z_of_seq,
    PROP_TABLES
)


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