import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path("/Users/velocityworks/IdeaProjects/wave-crispr-signal/validation/ontarget_gate/scripts/run_gate_v2.py")
SPEC = importlib.util.spec_from_file_location("run_gate_v2", MODULE_PATH)
gate = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = gate
SPEC.loader.exec_module(gate)


class GateV2Tests(unittest.TestCase):
    def test_assign_primary_split_is_deterministic(self):
        g = "TP53"
        s1 = gate.assign_primary_split(g)
        s2 = gate.assign_primary_split(g)
        self.assertEqual(s1, s2)
        self.assertIn(s1, {"dev_train", "dev_val", "primary_holdout"})

    def test_validate_split_integrity_detects_overlap(self):
        rows = [
            {"source": "s1", "split": "dev_train", "group_strength": "strong", "gene_or_target_group": "A"},
            {"source": "s1", "split": "dev_val", "group_strength": "strong", "gene_or_target_group": "A"},
            {"source": "s1", "split": "primary_holdout", "group_strength": "strong", "gene_or_target_group": "B"},
            {"source": "s2", "split": "external_holdout", "group_strength": "strong", "gene_or_target_group": "C"},
        ]
        out = gate.validate_split_integrity(rows)
        self.assertFalse(out["ok"])
        self.assertTrue(any(v.startswith("group_overlap:") for v in out["violations"]))

    def test_validate_split_integrity_passes_clean_case(self):
        rows = [
            {"source": "s1", "split": "dev_train", "group_strength": "strong", "gene_or_target_group": "A"},
            {"source": "s1", "split": "dev_val", "group_strength": "strong", "gene_or_target_group": "B"},
            {"source": "s1", "split": "primary_holdout", "group_strength": "strong", "gene_or_target_group": "C"},
            {"source": "s2", "split": "external_holdout", "group_strength": "strong", "gene_or_target_group": "D"},
        ]
        out = gate.validate_split_integrity(rows)
        self.assertTrue(out["ok"])

    def test_bootstrap_weighted_delta_returns_ci(self):
        rng = np.random.default_rng(7)
        y1 = rng.normal(size=100)
        y2 = rng.normal(size=100)
        m1 = y1 + rng.normal(scale=0.1, size=100)
        b1 = -y1 + rng.normal(scale=0.1, size=100)
        m2 = y2 + rng.normal(scale=0.1, size=100)
        b2 = -y2 + rng.normal(scale=0.1, size=100)
        lo, hi = gate.bootstrap_weighted_delta_spearman(
            [
                {"y": y1, "model": m1, "baseline_a": b1},
                {"y": y2, "model": m2, "baseline_a": b2},
            ],
            n_boot=200,
            seed=11,
        )
        self.assertLessEqual(lo, hi)
        self.assertGreater(hi, 0.0)

    def test_decide_outcome_go_path(self):
        dev = {
            "delta_model_minus_baseline_a_spearman": 0.05,
            "delta_model_minus_baseline_a_spearman_ci95": [0.01, 0.09],
            "delta_model_minus_baseline_b_spearman": 0.02,
        }
        primary = {
            "delta_model_minus_baseline_a_spearman": 0.04,
            "delta_model_minus_baseline_a_spearman_ci95": [0.01, 0.07],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        external = {
            "delta_model_minus_baseline_a_spearman": 0.02,
            "delta_model_minus_baseline_a_spearman_ci95": [0.001, 0.04],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        out = gate.decide_outcome(dev, primary, external, agg_delta=0.03, agg_ci=(0.01, 0.05))
        self.assertEqual(out["decision"], "GO")

    def test_decide_outcome_hard_fail_dev_is_no_go(self):
        dev = {
            "delta_model_minus_baseline_a_spearman": -0.2,
            "delta_model_minus_baseline_a_spearman_ci95": [-0.25, -0.02],
            "delta_model_minus_baseline_b_spearman": -0.1,
        }
        primary = {
            "delta_model_minus_baseline_a_spearman": 0.05,
            "delta_model_minus_baseline_a_spearman_ci95": [0.01, 0.08],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        external = {
            "delta_model_minus_baseline_a_spearman": 0.02,
            "delta_model_minus_baseline_a_spearman_ci95": [0.005, 0.03],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        out = gate.decide_outcome(dev, primary, external, agg_delta=0.03, agg_ci=(0.01, 0.05))
        self.assertEqual(out["decision"], "NO-GO")

    def test_decide_outcome_inconclusive_path(self):
        dev = {
            "delta_model_minus_baseline_a_spearman": 0.01,
            "delta_model_minus_baseline_a_spearman_ci95": [-0.005, 0.03],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        primary = {
            "delta_model_minus_baseline_a_spearman": 0.031,
            "delta_model_minus_baseline_a_spearman_ci95": [0.001, 0.05],
            "delta_model_minus_baseline_b_spearman": 0.01,
        }
        external = {
            "delta_model_minus_baseline_a_spearman": 0.005,
            "delta_model_minus_baseline_a_spearman_ci95": [-0.001, 0.02],
            "delta_model_minus_baseline_b_spearman": 0.005,
        }
        out = gate.decide_outcome(dev, primary, external, agg_delta=0.01, agg_ci=(-0.001, 0.03))
        self.assertEqual(out["decision"], "INCONCLUSIVE")


if __name__ == "__main__":
    unittest.main()
