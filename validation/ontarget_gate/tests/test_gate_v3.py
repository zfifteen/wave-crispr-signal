import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path("/Users/velocityworks/IdeaProjects/wave-crispr-signal")
MODULE_PATH = ROOT / "validation/ontarget_gate/scripts/run_gate_v3.py"
SPEC = importlib.util.spec_from_file_location("run_gate_v3", MODULE_PATH)
gate = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = gate
SPEC.loader.exec_module(gate)


class GateV3Tests(unittest.TestCase):
    def test_decide_inconclusive_when_precondition_fails(self):
        base = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.1],
            "delta_model_minus_baseline_c_spearman": 0.02,
        }
        out = gate.decide_outcome_v3(
            base,
            base,
            base,
            {
                "comparator_self_check_ok": False,
                "overlap_audit_ok": True,
                "holdout_min_n_ok": True,
            },
        )
        self.assertEqual(out["decision"], "INCONCLUSIVE")

    def test_decide_no_go_when_threshold_fails(self):
        dev = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.02,
        }
        primary = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.005,
        }
        external = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.02,
        }
        out = gate.decide_outcome_v3(
            dev,
            primary,
            external,
            {
                "comparator_self_check_ok": True,
                "overlap_audit_ok": True,
                "holdout_min_n_ok": True,
            },
        )
        self.assertEqual(out["decision"], "NO-GO")

    def test_decide_go_when_thresholds_pass(self):
        dev = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.02,
        }
        primary = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.011,
        }
        external = {
            "delta_model_minus_baseline_a_spearman_ci95": [-0.01, 0.2],
            "delta_model_minus_baseline_c_spearman": 0.012,
        }
        out = gate.decide_outcome_v3(
            dev,
            primary,
            external,
            {
                "comparator_self_check_ok": True,
                "overlap_audit_ok": True,
                "holdout_min_n_ok": True,
            },
        )
        self.assertEqual(out["decision"], "GO")

    def test_overlap_audit_training_manifest_unavailable(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            split = td / "split.csv"
            split.write_text("source,guide,guide_seq,label,gene_or_target_group,group_strength,group_method,split\n"
                             "x,g1,AAAAAAAAAAAAAAAAAAAA,0.1,G,strong,m,primary_holdout\n", encoding="utf-8")
            status = td / "status.json"
            status.write_text(json.dumps({"available": False}), encoding="utf-8")
            out = td / "out.json"

            cmd = [
                sys.executable,
                str(ROOT / "validation/ontarget_gate/scripts/audit_overlap.py"),
                "--split-manifest",
                str(split),
                "--training-manifest",
                str(td / "training_manifest.fasta"),
                "--training-status",
                str(status),
                "--output",
                str(out),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            obj = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(obj["message"], "training_manifest_unavailable")

    def test_comparator_self_check_fails_without_venv(self):
        from validation.ontarget_gate.comparators.crisprpred.adapter import CRISPRpredAdapter

        adapter = CRISPRpredAdapter()
        adapter.venv_python = Path("/nonexistent/comparator-python")
        result = adapter.self_check()
        self.assertFalse(result.ok)
        self.assertEqual(result.message, "comparator_venv_missing")


if __name__ == "__main__":
    unittest.main()
