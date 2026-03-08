"""Change-scoped CLI contract tests for preserved primary interfaces."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=120,
    )


def test_phase_weighted_scorecard_cli_score_contract() -> None:
    cmd = [
        sys.executable,
        "applications/phase_weighted_scorecard_cli.py",
        "score",
        "--guide",
        "GCTGCGGAGACCTGGAGAGA",
    ]
    result = _run(cmd)
    assert result.returncode == 0, result.stderr
    assert "Phase-Weighted CRISPR Scorecard" in result.stdout


def test_crispr_cli_help_contract() -> None:
    cmd = [sys.executable, "applications/crispr_cli.py", "--help"]
    result = _run(cmd)
    assert result.returncode == 0, result.stderr
    assert "CRISPR" in result.stdout


def test_genomic_disruption_api_help_contract() -> None:
    cmd = [sys.executable, "applications/genomic_disruption_api.py", "--help"]
    result = _run(cmd)
    assert result.returncode == 0, result.stderr
    assert "Genomic Disruption Analyzer" in result.stdout
