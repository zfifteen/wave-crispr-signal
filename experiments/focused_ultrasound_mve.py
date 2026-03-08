"""Compatibility wrapper for relocated legacy experiment module."""
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / "scripts") not in sys.path:
    sys.path.insert(0, str(repo_root / "scripts"))

from legacy.experiments.focused_ultrasound_mve import *  # noqa: F401,F403
from legacy.experiments.focused_ultrasound_mve import main as _legacy_main

if __name__ == "__main__":
    _legacy_main()
