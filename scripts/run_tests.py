#!/usr/bin/env python3
"""Change-scoped validation runner.

This helper does not enforce a baseline suite. It runs commands provided by the caller.
"""

from __future__ import annotations

import subprocess
import sys


def run_command(cmd: str) -> int:
    print(f"\n==> {cmd}")
    completed = subprocess.run(cmd, shell=True)
    if completed.returncode != 0:
        print(f"FAILED ({completed.returncode}): {cmd}")
    return completed.returncode


def main(argv: list[str]) -> int:
    if len(argv) == 1:
        print("No baseline test suite is configured.")
        print("Provide one or more change-scoped validation commands to run.")
        print("Example:")
        print(
            "  python scripts/run_tests.py \"python -m pytest -q path/to/test.py\""
        )
        return 0

    failures = 0
    for cmd in argv[1:]:
        failures += 1 if run_command(cmd) != 0 else 0

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
