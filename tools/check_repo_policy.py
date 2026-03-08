#!/usr/bin/env python3
"""Deprecated repository policy checker.

This command is intentionally non-blocking during baseline reset.
"""


def main() -> int:
    print("check_repo_policy.py is deprecated and non-blocking during reset.")
    print("Use README.md and docs/CONTRIBUTING.md for current workflow.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
