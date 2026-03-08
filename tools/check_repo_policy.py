#!/usr/bin/env python3
"""Repository policy checker for CRISPR-focused layout and hygiene."""

from __future__ import annotations

import fnmatch
import sys
from pathlib import Path


ALLOWED_TOP_LEVEL_DIRS = {
    ".git",
    ".github",
    ".idea",
    ".pytest_cache",
    ".venv314",
    "applications",
    "bin",
    "configs",
    "data",
    "docs",
    "experiments",
    "legacy",
    "modules",
    "notebooks",
    "proof_pack",
    "scripts",
    "static",
    "templates",
    "test_data",
    "tests",
    "tools",
    "wave_crispr_signal",
    "web_apps",
}

ALLOWED_TOP_LEVEL_FILES = {
    ".gitignore",
    ".pre-commit-config.yaml",
    "LICENSE",
    "Makefile",
    "README.md",
    "pyproject.toml",
    "repo_metadata.yaml",
    "requirements.txt",
    "targets.yaml",
}

FORBIDDEN_ROOT_PATTERNS = [
    "acf_scale_*.png",
    "null_dist_output.txt",
    "falsification_results.json",
    "FALSIFICATION_REPORT.md",
]

NON_CRISPR_KEYWORDS = [
    "ultrasound",
    "mri",
    "pain",
    "orcs",
    "z5d",
]

CORE_PATHS = [
    "wave_crispr_signal",
    "applications",
]


class PolicyError(Exception):
    pass


def check_root_layout(repo_root: Path) -> list[str]:
    issues: list[str] = []
    for p in repo_root.iterdir():
        name = p.name
        if p.is_dir():
            if name.startswith(".") and name in {".git", ".github", ".idea", ".pytest_cache", ".venv314"}:
                continue
            if name.startswith(".") and name not in ALLOWED_TOP_LEVEL_DIRS:
                issues.append(f"Unexpected top-level hidden directory: {name}")
            elif not name.startswith(".") and name not in ALLOWED_TOP_LEVEL_DIRS:
                issues.append(f"Unexpected top-level directory: {name}")
        else:
            if name.startswith("."):
                if name not in {".gitignore", ".pre-commit-config.yaml"}:
                    issues.append(f"Unexpected top-level hidden file: {name}")
            elif name not in ALLOWED_TOP_LEVEL_FILES:
                issues.append(f"Unexpected top-level file: {name}")
    return issues


def check_root_artifacts(repo_root: Path) -> list[str]:
    issues: list[str] = []
    for pattern in FORBIDDEN_ROOT_PATTERNS:
        for match in repo_root.glob(pattern):
            if match.is_file():
                issues.append(f"Forbidden generated artifact at root: {match.name}")
    return issues


def check_core_scope(repo_root: Path) -> list[str]:
    issues: list[str] = []
    for core in CORE_PATHS:
        root = repo_root / core
        if not root.exists():
            issues.append(f"Missing core path: {core}")
            continue
        for path in root.rglob("*.py"):
            rel = str(path.relative_to(repo_root)).lower()
            for keyword in NON_CRISPR_KEYWORDS:
                if keyword in rel:
                    issues.append(
                        f"Non-CRISPR keyword '{keyword}' found in core path: {path.relative_to(repo_root)}"
                    )
                    break
    return issues


def check_docs_index(repo_root: Path) -> list[str]:
    issues: list[str] = []
    docs_index = repo_root / "docs" / "INDEX.md"
    inventory = repo_root / "docs" / "INVENTORY_CLASSIFICATION.md"
    legacy_readme = repo_root / "legacy" / "README.md"
    readme = repo_root / "README.md"

    for required in [docs_index, inventory, legacy_readme, readme]:
        if not required.exists():
            issues.append(f"Missing required documentation file: {required.relative_to(repo_root)}")

    if docs_index.exists():
        text = docs_index.read_text(encoding="utf-8")
        if "INVENTORY_CLASSIFICATION.md" not in text:
            issues.append("docs/INDEX.md missing link/reference to INVENTORY_CLASSIFICATION.md")
        if "../legacy/README.md" not in text:
            issues.append("docs/INDEX.md missing link/reference to legacy/README.md")

    if readme.exists():
        text = readme.read_text(encoding="utf-8")
        if "docs/INDEX.md" not in text:
            issues.append("README.md missing link/reference to docs/INDEX.md")
        if "legacy/" not in text:
            issues.append("README.md missing legacy scope reference")

    return issues


def run(repo_root: Path) -> int:
    issues: list[str] = []
    issues.extend(check_root_layout(repo_root))
    issues.extend(check_root_artifacts(repo_root))
    issues.extend(check_core_scope(repo_root))
    issues.extend(check_docs_index(repo_root))

    if issues:
        print("Repository policy violations detected:")
        for item in issues:
            print(f"- {item}")
        return 1

    print("Repository policy checks passed.")
    return 0


def main() -> None:
    repo_root = Path(".").resolve()
    raise SystemExit(run(repo_root))


if __name__ == "__main__":
    main()
