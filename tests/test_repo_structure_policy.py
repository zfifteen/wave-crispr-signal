"""Structure and artifact hygiene tests for repository layout policy."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

ALLOWED_DIRS = {
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

ALLOWED_FILES = {
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

FORBIDDEN_ROOT_FILES = {
    "FALSIFICATION_REPORT.md",
    "falsification_results.json",
    "null_dist_output.txt",
}

FORBIDDEN_ROOT_GLOBS = [
    "acf_scale_*.png",
]


def test_allowed_top_level_layout_only():
    top_level = [p for p in REPO_ROOT.iterdir() if p.name != ".DS_Store"]

    unexpected_dirs = []
    unexpected_files = []

    for p in top_level:
        if p.is_dir():
            if p.name.startswith(".") and p.name in {".git", ".github", ".idea", ".pytest_cache", ".venv314"}:
                continue
            if p.name not in ALLOWED_DIRS:
                unexpected_dirs.append(p.name)
        else:
            if p.name.startswith("."):
                if p.name not in {".gitignore", ".pre-commit-config.yaml"}:
                    unexpected_files.append(p.name)
            elif p.name not in ALLOWED_FILES:
                unexpected_files.append(p.name)

    assert not unexpected_dirs, f"Unexpected top-level directories: {unexpected_dirs}"
    assert not unexpected_files, f"Unexpected top-level files: {unexpected_files}"


def test_no_generated_artifacts_in_repo_root():
    for name in FORBIDDEN_ROOT_FILES:
        assert not (REPO_ROOT / name).exists(), f"Generated artifact must not be in root: {name}"

    for pattern in FORBIDDEN_ROOT_GLOBS:
        matches = [p.name for p in REPO_ROOT.glob(pattern) if p.is_file()]
        assert not matches, f"Generated artifacts must not be in root for pattern '{pattern}': {matches}"


def test_docs_index_consistency():
    docs_index = REPO_ROOT / "docs" / "INDEX.md"
    inventory = REPO_ROOT / "docs" / "INVENTORY_CLASSIFICATION.md"
    legacy_readme = REPO_ROOT / "legacy" / "README.md"

    assert docs_index.exists(), "Missing docs/INDEX.md"
    assert inventory.exists(), "Missing docs/INVENTORY_CLASSIFICATION.md"
    assert legacy_readme.exists(), "Missing legacy/README.md"

    text = docs_index.read_text(encoding="utf-8")
    assert "INVENTORY_CLASSIFICATION.md" in text
    assert "../legacy/README.md" in text
