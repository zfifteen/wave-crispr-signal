# Repository Structure Policy

## Overview

This document establishes the formal repository structure and organization standards for the **wave-crispr-signal** project. This policy ensures consistency, maintainability, and discoverability of project components across scientific research, software development, and collaboration workflows.

## Project Classification

**Type**: Scientific Research Repository  
**Domain**: Computational Biology, Signal Processing, Mathematical Analysis  
**Language**: Python 3.12+  
**License**: MIT  

## Directory Structure Standards

### Root Level Organization

The repository follows a **flat scientific module** pattern with supporting infrastructure:

```
wave-crispr-signal/
├── .github/                    # GitHub configuration and policies
├── .grok/                      # AI assistant instructions  
├── applications/               # Domain-specific applications and tools
├── data/                       # Research datasets and reference data
├── docs/                       # Comprehensive documentation
├── notebooks/                  # Jupyter analysis notebooks
├── proof_pack/                 # Validation and verification tools
├── static/                     # Web interface static assets
├── templates/                  # Web interface templates
├── tests/                      # Test infrastructure (legacy location)
├── wave-ribonn/               # Specialized analysis subdirectory
├── *.py                       # Core scientific modules (root level)
├── *.md                       # Research papers and reports (root level)
└── *.csv, *.txt              # Immediate access datasets (root level)
```

### Detailed Directory Specifications

#### `.github/` - Repository Management
- **Purpose**: GitHub-specific configuration, workflows, and policies
- **Required Files**: 
  - `REPOSITORY_POLICY.md` (this file)
  - `copilot-instructions.yml`
- **Naming**: Kebab-case for files
- **Access**: Public configuration

#### `.grok/` - AI Assistant Configuration  
- **Purpose**: Instructions and guidelines for AI code assistants
- **Required Files**: `grok_instructions_*.md`
- **Naming**: Snake_case with descriptive suffixes
- **Access**: Project-internal

#### `applications/` - Domain Applications
- **Purpose**: CRISPR guide design tools, CLI interfaces, domain-specific applications
- **Structure**: Python package with `__init__.py`
- **Required Files**: 
  - `__init__.py` (package marker)
  - `crispr_*.py` (CRISPR-specific tools)
  - `wave_*.py` (signal processing tools)
- **Naming**: Snake_case with domain prefixes
- **Documentation**: Each major application must have corresponding docs in `docs/applications/`

#### `data/` - Research Datasets
- **Purpose**: Reference datasets, validation data, research inputs
- **Formats**: CSV (preferred), TXT, JSON
- **Naming**: Descriptive lowercase with underscores
- **Size Limit**: <10MB per file (larger files use external references)
- **Documentation**: Dataset descriptions in `README.md` within subdirectories

#### `docs/` - Documentation Hub
- **Purpose**: Comprehensive project documentation
- **Structure**: 
  ```
  docs/
  ├── applications/           # Application-specific docs
  ├── PR-*/                  # Pull request documentation
  ├── *.md                   # Core methodology docs
  ```
- **Naming**: UPPERCASE for major docs, kebab-case for files
- **Requirements**: Every major component must have documentation

#### `notebooks/` - Interactive Analysis
- **Purpose**: Jupyter notebooks for research, validation, examples
- **Naming**: Descriptive with version numbers where applicable
- **Requirements**: All notebooks must be runnable with current `requirements.txt`

#### `proof_pack/` - Validation Framework
- **Purpose**: Standalone validation tools, verification scripts, proof-of-concept demonstrations
- **Structure**: Self-contained with own `README.md`
- **Naming**: Snake_case for scripts, kebab-case for documentation
- **Requirements**: Must be executable independently

#### `tests/` - Test Infrastructure
- **Purpose**: Automated testing for core modules
- **Organization**: 
  - Root level: Legacy test files (`test_*.py`)
  - `tests/` subdirectory: Modern test structure
- **Requirements**: All core modules must have corresponding tests
- **Execution**: Via `run_tests.py` orchestrator

#### `static/` & `templates/` - Web Interface
- **Purpose**: Flask/web application assets and templates
- **Structure**: Standard web application layout
- **Naming**: Standard web conventions (CSS, JS, HTML)

### Core Module Organization (Root Level)

#### Scientific Modules
Core mathematical and computational modules reside at root level for direct import:

- `z_framework.py` - Core Z Framework implementation
- `topological_analysis.py` - Topological analysis functions  
- `invariant_features.py` - Invariant feature calculations
- `molecular_dynamics_framework.py` - MD simulation tools
- `*_hypothesis.py` - Hypothesis testing modules
- `*_validation.py` - Validation frameworks

#### Research Documents (Root Level)
Major research outputs and reports:

- `README.md` - Primary project documentation
- `*_WHITE_PAPER.md` - Research papers and academic outputs
- `*_REPORT.md` - Empirical findings and analysis reports
- `*_SUMMARY.md` - Issue resolution and methodology summaries
- `EXPERIMENTAL_*.md` - Experimental protocols and setups

#### Immediate Access Data (Root Level)
Frequently accessed datasets:

- `*.csv` - Tabular research data
- `*.txt` - Text-based datasets (e.g., `data/zeta.txt`)
- `requirements.txt` - Python dependencies

## File Naming Conventions

### Python Files
- **Modules**: `snake_case.py`
- **Domain Prefixes**: Use consistent prefixes (`crispr_`, `wave_`, `z5d_`, etc.)
- **Test Files**: `test_*.py`
- **Demo Files**: `demo_*.py`
- **Validation**: `*_validation.py`

### Documentation Files
- **Major Documents**: `UPPERCASE_WITH_UNDERSCORES.md`
- **Application Docs**: `TITLE_CASE.md`
- **Technical Docs**: `kebab-case.md` in subdirectories

### Data Files
- **Datasets**: `lowercase_with_underscores.csv`
- **Configuration**: `kebab-case.yml`, `snake_case.toml`

## Code Organization Standards

### Module Structure
```python
"""
Module docstring with Z Framework context
"""

# Standard library imports
import os
import sys

# Third-party imports  
import numpy as np
import matplotlib.pyplot as plt

# Local imports
from applications.crispr_guide_designer import CRISPRGuideDesigner
from z_framework import ZFrameworkCore
```

### Documentation Requirements
1. **Docstrings**: All public functions must have comprehensive docstrings
2. **Type Hints**: Required for new code (Python 3.12+ style)
3. **Examples**: Include usage examples in module docstrings
4. **Z Framework Context**: Reference Z Framework principles where applicable

### Testing Requirements
1. **Coverage**: All core modules must have test files
2. **Execution**: Tests must pass via `run_tests.py`
3. **Dependencies**: Tests must work with pinned `requirements.txt`
4. **Validation**: Include both unit tests and integration validation

## Dependency Management

### Requirements Management
- **File**: `requirements.txt` with exact version pinning
- **Format**: `package==exact.version.number`
- **Comments**: Include research context and version justification
- **Validation**: All requirements must support core functionality

### Scientific Dependencies
**Core Requirements**:
- `mpmath==1.3.0` (high-precision mathematics)
- `numpy==1.26.4` (numerical computing)
- `scipy==1.16.1` (scientific computing)
- `sympy==1.13.1` (symbolic mathematics)
- `matplotlib==3.10.5` (visualization)

**Domain-Specific**:
- `biopython==1.83` (biological sequence analysis)
- `scikit-learn==1.5.1` (machine learning)

## Quality Assurance

### Code Quality
1. **Linting**: Follow scientific Python conventions
2. **Testing**: Maintain >80% test coverage for core modules
3. **Documentation**: Every public API must be documented
4. **Validation**: Include empirical validation for mathematical claims

### Research Integrity
1. **Reproducibility**: All analyses must be reproducible with pinned dependencies
2. **Validation**: Include validation scripts for major claims
3. **Documentation**: Research methodologies must be documented
4. **Attribution**: Proper citation and attribution for external methods

## Continuous Integration

### Test Execution
- **Command**: `python run_tests.py`
- **Requirements**: All tests must pass before merge
- **Coverage**: Test both individual modules and integration

### Documentation Validation
- **Consistency**: Documentation must match implementation
- **Examples**: All code examples must be functional
- **Links**: Internal links must be valid

## Git Workflow Standards

### Branch Naming
- **Features**: `feature/descriptive-name`
- **Fixes**: `fix/issue-number-description`
- **Documentation**: `docs/component-update`

### Commit Messages
- **Format**: `Component: Brief description`
- **Examples**: 
  - `Z Framework: Add geometric resolution validation`
  - `CRISPR Tools: Improve guide scoring algorithm`
  - `Docs: Update repository policy standards`

### File Exclusions (`.gitignore`)
```
# Python
__pycache__/
*.py[cod]
*$py.class
*.so

# Research
*.png
*.jpg
*.pdf
results/
output/

# Environment
.env
.venv/
env/
venv/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db
```

## Enforcement and Compliance

### Automated Enforcement
- **CI/CD**: Repository structure validation in automated workflows
- **Linting**: Code style enforcement
- **Testing**: Comprehensive test suite execution

### Manual Review
- **Pull Requests**: All changes must maintain policy compliance
- **Documentation**: Updates must include corresponding documentation
- **Testing**: New features must include appropriate tests

### Policy Updates
- **Process**: Changes to this policy require issue discussion and approval
- **Notification**: Policy changes must be announced to all contributors
- **Versioning**: Policy changes should be versioned and tracked

## AI Assistant Integration

### Copilot Configuration
This policy is integrated with GitHub Copilot via `.github/copilot-instructions.yml` to ensure:

1. **Structure Adherence**: Automatic compliance with directory structure
2. **Naming Consistency**: Consistent file and variable naming
3. **Documentation Standards**: Proper docstring and comment generation
4. **Testing Requirements**: Automatic test file creation and updates

### Grok Configuration
Grok AI assistant instructions in `.grok/` provide:

1. **Scientific Context**: Z Framework methodology and principles
2. **Domain Knowledge**: CRISPR, signal processing, and mathematical analysis context
3. **Code Standards**: Project-specific coding patterns and conventions

---

**Policy Version**: 1.0  
**Last Updated**: 2025-01-20  
**Next Review**: 2025-04-20  

This policy ensures consistent, maintainable, and scientifically rigorous repository organization while supporting both human contributors and AI-assisted development workflows.