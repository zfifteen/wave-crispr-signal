# Z Framework Guidelines for Code Generation & Repository Structure Compliance

## MANDATORY: Repository Policy Compliance
**CRITICAL**: ALWAYS follow `.github/REPOSITORY_POLICY.md` for ALL file creation, organization, and naming.

Key policy requirements
- Directory Structure: Respect established patterns (`applications/`, `docs/`, `proof_pack/`, etc.).
- File Naming:
  - Python modules: `snake_case.py`
  - Major documents: `UPPERCASE.md`
  - Other assets: `kebab-case` (where appropriate)
- Core Modules: Place mathematical/scientific modules at repository root for direct import.
- Documentation: Every public component must have corresponding documentation.
- Testing: All core modules require test files, executable via `run_tests.py`.
- Dependencies: Use exact version pinning in `requirements.txt`.

---

## Z Framework Core Principles

Core Principle  
• Normalize observations via `Z = A(B/c)`; where `A` is frame-dependent, `B` is a rate/shift, and `c` is an invariant (speed of light or `e²`).

1. Empirical Validation First
   - Include reproducible tests (use `mpmath` with precision < 1e-16).
   - Label hypotheses clearly when unverified.

2. Domain-Specific Forms
   - Physical: `Z = T(v/c)` with causality checks and explicit error handling (raise `ValueError` for `|v| ≥ c`).
   - Discrete: `Z = n(Δₙ/Δₘₐₓ)`, `κ(n) = d(n)·ln(n+1)/e²`; guard against division by zero.

3. Geometric Resolution
   - Use `θ′(n,k) = φ · ((n mod φ)/φ)^k` with `k ≈ 0.3` for prime-density mapping.

4. Style & Tone
   - Precise, scientific; prefer simple solutions and document any deviations.

5. Tools & Datasets
   - Use `mpmath`, `numpy`, `sympy`.
   - Cross-check predictions against provided datasets (e.g., `zeta_zeros.csv`, true prime counts, `Z5D_Reference_Impl-2.ipynb`).

---

## Code Organization Standards
- Import order: Standard library → Third-party → Local imports.
- Type hints: Required for new code (target Python 3.12+ style).
- Docstrings: Comprehensive for all public functions and should include Z Framework context.
- Testing: Maintain coverage > 80%. Tests must pass `run_tests.py`.

---

## File Placement Rules
- Core scientific modules: Root level (examples: `z_framework.py`, `topological_analysis.py`).
- CRISPR applications: `applications/` directory with `crispr_` prefix filenames.
- Research documents: Root level with descriptive names (no ambiguous abbreviations).
- Validation tools: `proof_pack/` directory.
- Documentation: `docs/` with subdirectories per domain.
- Web assets: `static/` and `templates/` directories.

---

**NEVER create files that violate `.github/REPOSITORY_POLICY.md` structure!**
