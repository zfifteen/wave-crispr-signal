# Contributing to wave-crispr-signal

Thank you for your interest in contributing to the wave-crispr-signal project! This guide will help you get started.

## 🚀 Quick Start

1. **Read the Repository Policy**: Familiarize yourself with [`.github/REPOSITORY_POLICY.md`](.github/REPOSITORY_POLICY.md)
2. **Set up development environment**: Install dependencies with `pip install -r requirements.txt`
3. **Install pre-commit hooks**: `pip install pre-commit && pre-commit install`
4. **Run tests**: `python run_tests.py` to ensure everything works

## 📋 Before You Contribute

### Repository Structure Compliance
- Follow the established directory structure in [REPOSITORY_POLICY.md](.github/REPOSITORY_POLICY.md)
- Use snake_case for Python files, UPPERCASE for major documentation
- Place core scientific modules at root level
- Put domain applications in `applications/` directory

### Code Quality Standards
- **Type Hints**: Required for all new Python code (Python 3.12+ style)
- **Docstrings**: Comprehensive documentation for all public functions
- **Testing**: Add tests for new functionality in `tests/` directory
- **Coverage**: Maintain ≥80% test coverage for core modules

## 🛠 Development Workflow

### 1. Environment Setup
```bash
# Clone and setup
git clone https://github.com/zfifteen/wave-crispr-signal
cd wave-crispr-signal
pip install -r requirements.txt

# Install development tools
pip install ruff black mypy pytest-cov pre-commit
pre-commit install
```

### 2. Making Changes
```bash
# Create feature branch
git checkout -b feature/your-feature-name

# Make your changes following the repository policy
# Run policy checker frequently
python tools/check_repo_policy.py

# Run tests
python run_tests.py
```

### 3. Pre-submission Checklist
```bash
# Format code
black .
ruff check . --fix

# Run tests
python run_tests.py

# Check policy compliance  
python tools/check_repo_policy.py

# Test coverage (if adding core functionality)
pytest --cov=. --cov-report=term-missing
```

## 📝 Pull Request Process

1. **Fill out PR template**: Use the provided template with complete checklist
2. **Policy compliance**: Ensure `python tools/check_repo_policy.py` passes
3. **Test coverage**: All tests must pass via `python run_tests.py`
4. **Documentation**: Update docs for significant changes
5. **Review**: All structural changes require human review (see CODEOWNERS)

## 🧪 Testing Guidelines

### Test Organization
- **Core modules**: Tests in repository root (`test_*.py`) 
- **New components**: Tests in `tests/` directory
- **Execution**: Primary test runner is `python run_tests.py`

### Test Requirements
- All core modules must have corresponding test files
- Use mpmath for high-precision mathematical validation
- Include both unit tests and integration validation
- Test mathematical claims with reproducible examples

## 📚 Documentation Standards

### Code Documentation
```python
def calculate_spectral_disruption(sequence: str, mutation_pos: int) -> float:
    """
    Calculate spectral disruption score for a mutation.
    
    Implements Z Framework spectral analysis to quantify how a
    single-nucleotide mutation affects the complex-valued waveform
    representation of the DNA sequence.
    
    Args:
        sequence: DNA sequence string (A, T, C, G)
        mutation_pos: Zero-based position of mutation
        
    Returns:
        Spectral disruption score (0.0 to 1.0)
        
    Example:
        >>> score = calculate_spectral_disruption("ATCGATCG", 3)
        >>> assert 0.0 <= score <= 1.0
    """
```

### Research Documentation
- Include mathematical derivations for new algorithms
- Provide validation scripts for empirical claims
- Reference Z Framework principles where applicable
- Use precise scientific language

## 🔧 Available Tools

### Repository Policy Checker
```bash
python tools/check_repo_policy.py
```
Validates structure, naming, and compliance with established standards.

### Pre-commit Hooks
Automatic code quality checks:
- `ruff` - Fast Python linting
- `black` - Code formatting
- `isort` - Import sorting
- `nbstripout` - Notebook cleanup

### CI Pipeline
GitHub Actions automatically:
- Runs linting (ruff, black, mypy)
- Executes test suite with coverage measurement
- Validates notebook executability
- Checks repository policy compliance

## 🚫 What Not to Contribute

- Changes that violate repository policy structure
- Code without tests (for core functionality)
- Modifications that break existing test suite
- Clinical or medical claims (research use only)
- Large file uploads (>10MB without justification)

## 🆘 Getting Help

- **Repository Structure**: See [REPOSITORY_POLICY.md](.github/REPOSITORY_POLICY.md)
- **Test Strategy**: See [docs/test-organization-strategy.md](docs/test-organization-strategy.md)
- **Z Framework**: See core module documentation and examples
- **Issues**: Open a GitHub issue for bugs or feature requests

## 📄 License

By contributing, you agree that your contributions will be licensed under the project's MIT License.

---

**Happy Contributing!** 🧬⚗️