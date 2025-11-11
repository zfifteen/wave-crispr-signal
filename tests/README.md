# Tests

This directory contains unit and integration tests for all modules and tools in the wave-crispr-signal project.

## Test Structure

Tests are organized to mirror the module structure using pytest conventions:

- **test_z_framework.py** - Tests for core Z Framework functionality
- **test_topological_analysis.py** - Tests for topological analysis functions
- **test_invariant_features.py** - Tests for invariant feature calculations
- **test_molecular_dynamics.py** - Tests for molecular dynamics framework
- **test_dna_storage_hypothesis.py** - Tests for DNA storage hypothesis
- **test_falsification_experiments.py** - Tests for falsification experiments
- **test_pain_management.py** - Tests for pain management applications
- **test_geodesic_bridge.py** - Tests for geodesic bridge functionality
- **test_crispr_*.py** - Tests for CRISPR-related functionality

## Coverage Requirements

- **Minimum Coverage**: ≥80% for core modules
- **Test Types**: Both unit tests and integration tests
- **Dependencies**: All tests use pinned dependencies from `requirements.txt`

## Test Execution

### Using pytest (recommended)
```bash
pytest tests/
pytest tests/test_z_framework.py  # Run specific test
```

### Using the test runner script
```bash
python scripts/run_tests.py
```

### Coverage measurement
```bash
pytest --cov=modules tests/
```

## Writing Tests

- Follow pytest conventions
- Use descriptive test names
- Include both positive and negative test cases
- Test edge cases and error conditions
- Mock external dependencies when appropriate

## CI Integration

These tests are executed automatically in the CI/CD pipeline to ensure code quality and functionality.