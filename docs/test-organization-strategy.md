# Test Organization Strategy

## Current Structure

This repository maintains **dual test locations** for historical and practical reasons:

### Root-Level Tests (Legacy)
- `test_*.py` files in repository root
- Executed via `python run_tests.py` orchestrator
- Covers core framework modules: Z Framework, Invariant Features, Geodesic Bridge
- **Status**: Actively maintained, required for compatibility

### Tests Directory (Modern)
- `tests/` subdirectory for organized test suites
- Currently contains: `test_crispr_guide_designer.py`
- Follows pytest discovery conventions
- **Status**: Modern standard for new test files

## Execution Strategy

### Primary Test Runner
```bash
python run_tests.py  # Executes core framework tests
```

### Pytest Discovery
```bash
pytest tests/        # Runs organized test suites
pytest .             # Discovers all test files
```

### CI Integration
The policy enforcement workflow (`policy-check.yml`) uses `run_tests.py` as the authoritative test executor to ensure all critical framework tests pass.

## Migration Policy

**Current Decision**: Maintain dual structure during transition period

### For New Components
- Place new test files in `tests/` directory
- Follow `tests/test_component_name.py` naming pattern
- Use pytest conventions and fixtures

### For Existing Components
- Root-level tests remain in place for stability
- Critical path validation stays with `run_tests.py`
- Future consolidation considered based on project evolution

## Coverage Requirements

Both test locations contribute to overall coverage:
- **Minimum**: ≥80% coverage for core modules
- **Measurement**: Via pytest-cov during CI
- **Enforcement**: Policy check workflow validates coverage threshold

## Benefits of Dual Structure

1. **Backward Compatibility**: Existing workflows unaffected
2. **Gradual Migration**: New code uses modern patterns
3. **Tool Flexibility**: Both pytest and custom runners supported
4. **Discovery**: pytest automatically finds both locations

---

**Policy Compliance**: This dual structure is explicitly permitted by REPOSITORY_POLICY.md section "Test Infrastructure" which acknowledges both locations.