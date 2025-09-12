# Makefile for wave-crispr-signal repository

.PHONY: test smoke help clean install lint

# Default target
help:
	@echo "Available targets:"
	@echo "  test          - Run full test suite"
	@echo "  smoke         - Run smoke tests for CI"
	@echo "  mve-smoke     - Run focused ultrasound MVE smoke test"
	@echo "  install       - Install dependencies"
	@echo "  clean         - Clean up temporary files"
	@echo "  lint          - Run code linting"

# Install dependencies
install:
	pip install -r requirements.txt

# Run full test suite
test:
	python run_tests.py

# Run smoke tests for CI
smoke: mve-smoke test-z-framework-import mri-z5d-smoke
	@echo "✓ All smoke tests completed"

# Test Z Framework import compatibility
test-z-framework-import:
	@echo "Testing Z Framework import compatibility..."
	@python -c "from z_framework import ZFrameworkCalculator; calc = ZFrameworkCalculator(); print('✓ Z Framework import successful')"

# Run focused ultrasound MVE smoke test
mve-smoke:
	@echo "Running Focused Ultrasound MVE smoke test..."
	python tests/test_focused_ultrasound_mve.py --smoke

# Run MRI Z5D analysis smoke test
mri-z5d-smoke:
	@echo "Running MRI Z5D Analysis smoke test..."
	@python -c "import sys; sys.path.insert(0, '.'); from experiments.mri_z5d_analysis import Z5DGeodeskAnalyzer, generate_synthetic_mri_signals; analyzer = Z5DGeodeskAnalyzer(seed=42); signals = generate_synthetic_mri_signals(n_samples=3, signal_length=10, seed=42); results = [analyzer.analyze_signal_pattern(s, f'smoke_{i}') for i, s in enumerate(signals)]; print('✓ MRI Z5D Analysis smoke test passed')"

# Clean up temporary files
clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -exec rm -rf {} +
	rm -rf .pytest_cache/
	rm -rf results/focused_ultrasound_mve/run-*/
	rm -rf /tmp/focused_ultrasound_mve/

# Lint code (placeholder)
lint:
	@echo "Linting (placeholder - add specific linter when available)"

# Run focused ultrasound MVE experiment with default parameters
run-mve:
	python experiments/focused_ultrasound_mve.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--splits single \
		--domain discrete \
		--k-parameter 0.3 \
		--n-trials 1000 \
		--grid-size 100 \
		--output-dir results \
		--visualize

# Run quick MVE test with reduced parameters
run-mve-quick:
	python experiments/focused_ultrasound_mve.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--splits single \
		--domain discrete \
		--k-parameter 0.3 \
		--n-trials 100 \
		--grid-size 50 \
		--output-dir results

# Run MRI Z5D analysis with default parameters
run-mri-z5d:
	python experiments/mri_z5d_analysis.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--n-samples 100 \
		--signal-length 256 \
		--output-dir results \
		--k-parameter 0.04449

# Run quick MRI Z5D test with reduced parameters
run-mri-z5d-quick:
	python experiments/mri_z5d_analysis.py \
		--seed 42 \
		--bootstrap 100 \
		--permutation 100 \
		--n-samples 20 \
		--signal-length 64 \
		--output-dir results