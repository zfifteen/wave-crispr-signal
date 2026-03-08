# Makefile for wave-crispr-signal repository

.PHONY: help install test smoke clean lint run-mve run-mve-quick run-mri-z5d run-mri-z5d-quick run-mri-z5d-synthetic run-fus-enhancer run-fus-enhancer-quick

help:
	@echo "Available targets:"
	@echo "  install               - Install dependencies"
	@echo "  test                  - Print change-scoped validation guidance"
	@echo "  smoke                 - Placeholder (CI currently disabled)"
	@echo "  clean                 - Clean temporary files"
	@echo "  lint                  - Run linting"
	@echo "  run-mve*              - Run MVE experiments"
	@echo "  run-mri-z5d*          - Run MRI Z5D experiments"
	@echo "  run-fus-enhancer*     - Run FUS enhancer experiments"

install:
	pip install -r requirements.txt

test:
	@echo "No baseline test suite is configured."
	@echo "Run change-scoped validation commands for the behavior you changed."
	@echo "Example: python scripts/run_tests.py \"python -m pytest -q path/to/test.py\""

smoke:
	@echo "CI smoke target is intentionally disabled during reset."

clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -exec rm -rf {} +
	rm -rf .pytest_cache/
	rm -rf results/focused_ultrasound_mve/run-*/
	rm -rf /tmp/focused_ultrasound_mve/

lint:
	ruff check .

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

run-mri-z5d:
	python experiments/mri_z5d_analysis.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--signal-length 256 \
		--output-dir results \
		--k-parameter 0.04449 \
		--max-files-per-series 20

run-mri-z5d-quick:
	python experiments/mri_z5d_analysis.py \
		--seed 42 \
		--bootstrap 100 \
		--permutation 100 \
		--signal-length 64 \
		--output-dir results \
		--max-files-per-series 5

run-mri-z5d-synthetic:
	python experiments/mri_z5d_analysis.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--signal-length 256 \
		--output-dir results \
		--k-parameter 0.04449 \
		--use-synthetic

run-fus-enhancer:
	python fus_enhancer.py \
		--seed 42 \
		--bootstrap 1000 \
		--permutation 1000 \
		--n-trials 1000000 \
		--batch-size 100000 \
		--grid-size 100 \
		--k-parameter 0.3 \
		--output-dir results \
		--visualize

run-fus-enhancer-quick:
	python fus_enhancer.py \
		--seed 42 \
		--bootstrap 100 \
		--permutation 100 \
		--n-trials 10000 \
		--batch-size 1000 \
		--grid-size 50 \
		--output-dir results
