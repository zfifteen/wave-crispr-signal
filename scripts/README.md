# Scripts

This directory contains CLI entry points and automation scripts for the wave-crispr-signal project.

## Script Overview

- **app.py** - Flask web application entry point for the web interface
- **wave-crispr-signal.py** - Main CLI script for command-line usage  
- **run_tests.py** - Test suite runner and orchestrator

## Categories

### Web Interface
- **app.py** - Provides web-based access to wave-crispr-signal functionality

### Command Line Interface
- **wave-crispr-signal.py** - Main command-line interface for the project

### Development & Testing
- **run_tests.py** - Automated test execution and reporting

## Usage

### Web Interface
```bash
python scripts/app.py
```

### Command Line
```bash
python scripts/wave-crispr-signal.py [options]
```

### Test Execution
```bash
python scripts/run_tests.py
```

## Design Principles

- Scripts should serve as entry points, not contain core logic
- Import functionality from `modules/` and use `experiments/` and `tools/` as needed
- Each script should include proper argument parsing and help text
- Scripts should handle errors gracefully and provide useful feedback

## Configuration

Scripts may read configuration from the `configs/` directory when available.