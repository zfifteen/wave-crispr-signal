# Scripts

Operational scripts for maintenance and one-off workflows.

## Available Scripts

- `run_tests.py` - change-scoped validation helper

## Usage

```bash
# No args: prints reset guidance and exits 0
python scripts/run_tests.py

# Run specific targeted checks
python scripts/run_tests.py "python -m pytest -q path/to/test_file.py" "python applications/phase_weighted_scorecard_cli.py score --guide GCTGCGGAGACCTGGAGAGA"
```
