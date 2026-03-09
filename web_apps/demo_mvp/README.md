# Demo MVP

FastAPI demo app for sequence-level breathing/spectral analysis.

## Entry Point

- `web_apps/demo_mvp/app.py`

## Run

```bash
python3 web_apps/demo_mvp/app.py
```

Or:

```bash
uvicorn web_apps.demo_mvp.app:app --reload
```

## Routes

- `GET /` - UI page.
- `GET /health` - health status.
- `GET /info` - parameter and mode metadata.
- `POST /analyze` - guide/target analysis response with spectral peaks, rotation curve, and derived features.

## Analyze Request Fields

- `guide` (required)
- `target` (required)
- `rate_ratio` (optional, default `20.0`)
- `weight_mode` (optional, accepts `rate_ratio` or `nn_thermo`)
- `temperature_c` (optional)
- `mg_mM` (optional)
- `log_consent` (optional, default `false`)

Current behavior note:
- `/info` reports `nn_thermo` as not yet implemented.
- `/analyze` currently computes using the rate-ratio path.
