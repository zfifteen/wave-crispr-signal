# DNA Breathing Dynamics ΔΔG Validation — Quick Summary

## Conclusion: HYPOTHESIS NOT SUPPORTED ❌

Testing the hypothesis from [PR #53](https://github.com/zfifteen/dna-breathing-dynamics-encoding/pull/53):
**DNA breathing dynamics encoding does NOT show strong sensitivity to thermodynamic perturbations.**

## Key Results

- **Maximum Effect Size**: |Cohen's d| = 0.287 (need ≥ 0.5)
- **Statistical Significance**: Yes (q = 0.014)
- **Practical Significance**: No (effect too small)
- **Dose-Response Trends**: 75/75 monotonic (all p < 0.05)

## Interpretation

While **statistically significant** trends exist, the **effect sizes are too small** for practical use. The maximum effect (0.287) is 57.4% of the pre-registered threshold (0.5).

## Files

- **FINDINGS.md** — Full results report with evidence
- **README.md** — Experiment documentation
- **ddg_validation.py** — Complete validation pipeline
- **manifest.yml** — Experiment metadata

## Quick Test

```bash
# Smoke test (3 seconds)
python ddg_validation.py --input ../../data/doench2016.csv \
    --output /tmp/smoke --n-pairs 30 --seed 42 --smoke
```

## Quality Assurance

✅ Code review passed  
✅ Security scan passed (0 alerts)  
✅ Repository policy compliant  
✅ Fully reproducible (seed 42)
