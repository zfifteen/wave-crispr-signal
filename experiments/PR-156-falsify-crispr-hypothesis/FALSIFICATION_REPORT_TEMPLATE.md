# Falsification Report Template: PR-156 CRISPR Hypotheses

**Experiment ID**: PR-156-falsify-crispr-hypothesis  
**Report Date**: [YYYY-MM-DD]  
**Investigator**: [Your Name]  
**Repository**: https://github.com/zfifteen/wave-crispr-signal  

---

## Executive Summary

This report documents the results of falsification testing for two hypotheses related to the geodesic resolution function θ′(n,k) and curvature weighting κ(n) in CRISPR guide efficiency prediction.

**Key Findings**:
- Hypothesis 1: [FALSIFIED / NOT FALSIFIED]
- Hypothesis 2: [FALSIFIED / NOT FALSIFIED]

**Recommendation**: [Based on results, what actions should be taken?]

---

## Hypothesis 1: θ′(n,k) CRISPR Guide Efficiency Prediction

### Statement
The geodesic resolution function θ′(n,k) with k≈0.3, embedding golden ratio phasing, quantifies single-nucleotide mutation disruptions via Δentropy and sidelobe metrics, improving CRISPR guide efficiency prediction by ΔROC-AUC +0.047 over baselines like RuleSet3.

### Methods
- **Dataset**: [Specify dataset used: synthetic, Doench 2016, Kim 2025, etc.]
- **Sample size**: [N sequences]
- **Cross-validation**: [N-fold CV]
- **Bootstrap resamples**: [N resamples]
- **Seed**: [Random seed for reproducibility]

### Results

#### Performance Metrics

| Metric | Spectral Features | Baseline Features | Δ (Spectral - Baseline) |
|--------|-------------------|-------------------|-------------------------|
| Mean ROC-AUC | [X.XXX ± X.XXX] | [X.XXX ± X.XXX] | [X.XXX] |
| 95% CI | [X.XXX, X.XXX] | [X.XXX, X.XXX] | [X.XXX, X.XXX] |

#### Statistical Tests

- **Paired t-test**: t = [X.XXX], p = [X.XXX]
- **Bootstrap mean ΔROC-AUC**: [X.XXX ± X.XXX]
- **Bootstrap 95% CI**: [[X.XXX, X.XXX]]

#### Falsification Criteria

- [ ] ΔROC-AUC > 0
- [ ] p-value < 0.05 (statistically significant)
- [ ] 95% CI excludes zero

### Conclusion

**Hypothesis 1 is [FALSIFIED / NOT FALSIFIED]**

**Reasons**:
- [List specific reasons based on criteria]

**Interpretation**:
[Discuss what this means for the hypothesis and future work]

---

## Hypothesis 2: Z-Invariant Scoring Abstraction

### Statement
The sequence diversity-dependent curvature weighting κ(n) = d(n) · ln(n+1)/e² creates a Z-invariant scoring abstraction that robustly handles variable-length DNA sequences, revealing emergent periodicity patterns in off-target profiling.

### Methods
- **Length range**: [Min-Max nt, step size]
- **Sequences per length**: [N sequences]
- **Motif preservation**: [Yes/No]
- **Seed**: [Random seed for reproducibility]

### Results

#### Z-Score Invariance

**ANOVA Test**:
- F-statistic: [X.XXX]
- p-value: [X.XXX]
- Overall variance: [X.XXX]

**Mean Z-scores by length**:

| Length (nt) | Mean Z-score | Std Dev |
|-------------|--------------|---------|
| 20 | [X.XXX] | [X.XXX] |
| 30 | [X.XXX] | [X.XXX] |
| 40 | [X.XXX] | [X.XXX] |
| ... | ... | ... |

#### Periodicity Analysis

- **Max autocorrelation**: [X.XXX]
- **Periodicity detected**: [Yes/No]

**Autocorrelation values** (first 10 lags):
```
[X.XXX, X.XXX, X.XXX, ...]
```

#### Falsification Criteria

- [ ] Z-scores invariant across lengths (ANOVA p ≥ 0.05)
- [ ] Emergent periodicity detected (max autocorr > 0.1)

### Conclusion

**Hypothesis 2 is [FALSIFIED / NOT FALSIFIED]**

**Reasons**:
- [List specific reasons based on criteria]

**Interpretation**:
[Discuss what this means for the hypothesis and future work]

---

## Overall Conclusions

### Summary of Findings

| Hypothesis | Status | Confidence |
|------------|--------|------------|
| H1: θ′(n,k) improves efficiency prediction | [F/NF] | [High/Medium/Low] |
| H2: κ(n) creates Z-invariant abstraction | [F/NF] | [High/Medium/Low] |

### Implications

**If both hypotheses are NOT falsified**:
- The Z Framework's geometric resolution function and curvature weighting are supported by empirical evidence
- Consider advancing to real-world dataset validation (Doench 2016, Kim 2025)
- Prepare for comparison with state-of-the-art baselines (RuleSet3, DeepCRISPR)

**If Hypothesis 1 is falsified**:
- The spectral features derived from θ′(n,k) may not provide practical improvement over simple baselines
- Consider alternative k values, different phasing schemes, or additional spectral features
- Re-evaluate the claim of ΔROC-AUC +0.047

**If Hypothesis 2 is falsified**:
- The Z-invariant property may not hold across variable-length sequences
- Consider refining κ(n) formula or using different normalization schemes
- Investigate why invariance fails (length-dependent bias in encoding?)

**If both hypotheses are falsified**:
- The current formulation may need fundamental revision
- Return to theoretical foundations and mathematical derivations
- Consider whether the framework is applicable to CRISPR guide prediction at all

### Recommendations

1. **Short-term** (1-2 weeks):
   - [Specific actionable recommendations based on results]

2. **Medium-term** (1-2 months):
   - [Strategic recommendations for methodology refinement]

3. **Long-term** (3-6 months):
   - [Big-picture recommendations for research direction]

### Limitations

- **Dataset**: [Limitations of dataset used]
- **Sample size**: [Statistical power considerations]
- **Baseline comparison**: [RuleSet3 not yet integrated, using proxy features]
- **Synthetic data**: [If synthetic data was used, acknowledge limitations]

### Future Work

1. **Dataset integration**: Incorporate real CRISPR efficiency datasets
2. **Baseline implementation**: Integrate RuleSet3 for direct comparison
3. **Extended analysis**: Test additional k values, different φ periods
4. **Visualization**: Create ROC curves, Z-score distributions
5. **Off-target profiling**: Validate on off-target datasets

---

## Reproducibility

### Environment

- **Python version**: 3.12.x
- **Key dependencies**: See `requirements.txt`
- **Random seed**: [seed value]
- **Commit hash**: [git commit hash]

### Commands to Reproduce

```bash
# Hypothesis 1
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis1.py \
  --seed [seed] \
  --n-samples [N] \
  --n-bootstrap [N] \
  --n-folds [N] \
  --k-parameter [k] \
  --output-dir results/PR-156-falsify-hypothesis1

# Hypothesis 2
python experiments/PR-156-falsify-crispr-hypothesis/falsify_hypothesis2.py \
  --seed [seed] \
  --min-length [min] \
  --max-length [max] \
  --length-step [step] \
  --n-per-length [N] \
  --k-parameter [k] \
  --output-dir results/PR-156-falsify-hypothesis2

# Combined run
python experiments/PR-156-falsify-crispr-hypothesis/run_all_falsification_tests.py \
  --seed [seed] \
  [... other parameters ...]
```

### Data Files

Results saved to:
- `results/PR-156-falsify-hypothesis1/hypothesis1_results.json`
- `results/PR-156-falsify-hypothesis2/hypothesis2_results.json`
- `results/PR-156-falsify-combined/combined_report.json`

---

## Appendices

### Appendix A: Spectral Features (Hypothesis 1)

**Computed features**:
1. **Δentropy**: Change in spectral entropy between wild-type and mutant
2. **Δsidelobes**: Change in number of spectral sidelobes
3. **Frequency shift**: Shift in dominant frequency component

**Baseline features**:
1. **GC content**: Overall, first 5 nt, last 5 nt
2. **Sequence length**

### Appendix B: Z-Score Calculation (Hypothesis 2)

**Formula**:
```
Z = S(H / κ)
```

where:
- `H` = Spectral entropy from phase-weighted FFT
- `κ(n)` = `d(n) · ln(n+1) / e²` (curvature weight)
- `d(n)` = Sequence diversity (normalized Shannon entropy)
- `S(x)` = Sigmoid `1 / (1 + exp(-κ·x))`

### Appendix C: References

1. Doench et al. (2016). "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nature Biotechnology. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4745674/

2. Hsu et al. (2013). "DNA targeting specificity of RNA-guided Cas9 nucleases." Nature Biotechnology. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3773023/

3. Fourier methods in DNA analysis: https://pubmed.ncbi.nlm.nih.gov/16177379/

4. Golden ratio in biological modeling: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4345425/

---

**Report prepared by**: [Name]  
**Date**: [YYYY-MM-DD]  
**Contact**: [Email/GitHub]  
**License**: Research use only
