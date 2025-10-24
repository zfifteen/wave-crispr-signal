# Statistical Concerns and Corrections

## Critical Issue Identified

**Date**: 2025-01-10 (Post-PR creation)
**Reported by**: Repository reviewer
**Status**: ⚠️ **REQUIRES CORRECTION**

---

## Problem Statement

The current statistical analysis in `test_breathing_dynamics_encoding.py` (line 357) contains a **methodological flaw**:

```python
# CURRENT (INCORRECT):
t_stat, p_value = stats.ttest_ind(breath_scores, arb_means_list)
```

Where:
- `breath_scores`: n=100 (individual Z-scores, one per sequence)
- `arb_means_list`: n=10 (encoder-level means, one per arbitrary encoder)

This compares **different levels of aggregation** (individual scores vs. means), which **inflates the effect size** and may produce misleadingly low p-values.

---

## Why This Is Invalid

### 1. **Unequal Sample Sizes with Different Meanings**
- **Breathing group**: 100 independent sequence measurements
- **Arbitrary group**: 10 encoder-level aggregates (each summarizing 100 sequences)

The t-test treats these as comparable units, but they're not.

### 2. **Violation of Independence Assumption**
The 10 arbitrary values are **means**, not independent samples. Their variance is reduced by aggregation (Central Limit Theorem), making them appear more consistent than raw data.

### 3. **Inflated Effect Size**
Cohen's d uses pooled standard deviation:
```python
pooled_std = sqrt(
    ((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)
)
```

When comparing 100 individual values (higher variance) against 10 means (lower variance), the pooled std is **artificially low**, inflating Cohen's d.

### 4. **Assumption Violations**

**Normality**:
- 100 individual Z-scores: May or may not be normal
- 10 encoder means: By CLT, likely more normal
- But they're measuring different things!

**Equal variance**:
- Individual scores: Natural sampling variance
- Encoder means: Reduced variance due to averaging
- **Levene's test would likely fail**

---

## Correct Statistical Approaches

### **Option 1: Compare All Individual Scores** (RECOMMENDED)

```python
# Collect ALL individual scores for arbitrary encoders
all_arbitrary_scores = []
for trial in range(n_arbitrary_trials):
    arb_encoder = ArbitraryEncoder(seed=trial)
    arb_scores_gc = []

    for seq in sequences:
        z_gc = analyze_mutation_sensitivity(seq, arb_encoder, 'gc_affecting')
        arb_scores_gc.append(z_gc)

    all_arbitrary_scores.extend(arb_scores_gc)  # Collect ALL scores

# Now compare:
# breath_scores: n=100
# all_arbitrary_scores: n=1000 (100 sequences × 10 encoders)

t_stat, p_value = stats.ttest_ind(breath_scores, all_arbitrary_scores)
```

**Pros**:
- Correct comparison (both at sequence level)
- Respects independence
- Valid variance estimation

**Cons**:
- Assumes arbitrary encoder trials are independent (they are, different seeds)
- Larger arbitrary sample may have more power

---

### **Option 2: Hierarchical/Mixed-Effects Model** (MORE RIGOROUS)

Account for nested structure:
- **Level 1**: Sequences (n=100)
- **Level 2**: Encoders (1 breathing, 10 arbitrary)

```python
from statsmodels.formula.api import mixedlm

# Create dataframe with hierarchical structure
data = []
for seq_id, z_score in enumerate(breath_scores):
    data.append({'z_score': z_score, 'encoder_type': 'breathing',
                 'encoder_id': 0, 'seq_id': seq_id})

for encoder_id, result in enumerate(arbitrary_results):
    for seq_id, z_score in enumerate(result['gc_scores']):
        data.append({'z_score': z_score, 'encoder_type': 'arbitrary',
                     'encoder_id': encoder_id+1, 'seq_id': seq_id})

df = pd.DataFrame(data)

# Fit mixed model: z_score ~ encoder_type + (1|encoder_id) + (1|seq_id)
model = mixedlm("z_score ~ encoder_type", df, groups=df["encoder_id"])
result = model.fit()
```

**Pros**:
- Accounts for encoder-level and sequence-level variance
- Most statistically rigorous
- Properly handles nested structure

**Cons**:
- More complex
- Requires additional dependencies (statsmodels)

---

### **Option 3: Permutation Test** (NON-PARAMETRIC)

Avoid distributional assumptions entirely:

```python
# Observed difference
observed_diff = np.mean(breath_scores) - np.mean(all_arbitrary_scores)

# Permutation test
n_perms = 10000
perm_diffs = []

all_scores = np.concatenate([breath_scores, all_arbitrary_scores])

for _ in range(n_perms):
    np.random.shuffle(all_scores)
    perm_breath = all_scores[:len(breath_scores)]
    perm_arb = all_scores[len(breath_scores):]
    perm_diff = np.mean(perm_breath) - np.mean(perm_arb)
    perm_diffs.append(perm_diff)

# Two-tailed p-value
p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))
```

**Pros**:
- No distributional assumptions
- Robust to outliers
- Intuitive interpretation

**Cons**:
- Computationally intensive
- Doesn't provide effect size directly

---

## Re-Analysis with Correct Statistics

Running Option 1 (compare all individual scores):

### **Original (INCORRECT) Results:**
```
Breathing: n=100, mean=0.0801
Arbitrary: n=10 (MEANS!), mean=0.0570
Difference: +0.0231
Cohen's d: +4.130 (INFLATED!)
p-value: <0.000001
```

### **Corrected Results** (estimated):

If we properly compare 100 breathing scores against 1000 arbitrary scores (100 seqs × 10 encoders):

**Expected correction**:
- Cohen's d would be **smaller** (pooled std increases with more raw data)
- p-value might **increase** (though still likely significant)
- Effect might still be **real**, but **magnitude overstated**

**Need to re-run to get exact corrected values.**

---

## Impact on Conclusions

### **What Remains Valid:**
1. ✅ Breathing encoding shows **different** behavior than arbitrary
2. ✅ Selectivity (0.0000 for AT-affecting) is **real** - doesn't depend on stats
3. ✅ Directional finding (breathing > arbitrary for GC-affecting) likely **real**

### **What Is Questionable:**
1. ❌ **Magnitude** of Cohen's d = +4.130 is **inflated**
2. ❌ **"Very large effect"** claim may be **overstated**
3. ❌ Exact **p-value** is **not trustworthy**

### **Worst Case:**
- True Cohen's d might be +1.5 to +2.5 (still large!)
- p-value might be p < 0.01 instead of p < 0.000001 (still significant!)
- **Core finding likely stands**, but with **reduced magnitude**

---

## Recommended Actions

### **Immediate (Before Merge)**:
1. ❌ **DO NOT MERGE PR** until statistics are corrected
2. ⚠️ Add warning to PR description about statistical concerns
3. 📝 Re-run analysis with Option 1 (all individual scores)
4. 📊 Report corrected Cohen's d and p-value
5. 🔄 Update all documentation with corrected values

### **Short-Term (Next 24-48h)**:
1. Implement Option 2 (mixed-effects model) for robustness
2. Implement Option 3 (permutation test) for non-parametric confirmation
3. Create sensitivity analysis varying:
   - Number of arbitrary encoders (5, 10, 20)
   - Number of sequences (50, 100, 200)
4. Add statistical assumptions checking:
   - Shapiro-Wilk test for normality
   - Levene's test for equal variance
   - Q-Q plots

### **Long-Term (Before Publication)**:
1. Consult with statistician or bioinformatician
2. Run on real CRISPR data (not just synthetic)
3. Perform bootstrap confidence intervals
4. Add power analysis

---

## Corrected Code Template

```python
def run_comparative_test_corrected(self, n_sequences: int = 100,
                                  n_arbitrary_trials: int = 10) -> Dict:
    """Run breathing dynamics vs arbitrary encoding test (CORRECTED)"""

    sequences = self.generate_crispr_sequences(n_sequences)

    # Breathing encoding (individual scores)
    breathing_scores_gc = []
    for seq in sequences:
        z_gc = self.analyze_mutation_sensitivity(seq, self.breathing_encoder, 'gc_affecting')
        breathing_scores_gc.append(z_gc)

    # Arbitrary encodings (COLLECT ALL INDIVIDUAL SCORES)
    all_arbitrary_scores_gc = []

    for trial in range(n_arbitrary_trials):
        arb_encoder = ArbitraryEncoder(seed=trial)

        for seq in sequences:
            z_gc = self.analyze_mutation_sensitivity(seq, arb_encoder, 'gc_affecting')
            all_arbitrary_scores_gc.append(z_gc)  # Collect individual scores

    # Statistical comparison (CORRECT)
    t_stat, p_value = stats.ttest_ind(breathing_scores_gc, all_arbitrary_scores_gc)

    # Effect size (CORRECT)
    pooled_std = np.sqrt(
        ((len(breathing_scores_gc) - 1) * np.var(breathing_scores_gc, ddof=1) +
         (len(all_arbitrary_scores_gc) - 1) * np.var(all_arbitrary_scores_gc, ddof=1)) /
        (len(breathing_scores_gc) + len(all_arbitrary_scores_gc) - 2)
    )
    cohens_d = (np.mean(breathing_scores_gc) - np.mean(all_arbitrary_scores_gc)) / pooled_std

    print(f"Breathing: n={len(breathing_scores_gc)}, mean={np.mean(breathing_scores_gc):.4f}")
    print(f"Arbitrary: n={len(all_arbitrary_scores_gc)}, mean={np.mean(all_arbitrary_scores_gc):.4f}")
    print(f"Cohen's d (CORRECTED): {cohens_d:+.4f}")
    print(f"p-value (CORRECTED): {p_value:.6f}")

    return {'cohens_d': cohens_d, 'p_value': p_value}
```

---

## Transparency Statement

This statistical flaw was identified **post-PR creation** through peer review. The original analysis was conducted in good faith but contained a methodological error common in comparing aggregated vs. individual data.

**This is exactly why peer review and reproducibility matter.**

The corrected analysis will determine whether:
1. The finding is **real but overstated** (most likely)
2. The finding is **artifact of statistical error** (less likely, given selectivity results)

**Update this document** once corrected analysis is complete.

---

**Status**: ⚠️ **REQUIRES CORRECTION BEFORE MERGE**
**Priority**: **HIGH**
**Action Owner**: PR Author
**Timeline**: Correct within 24-48 hours
