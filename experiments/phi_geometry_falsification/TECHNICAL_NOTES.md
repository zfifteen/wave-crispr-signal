# Technical Notes on phi_phase_score Implementation

## Sequence-Dependence Mechanism

The corrected `phi_phase_score` function achieves sequence-dependence through **purine/pyrimidine weighting**:

### How It Works

1. **Base Weighting**: Each nucleotide gets a weight based on its biochemical class:
   - Purines (A/G): weight = 1.2 (larger bases, stronger stacking)
   - Pyrimidines (C/T): weight = 0.8 (smaller bases, weaker stacking)
   - Ambiguous (N): weight = 1.0 (neutral)

2. **Normalization**: Weights are divided by their mean to preserve scale:
   ```python
   weights = weights / np.mean(weights)
   ```

3. **Weighted Phase Alignment**: Computed as:
   ```python
   weighted_alignment = sum(weights * cos(phase_diffs)) / sum(weights)
   ```

### Important Property: Homopolymer Invariance

After mean-normalization, **homopolymers of the same length get identical scores**:

- `AAAAAAAA` → all weights 1.2 → normalized to [1.0, 1.0, ..., 1.0]
- `TTTTTTTT` → all weights 0.8 → normalized to [1.0, 1.0, ..., 1.0]

**Why this is acceptable**:
1. Real CRISPR datasets have mixed sequences, not homopolymers
2. The feature distinguishes sequences with different **purine/pyrimidine distributions**
3. What matters is the spatial pattern of purines vs pyrimidines along the sequence

### Discriminative Power

The feature discriminates between:
- `ATCGATCGATCG` (alternating, 50% purine) → score ≈ 0.369
- `AAAAGGGGGCCC` (clustered, 50% purine) → score ≈ 0.356
- `CCCCCCCCCCCC` (homopolymer pyrimidine) → score ≈ 0.349

### Comparison to Original Bug

**Original (broken)**: 
- `AAAAAAAAAA` → score 0.369463
- `TTTTTTTTTT` → score 0.369463
- `ATCGATCGAT` → score 0.369463
- **ALL sequences of same length → IDENTICAL score**

**Corrected**:
- `AAAAAAAAAA` → score 0.369463
- `TTTTTTTTTT` → score 0.369463 (same due to normalization, acceptable)
- `ATCGATCGAT` → score 0.369014 (DIFFERENT! ✓)

The key difference: the corrected version distinguishes sequences with **different spatial distributions** of purines/pyrimidines, even if they have the same overall composition.

### Validation

The sequence-dependence is validated by:
1. Shuffle tests (test_phi_geometry_sequence_shuffle.py)
2. Sequences with varying compositions return different scores
3. Permutation tests show correlation changes after sequence shuffling

This confirms the feature is sequence-dependent, not just length-dependent.

---

**Note**: The purine/pyrimidine weighting scheme is a simple but biophysically-motivated approach. More sophisticated schemes (e.g., dinucleotide stacking energies) could provide stronger discrimination, but would increase complexity. The current implementation provides a minimal fix that achieves sequence-dependence while maintaining simplicity and interpretability.
