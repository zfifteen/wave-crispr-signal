# Z Framework — Copilot Repository Instructions (wave-crispr-signal)

> Purpose: Steer Copilot so it **only** creates experiment designs and code that pass our scientific gates and repo policy.

---

## 0) Repository Policy (MANDATORY)
**Always** follow `.github/REPOSITORY_POLICY.md` for file layout, naming, tests, and version pinning.

- **Dirs**: keep to established structure (`applications/`, `docs/`, `proof_pack/`, `experiments/`, `src/`, etc.).
- **Names**: Python `snake_case.py`, major docs `UPPERCASE.md`, other assets `kebab-case`.
- **Docs**: every public component must have documentation.
- **Tests**: create tests for each core module; runnable in CI.
- **Deps**: exact pins in `requirements.txt`.

---

## 1) Absolute Scientific Gates (no exceptions)
- **Human DNA only**: Use **human nucleotide FASTA** (GRCh38/hg38 or curated CRISPR sets).  
  For **DNA sequences**: Validate **A/C/G/T/N only** (case-insensitive). **Reject** `U` and IUPAC ambiguity codes beyond `N`.  
  For **RNA sequences** (e.g., guide RNA, sgRNA): Validate **A/C/G/U/N only** (case-insensitive). **Reject** `T` and IUPAC ambiguity codes beyond `N`.
- **No fabrication**: Never derive DNA from protein/RNA or fabricate nucleotides.
- **Fail-fast validation**: Start every pipeline with nucleotide-only checks; raise clear `ValueError` on violation.
- **Z invariants (domain-correct)**
    - **Discrete/biological (DEFAULT)**: `Z = A(B / e^2)`; guard divide-by-zero; document A and B.
    <!-- Physical domain constraint removed: not applicable to CRISPR/biological analysis. -->
    <!-- Only the discrete/biological domain is supported. -->
- **Geometric resolution**: θ′(n,k) = φ·((n mod φ)/φ)^k (where φ is the geometric period, e.g., φ = 21 for 21-nt guides; see docs/Z_FRAMEWORK.md, section "Geometric resolution") with default `k ≈ 0.3`. Document any deviation.

---

## 2) Dataset Provenance Gates
- Record **dataset name, version, URL, license, taxonomy**, and **SHA256** of local file(s).
- **Human filter required**: Ensure **Homo sapiens**; print the dataset version at runtime.
- Default dataset allowed: **BioGRID-ORCS Homo sapiens v1.1.17** (with citation in docs).

---

## 3) Statistical Validity & Leakage Gates
- **Pre-register endpoints**:
    - Pearson **r** with **95% bootstrap CI** (≥1,000 resamples).
    - **Partial r** controlling for **GC%**, **guide length**, **guide position**.
    - **Effect size** (Cohen’s d with CI) for defined contrasts.
- **Null model**: ≥1,000× **permutation/shuffle** for empirical p-values.
- **Leakage control**: **split-by-screen** and **split-by-gene**; no entity appears in both train/eval.
- **Multiple comparisons**: apply **Benjamini–Hochberg FDR** when comparing ≥3 metrics.
- **Power/sample size**: include a brief justification (or cite prior power analysis).

---

## 4) Reproducibility & Environment Gates
- **CLI contract** (required flags):  
  `--seed`, `--bootstrap`, `--permutation`, `--splits`, `--domain`, `--pam`
- **Persist metadata** to results: **seed**, **git commit**, **dataset name/version**, **SHA256**, **runtime**, and **pip freeze** as `results/env.txt`.
- **Precision**: use `mpmath` with `mp.dps = 50` where high precision is required.
- **Pinned env**: Python 3.12.*, exact `requirements.txt` pins.

---

## 5) CI & Layout Gates
- Provide a **smoke** dataset + golden outputs so CI completes **<5s**.
- CI runs: `pytest -q` (validators, domain guards, stats wrappers) **and** `make smoke`.
- Directory & artifacts:
    - `experiments/<exp_id>/README.md` (copy-exact commands, timing, endpoints).
    - `experiments/<exp_id>/manifest.yml` (see template below).
    - `results/<exp_id>/run-YYYYMMDD-HHMMSS/` → `results.json`, `results.csv`, `env.txt`, `log.txt`.

---

## 6) Licensing & Citation Gates
- Include proper **citation** and **license** for BioGRID-ORCS (v1.1.17) and any other sources.
- Print dataset **version** at runtime.

---

## 7) Refusal Policy (what Copilot must do if a gate can’t be met)
If any gate cannot be satisfied, respond with:

