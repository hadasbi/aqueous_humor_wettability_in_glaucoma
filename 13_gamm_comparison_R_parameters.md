# Parameter Choices in `13_gamm_comparison_R.Rmd`

Comprehensive catalogue of every analytical parameter, the rationale for the
current choice, and what alternatives could have been used.

---

## Model Structure

### 1. Time smooth basis dimension — `k = 20`

- **Why chosen**: Each case contributes ~41 time points over ~10 seconds.
  `k = 20` provides sufficient flexibility to capture the non-linear growth
  trajectory (rapid initial rise followed by plateau) while ML penalisation
  shrinks away unused complexity. Setting k to roughly half the number of
  unique time points is a common rule of thumb that balances resolution
  against overfitting.
- **Alternatives**:
  - `k = 15`: Slightly more conservative; adequate for smoother curves.
  - `k = 30`: Higher flexibility; safe with penalisation but slower.
- **Diagnostic**: `gam.check()` tests whether k is too low (significant
  p-value on the basis dimension test indicates more flexibility is needed).

### 2. Subject random smooths — `s(time_sec, case, bs = "fs", k = 10)`

- **Why chosen**: `bs = "fs"` (factor smooth interaction) gives each subject
  its own smooth deviation from the population curve, shrunk toward zero by
  a shared smoothing penalty. This accounts for the repeated-measures design
  (multiple time points per subject) and captures both random intercepts and
  random non-linear trajectory differences. `k = 10` is deliberately lower
  than the population smooth (`k = 20`) because individual deviations from
  the group mean are expected to be simpler than the overall trajectory.
  With 23 subjects, this produces 23 × 10 = 230 random-effect basis
  functions — manageable for the dataset.
- **Alternatives**:
  - `k = 15`: More flexible individual deviations at the cost of model rank
    (23 × 15 = 345 basis functions for random effects alone).
  - `k = 5`: Very constrained; may miss genuine inter-subject variability.
  - `bs = "re"`: Random intercepts only (no random slopes or curvature) —
    too restrictive for trajectories that may differ in shape across subjects.

### 3. Spline basis type — `bs = "cr"` (cubic regression splines)

- **Why chosen**: Computationally efficient (tridiagonal penalty matrix),
  natural boundary behaviour (linear extrapolation beyond the data range),
  and well-behaved derivatives under penalty-based smoothing. A sensible
  default for univariate time smooths.
- **Alternatives**: `"tp"` (thin plate), `"ps"` (P-splines), `"bs"`
  (B-splines). Minimal practical difference for univariate smooths.

---

## Estimation and Comparison

### 4. Smoothing parameter estimation — `method = "ML"`

- **Why chosen**: Maximum likelihood is the standard criterion for
  likelihood-based model comparison (AIC and LRT). Both AIC and the
  likelihood ratio test statistic require the full (not restricted)
  log-likelihood. Using the same estimation method for both models
  ensures their likelihoods are directly comparable.
- **Caveat**: ML can slightly underestimate variance components compared
  to REML, which may affect smoothing parameter estimates. In practice
  the difference is negligible for this dataset.
- **Alternatives**:
  - `method = "REML"`: Better for final inference and smoothing parameter
    estimation, but the REML criterion changes when the smooth structure
    changes. Some authors argue REML comparison is valid when parametric
    fixed effects are identical (as in our case); others recommend ML
    for any model comparison. ML is the conservative choice.

### 5. Likelihood Ratio Test — `anova(m1, m2, test = "Chisq")` and `test = "F"`

- **Why chosen**: The LRT is the standard formal test for comparing nested
  models. Model 1 (Parallel) is nested inside Model 2 (Divergent) because
  constraining all group-specific smooths to be identical recovers Model 1.
- **Test statistic**: $D = -2 \ln(L_1 / L_2)$, approximately $\chi^2$
  distributed with degrees of freedom equal to the difference in effective
  model degrees of freedom.
- **Both Chi-squared and F-test are reported**: The Chi-squared test is
  the classical LRT. The F-test is more conservative and is recommended
  for Gaussian family models where the dispersion parameter (residual
  variance) is estimated rather than known. Reporting both provides
  robustness: if they agree, the conclusion is firm; if they disagree,
  the more conservative F-test should be trusted.
- **Caveat**: The $\chi^2$ approximation is approximate for penalised
  smooths because the effective degrees of freedom are not integers. mgcv
  handles this via the approach of Wood (2013).
- **Alternative**:
  - **Parametric bootstrap LRT**: Simulate data under Model 1, refit both
    models, build the null distribution of D empirically. More accurate
    but computationally expensive.

### 6. AIC interpretation — ΔAIC thresholds

- **Thresholds used**: ΔAIC > 2 = meaningful difference; ΔAIC > 10 =
  strong evidence (Burnham & Anderson, 2002).
- **Why chosen**: Standard thresholds in the ecological and biomedical
  literature. Widely recognised and easy to communicate.
- **Alternatives**:
  - **BIC** (`BIC()`): Penalises complexity more heavily; favours simpler
    models with large n.
  - **AICc** (corrected AIC): Better for small samples. For n ≈ 933 and
    moderate k, the correction is negligible.

---

## Visualisation

### 7. Prediction grid — 200 points from min to max time

- **Why chosen**: Dense enough for smooth curves and accurate
  visualisation without being computationally wasteful.
- **Alternatives**: 100 (coarser) or 500 (unnecessarily fine) points.

### 8. Excluding random effects — `exclude = "s(time_sec,case)"`

- **Why chosen**: Population-level predictions represent the average group
  trajectory, not any specific subject. Excluding the random smooth
  removes subject-specific deviations.
- **Alternative**: Include random effects for individual-level predicted
  curves (useful for diagnostic overlays but not for group comparison).

### 9. Group difference SE — Exact via linear predictor matrix

- **Why chosen**: Uses `Xp %*% Vcov %*% t(Xp)` from Wood (2017, §7.2.5)
  to properly account for shared model components when computing
  `SE(group_A − group_B)`.
- **Alternative**: Independence approximation `sqrt(SE_A² + SE_B²)` which
  overestimates SE (conservative) because the groups share the random
  effects term.

### 10. Random seed — `set.seed(42)`

- **Why chosen**: Ensures reproducibility. The value 42 is arbitrary.
  The main results (AIC, LRT) are deterministic given the data and do not
  depend on the seed, but `gam.check()` simulation-based tests do.

---

## Summary Table

| # | Parameter | Current Value | Key Alternatives |
|---|-----------|---------------|------------------|
| 1 | Time smooth k | `k = 20` | 15, 30 |
| 2 | Random smooth bs/k | `bs = "fs"`, `k = 10` | `"re"`, k = 5 or 15 |
| 3 | Spline basis | `bs = "cr"` | `"tp"`, `"ps"` |
| 4 | Estimation method | `ML` | `REML` |
| 5 | Nested model test | LRT (Chi-sq + F-test) | Parametric bootstrap |
| 6 | AIC thresholds | ΔAIC > 2 / > 10 | BIC, AICc |
| 7 | Prediction grid | 200 points | 100, 500 |
| 8 | Population predictions | Exclude RE | Include RE |
| 9 | Difference SE | Exact (Xp-based) | Independence approx. |
| 10 | Random seed | 42 | Any integer |
