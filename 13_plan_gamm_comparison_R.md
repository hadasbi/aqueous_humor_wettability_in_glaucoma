# GAMM Model Comparison: Parallel vs Divergent Group Effects

## Context and Motivation

When comparing liquid spreading trajectories across clinical groups, two
distinct questions arise:

1. **Do groups differ in their baseline level** (vertical shift)?
2. **Do groups differ in their non-linear spreading dynamics** (shape)?

A single model that includes both a global smooth and group-specific smooths
conflates these questions and can suffer from concurvity (the group-specific
smooths may compete with the global smooth for the same variance). This
notebook separates the questions cleanly by fitting two nested models and
formally testing whether the more complex one is justified:

- **Model 1 (Parallel)**: One shared non-linear time trend; groups differ
  only by an intercept shift.
- **Model 2 (Divergent)**: Each group has its own non-linear time trend.

If the Parallel model fits just as well, the groups spread at the same rate
and only differ in baseline area. If the Divergent model is significantly
better, the spreading dynamics themselves differ between groups.

## Data

- Input: [results/eye_liquid_spread_measurements.csv](results/eye_liquid_spread_measurements.csv)
- Columns used: `case`, `time_sec`, `area_cm2`, `group` (CONTROL / GLAUCOMA + URAMOX / GLAUCOMA - URAMOX)
- 23 cases (each a different person), ~933 rows, ~41 frames per case at 0.25s intervals over ~10s
- Two grouping variables are tested:
  - **2-group** (`condition`): CONTROL vs GLAUCOMA (pooled)
  - **3-group** (`group`): CONTROL vs GLAUCOMA − URAMOX vs GLAUCOMA + URAMOX

## Analysis Approach

### Model Specifications

**Model 1 — Parallel** (shared shape, group intercept shift):

```r
area_cm2 ~ group + s(time_sec, k=20, bs="cr") + s(time_sec, case, bs="fs", k=10)
```

- `group`: parametric intercept shift per group
- `s(time_sec)`: one master non-linear trajectory shared by all groups
- `s(time_sec, case, bs="fs")`: per-subject random smooth deviations

**Model 2 — Divergent** (group-specific shapes):

```r
area_cm2 ~ group + s(time_sec, by=group, k=20, bs="cr") + s(time_sec, case, bs="fs", k=10)
```

- `group`: parametric intercept shift per group
- `s(time_sec, by=group)`: separate non-linear trajectory per group
- `s(time_sec, case, bs="fs")`: per-subject random smooth deviations
- No redundant global smooth — each group smooth captures the full shape

### Why These Two Models Are Nested

Model 1 is a special case of Model 2: if all group-specific smooths in Model 2
are constrained to be identical, the result is Model 1. This nesting relationship
is what makes the Likelihood Ratio Test valid.

### Comparison Methods

1. **Likelihood Ratio Test (LRT)** via `anova(m_parallel, m_divergent)`:
   Formally tests whether the extra complexity of group-specific shapes is
   statistically justified. Both Chi-squared and F-test versions are reported;
   the F-test is more conservative and recommended for Gaussian models with
   estimated dispersion.

2. **Akaike Information Criterion (AIC)**: Balances goodness-of-fit against
   complexity. Lower AIC = better. ΔAIC > 2 is meaningful; > 10 is strong
   evidence (Burnham & Anderson, 2002).

Both models are fitted with `method = "ML"` (maximum likelihood) because
likelihood-based comparison criteria (AIC, LRT) require the full likelihood,
not the restricted likelihood from REML.

### Visualisation

- Population-level predicted curves from both models (random effects excluded)
- Group difference curves with 95% CIs computed via the exact linear predictor
  matrix approach (Wood 2017, §7.2.5)

## Key R Packages

- `mgcv`: GAM fitting with ML/REML, factor smooth interactions (`bs="fs"`)
- `gratia`: Modern ggplot2-based visualization of mgcv models
- `tidyverse`: Data wrangling and plotting

## Rmd Section Outline

1. **Setup**: Load packages, read CSV, prepare factors (3-group and 2-group)
2. **Two-Group Comparison** (CONTROL vs GLAUCOMA):
   - Fit Parallel and Divergent models
   - AIC comparison, LRT (Chi-squared + F-test)
   - Predicted curves from both models
   - Group difference curve with CI (Divergent model)
3. **Three-Group Comparison** (CONTROL vs GLAUCOMA−URAMOX vs GLAUCOMA+URAMOX):
   - Fit Parallel and Divergent models
   - AIC comparison, LRT (Chi-squared + F-test)
   - Predicted curves from both models
   - Group difference curves with CI (Divergent model)
4. **Diagnostics**: `gam.check()` for both models
5. **Session Info**
