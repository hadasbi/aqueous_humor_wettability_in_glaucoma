# Aqueous Humor Wettability in Glaucoma patients

Analysis code for the paper on eye-surface liquid spreading dynamics, comparing spreading area over time between clinical groups.

## Groups

- **CONTROL** (n = 8) — healthy control eyes
- **GLAUCOMA + URAMOX** (n = 10) — glaucoma, with oral acetazolamide
- **GLAUCOMA − URAMOX** (n = 5) — glaucoma, without oral acetazolamide

For 2-group analyses, both glaucoma groups are pooled into **GLAUCOMA** (n = 15). Cases are identified by number only (`Case#1` … `Case#23`); no patient-identifying information is included.

## Data

The de-identified measurements used by all analysis scripts are provided in:

```
results/eye_liquid_spread_measurements.csv
```

Columns: `case, slide, time_sec, area_cm2, red_pixel_count, velocity_cm2_per_sec, group, condition`.

`01_data_processing.ipynb` is the script that regenerates this CSV from the raw image archives. All other scripts read the CSV directly.

## Scripts and the figures they produce

Run order: `01` first (only if regenerating the CSV), then any analysis script.

| Figure | Script | Description |
|--------|--------|-------------|
| — | `01_data_processing.ipynb` | Extracts red-pixel spread area (calibrated to cm²) from each frame and writes `results/eye_liquid_spread_measurements.csv`. Produces no paper figure. |
| **Figure 1** | `15_spaghetti_plots_R.Rmd` | Raw individual per-case area-vs-time trajectories, faceted by 2 groups (control vs glaucoma). |
| **Figure 2** | `02_control_vs_glaucoma.ipynb` | Box-and-swarm plots of initial, final, and change (Δ) in area — CONTROL vs GLAUCOMA. |
| **Figure 3** | `04a_pointwise_area_control_vs_glaucoma.ipynb` | Pointwise area comparison (CONTROL vs GLAUCOMA), median + IQR over time with per-timepoint Mann–Whitney U tests and Benjamini–Hochberg-adjusted p-values. |
| **Figure 4** | `15_spaghetti_plots_R.Rmd` | Raw individual per-case area-vs-time trajectories, faceted by 3 groups. |
| **Figure 5** | `03_control_vs_glaucoma_uramox.ipynb` | Box-and-swarm plots of initial, final, and Δ area — 3 groups. |
| **Figure 6** | `05a_pointwise_area_control_vs_glaucoma_uramox.ipynb` | Pointwise pairwise area comparison across the 3 groups, median + IQR with BH-adjusted p-values. |
| **Figure 7** | `13_gamm_comparison_R.Rmd` | Population-level GAMM predicted curves, Parallel vs Divergent models — 2-group (Supplementary). |
| **Figure 8** | `13_gamm_comparison_R.Rmd` | Population-level GAMM predicted curves, Parallel vs Divergent models — 3-group (Supplementary). |

*(Figures 7 and 8 are presented in the paper's Supplementary Material but are numbered continuously as main figures. There are 8 figures total; Figures 1–6 are main-text figures.)*

`13_plan_gamm_comparison_R.md` and `13_gamm_comparison_R_parameters.md` document the modeling rationale and parameter choices for `13`.

## Requirements

**Python** (`01`, `02`, `03`, `04a`, `05a`): `numpy`, `pandas`, `matplotlib`, `scipy`, `statsmodels`, `pillow`, `openpyxl`

**R** (`13`, `15`): `tidyverse`, `mgcv`, `gratia`

Render an R script to HTML with:

```powershell
Rscript -e "rmarkdown::render('13_gamm_comparison_R.Rmd')"
```
