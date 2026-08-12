# CAN-MD Database: PCA Analysis and Factor Diagnostics

A partial replication and extension of **Fortin-Gagnon, Leroux, Stevanovic, and Surprenant (2022)**, "A Large Canadian Database for Macroeconomic Analysis", *Canadian Journal of Economics*, 55(4): 1799–1833. <https://doi.org/10.1111/caje.12618>

> **Status:** Work in progress. Figures 2, 3, 4 and Table 1 complete. Recession prediction models (Figures 6, 7, Table 2) in progress.

------------------------------------------------------------------------

## Overview

This repository replicates the factor analysis component of the CAN-MD paper using MATLAB. The CAN-MD database is a large monthly macroeconomic dataset for Canada covering 1981M01–2019M12 with 387 series spanning real activity, labour markets, trade, financial conditions, prices, and housing — both at the national and provincial levels.

The analysis proceeds in four scripts, each handling a distinct stage of the pipeline:

| Script | Dataset | Purpose | Outputs |
|----|----|----|----|
| `a_canmad_pca.m` | `TR_CAN_MD.csv` | Data prep, cleaning, standardisation, PCA | `pca_data.mat` |
| `b_canmad_viz.m` | `pca_data.mat` | Factor visualisation | Figures 2, 4 |
| `c_canmad_factor_select.m` | `balanced_can_md.csv` | Factor selection criteria, recursive estimation | Table 1, Figure 3 |
| `d_canmad_recession.m` | TBD | Recession prediction (probit, lasso) | Figures 6, 7, Table 2 |

------------------------------------------------------------------------

## Data

Two datasets from the CAN-MD replication package are used, reflecting the paper's own design:

**`TR_CAN_MD.csv`** — the stationarity-transformed raw panel (March 2020 vintage). Used for scripts `a_` and `b_` to examine individual series and their factor loadings. Contains missing values for some series (handled by dropping or trimming).

**`balanced_can_md.csv`** — the EM-algorithm balanced panel (January 2020 vintage), starting 1981M01 with no missing values. Used for script `c_` for factor selection criteria, which require a complete matrix.

Both files are from:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "Replication Data and Code for: A Large Canadian Database for Macroeconomic Analysis", Borealis, V1. <https://doi.org/10.5683/SP3/VKTBBO>

Data are publicly available at the above link and are not redistributed in this repository.

------------------------------------------------------------------------

## Sample and Variable Selection

- **Estimation window:** 1992M01–2019M12 (scripts `a_`, `b_`); 1981M01–2019M12 (script `c_`)
- **Series:** 115 aggregate Canadian series retained after dropping 271 provincial series (identified via province-code suffixes: `_NF`, `_PEI`, `_NS`, `_NB`, `_QC`, `_ONT`, `_MAN`, `_SAS`, `_ALB`, `_BC`)
- **`FIN_new` dropped** in script `a_` due to 62 missing values at the start of the 1992 window; retained in `c_` via the EM-balanced panel

------------------------------------------------------------------------

## Results

### Table 1 — Number of factors estimated (Canada panel)

| Criterion                 | Estimate |
|---------------------------|----------|
| BN02 (Bai-Ng 2002)        | 6        |
| ABC (Alessi et al. 2010)  | 6        |
| ON (Onatski 2010)         | 0        |
| AH (Ahn-Horenstein 2013)  | 2        |
| HL (Hallin-Liska 2007)    | 4        |
| BN07 (Bai-Ng 2007)        | 3        |
| AW (Amengual-Watson 2007) | 4        |

All seven estimates match the paper's reported values exactly.

### Figure 2 — Eigenvalues and cumulative variance explained

Scree plot and trace for the first 20 principal components. First component explains \~9.9% of variance; first 10 components cumulate to \~45.7%. Slight difference from paper due to longer sample window (1992 vs 1981 start).

### Figure 3 — Recursive factor selection (1991M01–2019M12)

Expanding window estimates of factor count by BN02 and HL criteria. BN02 (solid black) shows a stable step function rising from 4 to 6 over time. HL (blue dashed) fluctuates between 2 and 4, with upward drift post-2008 — consistent with the paper's finding that the number of detectable factors increases as sample size grows.

### Figure 4 — Factor time series and highest-loading series

Four-panel plot of the first four factor scores (black, left axis) alongside the series with the highest absolute loading on each factor (red, right axis). Factors identified as: (1) real activity — GDP, (2) external sector — USD/CAD exchange rate, (3) monetary conditions — 3-month commercial paper rate, (4) financial conditions — 5–10Y government bond spread. Minor differences from paper's factor labelling reflect sample window and exclusion of `FIN_new`.

------------------------------------------------------------------------

## Repository Structure

```         
canmd-gagnon-et-al-2022-pca/
├── Code/
│   ├── a_canmad_pca.m
│   ├── b_canmad_viz.m
│   ├── c_canmad_factor_select.m
│   └── d_canmad_recession.m          ← in progress
├── Input/
│   ├── TR_CAN_MD.csv                 ← not redistributed; download from Borealis
│   └── balanced_can_md.csv           ← not redistributed; download from Borealis
├── Output/
│   ├── Figure2_scree.pdf
│   ├── Figure3_recursive.pdf
│   ├── Figure4_factors.pdf
│   └── Table1_factors.csv
├── Tools/
│   └── EstimationNumberFactors/      ← toolbox from replication package
├── LICENSE
└── README.md
```

------------------------------------------------------------------------

## Dependencies

**MATLAB** (tested on R2020b+) with the following toolboxes: - Statistics and Machine Learning Toolbox (`pca`, `zscore`)

**Estimation toolbox** from the CAN-MD replication package (`Tools/EstimationNumberFactors/`), which provides `select_numfac`, `nbplog`, `HLestimate`, `Ahn_Horenstein`, and associated Gauss utility functions (`cols`, `rows`, `lagn`, `trimr`, `standard`). Download the full replication package from <https://doi.org/10.5683/SP3/VKTBBO> and copy the `40_tools/Datasets/` contents into `Tools/EstimationNumberFactors/`.

------------------------------------------------------------------------

## Notes on Replication

This is a **practice replication** built for learning purposes, not a certified exact replication. Key differences from the paper:

- Scripts `a_` and `b_` use the 1992M01 sample start (post inflation-targeting adoption) rather than the paper's 1981M01, since the focus is on the modern monetary policy regime
- `FIN_new` is excluded from the PCA in `a_`/`b_` due to coverage gaps in the 1992-start window; it is retained in `c_` via the EM-balanced panel
- The paper's Figure 4 is produced in R (`update_heatmaps.R`); our Figure 4 is an independent MATLAB implementation

All code is written from scratch with genuine understanding. The `Tools/` folder contains third-party functions from the authors' replication package used under CC BY-NC 4.0.

------------------------------------------------------------------------

## Citation

If you use this repository, please cite the original paper:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "A Large Canadian Database for Macroeconomic Analysis", *Canadian Journal of Economics*, 55(4): 1799–1833. <https://doi.org/10.1111/caje.12618>

And the replication data:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "Replication Data and Code for: A Large Canadian Database for Macroeconomic Analysis", Borealis, V1. <https://doi.org/10.5683/SP3/VKTBBO>

------------------------------------------------------------------------

## Author

Amal Varghese \| MA Economics, Yonsei University \| [github.com/amalzmv-013](https://github.com/amalzmv-013)
