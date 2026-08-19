# CAN-MD Database: PCA Analysis and Factor Diagnostics

A partial replication of **Fortin-Gagnon, Leroux, Stevanovic, and Surprenant (2022)**, "A Large Canadian Database for Macroeconomic Analysis", *Canadian Journal of Economics*, 55(4): 1799–1833. https://doi.org/10.1111/caje.12618

> **Scope:** This repository replicates the factor analysis component of the paper (Figures 2, 3, 4 and Table 1). The recession prediction exercise (Figures 6, 7, Table 2) is deferred — the paper's implementation involves BIC-optimal factor combination search, EM-algorithm factor extraction, and custom pseudo-R² computation for lasso models that go beyond the scope of this exercise. The primary purpose of this repository is to demonstrate MATLAB-based PCA analysis of a large macroeconomic panel.

---

## Overview

The CAN-MD database is a large monthly macroeconomic dataset for Canada covering 1981M01–2019M12 with 387 series spanning real activity, labour markets, trade, financial conditions, prices, and housing — at both the national and provincial levels.

The analysis proceeds in three scripts, each handling a distinct stage of the pipeline:

| Script | Dataset | Purpose | Outputs |
|---|---|---|---|
| `a_canmad_pca.m` | `TR_CAN_MD.csv` | Data prep, cleaning, standardisation, PCA | `pca_data.mat` |
| `b_canmad_viz.m` | `pca_data.mat` | Factor visualisation | Figures 2, 4 |
| `c_canmad_factor_select.m` | `balanced_can_md.csv` | Factor selection criteria, recursive estimation | Table 1, Figure 3 |

---

## Data

Two datasets from the CAN-MD replication package are used, reflecting the paper's own design:

**`TR_CAN_MD.csv`** — the stationarity-transformed raw panel (March 2020 vintage). Used for scripts `a_` and `b_` to examine individual series and their factor loadings.

**`balanced_can_md.csv`** — the EM-algorithm balanced panel (January 2020 vintage), 1981M01–2019M12 with no missing values. Used for script `c_` for factor selection criteria, which require a complete matrix.

Both files are publicly available from:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "Replication Data and Code for: A Large Canadian Database for Macroeconomic Analysis", Borealis, V1. https://doi.org/10.5683/SP3/VKTBBO

Data are not redistributed in this repository. Download from the above link and place in `Input/`.

---

## Sample and Variable Selection

**Scripts `a_` and `b_`** use 1992M01–2019M12 (336 months). The 1992M01 start reflects the post inflation-targeting adoption period — the Bank of Canada formally adopted inflation targeting in February 1991, with 1992M01 being the conventional clean start date in Canadian macro empirical work. This is a deliberate deviation from the paper's 1981M01 start, motivated by interest in the modern monetary policy regime.

**Script `c_`** uses the full 1981M01–2019M12 sample (468 months) from the balanced panel — consistent with the paper's factor selection exercise.

**Series selection:** 115 aggregate Canadian series retained after dropping 271 provincial series identified via province-code suffixes (`_NF`, `_PEI`, `_NS`, `_NB`, `_QC`, `_ONT`, `_MAN`, `_SAS`, `_ALB`, `_BC`). `FIN_new` dropped in scripts `a_`/`b_` due to 62 missing values at the 1992 start; retained in `c_` via the EM-balanced panel.

---

## Results

### Table 1 — Number of factors (Canada panel)

| Criterion | Paper | This replication |
|---|---|---|
| BN02 (Bai-Ng 2002) | 6 | 6 ✓ |
| ABC (Alessi et al. 2010) | 6 | 6 ✓ |
| ON (Onatski 2010) | 0 | 0 ✓ |
| AH (Ahn-Horenstein 2013) | 2 | 2 ✓ |
| HL (Hallin-Liska 2007) | 4 | 4 ✓ |
| BN07 (Bai-Ng 2007) | 3 | 3 ✓ |
| AW (Amengual-Watson 2007) | 4 | 4 ✓ |

Perfect 7/7 match on the Canada-only panel.

### Figure 2 — Eigenvalues and cumulative variance explained
Scree plot and trace for the first 20 principal components. First component explains ~9.9% of variance; first 10 components cumulate to ~45.7%. Slightly lower than the paper due to the shorter 1992-start window — fewer observations capture less of the low-frequency variance driven by the 1981–82 and 1990–92 recessions.

### Figure 3 — Recursive factor selection (1991M01–2019M12)
Expanding window estimates of factor count by BN02 (solid black) and HL (blue dashed) criteria. BN02 shows a stable step function rising from 4 to 6 over the sample period. HL fluctuates between 2 and 4, with upward drift post-2008 — consistent with the paper's finding that the number of detectable factors grows with sample size. Minor level differences from the paper reflect the different `TR_CAN_MD.csv` vintage and sample start.

### Figure 4 — Factor time series and highest-loading series
Four-panel dual-axis plot of the first four factor scores (black) alongside the series with the highest absolute loading on each factor (red). Factors identified as: (1) real activity — GDP, (2) external sector — USD/CAD exchange rate, (3) monetary conditions — 3-month commercial paper rate, (4) financial conditions — 5–10Y government bond spread. F2 and F3 differ from the paper's labelling due to the exclusion of `FIN_new` and the shorter 1992-start window, which reshuffles which series loads highest on each factor.

---

## Notes on Deviations from the Paper

| Aspect | Paper | This replication | Reason |
|---|---|---|---|
| Sample start | 1981M01 | 1992M01 (a_, b_) | Post inflation-targeting regime focus |
| Factor extraction | EM algorithm | Standard PCA | Pedagogical simplicity; EM-balanced panel used for c_ |
| `FIN_new` | Retained | Dropped (a_, b_) | Coverage gap in 1992-start window |
| Recession prediction | Full pipeline | Deferred | Complexity of BIC-optimal factor search and lasso pseudo-R² |
| Figure 4 | R (`update_heatmaps.R`) | MATLAB (b_) | Independent implementation |

---

## Repository Structure

```
canmd-gagnon-et-al-2022-pca/
├── Code/
│   ├── a_canmad_pca.m                ← data prep and PCA
│   ├── b_canmad_viz.m                ← Figures 2 and 4
│   └── c_canmad_factor_select.m      ← Table 1 and Figure 3
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
├── canmd_analysis.qmd                ← Quarto analysis document
├── LICENSE
└── README.md
```

---

## Dependencies

**MATLAB** (tested on R2020b+):
- Statistics and Machine Learning Toolbox (`pca`, `zscore`)
- `EstimationNumberFactors` toolbox from the CAN-MD replication package — provides `select_numfac`, `nbplog`, `HLestimate`, `Ahn_Horenstein`, and Gauss utility functions (`cols`, `rows`, `lagn`, `trimr`, `standard`). Download from https://doi.org/10.5683/SP3/VKTBBO and copy `40_tools/Datasets/` contents into `Tools/EstimationNumberFactors/`.

**Known toolbox issue:** `HLestimate` internally calls `rand('state', 0)`, switching MATLAB to a legacy random number generator. In the recursive estimation loop (Figure 3), `HLestimate` must run *before* `nbplog` within each iteration to avoid corruption of the BN02 estimates. See `c_canmad_factor_select.m` comments for details.

---

## Citation

If referencing this repository, please cite the original paper:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "A Large Canadian Database for Macroeconomic Analysis", *Canadian Journal of Economics*, 55(4): 1799–1833. https://doi.org/10.1111/caje.12618

And the replication data:

> Fortin-Gagnon, O., Leroux, M., Stevanovic, D. and Surprenant, S. (2022), "Replication Data and Code for: A Large Canadian Database for Macroeconomic Analysis", Borealis, V1. https://doi.org/10.5683/SP3/VKTBBO

---

## Author

Amal Varghese | MA Economics, Yonsei University | [github.com/amalzmv-013](https://github.com/amalzmv-013)
