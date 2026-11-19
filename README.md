# Cross-cancer heterogeneity in the prognostic effects of immune microenvironment signatures

A pan-cancer random-effects meta-analysis with meta-regression attribution.

## Overview

Immune microenvironment signatures (CD8 T cells, tertiary lymphoid structures, cancer-associated fibroblasts, tumour-associated neutrophils, type-2 conventional dendritic cells, IFN-γ pathway, and proliferation) show strikingly inconsistent prognostic effects across cancer types. Using TCGA PanCancer expression and survival data from 31 cancer types (8,924 primary tumours), we quantified this heterogeneity and attributed it with univariate meta-regression.

**Key findings**

- Between-cancer heterogeneity was marked (I² = 51–86%), including direction reversals (the CD8 T-cell signature is a risk factor in glioblastoma but protective in ovarian cancer).
- Proliferation was the dominant moderator, explaining 27.6–37.5% of the between-cancer variance for CD8, CAF and IFN-γ signatures (permutation p < 0.05), consistent with stoichiometric dilution of immune signal by proliferating malignant cells.
- The tumour-associated neutrophil signature was driven by non-tumoral (stromal/immune) content rather than germ-layer origin.
- Tumour purity and tumour mutational burden contributed little, and the cDC2 signature had no significant moderator.

## Repository contents

- `scripts/` — all analysis code (R and Python). Entry points: `run_specB.R` (signature scoring + per-cancer Cox) and `run_meta_v5.R` (meta-analysis + meta-regression attribution), with sensitivity analyses (`run_sens_*.R`) and figure scripts (`make_figures.R`, `make_fig1.py`).
- `results/` — analysis outputs: `specB/` (main results), `sens/` (sensitivity analyses), and per-cancer / meta tables at the top level.

## Data availability

Expression data (Toil-reprocessed TCGA PanCancer; Vivian et al. 2017, *Nat Biotechnol*) and survival data (Liu et al. 2018, *Cell*) are controlled-access and are **not** included in this repository. Request access via the NCI Genomic Data Commons (GDC) / dbGaP.

## Reproducing the analysis

- **R** 4.3.3 with `metafor` 5.2.1, `survival` 3.5.8, `ggplot2` 4.0.4.
- **Python** 3.x with `pypdf`, `requests` (gene-set / annotation utilities).
- Run `Rscript scripts/run_specB.R`, then `Rscript scripts/run_meta_v5.R`.

## Preregistration

This study was preregistered on the Open Science Framework: https://osf.io/39htc

## Status

Submitted to *npj Precision Oncology*.
