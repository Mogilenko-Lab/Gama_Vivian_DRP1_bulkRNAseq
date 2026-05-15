# 03_Results/02_Analysis/DE_results — per-contrast differential expression tables (Supplementary Table S1_DE)

## Overview

This directory contains gene-level differential expression statistics for all 9 contrasts from the DRP1 mutation study. Each file is a direct extract from the `contrast_tables.rds` checkpoint and represents the primary DE output backing **Supplementary Table S1_DE**. Files are produced by limma-voom (TMM normalization → voom log2-CPM → `voomLmFit` → `eBayes(robust=TRUE)` → Benjamini-Hochberg FDR). All 9 files share 21,122 rows (genes passing `filterByExpr()`).

**Generating scripts:** `02_Analysis/1.1.main_pipeline.R` (full pipeline) or `02_Analysis/1.2.generate_contrast_tables.R` (standalone regeneration from checkpoint).

---

## File Inventory

| File | Contrast | Biological Question |
|------|----------|---------------------|
| `G32A_vs_Ctrl_D35_DE_results.csv` | G32A vs Control at Day 35 | Early G32A mutation effect |
| `G32A_vs_Ctrl_D65_DE_results.csv` | G32A vs Control at Day 65 | Late G32A mutation effect |
| `R403C_vs_Ctrl_D35_DE_results.csv` | R403C vs Control at Day 35 | Early R403C mutation effect |
| `R403C_vs_Ctrl_D65_DE_results.csv` | R403C vs Control at Day 65 | Late R403C mutation effect |
| `Time_Ctrl_DE_results.csv` | D65 − D35 in Controls | Normal maturation trajectory |
| `Time_G32A_DE_results.csv` | D65 − D35 in G32A | G32A maturation trajectory |
| `Time_R403C_DE_results.csv` | D65 − D35 in R403C | R403C maturation trajectory |
| `Maturation_G32A_specific_DE_results.csv` | (D65_G32A − D35_G32A) − (D65_Ctrl − D35_Ctrl) | G32A-specific trajectory deviation (TrajDev) |
| `Maturation_R403C_specific_DE_results.csv` | (D65_R403C − D35_R403C) − (D65_Ctrl − D35_Ctrl) | R403C-specific trajectory deviation (TrajDev) |

All files regenerated from checkpoint `contrast_tables.rds` on 2026-01-13; consistent 21,122-gene universe.

---

## Column Dictionary

| Column | Description | Units / Range |
|--------|-------------|---------------|
| *(first unnamed column)* | HGNC gene symbol (row identifier) | character |
| `logFC` | Log2 fold-change (positive = higher in first condition) | log2 ratio; typically −10 to +10 |
| `AveExpr` | Average log2-CPM expression across all samples | log2-CPM; typically −5 to +15 |
| `t` | Moderated t-statistic from eBayes (used as GSEA ranking metric) | dimensionless; typically −15 to +15 |
| `P.Value` | Raw p-value from moderated t-test | 0–1 |
| `adj.P.Val` | Benjamini-Hochberg FDR-adjusted p-value; primary significance threshold (< 0.05) | 0–1 |
| `B` | B-statistic: log-odds that the gene is differentially expressed; B > 0 means more likely DE than not | dimensionless; typically −10 to +20 |

---

## Statistical Methods

### Differential Expression Pipeline

```
1. TMM normalization (edgeR)
       ↓
2. voom log2-CPM transformation with precision weights
       ↓
3. Linear model fitting (voomLmFit)
       ↓
4. Contrast estimation (contrasts.fit)
       ↓
5. Empirical Bayes moderation (eBayes, robust=TRUE)
       ↓
6. Multiple testing correction (Benjamini-Hochberg)
```

### Key Parameters

| Parameter | Value |
|-----------|-------|
| FDR threshold | 0.05 (BH-adjusted; primary significance cutoff) |
| Volcano FC threshold | \|log2FC\| ≥ 1 (2-fold; visualization only) |
| eBayes moderation | `robust = TRUE` |
| Sort order | t-statistic (descending) |

---

## Contrast Definitions (quick reference)

| Contrast | Formula |
|----------|---------|
| `G32A_vs_Ctrl_D35` | D35_G32A − D35_Control |
| `R403C_vs_Ctrl_D35` | D35_R403C − D35_Control |
| `G32A_vs_Ctrl_D65` | D65_G32A − D65_Control |
| `R403C_vs_Ctrl_D65` | D65_R403C − D65_Control |
| `Time_Ctrl` | D65_Control − D35_Control |
| `Time_G32A` | D65_G32A − D35_G32A |
| `Time_R403C` | D65_R403C − D35_R403C |
| `Maturation_G32A_specific` | (D65_G32A − D35_G32A) − (D65_Control − D35_Control) |
| `Maturation_R403C_specific` | (D65_R403C − D35_R403C) − (D65_Control − D35_Control) |

The two `Maturation_*_specific` contrasts are the **TrajDev** interaction terms used to classify GSEA patterns as Compensation, Sign_reversal, etc.

---

## Source Checkpoints

These CSV files are derived from `checkpoints/contrast_tables.rds` (all 9 contrasts as a single R named list) and `checkpoints/model_objects.rds` (full limma model). Standalone regeneration: `Rscript 02_Analysis/1.2.generate_contrast_tables.R`.
