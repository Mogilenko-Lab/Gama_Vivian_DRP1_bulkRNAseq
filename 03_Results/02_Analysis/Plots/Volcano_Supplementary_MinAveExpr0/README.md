# Volcano_Supplementary_MinAveExpr0 — AveExpr > 0 Filtered Volcano Plots (Supp Fig S2E, S2F)

This directory contains volcano plots and AveExpr-significance scatter plots generated with the minimum average expression threshold set to 0 (no expression floor). This set backs **Supp Fig S2E** (vertical volcano panels) and **Supp Fig S2F** (scatter_aveexpr_significance/ panels) of the manuscript.

The key difference from `../Volcano/` is that genes with AveExpr ≤ 0 (very lowly expressed, often below detection) are **not filtered out** before differential expression testing. These plots serve as a sensitivity analysis to demonstrate that the main findings are not driven by the standard expression filter.

## Directory structure

```
Volcano_Supplementary_MinAveExpr0/
├── vertical_fdr/               # Vertical volcanos, FDR threshold, minAveExpr=0
├── vertical_p/                 # Vertical volcanos, p-value threshold, minAveExpr=0
├── vertical_fdr_calcium/       # Same + calcium gene highlighting, FDR
├── vertical_p_calcium/         # Same + calcium gene highlighting, p-value
└── scatter_aveexpr_significance/ # AveExpr vs significance scatter (Supp Fig S2F)
```

## File inventory

### vertical_fdr/, vertical_p/, vertical_fdr_calcium/, vertical_p_calcium/

Each subdirectory contains 15 files (9 per-contrast + 6 composite) following the same convention as `../Volcano/`. All filenames carry the `_minAveExpr0` suffix, e.g., `G32A_vs_Ctrl_D35_vertical_minAveExpr0.pdf/.png`.

Composite files per subdirectory: `all_contrasts_vertical_minAveExpr0.pdf`, `all_disease_vs_control_vertical_minAveExpr0.pdf`, `D35_comparisons_vertical_minAveExpr0.pdf`, `D65_comparisons_vertical_minAveExpr0.pdf`, `time_effects_vertical_minAveExpr0.pdf`, `maturation_effects_vertical_minAveExpr0.pdf`.

### scatter_aveexpr_significance/ (Supp Fig S2F)

| File | Description |
|---|---|
| `aveexpr_scatter_all_mutation_contrasts_2x2.pdf/.png` | 2×2 grid showing AveExpr vs −log10(FDR) for all 4 mutation-vs-control contrasts |
| `aveexpr_scatter_G32A_vs_Ctrl_D35.pdf/.png` | AveExpr vs significance scatter for G32A vs Control at D35 |
| `aveexpr_scatter_G32A_vs_Ctrl_D65.pdf/.png` | AveExpr vs significance scatter for G32A vs Control at D65 |
| `aveexpr_scatter_R403C_vs_Ctrl_D35.pdf/.png` | AveExpr vs significance scatter for R403C vs Control at D35 |
| `aveexpr_scatter_R403C_vs_Ctrl_D65.pdf/.png` | AveExpr vs significance scatter for R403C vs Control at D65 |

## Generating scripts

- `01_Scripts/R_scripts/generate_vertical_volcanos.R` — vertical minAveExpr0 volcano panels
- `01_Scripts/R_scripts/generate_fdr_raster_volcanos.R` — rasterized FDR volcano variants

## Reading guide

**vertical_fdr/ etc.**: Identical interpretation to `../Volcano/vertical_fdr/` but includes genes with AveExpr ≤ 0. Expect more points at high −log10(p) with very small fold-changes — these are typically lowly expressed genes where noise inflates apparent significance. The overall distribution of significant genes in the high-AveExpr range should be nearly identical to the standard filtered set.

**scatter_aveexpr_significance/ (Supp Fig S2F)**: X-axis = AveExpr (mean log2CPM across all samples); Y-axis = −log10(FDR). Color typically encodes upregulated (orange) vs downregulated (blue) vs not significant (gray). The scatter reveals whether significant genes cluster at a specific expression range. The standard analysis filter (AveExpr > threshold) removes the low-expression tail visible in these plots. Genes with AveExpr ≤ 0 should be interpreted cautiously as they are prone to technical noise and disproportionately large fold-changes.

## Manuscript supplementary caption (Supp Fig S2F, `aveexpr_scatter_all_mutation_contrasts_2x2.pdf`)

**Average expression vs differential expression significance, no-AveExpr-filter sensitivity panel.** Two-by-two grid of per-contrast scatter plots for the four mutation-vs-control contrasts (G32A and R403C × D35 and D65). X-axis = AveExpr (mean log2-CPM across all 25 samples after TMM normalisation and voom transformation); Y-axis = −log10(BH-adjusted FDR) from the limma-voom moderated t-test. Each point is one gene; colour distinguishes upregulated (orange) vs downregulated (blue) vs non-significant (grey). No AveExpr filter is applied; the standard analysis floor (applied in `02_Analysis/1.1.main_pipeline.R`) would remove the low-AveExpr tail visible on the left of each panel. The panel demonstrates that significant genes are present across the full AveExpr range and are not artefacts of the filter; conversely, the most extreme −log10(FDR) values in the low-AveExpr tail correspond to genes near or below detection and should not be over-interpreted.

---

**Last Updated**: 2026-05-15
**Generating scripts**: `01_Scripts/R_scripts/generate_vertical_volcanos.R`, `01_Scripts/R_scripts/generate_fdr_raster_volcanos.R`
