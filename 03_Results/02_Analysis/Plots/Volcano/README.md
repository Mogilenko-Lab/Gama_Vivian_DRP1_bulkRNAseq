# Volcano — Differential Expression Volcano Plots

This directory contains volcano plots and related diagnostic plots for all 9 experimental contrasts, generated in multiple formats and significance thresholds. These are the standard analysis-set volcano plots; the MinAveExpr=0 sensitivity set is in `../Volcano_Supplementary_MinAveExpr0/`.

## Directory structure

```
Volcano/
├── p/                        # Standard horizontal, p-value threshold (p < 0.05)
├── fdr/                      # Standard horizontal, FDR threshold (FDR < 0.10)
├── fdr_raster/               # Rasterized FDR horizontal volcanoes
├── vertical_fdr/             # Vertical layout, FDR threshold
├── vertical_p/               # Vertical layout, p-value threshold
├── vertical_fdr_calcium/     # Vertical, FDR threshold, calcium genes highlighted
├── vertical_p_calcium/       # Vertical, p-value threshold, calcium genes highlighted
├── MD/                       # Mean-Difference (MA) plots
└── FC-B/                     # Fold-Change vs B-statistic plots
```

## File inventory by subdirectory

### p/ and fdr/

9 files each, one per contrast: `{contrast}_standard.pdf`

Contrasts: G32A_vs_Ctrl_D35, G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D35, R403C_vs_Ctrl_D65, Time_Ctrl, Time_G32A, Time_R403C, Maturation_G32A_specific, Maturation_R403C_specific

### vertical_fdr/ and vertical_p/ (15 files each)

9 per-contrast files (`{contrast}_vertical.pdf`) plus 6 composite multi-panel files:

| File | Contents |
|---|---|
| `all_contrasts_vertical.pdf` | All 9 contrasts in one figure |
| `all_disease_vs_control_vertical.pdf` | 4 mutation-vs-control panels |
| `D35_comparisons_vertical.pdf` | G32A and R403C vs Control at D35 |
| `D65_comparisons_vertical.pdf` | G32A and R403C vs Control at D65 |
| `time_effects_vertical.pdf` | Ctrl, G32A, and R403C maturation (3 panels) |
| `maturation_effects_vertical.pdf` | Interaction effects: G32A- and R403C-specific (2 panels) |

### vertical_fdr_calcium/ and vertical_p_calcium/ (15 files each)

Same structure as vertical_fdr/vertical_p/ but with calcium signaling genes highlighted in red and always labeled (NNAT, CACNG3, CACNA1S, ATP2A1, RYR1, MYLK3, VDR, STIM1, STIM2, ORAI1, CALB1, CALR, PNPO).

### MD/ (9 files): `{contrast}_MDplot.pdf`

Mean-Difference (MA) plots: x-axis = average log2 expression (AveExpr), y-axis = log2FC.

### FC-B/ (9 files): `{contrast}_FC_vs_B.pdf`

Fold-Change vs B-statistic plots: x-axis = log2FC, y-axis = B-statistic (log-odds of DE).

## Generating scripts

- `02_Analysis/1.1.main_pipeline.R` — generates all standard volcanic, MD, and FC-B plots
- `01_Scripts/R_scripts/generate_vertical_volcanos.R` — generates vertical and composite multi-panel figures

## Reading guide

**Standard volcano (p/, fdr/)**: X-axis = log2FC; Y-axis = −log10(p-value). Points are colored red (significant upregulated), blue (significant downregulated), or gray (not significant). Horizontal dashed line = significance threshold; vertical dashed lines = fold-change cutoff (|log2FC| = 1). Top genes by significance are labeled.

**Vertical layout**: Same axes as standard but rotated 90° — log2FC on the y-axis, −log10(p) on the x-axis. This layout makes it easier to compare gene labels across multiple contrasts placed side-by-side.

**Significance thresholds**: `p/` variants color genes at raw p < 0.05; `fdr/` variants use FDR < 0.10. Both display −log10(raw p-value) on the y-axis for visual resolution.

**Calcium-highlighted variants**: All calcium genes listed above are always labeled and shown in red regardless of significance. Useful for tracking calcium-related genes (especially NNAT, which is consistently strongly downregulated across all mutation-vs-control contrasts, log2FC ≈ −4).

**MD plots**: Used to assess mean-variance relationship and identify expression-dependent biases. Genes at very low AveExpr are more prone to noise; the AveExpr filter in the main pipeline removes these.

**FC-B plots**: B > 0 means a gene is more likely DE than not-DE. These complement volcano plots by combining effect size and significance into a single score.

## Manuscript figure mapping

| File family | Manuscript role |
|---|---|
| `vertical_fdr/all_contrasts_vertical.pdf`, `vertical_p/all_contrasts_vertical.pdf` | Standard-filter vertical volcano composite (Supp Fig S2E) |
| `vertical_fdr_calcium/*`, `vertical_p_calcium/*` | Calcium-gene highlighted volcanos (extended supplementary) |
| `fdr_raster/*` | Rasterised drop-in replacement for the per-contrast FDR volcanos (production-print only) |
| `MD/*`, `FC-B/*` | Diagnostic only, not in the manuscript |

## Manuscript supplementary caption (Supp Fig S2E, standard-filter set)

**Differential expression landscape across the nine experimental contrasts in DRP1 mutant cortical neurons.** Vertical volcano panels: each subpanel corresponds to one experimental contrast. Y-axis = log2 fold-change (limma-voom coefficient); x-axis = −log10(p-value) from the moderated t-statistic. Points are coloured red (significant, upregulated) or blue (significant, downregulated) at the threshold indicated by the panel filename (`fdr` set = BH-corrected FDR < 0.10; `p` set = raw p < 0.05). The fold-change cutoff |log2FC| = 1 is shown as vertical dashed lines. Up to 20 top genes by significance are labelled per panel. Time-effect contrasts (Time_Ctrl, Time_G32A, Time_R403C) yield several thousand DEGs each as expected for developmental maturation; mutation-vs-control and interaction contrasts (Maturation_G32A_specific, Maturation_R403C_specific) yield smaller, more targeted DEG sets. Standard analysis filter (AveExpr threshold applied upstream); the matched no-filter sensitivity panel is in `../Volcano_Supplementary_MinAveExpr0/`.

---

**Last Updated**: 2026-05-15
**Generating scripts**: `02_Analysis/1.1.main_pipeline.R`, `01_Scripts/R_scripts/generate_vertical_volcanos.R`
