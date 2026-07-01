# 03_Results/02_Analysis/Sensitivity_Analysis — 81-combination threshold robustness sweep

## Overview

This directory holds outputs from the threshold-sensitivity analysis that validates the 8-pattern trajectory taxonomy across a systematic grid of classifier parameters. The analysis re-runs the `classify_pattern()` function (defined in `01_Scripts/Python/pattern_definitions.py`) over an 81-combination grid formed by varying three threshold families: NES effect size (`NES_EFFECT`: 3 levels), NES strong-signal cutoff (`NES_STRONG`: 3 levels), and improvement/worsening ratio fractions (`IMPROVEMENT_RATIO` × `WORSENING_RATIO`: 3 × 3 = 9 combinations). The key finding is that Compensation is the dominant classifiable pattern in both mutations across all 81 combinations (comp_exceeds_passive 81/81 for G32A, 81/81 for R403C), while Progressive remains rare (0% of classifiable) throughout the grid. The companion digests in `Supplementary/` summarise this grid for direct paper citation. The 12,221-universe version (all GSEA pathways) is kept alongside the manuscript-aligned 5,267-universe version; see `Supplementary/README.md` for a cross-reference table of the two denominators.

## File Inventory

| File | Description |
|------|-------------|
| `sensitivity_results.csv` | Full 81-combination grid (162 data rows: 81 combos × 2 mutations); per-pattern counts and percentages with both `n_pathways` and `comp_of_classifiable_*` columns. 12,221-universe denominator. |
| `claim_stability.csv` | Per-combination stability flags (`comp_exceeds_passive_G32A/R403C`, `R403C_more_compensation`, `progressive_rare_*`, etc.) suitable for counting "how many of 81 combinations support each key claim". |
| `sensitivity_analysis_summary.pdf` | Multi-panel summary figure showing pattern-count stability across the grid. |
| `sensitivity_analysis_summary.png` | Raster version of the same summary figure. |
| `sensitivity_heatmap_G32A.pdf` | Heatmap of Compensation % (of classifiable) across the 81 threshold combinations for G32A. |
| `sensitivity_heatmap_G32A.png` | Raster version. |
| `sensitivity_heatmap_R403C.pdf` | Heatmap of Compensation % (of classifiable) across the 81 threshold combinations for R403C. |
| `sensitivity_heatmap_R403C.png` | Raster version. |

## Generating Script

**`/workspaces/Gama_Vivian_DRP1_bulkRNAseq/02_Analysis/revision/supplements/Supp4.sensitivity_analysis.py`**

This script reads `03_Results/02_Analysis/master_gsea_table.csv`, sweeps the 81-combination threshold grid, writes `sensitivity_results.csv` and `claim_stability.csv`, and produces the PDF/PNG figures. It also writes the 12,221-universe digest to `Supplementary/6a_sensitivity_stability_digest.tsv`. Deterministic; run in 1–3 minutes.

The 5,267-universe digest (paper-aligned; `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`) is produced by the separate script `/workspaces/Gama_Vivian_DRP1_bulkRNAseq/02_Analysis/revision/supplements/6a.sensitivity_sig_universe.py`.

## Column Dictionary

### sensitivity_results.csv

| Column | Description |
|--------|-------------|
| `combination_id` | Human-readable label (e.g., `NES_eff=0.4_strong=0.8_imp=0.6_wors=1.25`) |
| `NES_EFFECT` | Minimum |NES| to count as an "effect" (threshold level 1 of 3) |
| `NES_STRONG` | Minimum |NES| to count as a "strong" enrichment (threshold level 2 of 3) |
| `IMPROVEMENT_RATIO` | Fraction of Early signal remaining at Late that constitutes improvement |
| `WORSENING_RATIO` | Fold-amplification factor at Late that constitutes worsening |
| `mutation` | `G32A` or `R403C` |
| `n_pathways` | Total pathways entering classification (12,221 universe) |
| `Compensation` … `Complex` | Absolute count assigned to each of the 8 patterns |
| `*_pct` columns | Pattern count as % of `n_pathways` |
| `Passive_total` | Sum of Natural_improvement + Natural_worsening |
| `Insufficient_data` | Pathways excluded due to missing trajectory stages |
| `is_default` | `True` for the manuscript default threshold combination |

### claim_stability.csv

| Column | Description |
|--------|-------------|
| `comp_exceeds_passive_G32A/R403C` | Boolean: Compensation count > Passive total for that mutation |
| `comp_dominates_classifiable_G32A/R403C` | Boolean: Compensation is the plurality pattern among classifiable pathways |
| `R403C_more_compensation` | Boolean: R403C Compensation count > G32A Compensation count |
| `progressive_rare_G32A/R403C` | Boolean: Progressive < 1% of classifiable |
| `comp_of_classifiable_G32A/R403C` | Compensation as % of classifiable (excl. Insufficient_data) |

## Interpretation

The stability heatmaps show how Compensation % of classifiable varies as NES thresholds sweep across the grid. A narrow range (e.g., 52–58% for G32A across 81 combinations) indicates the finding is not an artefact of threshold choice. Values in `claim_stability.csv` can be summarised as "claim holds in X of 81 combinations." For the manuscript-aligned numbers (5,267-universe denominator), use `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv` rather than the files in this directory.

## Caption

**Supplementary Table X | Stability of key pathway-classification claims across an 81-combination NES-threshold sweep.** Robustness of the GSEA-based pathway classifier underlying the main-text compensation analysis (and Supplementary Fig. 6a) was assessed by re-running the classifier over a four-dimensional grid of decision thresholds — effect-size cut-off (NES_EFFECT ∈ {0.4, 0.5, 0.6}), strong-effect cut-off (NES_STRONG ∈ {0.8, 1.0, 1.2}), rescue-ratio cut-off (IMPROVEMENT_RATIO ∈ {0.6, 0.7, 0.8}) and worsening-ratio cut-off (WORSENING_RATIO ∈ {1.25, 1.30, 1.40}) — yielding 81 threshold combinations per mutant. Each row corresponds to one combination and reports, separately for the G32A and R403C mutants: (i) the number and percentage of pathways classified as compensatory (*comp_count_\**, *comp_pct_\**); (ii) the number of pathways receiving any non-neutral classification (*classifiable_\**) together with the fraction of that denominator accounted for by compensation (*comp_of_classifiable_\**); and (iii) the percentage of pathways classified as progressive-only-rare (*prog_pct_\**). Five Boolean indicators summarise whether each combination supports the manuscript's headline claims: compensation is the plurality pattern in the classifiable pool of each mutant (*comp_dominates_classifiable_\**); R403C exhibits greater compensation than G32A (*R403C_more_compensation*; raw gap, *compensation_diff*); the progressive pattern remains rare (<1% of classifiable; *progressive_rare_\**); and the compensatory fraction exceeds the passive (untreated) baseline in each mutant (*comp_exceeds_passive_\**). The threshold combination used in the main text is flagged by *is_default = True*. 
