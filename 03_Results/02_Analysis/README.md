# 03_Results/02_Analysis — DRP1 bulk RNA-seq analysis results

## Overview

This directory is the root output tree for the complete differential-expression and pathway-enrichment analysis of iPSC-derived cortical neurons carrying DRP1 mutations (G32A, GTPase domain; R403C, stalk domain) against isogenic controls at two developmental timepoints (D35, D65). It backs every main-text figure and supplementary table cited in the manuscript:

- **Supplementary Table S1_DE** — per-contrast DE statistics in `DE_results/` (9 contrasts × 21,122 genes).
- **Supplementary Table S2** — `master_gsea_table.csv`, the 12,221-pathway × 9-contrast GSEA universe with trajectory pattern labels.
- **Supplementary Table S3** — `Tables/ribosome_compartment_summary.csv`, three-compartment ribosome Jaccard comparison.
- **Supplementary Table S4** — `Tables/per_database_pattern_summary.csv`, per-collection pattern frequencies.
- **6a sensitivity digest** — `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`, the 81-combination threshold-robustness summary.

The core narrative — the synaptic-ribosome translation paradox (ribosome biogenesis UP / synaptic and cytoplasmic translation DOWN) and Compensation as the dominant active pattern in 81/81 threshold combinations — is traceable from files in this tree. GSVA trajectory driver labels use `mutant_driven` / `ctrl_driven` / `both_moving` / `neither_moving` at arm-level |Δ| ≥ 0.10; structural-pool modules are `both_moving`, biogenesis modules are `ctrl_driven`. The full pattern taxonomy is canonicalised in `docs/PATTERN_CLASSIFICATION.md` and implemented in `01_Scripts/Python/pattern_definitions.py`.

## Experimental design

| Factor | Levels | Notes |
|--------|--------|-------|
| Genotype | Ctrl, G32A, R403C | DRP1 mutation status (G32A = GTPase domain, R403C = stalk domain) |
| Timepoint | D35, D65 | Days of cortical-neuron maturation |
| Replicates | 25 samples total, n = 3–6 per group | Unbalanced; smallest cell from `metadata.csv` |

Replicate counts per group are read directly from `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv`

---

## Directory map

Every immediate child of `03_Results/02_Analysis/`. Drill-down READMEs in each subdirectory.

| Path | Purpose | Backs |
|------|---------|-------|
| `DE_results/` | Per-contrast DE statistics — 9 CSVs, 21,122 genes each | Supp Table S1_DE |
| `Tables/` | Aggregate pathway summary tables (ribosome compartments, per-database patterns, shared DEGs) | Supp Tables S3, S4 |
| `Supplementary/` | Reviewer-response digests (5a AveExpr filter, 5b intersections, 6a sensitivity, 6c ribosome compartments) | Reviewer Q5/Q6 panels |
| `Summary/` | `DE_summary.csv` — up/down DEG counts per contrast at FDR 0.05 | Fig 1 DEG bars |
| `Calcium_genes/` | Curated calcium-gene DE table, volcano panels, boxplots, heatmap | Calcium-gene supplementary panel |
| `Sensitivity_Analysis/` | 81-combination threshold-robustness sweep (heatmaps + `claim_stability.csv`) | Supp Fig 6a |
| `Python_exports/` | Intermediate R→Python GSEA exports consumed by `1.5.create_master_pathway_table.py` | None directly (intermediate) |
| `Verification_reports/` | QC reports, ribosome / SynGO / calcium gene-membership audits | Reviewer verification |
| `checkpoints/` | Cached R objects (RDS, ~147 MB) — `model_objects.rds`, `all_gsea_results.rds`, `gsva_all_pathways.rds`, etc. | Reproducibility |
| `Plots/` | All visualisations across 17 subdirectories (see `Plots/README.md`) | Main + supplementary figures |
| (top-level CSVs) | See "Top-level data files" below | See per-file row |

Sibling READMEs to drill into: [`Plots/README.md`](Plots/README.md), [`Tables/README.md`](Tables/README.md), [`DE_results/README.md`](DE_results/README.md), [`Summary/README.md`](Summary/README.md), [`Supplementary/README.md`](Supplementary/README.md), [`Sensitivity_Analysis/README.md`](Sensitivity_Analysis/README.md), [`Calcium_genes/README.md`](Calcium_genes/README.md), [`Verification_reports/README.md`](Verification_reports/README.md), [`Python_exports/README.md`](Python_exports/README.md), [`checkpoints/README.md`](checkpoints/README.md).

Note: the `Ribosome_paradox/` plot directory referenced in earlier revisions has been removed along with `2.1.viz_ribosome_paradox.R`; the corresponding panels are now produced by `2.3.viz_synaptic_ribosomes.R`, `2.4.viz_critical_period_trajectories_gsva.R`, and the publication-figure Python scripts in `02_Analysis/3.x`.

---

## Top-level data files

All files sit directly in `03_Results/02_Analysis/` and are referenced by either the manuscript or downstream visualisation scripts.

| File | Caption | Rows | Generator | Key columns |
|------|---------|------|-----------|-------------|
| `master_gsea_table.csv` | Canonical GSEA universe: every pathway × contrast combination with NES, FDR, and trajectory-pattern labels (Supp Table S2 source) | 109,989 | `02_Analysis/1.5.create_master_pathway_table.py` | `pathway_id`, `database`, `Description`, `contrast`, `category`, `NES`, `pvalue`, `p.adjust`, `setSize`, `Pattern_G32A`, `Pattern_R403C`, `Confidence_*`, `Super_Category_*`, `NES_{Early,TrajDev,Late}_{G32A,R403C}`, `ever_significant`, `ever_significant_trajectory` |
| `master_gsva_all_table.csv` | Comprehensive GSVA scores: every pathway (10–500 gene filter) × genotype × timepoint, with driver labels per arm | 73,326 | `02_Analysis/1.6.gsva_analysis.R` | `pathway_id`, `database`, `pathway_name`, `Genotype`, `Day`, `Mean_GSVA`, `SD_GSVA`, `SE_GSVA`, `Expression_vs_CtrlD35`, `Divergence_vs_Ctrl`, `delta_{Ctrl,G32A,R403C}`, `Driver_G32A`, `Driver_R403C`, `t_statistic`, `p_value`, `p_adjusted`, `significant` |
| `master_gsva_focused_table.csv` | GSVA scores for the 7 curated trajectory modules (long format, per-group rows) | 42 | `02_Analysis/1.7.create_master_gsva_table.R` | `Module`, `Display_Name`, `Panel_ID`, `Source_Database`, `N_genes`, `Genotype`, `Day`, `Mean_GSVA`, `SD_GSVA`, `Expression_vs_CtrlD35`, `Divergence_vs_Ctrl`, `t_statistic`, `p_adjusted`, `significant` |
| `gsva_pattern_summary.csv` | One row per focused module summarising trajectory deviation and pattern classification for both mutations | 7 | `02_Analysis/1.7.create_master_gsva_table.R` | `Module`, `Display_Name`, `Source_Database`, `Expression_vs_CtrlD35_*`, `Divergence_vs_Ctrl_*`, `TrajDev_{G32A,R403C}`, `Pattern_{G32A,R403C}`, `Confidence_*`, `Super_Category_*`, `Change_Consistency` |
| `gsva_statistics_summary.txt` | Plain-text dump of focused-GSVA modules, gene counts, and t-test counts | — | `02_Analysis/1.7.create_master_gsva_table.R` | (free text) |
| `DE_threshold_summary.csv` | DEG counts per contrast at FDR 0.5 / 0.1 / 0.05 / 0.01 | 9 | `02_Analysis/revision/supplements/Supp3.analyze_de_thresholds.R` | `Contrast`, `FDR_0.5`, `FDR_0.1`, `FDR_0.05`, `FDR_0.01` |
| `Top10_genes_by_contrast.csv` | Top 10 DE genes (ranked by t-statistic) for each of the 9 contrasts | 90 | `02_Analysis/revision/supplements/Supp3.analyze_de_thresholds.R` | `Contrast`, `Gene`, `logFC`, `AveExpr`, `P.Value`, `adj.P.Val`, `B` |
| `replicate_level_gsva_long.csv` | Per-sample GSVA scores in long format for every pathway in the 5,267-universe (Supp Fig 10 source) | 305,525 | `02_Analysis/X.export_gsva_long.R` | `pathway_id`, `sample_id`, `genotype`, `day`, `gsva_score` |

---

## Generating-script map

All scripts live under `02_Analysis/`. The 1.x stage runs the heavy compute; 2.x/3.x are visualisation only (no recomputation).

| Output category | Script | Notes |
|-----------------|--------|-------|
| DE statistics (`DE_results/*.csv`, `checkpoints/{de_results,model_objects,fit_object,voom_object,DGE_object}.rds`, `Plots/Volcano/*`, `Plots/GSEA/*`) | `1.1.main_pipeline.R` | Master entrypoint; runs limma-voom + edgeR + MSigDB GSEA; 30–60 min cold, instant with checkpoints |
| Per-contrast DE tables (regeneration only) | `1.2.generate_contrast_tables.R` | Reads `de_results.rds`; writes CSVs into `DE_results/` |
| MitoCarta GSEA | `1.3.add_mitocarta.R` | Augments `checkpoints/mitocarta_gsea_results.rds`; <5 min |
| R→Python GSEA export (`Python_exports/`) | `1.4.export_gsea_for_python.R` | Feeds `1.5` |
| `master_gsea_table.csv` | `1.5.create_master_pathway_table.py` | Joins GSEA stats with trajectory pattern labels from `01_Scripts/Python/pattern_definitions.py` |
| `master_gsva_all_table.csv` | `1.6.gsva_analysis.R` | 15–40 min; caches to `checkpoints/gsva_all_pathways.rds` |
| `master_gsva_focused_table.csv`, `gsva_pattern_summary.csv`, `gsva_statistics_summary.txt` | `1.7.create_master_gsva_table.R` | <1 min |
| SynGO ribosome gene extracts (`Verification_reports/syngo_*`) | `1.8.extract_syngo_ribosome_genes.R` | <1 min |
| Mito-translation cascade plots (`Plots/Mito_translation_cascade/`) | `2.2.viz_mito_translation_cascade.R` | Reads checkpoints; no recompute |
| Synaptic-ribosome panel (`Plots/Synaptic_ribosomes/`) | `2.3.viz_synaptic_ribosomes.R` | Reads `master_gsva_focused_table.csv` + GSEA checkpoints |
| Critical-period trajectory panels (`Plots/Critical_period_trajectories/`) | `2.4.viz_critical_period_trajectories_gsva.R` | GSVA-driven; uses driver labels |
| Complex V (`Plots/Complex_V_analysis/`) | `2.5.viz_complex_v_analysis.R` | OXPHOS sub-panel |
| Calcium-gene panels (`Calcium_genes/`) | `2.6.viz_calcium_genes.R` | Volcanoes, boxplots, heatmap |
| Publication figures (heatmaps, dotplots, upset) | `3.1.publication_figures.py`, `3.2.publication_figures_dotplot.py`, `3.3.ribosome_upset_plot.py`, `3.4.pattern_summary_normalized.py` | All Python; consume `master_gsea_table.csv` |
| Trajectory flow / alluvial / bump / chord | `3.5.viz_trajectory_flow.py`, `3.6.viz_alluvial_ggalluvial.R`, `3.7.viz_bump_chart.py`, `3.7.viz_chord_diagrams.py`, `3.8.viz_interactive_bump*.py`, `3.9.viz_pooled_dotplots.R` | See `Plots/README.md` |
| SynGO running-sum normalised | `3.10.viz_syngo_running_sum_normalised.R` | Reads `syngo_gsea_results.rds` |
| `replicate_level_gsva_long.csv` | `X.export_gsva_long.R` | Standalone exporter; reads `master_gsea_table.csv` for pathway universe + `gsva_all_pathways.rds` |
| `DE_threshold_summary.csv`, `Top10_genes_by_contrast.csv` | `revision/supplements/Supp3.analyze_de_thresholds.R` | Reviewer-response digest |
| Sensitivity sweep (`Sensitivity_Analysis/`, `Supplementary/6a_*`) | `revision/supplements/6a.*.py`, `Supp4.sensitivity_analysis.py` | 81-combination FDR × |NES| × min-set-size grid |
| Supplements 5a/5b/6b/6c/7–10 | `revision/supplements/Supp{1..10}.*`, `5a.*`, `5b.*`, `6b_*`, `6c.*` | Reviewer-specific panels |
| Interactive bump dashboard | `01_Scripts/Python/bump_dashboard/` + `3.8.viz_interactive_bump_dashboard.py` | Modular Python package; see `01_Scripts/Python/README.md` |

`SCRIPTS.md` inside `02_Analysis/` carries the unabridged inventory.

---

## Statistical methods

### Differential expression (`1.1.main_pipeline.R`)

| Parameter | Value | Usage |
|-----------|-------|-------|
| Normalisation | TMM (edgeR `calcNormFactors`) | Library-size correction |
| Filtering | `edgeR::filterByExpr()` | Replaces manual count cutoffs |
| Transformation | voom log2-CPM with precision weights | `limma::voomLmFit` |
| Empirical Bayes | `eBayes(robust = TRUE)` | Moderated t |
| Multiple testing | Benjamini–Hochberg | Within-contrast |
| DEG cutoff | FDR < 0.05 (no FC threshold) | Bar charts, UpSet plots |
| Volcano colouring (fdr mode) | FDR ≤ 0.1 | Highlight only |
| Volcano colouring (p mode) | raw p ≤ 0.05, |log2FC| ≥ 2 | Alternate display |

### Contrasts (9 total)

Mutation effects (4):
- `G32A_vs_Ctrl_D35`, `R403C_vs_Ctrl_D35`, `G32A_vs_Ctrl_D65`, `R403C_vs_Ctrl_D65`.

Maturation effects (3):
- `Time_Ctrl` (D65 − D35 within control), `Time_G32A`, `Time_R403C`.

Maturation-specific interactions (2, difference-in-difference):
- `Maturation_G32A_specific` = (D65_G32A − D35_G32A) − (D65_Ctrl − D35_Ctrl)
- `Maturation_R403C_specific` = (D65_R403C − D35_R403C) − (D65_Ctrl − D35_Ctrl)

### GSEA (`1.1.main_pipeline.R`, `1.3.add_mitocarta.R`)

fgsea v1.26+ via clusterProfiler, ranked by limma t-statistic, 10,000 permutations, BH FDR, HGNC→Entrez via `org.Hs.eg.db`. Twelve databases: Hallmark (50), KEGG (186), Reactome (1,615), GO:BP (7,658), GO:CC (1,006), GO:MF (1,738), WikiPathways (664), Canonical (2,922), CGP (3,358), TF (1,137), SynGO v1.1 (~300), MitoCarta 3.0 (149).

### GSVA (`1.6.gsva_analysis.R`, `1.7.create_master_gsva_table.R`)

GSVA v1.48+ Gaussian kernel; single-sample enrichment scores feeding per-genotype × per-timepoint summaries. Pathway-size filter 10–500 genes. Driver labels are assigned per arm with |Δ| ≥ 0.10 threshold (see `1.6.gsva_analysis.R` and `01_Scripts/Python/patterns.py`). The pathway universe is intentionally aligned to `master_gsea_table.csv` so GSVA-derived trajectories can be cross-referenced row-for-row.

### Trajectory pattern taxonomy

Eight mutually exclusive patterns over Early → TrajDev → Late dynamics. Active patterns require significant trajectory deviation; passive patterns are developmental buffering without active compensation.

| Pattern | Class | Criterion |
|---------|-------|-----------|
| Compensation | Active | Early defect + TrajDev opposes + Late improved |
| Sign_reversal | Active | Early defect + TrajDev opposes + Late opposite sign |
| Progressive | Active | Early defect + TrajDev amplifies + Late worsened |
| Natural_improvement | Passive | Early defect + no TrajDev + Late improved |
| Natural_worsening | Passive | Early defect + no TrajDev + Late worsened |
| Late_onset | — | No early defect, late defect emerges |
| Transient | — | Early defect resolved by late |
| Complex | — | Non-linear or multiphasic |

Super-categories used in the main text: `Active_Compensation`, `Active_Reversal`, `Active_Progression`, `Passive`, `Late_onset`, `Other`. Canonical spec: `docs/PATTERN_CLASSIFICATION.md`; implementation: `01_Scripts/Python/pattern_definitions.py`.

---

## How to regenerate

The pipeline is checkpoint-cached: a clean re-run from raw counts takes ~30–60 min for DE+GSEA plus 15–40 min for the comprehensive GSVA pass; warm re-runs are seconds-minutes. To force a full recompute, set `config$force_recompute <- TRUE` at the top of `1.1.main_pipeline.R`.

```bash
# Phase 1 — core analysis (DE + MSigDB GSEA; produces all checkpoints)
Rscript 02_Analysis/1.1.main_pipeline.R

# Phase 1 supplements (MitoCarta GSEA, R->Python export)
Rscript 02_Analysis/1.3.add_mitocarta.R
Rscript 02_Analysis/1.4.export_gsea_for_python.R

# Master tables
python3 02_Analysis/1.5.create_master_pathway_table.py     # master_gsea_table.csv
Rscript 02_Analysis/1.6.gsva_analysis.R                    # master_gsva_all_table.csv (slow)
Rscript 02_Analysis/1.7.create_master_gsva_table.R         # focused GSVA + pattern summary
Rscript 02_Analysis/X.export_gsva_long.R                   # replicate_level_gsva_long.csv

# Reviewer-response digests (DE thresholds, sensitivity, ribosome compartments)
Rscript 02_Analysis/revision/supplements/Supp3.analyze_de_thresholds.R
python3 02_Analysis/revision/supplements/Supp4.sensitivity_analysis.py
python3 02_Analysis/revision/supplements/6a.per_database_pattern_dual_denominator.py
python3 02_Analysis/revision/supplements/6c.ribosome_compartment_summary.py

# Visualisation (no recompute; reads checkpoints + master CSVs)
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.2.publication_figures_dotplot.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.7.viz_chord_diagrams.py
```

Caching pattern used throughout:

```r
result <- load_or_compute(
  checkpoint_file = "03_Results/02_Analysis/checkpoints/<name>.rds",
  compute_fn      = function() { <expensive_computation>() },
  force_recompute = config$force_recompute,
  description     = "<what this computes>"
)
```

---

## Related documentation

- `AGENTS.md` (repo root) — agent + reproducibility rules.
- `docs/PATTERN_CLASSIFICATION.md` — canonical pattern taxonomy spec.
- `docs/DATA_LINEAGE.md` — input → checkpoint → output flow.
- `02_Analysis/SCRIPTS.md` — full script inventory.
- `Plots/README.md` — figure-by-figure documentation.
- `DE_results/README.md`, `Tables/README.md`, `Summary/README.md`, `Supplementary/README.md`, `Sensitivity_Analysis/README.md`, `Calcium_genes/README.md`, `Verification_reports/README.md`, `Python_exports/README.md`, `checkpoints/README.md` — sibling drill-downs.
