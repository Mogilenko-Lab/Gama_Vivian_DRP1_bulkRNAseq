# Data Lineage: DRP1 Bulk RNA-seq Analysis Pipeline

**Last Updated**: 2026-05-14
**Purpose**: Documents the complete flow of raw inputs to manuscript figures, tables, and
supplementary artefacts. This file covers lineage edges only; biological interpretation of
patterns lives in `docs/PATTERN_CLASSIFICATION.md`.

---

## 1. Top-Level Orientation

This is an end-to-end data lineage document for a bulk RNA-seq differential expression study
of DRP1 mutations (G32A, R403C) in iPSC-derived cortical neurons at two developmental
timepoints (Day 35 and Day 65). The pipeline is entirely scripted and version-controlled.
Starting data are gene-level count matrices (FASTQ → STAR → featureCounts performed outside
this repo). All downstream analysis steps live under `02_Analysis/` and produce artefacts in
`03_Results/`.

---

## 2. Inputs (`00_Data/` and `03_Results/01_Preprocessing/`)

### Raw count data (preprocessing outputs — not regenerated in this repo)

| Input | Path | Description |
|-------|------|-------------|
| Counts matrix | `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt` | Gene × sample integer count matrix from featureCounts (Subread 2.0.2); 25 samples, GRCh38.p14 |
| Sample metadata | `03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv` | Genotype (Ctrl / G32A / R403C), timepoint (D35 / D65), replicate IDs |

### Reference databases

| Input | Path | Description |
|-------|------|-------------|
| MitoCarta 3.0 GMT | `00_Data/MitoCarta_3.0/MitoPathways3.0.gmx` | Mitochondrial pathway gene sets |
| SynGO bulk (2023-12-01) | `00_Data/SynGO_bulk_20231201/syngo_annotations.xlsx`, `syngo_ontologies.xlsx`, `syngo_genes.xlsx` | Synaptic gene ontology (CC namespace) |
| MSigDB (via msigdbr) | fetched at runtime | Hallmark (H), GO:BP (C5:BP), GO:CC (C5:CC), GO:MF (C5:MF), Reactome (C2:CP:REACTOME), KEGG (C2:CP:KEGG), WikiPathways (C2:CP:WIKIPATHWAYS), CGP (C2:CGP), TF targets (C3:TFT) — 9 collections totalling ~12,000 gene sets together with MitoCarta and SynGO |

---

## 3. Pipeline Stages

### Stage 1 — Preprocessing (upstream of this repo)

```
FASTQ files
  → Trimmomatic v0.39 (adapter/quality trimming)
  → STAR 2.7.11b two-pass alignment → GRCh38.p14
  → featureCounts 2.0.2 (paired-end, non-stranded, gene_id summarisation)
  → sorted_counts_matrix.txt   [03_Results/01_Preprocessing/04_FeatureCounts/]
```

Post-condition: counts matrix and metadata CSV available at the paths in §2.

---

### Stage 2 — Differential Expression (`1.1.main_pipeline.R`)

**Script**: `02_Analysis/1.1.main_pipeline.R`

**Inputs**:
- `sorted_counts_matrix.txt` (counts matrix)
- `metadata.csv`
- `01_Scripts/RNAseq-toolkit/` (git submodule; provides `run_gsea_hsmm()`, `create_standard_volcano()`, and other helpers)

**Processing**:
- `filterByExpr()` retains 21,122 genes (min.count = 10, min.prop = 0.7)
- TMM normalisation via edgeR
- Factorial design matrix (genotype × timepoint), fit with `voomLmFit()` (robust eBayes, no sample weights)
- Nine contrasts: four mutation-vs-control (Early D35, Late D65 for each mutation), three maturation (Time_Ctrl, Time_G32A, Time_R403C), two interaction/TrajDev (Maturation_G32A_specific, Maturation_R403C_specific)
- GSEA across 9 MSigDB collections + SynGO + MitoCarta (12 collections; 100,000 permutations; clusterProfiler/fgsea; gene-set size 10–500; FDR BH)

**Outputs** (all under `03_Results/02_Analysis/`):

| Output | Path | Notes |
|--------|------|-------|
| DE tables (9 contrasts) | `DE_results/<contrast>_DE_results.csv` | limma topTable, all genes |
| Model checkpoint | `checkpoints/model_objects.rds` | fit, voom object, DGE, contrasts |
| Contrast tables checkpoint | `checkpoints/contrast_tables.rds` | topTable results for all 9 contrasts |
| MSigDB GSEA checkpoint | `checkpoints/all_gsea_results.rds` | nested list: contrast → database → gseaResult |
| SynGO GSEA checkpoint | `checkpoints/syngo_gsea_results.rds` | contrast → gseaResult |
| QC variables checkpoint | `checkpoints/qc_variables.rds` | logCPM matrix, sample annotation, ordered_samples |
| Gene intersection checkpoint | `checkpoints/gene_intersections.rds` | baseline9, mat38 shared-DEG lists |
| Standard volcano plots | `Plots/Volcano/p/` and `Plots/Volcano/fdr/` | One PDF per contrast, two threshold modes |
| Vertical volcano plots | `Plots/Volcano/` (vertical) | With calcium gene highlights |
| MD plots | `Plots/Volcano/MD/` | Mean-difference plots |
| FC-vs-B plots | `Plots/Volcano/FC-B/` | log2FC vs B-statistic |
| DEG barplot | `Plots/General/` | Up/down per contrast |
| Sample QC plots | `Plots/General/sample_correlation_heatmap_ordered.pdf`, `PCA_plot.pdf` | |
| GSEA per-contrast plots | `Plots/GSEA/<contrast>/<database>/` | Dotplots, barplots, running-sum plots |

---

### Stage 3 — MitoCarta GSEA (`1.3.add_mitocarta.R`)

**Script**: `02_Analysis/1.3.add_mitocarta.R`

**Inputs**:
- `checkpoints/contrast_tables.rds`
- `00_Data/MitoCarta_3.0/MitoPathways3.0.gmx`

**Outputs**:
- `checkpoints/mitocarta_gsea_results.rds`

---

### Stage 4 — Export GSEA to Python CSVs (`1.4.export_gsea_for_python.R`)

**Script**: `02_Analysis/1.4.export_gsea_for_python.R`

**Inputs**:
- `checkpoints/all_gsea_results.rds`
- `checkpoints/syngo_gsea_results.rds`
- `checkpoints/mitocarta_gsea_results.rds`

**Processing**:
- Extracts and combines all 12 collections × 9 contrasts into long format
- Renames contrasts to trajectory framework names (Early/TrajDev/Late for each mutation)
- Adds `ever_significant` and `ever_significant_trajectory` flags
- Pivots to wide format (one row per pathway, NES and p.adjust columns per contrast)

**Outputs** (under `03_Results/02_Analysis/Python_exports/`):

| File | Description |
|------|-------------|
| `gsea_results_long.csv` | ~31 MB; one row per pathway × contrast × database |
| `gsea_results_wide.csv` | ~5 MB; one row per pathway; NES and p.adjust as wide columns |

---

### Stage 5 — GSVA per-replicate scoring (`1.6.gsva_analysis.R`)

**Script**: `02_Analysis/1.6.gsva_analysis.R`

**Inputs**:
- `checkpoints/qc_variables.rds` (logCPM matrix, sample annotation)
- `checkpoints/all_gsea_results.rds`, `syngo_gsea_results.rds`, `mitocarta_gsea_results.rds` (gene-set universe; loaded from checkpoints to guarantee pathway_id identity with GSEA)
- `00_Data/SynGO_bulk_20231201/`, `00_Data/MitoCarta_3.0/MitoPathways3.0.gmx`

**Processing**:
- GSVA (Bioconductor GSVA 2.2.0) with `gsvaParam(kcdf = "Gaussian", minSize = 10, maxSize = 500)`
- ~14,500 gene sets from same collections as GSEA
- Group-level Welch t-tests (mutant vs Ctrl within each timepoint; BH FDR correction)
- GSVA driver classification per pathway × mutation: Δ_arm = median(D65) − median(D35) for each genotype arm; arms with |Δ| ≥ 0.10 classified as 'moving'; labels: `mutant_driven`, `ctrl_driven`, `both_moving`, `neither_moving`

**Outputs** (under `03_Results/02_Analysis/`):

| Output | Path | Notes |
|--------|------|-------|
| GSVA scores checkpoint | `checkpoints/gsva_all_pathways.rds` | Per-sample score matrix + driver labels |
| Master GSVA table (all) | `master_gsva_all_table.csv` | ~87 K rows; pathway × genotype × timepoint |

---

### Stage 6 — GSVA long-format export (`X.export_gsva_long.R`)

**Script**: `02_Analysis/X.export_gsva_long.R`

**Inputs**:
- `checkpoints/gsva_all_pathways.rds`
- `checkpoints/qc_variables.rds`

**Outputs**:
- `03_Results/02_Analysis/replicate_level_gsva_long.csv` — per-sample GSVA scores in long format (pathway_id, sample_id, genotype, day, gsva_score); consumed by the interactive dashboard

---

### Stage 7 — Master GSEA table + pattern classification (`1.5.create_master_pathway_table.py`)

**Script**: `02_Analysis/1.5.create_master_pathway_table.py`

**Inputs**:
- `Python_exports/gsea_results_long.csv`
- `Python_exports/gsea_results_wide.csv`
- `01_Scripts/Python/pattern_definitions.py` (canonical pattern classifier)
- `01_Scripts/Python/config.py`

**Processing**:
- Merges long-format GSEA statistics with wide-format NES/p.adjust columns
- Calls `add_pattern_classification()` → assigns one of 8 trajectory patterns per pathway per mutation (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient, Complex)
- Classification thresholds: significance FDR < 0.05, |NES_effect| > 0.5; improvement ratio < 0.7; worsening ratio > 1.3; strong effect |NES| > 1.0
- Calls `add_super_category_columns()` → coarser grouping for main-figure summaries
- Tags `ever_significant` universe (5,267 pathways; FDR < 0.05 in at least one of 9 contrasts)

**Outputs**:
- `03_Results/02_Analysis/master_gsea_table.csv` — ~109,989 rows; primary downstream source for all Python visualizations (Supplementary Table S2)

---

### Stage 8 — Focused GSVA master tables (`1.7.create_master_gsva_table.R`)

**Script**: `02_Analysis/1.7.create_master_gsva_table.R`

**Inputs**:
- `checkpoints/gsva_module_scores.rds` (focused 7-module GSVA scores)
- `checkpoints/qc_variables.rds`

**Note**: This script operates on a focused 7-module subset separate from the comprehensive
GSVA run in Stage 5. The comprehensive all-pathway GSVA table is produced by `1.6.gsva_analysis.R`.

**Outputs**:
- `03_Results/02_Analysis/master_gsva_focused_table.csv` — 42 rows (7 modules × 6 groups)
- `03_Results/02_Analysis/gsva_pattern_summary.csv` — pattern classifications for 7 modules (wide)
- `03_Results/02_Analysis/gsva_statistics_summary.txt` — summary statistics

---

### Stage 9 — Supplementary sensitivity analysis

**Script**: Inline within the revision analysis suite (previously `02_Analysis/Supp4.sensitivity_analysis.py` — deleted in this revision; logic preserved in manuscript description).

**Processing**: 81-combination grid over thresholds (NES_EFFECT ∈ {0.4, 0.5, 0.6}; NES_STRONG ∈ {0.8, 1.0, 1.2}; IMPROVEMENT_RATIO ∈ {0.6, 0.7, 0.8}; WORSENING_RATIO ∈ {1.25, 1.3, 1.4}) applied to the 5,267 ever-significant pathway universe.

**Outputs**:
- `03_Results/02_Analysis/Supplementary/6a_sensitivity_stability_digest_5267universe.tsv` (Supplementary Data digest referenced in manuscript as "Supplementary Fig. S8 data")

---

### Stage 10 — AveExpr > 0 supplementary volcano diagnostics

**Script**: `02_Analysis/1.1.main_pipeline.R` (AveExpr-filtered volcano section)

**Outputs**:
- `03_Results/02_Analysis/Plots/Volcano_Supplementary_MinAveExpr0/vertical_fdr_calcium/` (Supp Fig S2E)
- `03_Results/02_Analysis/Plots/Volcano_Supplementary_MinAveExpr0/scatter_aveexpr_significance/` (Supp Fig S2F)
- `03_Results/02_Analysis/5a_aveexpr_*.csv` — AveExpr diagnostic tables

---

### Stage 11 — Python publication visualizations

#### 11a. Pattern summary normalized (Fig 5A)
**Script**: `02_Analysis/3.4.pattern_summary_normalized.py`
- Input: `master_gsea_table.csv`
- Output: `03_Results/02_Analysis/Plots/Pattern_Summary_Normalized/pattern_summary_normalized.pdf`

#### 11b. Trajectory flow / bump charts (Figs 5B, 5C)
**Script**: `02_Analysis/3.5.viz_trajectory_flow.py` and `02_Analysis/3.7.viz_bump_chart.py`
- Inputs: `master_gsea_table.csv`, `Python_exports/gsea_results_wide.csv`
- Outputs:
  - `03_Results/02_Analysis/Plots/Trajectory_Flow/bump_curved_nes_significant.pdf` (Fig 5B; 4,142 pathways)
  - `03_Results/02_Analysis/Plots/Trajectory_Flow/bump_focused_FINAL_paper_combined.pdf` (Fig 5C; 104 MitoCarta + SynGO pathways)
- Module: `01_Scripts/Python/viz_bump_charts.py`

#### 11c. Semantic pathway dotplot (Figs 4 / 6A)
**Script**: `02_Analysis/3.2.publication_figures_dotplot.py`
- Input: `master_gsea_table.csv` (via `01_Scripts/Python/data_loader.py`)
- Output: `03_Results/02_Analysis/Plots/Publication_Figures_Dotplot/Fig4_Semantic_Pathway_Overview_dotplot.pdf`

#### 11d. Per-database pattern heatmap (Supp Fig S7)
**Script**: `02_Analysis/3.1.publication_figures.py` or revision scripts
- Input: `master_gsea_table.csv`, `Tables/per_database_pattern_summary.csv`
- Output: `03_Results/02_Analysis/Plots/Supplementary_6b/per_database_pattern_heatmap.pdf`

#### 11e. Focused panel sensitivity classifications (Supp Fig S8)
- Output: `03_Results/02_Analysis/Plots/Supplementary_8/focused_panel_classifications.pdf`

#### 11f. Cross-compartment ribosome trajectory (Supp Fig S9)
- Output: `03_Results/02_Analysis/Plots/Supplementary_9/cross_compartment_ribosome_trajectory.pdf`
- Output: `03_Results/02_Analysis/Plots/Supplementary_9/geometric_scatter_both_mutations.pdf`

#### 11g. Ribosome UpSet + Euler set (Supp Fig S6)
**Script**: `02_Analysis/3.3.ribosome_upset_plot.py`
- Input: `master_gsea_table.csv`
- Output: `03_Results/02_Analysis/Plots/` (UpSet + Euler figures)

#### 11h. Chord diagrams (Fig 5E)
**Script**: `02_Analysis/3.7.viz_chord_diagrams.py`
- Input: `master_gsea_table.csv`, `Python_exports/gsea_results_wide.csv`
- Output: `03_Results/02_Analysis/Plots/Chord_Diagrams/chord_diagram_G32A.pdf`

#### 11i. Per-database pattern summary table (Supp Table S4)
- Output: `03_Results/02_Analysis/Tables/per_database_pattern_summary.csv`

#### 11j. Ribosome compartment summary table (Supp Table S3)
- Output: `03_Results/02_Analysis/Tables/ribosome_compartment_summary.csv`
  - Key statistics: Jaccard 0.561 synaptic ↔ cytoplasmic; 0 synaptic ↔ mitochondrial

#### 11k. Interactive bump dashboard (Supplementary Data File 1)
**Script**: `02_Analysis/3.8.viz_interactive_bump_dashboard.py` (entry point)
**Module**: `01_Scripts/Python/bump_dashboard/` (`DashboardPipeline` class)
- Inputs: `master_gsea_table.csv`, `replicate_level_gsva_long.csv`, `master_gsva_all_table.csv`
- Output: `03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html`
  (= Supplementary Data File 1; embeds 5,323 pathways with per-replicate GSVA panels and driver-label sidebar filters)

---

### Stage 12 — R publication visualizations (consume RDS checkpoints directly)

| Script | Inputs | Primary Outputs |
|--------|--------|-----------------|
| `02_Analysis/2.2.viz_mito_translation_cascade.R` | `checkpoints/contrast_tables.rds`, `checkpoints/all_gsea_results.rds` | `Plots/Mito_translation_cascade/Mechanistic_Cascade_Heatmap.pdf` (Fig 6B) |
| `02_Analysis/2.3.viz_synaptic_ribosomes.R` | `checkpoints/contrast_tables.rds`, `checkpoints/syngo_gsea_results.rds` | `Plots/Synaptic_ribosomes/Panel_C_Expression_Heatmap.pdf` (Fig 5F); GSEA running-sum PDFs for SynGO (Fig 5D) |
| `02_Analysis/2.4.viz_critical_period_trajectories_gsva.R` | `checkpoints/gsva_all_pathways.rds`, `checkpoints/qc_variables.rds` | `Plots/Supplementary_10/replicate_level_gsva.pdf` (Supp Fig S10); `Plots/Critical_period_trajectories/gsva/` |
| `02_Analysis/2.6.viz_calcium_genes.R` | `checkpoints/contrast_tables.rds`, `checkpoints/qc_variables.rds` | `Calcium_genes/calcium_genes_expression_heatmap.pdf` (Fig 6C); `Calcium_genes/calcium_genes_boxplots.pdf` (Supp Fig S2D); per-contrast calcium volcano PDFs |
| `02_Analysis/2.1.viz_ribosome_paradox.R` | `checkpoints/contrast_tables.rds`, GSEA checkpoints | Ribosome paradox figures |
| `02_Analysis/2.5.viz_complex_v_analysis.R` | GSEA checkpoints | Complex V / ATP synthase subunit figures |
| `02_Analysis/3.6.viz_alluvial_ggalluvial.R` | `master_gsea_table.csv` | Alluvial diagrams |
| `02_Analysis/3.9.viz_pooled_dotplots.R` | `master_gsea_table.csv`, GSEA checkpoints | Cross-database pooled dotplots |

---

## 4. Manuscript Artefact Map

All paths are relative to `03_Results/02_Analysis/`.

> **Note — internal navigation, not the canonical figure list.** This table maps
> analysis outputs to the scripts that produce them, for computational
> reproducibility. The figure/panel numbers reflect an internal working layout;
> after multiple manuscript-revision iterations they may no longer match the
> final published figure numbering, and some rows point to intermediate artefacts
> rather than final paper figures. Treat the published manuscript as the source
> of truth for figure identity — use this table to locate the code and outputs
> behind each analysis, not to enumerate the final figures.

| Figure / Table | Source File Path | Producing Script |
|----------------|-----------------|------------------|
| **Fig 4 (= Fig 6A)** — Semantic pathway dotplot | `Plots/Publication_Figures_Dotplot/Fig4_Semantic_Pathway_Overview_dotplot.pdf` | `02_Analysis/3.2.publication_figures_dotplot.py` |
| **Fig 4A** — Volcano plots | `Plots/Volcano/` | `02_Analysis/1.1.main_pipeline.R` |
| **Fig 4C** — Euler diagram | Euler set in `Plots/` (Supp Fig S6 holds UpSet) | `02_Analysis/3.3.ribosome_upset_plot.py` |
| **Fig 5A** — Pattern distribution | `Plots/Pattern_Summary_Normalized/pattern_summary_normalized.pdf` | `02_Analysis/3.4.pattern_summary_normalized.py` |
| **Fig 5B** — Bump chart (4,142 pathways) | `Plots/Trajectory_Flow/bump_curved_nes_significant.pdf` | `02_Analysis/3.5.viz_trajectory_flow.py` / `3.7.viz_bump_chart.py` |
| **Fig 5C** — Focused bump (MitoCarta + SynGO) | `Plots/Trajectory_Flow/bump_focused_FINAL_paper_combined.pdf` | `02_Analysis/3.5.viz_trajectory_flow.py` |
| **Fig 5D** — SynGO running-sum plots | `Plots/GSEA/{G32A,R403C}_vs_Ctrl_{D35,D65}/SynGO/*_SynGO_running_sum.pdf` | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| **Fig 5E** — Chord diagram (G32A) | `Plots/Chord_Diagrams/chord_diagram_G32A.pdf` | `02_Analysis/3.7.viz_chord_diagrams.py` |
| **Fig 5F** — Synaptic ribosome heatmap | `Plots/Synaptic_ribosomes/Panel_C_Expression_Heatmap.pdf` | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| **Fig 6A** (= Fig 4) | See Fig 4 row above | — |
| **Fig 6B** — Mito translation cascade heatmap | `Plots/Mito_translation_cascade/Mechanistic_Cascade_Heatmap.pdf` | `02_Analysis/2.2.viz_mito_translation_cascade.R` |
| **Fig 6C** — Calcium gene heatmap | `Calcium_genes/calcium_genes_expression_heatmap.pdf` | `02_Analysis/2.6.viz_calcium_genes.R` |
| **Supp Fig S2D** — Calcium boxplots | `Calcium_genes/calcium_genes_boxplots.pdf` | `02_Analysis/2.6.viz_calcium_genes.R` |
| **Supp Fig S2E** — Calcium volcano (AveExpr > 0, vertical, FDR) | `Plots/Volcano_Supplementary_MinAveExpr0/vertical_fdr_calcium/` | `02_Analysis/1.1.main_pipeline.R` |
| **Supp Fig S2F** — AveExpr scatter | `Plots/Volcano_Supplementary_MinAveExpr0/scatter_aveexpr_significance/` | `02_Analysis/1.1.main_pipeline.R` |
| **Supp Fig S6** — UpSet + Euler (ribosome overlap) | `Plots/` (UpSet + Euler set) | `02_Analysis/3.3.ribosome_upset_plot.py` |
| **Supp Fig S7** — Per-database pattern heatmap | `Plots/Supplementary_6b/per_database_pattern_heatmap.pdf` | `02_Analysis/3.1.publication_figures.py` |
| **Supp Fig S8** — Sensitivity focused panel | `Plots/Supplementary_8/focused_panel_classifications.pdf` | Revision sensitivity script |
| **Supp Fig S9** — Cross-compartment ribosome trajectory | `Plots/Supplementary_9/cross_compartment_ribosome_trajectory.pdf` and `geometric_scatter_both_mutations.pdf` | `02_Analysis/3.1.publication_figures.py` |
| **Supp Fig S10** — Replicate-level GSVA | `Plots/Supplementary_10/replicate_level_gsva.pdf` | `02_Analysis/2.4.viz_critical_period_trajectories_gsva.R` |
| **Supp Table S1_DE** | `DE_results/<contrast>_DE_results.csv` | `02_Analysis/1.1.main_pipeline.R` |
| **Supp Table S2** — Master GSEA table | `Tables/master_gsea_table.csv` | `02_Analysis/1.5.create_master_pathway_table.py` |
| **Supp Table S3** — Ribosome compartment summary | `Tables/ribosome_compartment_summary.csv` | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| **Supp Table S4** — Per-database pattern summary | `Tables/per_database_pattern_summary.csv` | `02_Analysis/3.1.publication_figures.py` or `3.4.pattern_summary_normalized.py` |
| **Supplementary Data File 1** — Interactive dashboard | `Plots/Trajectory_Flow/interactive_bump_dashboard.html` | `02_Analysis/3.8.viz_interactive_bump_dashboard.py` |
| **Supplementary Sensitivity digest** | `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv` | Sensitivity analysis script |

---

## 5. Data Loader Module

`01_Scripts/Python/data_loader.py` provides centralised access to the two Python-export files:

```python
from Python.data_loader import load_classified_pathways
df = load_classified_pathways()
# Loads master_gsea_table.csv (NES, patterns, classifications)
# Merges p.adjust columns from gsea_results_wide.csv (trajectory p-values)
```

This is the entry point for all `3.x` Python visualization scripts.

---

## 6. Dependency Graph (condensed)

```
sorted_counts_matrix.txt + metadata.csv
  └─► 1.1.main_pipeline.R
        ├─► checkpoints/contrast_tables.rds + model_objects.rds + qc_variables.rds
        │     ├─► 1.3.add_mitocarta.R → checkpoints/mitocarta_gsea_results.rds
        │     ├─► 1.4.export_gsea_for_python.R → Python_exports/gsea_results_{long,wide}.csv
        │     │     └─► 1.5.create_master_pathway_table.py → master_gsea_table.csv (Supp T S2)
        │     │           └─► 3.1–3.5, 3.7–3.9 Python viz → Plots/
        │     │                 └─► 3.8 (dashboard) also needs replicate_level_gsva_long.csv
        │     ├─► 1.6.gsva_analysis.R → checkpoints/gsva_all_pathways.rds
        │     │     ├─► master_gsva_all_table.csv
        │     │     └─► X.export_gsva_long.R → replicate_level_gsva_long.csv
        │     ├─► 1.7.create_master_gsva_table.R → master_gsva_{focused,all}_table.csv + gsva_pattern_summary.csv
        │     ├─► 2.2.viz_mito_translation_cascade.R → Fig 6B
        │     ├─► 2.3.viz_synaptic_ribosomes.R → Fig 5D, 5F, Supp T S3
        │     ├─► 2.4.viz_critical_period_trajectories_gsva.R → Supp Fig S10
        │     └─► 2.6.viz_calcium_genes.R → Fig 6C, Supp Fig S2D/S2E/S2F
        ├─► checkpoints/all_gsea_results.rds + syngo_gsea_results.rds
              └─► (same downstream as above)
```

---

## 7. Deprecated / Removed in This Revision

The following scripts were deleted in the current commit window (per `git status`). They must
not be cited in methods or lineage documentation going forward.

| Deleted Script | Former Purpose |
|----------------|----------------|
| `02_Analysis/Supp1.verify_enrichments.R` | Enrichment verification QC |
| `02_Analysis/Supp2.diagnose_calcium_genes.R` | Calcium gene expression diagnostics |
| `02_Analysis/Supp3.analyze_de_thresholds.R` | DE threshold sensitivity |
| `02_Analysis/Supp4.sensitivity_analysis.py` | 81-combination pattern sensitivity grid |
| `02_Analysis/Supp5.prepare_explorer_data.py` | Explorer data preparation |
| `02_Analysis/Supp6.app_bump_chart_explorer.py` | Interactive Shiny-style explorer app |
| `02_Analysis/generate_DE_counts_FDR_0.1.R` | FDR 0.10 DE count generation |
| `02_Analysis/regenerate_de_tables.R` | DE table regeneration helper |
| `02_Analysis/regenerate_gsea_plots.R` | GSEA plot regeneration helper |

---

## 8. Change Log

| Date | Change |
|------|--------|
| 2026-05-14 | Full rewrite: added pipeline stages 1–12, manuscript artefact map, GSVA driver classification, interactive dashboard lineage, sensitivity analysis, deprecated-script table; aligned to current submission methods and git status |
| 2025-12-04 | Removed unused intermediate files; renamed from SCRIPT_DEPENDENCY_ANALYSIS.md |
| 2025-12-04 | Comprehensive audit of data flow |
| 2025-12-01 | Initial investigation of 1.4 vs 1.5 relationship |
