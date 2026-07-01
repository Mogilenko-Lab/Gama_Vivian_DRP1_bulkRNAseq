# CLAUDE.md

General guidance to Claude Code when working with code in this repository.

## Project Overview

This is a bulk RNA-seq differential expression analysis of DRP1 mutations (G32A and R403C) in iPSC-derived cortical neurons across two maturation timepoints (Day 35 and Day 65). The repository contains a complete downstream analysis pipeline starting from count matrices.

**Key Context:**
- **Frozen analysis:** Git tag `v1.0` (commit `d6ec164`) represents the manuscript submission version
- **Container:** Docker-based dev environment (`scdock-r-dev:v0.5.1`) via VS Code Dev Containers
- **Data scope:** Counts matrices and metadata only (no raw FASTQ files)
- **Reference genome:** GRCh38.p14
- **Path handling:** All paths use `here::here()` for workspace-relative references

## Essential Commands

### Running the Analysis

**Main differential expression pipeline:**
```bash
Rscript 02_Analysis/1.1.main_pipeline.R
```
- Uses checkpoint caching to avoid recomputation
- Set `config$force_recompute = TRUE` to ignore cached results
- Outputs to `03_Results/02_Analysis/`

**Generate contrast-specific result tables:**
```bash
Rscript 02_Analysis/1.2.generate_contrast_tables.R
```

**Add MitoCarta pathway annotations:**
```bash
Rscript 02_Analysis/1.3.add_mitocarta.R
```

**Export GSEA data for Python analysis:**
```bash
Rscript 02_Analysis/1.4.export_gsea_for_python.R
```

**Create comprehensive master tables:**
```bash
python3 02_Analysis/1.5.create_master_pathway_table.py     # Master GSEA table (all databases)
Rscript 02_Analysis/1.7.create_master_gsva_table.R       # Master GSVA tables (focused + all)
# Note: 1.7 produces both focused (7 modules) and comprehensive (all pathways) tables
```

**Publication figure generation (Python):**
```bash
python3 02_Analysis/3.1.publication_figures.py          # Main publication figures
python3 02_Analysis/3.2.publication_figures_dotplot.py  # Dotplot version
python3 02_Analysis/3.3.ribosome_upset_plot.py          # Ribosome overlap UpSet plot
python3 02_Analysis/3.4.pattern_summary_normalized.py    # Normalized pattern visualizations
python3 02_Analysis/3.5.viz_trajectory_flow.py          # Trajectory flow (alluvial diagrams)
```

### Development Setup

**Initialize git submodules:**
```bash
git submodule update --init --recursive
```
- Required for `RNAseq-toolkit` helper functions
- Submodule location: `01_Scripts/RNAseq-toolkit`

**Verify environment:**
```r
library(here)
here::here()  # Should show /workspaces/Gama_Vivian_DRP1_bulkRNAseq

# Test data access
counts_file <- here::here("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt")
file.exists(counts_file)  # Should be TRUE
```

### Container Management

**Launch container:**
- Open folder in VS Code
- Command Palette → "Dev Containers: Reopen in Container"

**Access R environments:**
- `radian` - Interactive R console (recommended)
- `r-base` or `R` - Standard R session

## Architecture

### Repository Structure

```
00_Data/                        # Reference databases (SynGO, MitoCarta)
01_Scripts/
  ├── RNAseq-toolkit/          # Git submodule with GSEA & DE helpers
  └── R_scripts/               # Project-specific helper functions
02_Analysis/                    # Analysis pipeline scripts
  ├── 1.1.main_pipeline.R       # Core DE analysis with checkpointing
  ├── 1.2.generate_contrast_tables.R
  ├── 1.3.add_mitocarta.R
  ├── 1.4.export_gsea_for_python.R
  ├── 1.5.create_master_pathway_table.py   # Master GSEA table
  ├── 1.6.gsva_analysis.R                  # GSVA score computation
  ├── 1.7.create_master_gsva_table.R       # Master GSVA tables (focused + all)
  ├── viz_*.R                  # Visualization scripts
  └── *.py                     # Python analysis scripts
03_Results/
  ├── 01_Preprocessing/        # Counts matrices and metadata
  └── 02_Analysis/             # DE results, GSEA, plots
      ├── checkpoints/         # Cached computation results (.rds)
      ├── DE_results/          # Differential expression tables
      ├── Python_exports/      # Exported CSV files for Python analysis
      ├── Plots/               # All visualizations
      ├── Verification_reports/
      ├── master_gsea_table.csv            # Comprehensive GSEA results (109K rows)
      ├── master_gsva_all_table.csv        # Comprehensive GSVA all pathways (87K rows)
      ├── master_gsva_focused_table.csv    # Focused GSVA 7 modules (42 rows)
      ├── gsva_pattern_summary.csv         # GSVA pattern classifications
      └── gsva_statistics_summary.txt      # GSVA summary statistics
```

### Core Analysis Pipeline

The analysis follows a multi-stage pipeline:

1. **Data Loading & Normalization** (`1.1.main_pipeline.R`)
   - Loads counts matrix from `03_Results/01_Preprocessing/04_FeatureCounts/`
   - TMM normalization via edgeR
   - Limma-voom for differential expression

2. **Experimental Design**
   - Factorial design: genotype × timepoint
   - 3 genotypes: Ctrl, G32A (GTPase domain), R403C (stalk domain)
   - 2 timepoints: Day 35 (early), Day 65 (mature)

3. **Contrast Definitions**
   - **Mutation effects:** `G32A_vs_Ctrl_D35`, `R403C_vs_Ctrl_D35`, `G32A_vs_Ctrl_D65`, `R403C_vs_Ctrl_D65`
   - **Maturation effects:** `Time_Ctrl`, `Time_G32A`, `Time_R403C`
   - **Interactions:** `Maturation_G32A_specific`, `Maturation_R403C_specific`

4. **Gene Set Enrichment Analysis (GSEA)**
   - MSigDB collections: Hallmark, KEGG, Reactome, GO:BP/CC/MF, CGP, TF, WikiPathways
   - Custom databases: SynGO (synaptic ontology), MitoCarta (mitochondrial pathways)
   - Results cached in `checkpoints/gsea_*.rds`

5. **Visualization**
   - Volcano plots with calcium gene highlighting
   - GSEA dotplots, barplots, running sum plots
   - Trajectory heatmaps (Early → TrajDev → Late framework)

### Checkpoint System

The pipeline uses checkpointing to avoid expensive recomputation:

```r
# Example checkpoint usage
result <- load_or_compute(
  checkpoint_file = here::here(checkpoint_dir, "gsea_hallmark.rds"),
  compute_fn = function() { run_gsea(...) },
  force_recompute = config$force_recompute,
  description = "Hallmark GSEA"
)
```

**Checkpoint locations:** `03_Results/02_Analysis/checkpoints/*.rds`

**To force recomputation:**
- Edit `config$force_recompute = TRUE` in `1.1.main_pipeline.R`
- Or manually delete checkpoint files

### RNAseq-toolkit Submodule

The `01_Scripts/RNAseq-toolkit` submodule provides:

- **GSEA processing:** `run_gsea()`, `get_significant_pathways()`, `calculate_pathway_scores()`
- **GSEA plotting:** `gsea_dotplot()`, `gsea_barplot()`, `gsea_running_sum_plot()`
- **DE plotting:** `create_volcano_plot()`, `analyze_pathway_volcano()`
- **Custom theme:** `custom_minimal_theme_with_grid()`

**Sourcing helpers in analysis scripts:**
```r
config <- list(helper_root = "01_Scripts/RNAseq-toolkit")
source(here::here(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_dotplot.R"))
```

### Key Analysis Patterns

**1. Loading DE results:**
```r
de_results <- read.csv(here::here("03_Results/02_Analysis/DE_results/G32A_vs_Ctrl_D35_DE_results.csv"))
```

**2. Running GSEA on contrast:**
```r
library(clusterProfiler)
library(msigdbr)
source(here::here("01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R"))

gsea_result <- run_gsea(
  de_results = de_results,
  rank_metric = "t",
  species = "Homo sapiens",
  msigdb_collection = "H",
  pvalueCutoff = 0.05
)
```

**3. Creating volcano plots with gene highlighting:**
```r
source(here::here("01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R"))

# FDR mode (default): p_cutoff interpreted as FDR threshold
create_standard_volcano(
  de_results,
  highlight_gene = config$calcium_genes,
  p_cutoff = 0.1,           # FDR threshold when decision_by = "fdr"
  decision_by = "fdr",      # Default; uses adj.P.Val for coloring
  fc_cutoff = config$fc_cutoff
)

# Raw p-value mode: p_cutoff interpreted as raw p threshold
create_standard_volcano(
  de_results,
  highlight_gene = config$calcium_genes,
  p_cutoff = config$p_cutoff,  # Raw p threshold when decision_by = "p"
  decision_by = "p",
  fc_cutoff = config$fc_cutoff
)
```

**4. Trajectory framework mapping:**
```
Early: D35 mutation effects (G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35)
TrajDev: Mutation-specific maturation changes (Maturation_G32A_specific, Maturation_R403C_specific)
Late: D65 mutation effects (G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D65)
```

### Pattern Classification System

**Canonical source:** `01_Scripts/Python/pattern_definitions.py`
**Full documentation:** `docs/PATTERN_CLASSIFICATION.md`

Pathways are classified into 8 mutually exclusive patterns based on trajectory dynamics:

| Pattern | Active? | Description |
|---------|---------|-------------|
| **Compensation** | Active | Early defect + TrajDev opposes + Late improved |
| **Sign_reversal** | Active | TrajDev opposes + sign flipped between Early/Late |
| **Progressive** | Active | Early defect + TrajDev amplifies + Late worsened |
| **Natural_worsening** | Passive | Early defect + TrajDev NS + Late worsened |
| **Natural_improvement** | Passive | Early defect + TrajDev NS + Late improved |
| **Late_onset** | - | No Early defect + Late defect emerges |
| **Transient** | - | Strong Early defect + Late resolved |
| **Complex** | - | Inconsistent or multiphasic |

**Key thresholds:**
- Significance: p.adjust < 0.05 (High), < 0.10 (Medium)
- Effect size: |NES| > 0.5 (minimum), |NES| > 1.0 (strong)
- Improvement: |Late|/|Early| < 0.7 (≥30% reduction)
- Worsening: |Late|/|Early| > 1.3 (≥30% increase)

**Active vs Passive distinction:**
- **Active** (Compensation, Sign_reversal, Progressive): Requires significant TrajDev (p < 0.05, |NES| > 0.5)
- **Passive** (Natural_improvement, Natural_worsening): TrajDev not significant

**Super-categories (simplified for main figures/text):**

| Super-Category | Includes |
|----------------|----------|
| Active_Compensation | Compensation |
| Active_Reversal | Sign_reversal |
| Active_Progression | Progressive |
| Passive | Natural_improvement, Natural_worsening |
| Late_onset | Late_onset |
| Other | Transient, Complex |

**Complex subtypes** (for detailed analysis):
- Complex_opposing: TrajDev opposes Early but magnitude insufficient
- Complex_amplifying: TrajDev amplifies Early but magnitude insufficient
- Complex_multiphasic: Inconsistent directionality across trajectory stages

**Usage:**
```python
from Python.pattern_definitions import classify_pattern, add_pattern_classification, add_super_category_columns

# Single pathway
pattern, confidence = classify_pattern(early_nes, early_padj, trajdev_nes, trajdev_padj, late_nes, late_padj)

# Batch classification
df = add_pattern_classification(gsea_wide_df, mutations=['G32A', 'R403C'])
```

### Pattern System Synchronization Checklist

**CRITICAL:** The pattern classification system is implemented in both Python and R. When modifying the pattern system, you MUST maintain synchronization between implementations.

**Files that MUST stay aligned:**

1. **Python (canonical source):**
   - `01_Scripts/Python/pattern_definitions.py` (lines 47-60: thresholds, lines 276-411: classification logic)
   - `01_Scripts/Python/pattern_definitions.py` (lines 80-90: SUPER_CATEGORY_MAP)

2. **R (GSVA implementation):**
   - `02_Analysis/1.7.create_master_gsva_table.R` (lines 348-354: thresholds, lines 379-486: classification logic)
   - `02_Analysis/1.7.create_master_gsva_table.R` (lines 527-543: super_category_map)

3. **Documentation:**
   - `docs/PATTERN_CLASSIFICATION.md` (primary reference)
   - `CLAUDE.md` (this file, pattern summary section)

**When adding or modifying patterns:**

✅ **Required updates:**
1. Update pattern definition in `pattern_definitions.py` (add to PATTERN_DEFINITIONS dict)
2. Add pattern to appropriate constant lists (MEANINGFUL_PATTERNS, ACTIVE_PATTERNS, etc.)
3. Update classification logic in `pattern_definitions.py` (classify_pattern function)
4. Mirror changes in `1.7.create_master_gsva_table.R` (classify_gsva_pattern function)
5. Update super-category mappings in both files (SUPER_CATEGORY_MAP)
6. Update `docs/PATTERN_CLASSIFICATION.md` with full pattern specification
7. Update this file (CLAUDE.md) with summary

✅ **Verification steps:**
1. Regenerate master tables:
   ```bash
   python3 02_Analysis/1.5.create_master_pathway_table.py
   Rscript 02_Analysis/1.7.create_master_gsva_table.R
   ```
2. Compare pattern distributions in both outputs (should be consistent)
3. Verify threshold values match exactly:
   - Python: `pattern_definitions.py` lines 47-60
   - R: `1.7.create_master_gsva_table.R` lines 348-354
4. Check pattern classification order is identical in both implementations

✅ **Key differences to maintain:**
- **TrajDev significance:** Python uses p-values (p.adjust < 0.05), R uses magnitude threshold (|TrajDev| > GSVA_EFFECT)
- **Reason:** GSVA TrajDev is calculated (difference of differences), not statistically tested
- **Both approaches are valid** and documented in respective implementation comments

⚠️ **Common mistakes:**
- Adding pattern to Python but forgetting R implementation
- Changing thresholds in one file but not the other
- Updating super-category map in one place only
- Forgetting to update documentation after code changes

📝 **Testing alignment (future):**
- Create test suite with example trajectories
- Verify both R and Python classify identically
- Run on every pattern system change

## Important Configuration

### Analysis Parameters

Located in `02_Analysis/1.1.main_pipeline.R`:
```r
config <- list(
  counts_file   = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt",
  metadata_file = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv",
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",
  fdr_cutoff    = 0.05,      # FDR threshold for DEG classification (decideTests with BH)
  p_cutoff      = 0.05,      # Raw p-value threshold for volcano plots (p mode)
  fc_cutoff     = 2,         # |log2FC| >= 2 (4-fold change) for volcano visualization
  calcium_genes = c("NNAT","CACNG3","CACNA1S","ATP2A1","RYR1","MYLK3","VDR","STIM1","STIM2","ORAI1_1","CALB1","CALR","PNPO"),
  syngo_dir     = "00_Data/SynGO_bulk_20231201",
  syngo_ns      = "CC",
  mitocarta_file = "00_Data/MitoCarta_3.0/MitoPathways3.0.gmx",
  force_recompute = FALSE
)
```

### Species and Gene ID Handling

- **Species:** Homo sapiens (human)
- **Gene identifiers:** HGNC gene symbols
- **Annotation:** org.Hs.eg.db for Entrez ID conversion

### Critical SynGO Integration

SynGO (Synapse Gene Ontology) requires special handling:

```r
source(here::here("01_Scripts/R_scripts/run_syngo_gsea.R"))

syngo_result <- run_syngo_gsea(
  ranked_genes = ranked_gene_list,
  syngo_dir = config$syngo_dir,
  namespace = config$syngo_ns,  # "CC" for cellular component
  pvalueCutoff = 0.05
)
```

**SynGO data files:**
- `00_Data/SynGO_bulk_20231201/syngo_genes.xlsx`
- `00_Data/SynGO_bulk_20231201/syngo_ontologies.xlsx`
- `00_Data/SynGO_bulk_20231201/syngo_annotations.xlsx`

## Common Tasks

### Adding a New Contrast

1. Define contrast in design matrix (edit `1.1.main_pipeline.R`)
2. Add to contrasts vector
3. Run through GSEA pipeline
4. Generate visualizations

### Modifying GSEA Databases

Edit the database list in `1.1.main_pipeline.R`:
```r
gsea_databases <- c("H", "C2", "C5", "C6")  # Add/remove MSigDB collections
```

### Creating Custom Visualizations

Use existing viz scripts as templates:
- `2.1.viz_ribosome_paradox.R` - Translation crisis visualization
- `2.2.viz_mito_translation_cascade.R` - Pathway-specific deep dives
- `2.3.viz_synaptic_ribosomes.R` - Compartment-specific ribosome analysis
- `2.4.viz_critical_period_trajectories_gsva.R` - Temporal trajectory patterns with GSVA
- `viz_developmental_framework.R` - Developmental pattern analysis
- `2.5.viz_complex_v_analysis.R` - ATP synthase (Complex V) deep-dive
- `2.6.viz_calcium_genes.R` - Calcium signaling gene analysis
- `viz_pooled_dotplots.R` - Cross-database GSEA dotplots

### Creating Master Summary Tables

Three comprehensive master tables provide complete analysis results:

**Master GSEA Table:**
```bash
python3 02_Analysis/1.5.create_master_pathway_table.py
```
- **Output**: `master_gsea_table.csv` (109,989 rows)
- **Contents**: All GSEA pathways across all contrasts with pattern classifications
- **Use for**: Comprehensive pathway analysis, pattern distribution queries, cross-database comparisons

**Master GSVA Tables (Focused + All):**
```bash
Rscript 02_Analysis/1.7.create_master_gsva_table.R
```
- **Outputs**:
  - `master_gsva_focused_table.csv`: GSVA scores in long format (42 rows: 7 modules × 6 groups)
  - `gsva_pattern_summary.csv`: Pattern classifications in wide format (7 modules)
  - `gsva_statistics_summary.txt`: Summary statistics and usage guide
- **Contents**: GSVA enrichment scores for 7 key trajectory modules with statistical tests
- **Use for**: Module-level trajectory analysis, pattern classification, group comparisons

The `1.7.create_master_gsva_table.R` script also produces `master_gsva_all_table.csv` (87K rows) containing GSVA scores for all ~14,500 pathways. This is used by the interactive explorer.

See `03_Results/02_Analysis/README.md` for complete documentation.

### Regenerating All Figures

```bash
# From main pipeline (slow, uses checkpoints)
Rscript 02_Analysis/1.1.main_pipeline.R

# Create master tables (optional, for comprehensive queries)
python3 02_Analysis/1.5.create_master_pathway_table.py
Rscript 02_Analysis/1.7.create_master_gsva_table.R  # Produces focused + all tables

# Individual visualizations (fast, uses cached results)
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R  # Classical alluvial diagram

# Python visualizations (run after R exports)
python3 02_Analysis/3.1.publication_figures.py
python3 02_Analysis/3.4.pattern_summary_normalized.py
python3 02_Analysis/3.5.viz_trajectory_flow.py
```

## Critical Implementation Details

### Path Handling

**Always use `here::here()` for file paths:**
```r
# Good
counts <- read.table(here::here("03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt"))

# Bad - will break if working directory changes
counts <- read.table("/workspaces/GVDRP1/03_Results/...")
```

### Sourcing Helper Functions

Helper sourcing uses a safe wrapper:
```r
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found → ", path)
  }
}

source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
```

### Checkpoint Management

Checkpoints prevent recomputation but can become stale:

**When to delete checkpoints:**
- Changed input data (counts matrix or metadata)
- Modified GSEA parameters
- Updated contrast definitions
- Debugging unexpected results

**Delete specific checkpoints:**
```bash
rm 03_Results/02_Analysis/checkpoints/gsea_hallmark.rds
```

**Delete all checkpoints:**
```bash
rm 03_Results/02_Analysis/checkpoints/*.rds
```

### GSEA Result Storage

GSEA results are stored as nested lists:
```r
# Structure
gsea_results <- list(
  contrast_name = list(
    database_name = gseaResult_object
  )
)

# Access
hallmark_result <- gsea_results[["G32A_vs_Ctrl_D35"]][["hallmark"]]
```

### Handling Missing Pathways

Some databases may have no significant enrichments:
```r
# Check before plotting
if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
  plot <- gsea_dotplot(gsea_result)
} else {
  message("No significant pathways for ", contrast, " in ", database)
}
```

## Migration Context

This repository was migrated from an original workstation on 2025-11-19. Key changes:

**Removed dependencies:**
- Raw FASTQ files (preprocessing complete)
- Reference genome mounts (alignment complete)
- Bash preprocessing scripts (no longer needed)

**Added components:**
- Docker Compose configuration
- Environment variable support (`.env` file)
- Path portability improvements
- Enhanced documentation

**See also:**
- `MIGRATION.md` - Complete migration documentation
- `SETUP.md` - Container setup instructions
- `CHANGELOG.md` - Tracking analysis updates

## Key Findings Summary

**Maturation paradox:**
- Synaptic ribosome programs downregulated despite increased ribosome biogenesis
- Suggests ATP/translation bottleneck at synapses

**Mitochondrial compensation:**
- 4 pathways show shared compensation in both mutations:
  1. Mitochondrial ribosome
  2. Mitochondrial ribosome assembly
  3. mtDNA maintenance
  4. OXPHOS

**Mutation-specific patterns:**
- G32A: Strong synaptic pathway compensation
- R403C: Broader compensation profile (100 vs 74 pathways)

**Trajectory framework:**
- Early (D35): Initial mutation effects
- TrajDev: Developmental trajectory divergence
- Late (D65): Mature neuron state

## Troubleshooting

### "Cannot find RNAseq-toolkit module"
```bash
git submodule update --init --recursive
```

### "Permission denied" in container
Check `.devcontainer/.env` has correct UID/GID:
```bash
id -u  # Get your UID
id -g  # Get your GID
# Update .devcontainer/.env accordingly
```

### Checkpoint loading fails
Delete and regenerate:
```bash
rm 03_Results/02_Analysis/checkpoints/problematic_checkpoint.rds
Rscript 02_Analysis/1.1.main_pipeline.R
```

### GSEA fails with gene ID errors
Verify gene symbols in DE results match Entrez conversion:
```r
library(org.Hs.eg.db)
genes <- de_results$gene
entrez <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
table(is.na(entrez))  # Check conversion success rate
```

### Plot labels unreadable
Increase figure size or reduce number of pathways:
```r
# In plotting functions
ggsave("output.pdf", width = 12, height = 10)  # Increase dimensions

# Or filter to top N pathways
top_pathways <- head(pathways, 20)
```

## Additional Documentation

- `README.md` - Project overview and quick facts
- `MIGRATION.md` - Migration documentation
- `SETUP.md` - Container setup guide
- `CHANGELOG.md` - Analysis change tracking
- `ISSUES.md` - Known issues and limitations
- `VOLCANO_IMPROVEMENTS_SUMMARY.md` - Volcano plot refinements
- `SESSION_SUMMARY*.md` - Session handoff notes
- `02_Analysis/README_trajectory_analysis.md` - Trajectory framework documentation
- `03_Results/02_Analysis/Plots/Cross_database_validation/KEY_FINDINGS_SUMMARY.md` - Biological interpretation
