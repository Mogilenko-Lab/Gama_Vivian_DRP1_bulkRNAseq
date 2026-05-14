###############################################################################
##  Comprehensive GSVA Analysis - All Databases                              ##
###############################################################################
##  PURPOSE: Run GSVA on ALL pathways tested in GSEA analysis to create      ##
##           a comprehensive master GSVA table for exploration               ##
##                                                                            ##
##  SCOPE:   - All MSigDB collections (Hallmark, KEGG, Reactome, GO, etc.)   ##
##           - MitoCarta pathways                                            ##
##           - SynGO pathways                                                ##
##           - ~14,500 pathways total (post size-filter)                     ##
##                                                                            ##
##  OUTPUTS: - master_gsva_all_table.csv (comprehensive GSVA table)          ##
##           - gsva_all_pathways_checkpoint.rds (for reuse)                  ##
##                                                                            ##
##  GENE-SET SOURCE OF TRUTH                                                  ##
##  -------------------------                                                ##
##  Gene sets are loaded from the GSEA result checkpoints, NOT from msigdbr ##
##  or upstream Excel/GMT files. Reason: msigdbr changed its API and split  ##
##  KEGG between version 7.x and 10.x (and SynGO's GO-ID-vs-name divergence ##
##  is another source of mismatch). Loading from the GSEA checkpoints       ##
##  guarantees that the GSVA pathway_id universe equals the GSEA universe   ##
##  by construction, so dashboard click-through modals can always join GSVA ##
##  scores onto every dashboard record. See freeze_requirements.R for the   ##
##  full rationale. If the checkpoints are missing, fall back to            ##
##  re-running the GSEA pipeline first; do NOT swap in a fresh msigdbr      ##
##  load here.                                                              ##
###############################################################################

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(GSVA)
library(msigdbr)

message("=" |> rep(78) |> paste(collapse = ""))
message("Comprehensive GSVA Analysis - All Databases")
message("=" |> rep(78) |> paste(collapse = ""))

# ============================================================================ #
# 0. Configuration                                                             #
# ============================================================================ #

config <- list(
  # Paths
  checkpoint_dir = here("03_Results/02_Analysis/checkpoints"),
  out_root = here("03_Results/02_Analysis"),
  syngo_dir = here("00_Data/SynGO_bulk_20231201"),
  mitocarta_file = here("00_Data/MitoCarta_3.0/MitoPathways3.0.gmx"),

  # GSVA parameters
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500,  # Cap at 500 to avoid very large gene sets

  # Species
  species = "Homo sapiens",

  # Force recompute
  force_recompute_gsva = TRUE
)

# ============================================================================ #
# 1. Load Expression Data                                                     #
# ============================================================================ #

message("\n📂 Loading expression data and metadata...")

qc_vars <- readRDS(file.path(config$checkpoint_dir, "qc_variables.rds"))
logCPM <- qc_vars$logCPM
annot <- qc_vars$annot

# Prepare sample metadata
annot_df <- annot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Genotype = gsub("Control", "Ctrl", genotype),
    Genotype = factor(Genotype, levels = c("Ctrl", "G32A", "R403C")),
    Day = ifelse(days == "D35", 35, 65),
    Timepoint = days
  )

message(sprintf("  ✓ Expression matrix: %d genes × %d samples",
                nrow(logCPM), ncol(logCPM)))

# ============================================================================ #
# 2. Load All Gene Sets                                                       #
# ============================================================================ #
#
# The dashboard's source of truth is the GSEA universe. Earlier versions of
# this script independently loaded MSigDB / MitoCarta / SynGO from upstream
# files; that path diverged from the GSEA universe whenever the env's msigdbr
# version differed from the one used to produce the GSEA checkpoint (e.g.,
# msigdbr 7.x has only "CP:KEGG" while msigdbr 10.x has "CP:KEGG_LEGACY" +
# "CP:KEGG_MEDICUS", and SynGO IDs differ between term names and GO IDs).
#
# To guarantee the GSVA pathway_id set matches master_gsea_table.csv exactly,
# we reuse the gene-set lists embedded in the GSEA result checkpoints, which
# were produced by the same upstream pipeline that built the master table.
# Each `gseaResult` object carries an `@geneSets` slot — a named list of
# character vectors with the canonical pathway IDs as names.

message("\n🧬 Loading gene sets from GSEA checkpoints (canonical universe)...")

all_gene_sets <- list()
gene_set_metadata <- data.frame()

gsea_main_ck     <- readRDS(file.path(config$checkpoint_dir, "all_gsea_results.rds"))
mitocarta_ck     <- readRDS(file.path(config$checkpoint_dir, "mitocarta_gsea_results.rds"))
syngo_ck         <- readRDS(file.path(config$checkpoint_dir, "syngo_gsea_results.rds"))

# Each GSEA checkpoint is a per-contrast list of gseaResult objects; we just
# need the gene-set library, which is shared across contrasts.
.first_gsea_genesets <- function(ck) {
  for (i in seq_along(ck)) {
    obj <- ck[[i]]
    if (inherits(obj, "gseaResult") && length(obj@geneSets) > 0) return(obj@geneSets)
  }
  NULL
}

# --- MSigDB collections (10 DBs, multi-contrast checkpoint) ---
msigdb_dbs <- c("hallmark", "canon", "gobp", "gomf", "gocc",
                "kegg", "reactome", "wiki", "cgp", "tf")

contrast_name <- names(gsea_main_ck)[1]
message(sprintf("  Source contrast: %s", contrast_name))

for (db_name in msigdb_dbs) {
  gs_obj <- gsea_main_ck[[contrast_name]][[db_name]]
  if (is.null(gs_obj) || !inherits(gs_obj, "gseaResult")) {
    message(sprintf("    ⚠ %s — no gseaResult found in checkpoint; skipping", db_name))
    next
  }
  gene_sets <- gs_obj@geneSets
  all_gene_sets[[db_name]] <- gene_sets

  gene_set_metadata <- rbind(gene_set_metadata, data.frame(
    pathway_id = paste0(db_name, "::", names(gene_sets)),
    database = db_name,
    pathway_name = names(gene_sets),
    n_genes_total = vapply(gene_sets, length, integer(1)),
    stringsAsFactors = FALSE
  ))
  message(sprintf("    ✓ %-9s %5d sets from checkpoint", db_name, length(gene_sets)))
}

# --- MitoCarta ---
mc_sets <- .first_gsea_genesets(mitocarta_ck)
if (!is.null(mc_sets)) {
  all_gene_sets[["MitoCarta"]] <- mc_sets
  gene_set_metadata <- rbind(gene_set_metadata, data.frame(
    pathway_id = paste0("MitoCarta::", names(mc_sets)),
    database = "MitoCarta",
    pathway_name = names(mc_sets),
    n_genes_total = vapply(mc_sets, length, integer(1)),
    stringsAsFactors = FALSE
  ))
  message(sprintf("    ✓ %-9s %5d sets from checkpoint", "MitoCarta", length(mc_sets)))
}

# --- SynGO ---
syngo_sets <- .first_gsea_genesets(syngo_ck)
if (!is.null(syngo_sets)) {
  all_gene_sets[["SynGO"]] <- syngo_sets
  gene_set_metadata <- rbind(gene_set_metadata, data.frame(
    pathway_id = paste0("SynGO::", names(syngo_sets)),
    database = "SynGO",
    pathway_name = names(syngo_sets),
    n_genes_total = vapply(syngo_sets, length, integer(1)),
    stringsAsFactors = FALSE
  ))
  message(sprintf("    ✓ %-9s %5d sets from checkpoint", "SynGO", length(syngo_sets)))
}

# --- Summary ---

total_gene_sets <- sum(sapply(all_gene_sets, length))
message(sprintf("\n✓ Total gene sets loaded: %d from %d databases",
                total_gene_sets, length(all_gene_sets)))

message("\n  Per-database sets loaded:")
for (db in names(all_gene_sets)) {
  message(sprintf("    %-12s %5d sets", db, length(all_gene_sets[[db]])))
}

# ============================================================================ #
# 3. Filter Gene Sets by Size and Expression                                  #
# ============================================================================ #

message("\n🔍 Filtering gene sets by size and gene availability...")

# Flatten all gene sets with database prefix
all_gene_sets_flat <- list()
for (db in names(all_gene_sets)) {
  for (gs_name in names(all_gene_sets[[db]])) {
    pathway_id <- paste0(db, "::", gs_name)
    all_gene_sets_flat[[pathway_id]] <- all_gene_sets[[db]][[gs_name]]
  }
}

# Filter to genes in expression matrix
all_gene_sets_filtered <- lapply(all_gene_sets_flat, function(genes) {
  genes[genes %in% rownames(logCPM)]
})

# Filter by size
gene_set_sizes <- sapply(all_gene_sets_filtered, length)
valid_sets <- gene_set_sizes >= config$minSize & gene_set_sizes <= config$maxSize

all_gene_sets_filtered <- all_gene_sets_filtered[valid_sets]

# Update metadata
gene_set_metadata <- gene_set_metadata %>%
  mutate(
    n_genes_in_expression = sapply(pathway_id, function(pid) {
      if (pid %in% names(all_gene_sets_filtered)) {
        length(all_gene_sets_filtered[[pid]])
      } else {
        0
      }
    }),
    passed_filter = pathway_id %in% names(all_gene_sets_filtered)
  )

message(sprintf("  ✓ Valid gene sets: %d/%d (size %d-%d, genes in expression matrix)",
                sum(valid_sets), length(all_gene_sets_flat),
                config$minSize, config$maxSize))

message("\n  Per-database passed-filter counts:")
passed_by_db <- gene_set_metadata %>%
  group_by(database) %>%
  summarise(loaded = n(), passed = sum(passed_filter), .groups = "drop")
for (i in seq_len(nrow(passed_by_db))) {
  message(sprintf("    %-12s %5d / %5d passed",
                  passed_by_db$database[i],
                  passed_by_db$passed[i],
                  passed_by_db$loaded[i]))
}

# ============================================================================ #
# 4. Run GSVA                                                                  #
# ============================================================================ #

message("\n🔬 Running GSVA on all pathways...")

gsva_checkpoint_file <- file.path(config$checkpoint_dir, "gsva_all_pathways.rds")

if (!config$force_recompute_gsva && file.exists(gsva_checkpoint_file)) {
  message("  Loading cached GSVA results...")
  gsva_checkpoint <- readRDS(gsva_checkpoint_file)
  gsva_scores <- gsva_checkpoint$scores

} else {
  message(sprintf("  Computing GSVA for %d pathways...", length(all_gene_sets_filtered)))
  message("  This may take 10-30 minutes depending on system resources...")

  # Create GSVA parameter object
  gsva_param <- gsvaParam(
    exprData = as.matrix(logCPM),
    geneSets = all_gene_sets_filtered,
    kcdf = config$kcdf,
    minSize = config$minSize,
    maxSize = config$maxSize
  )

  # Run GSVA with progress
  start_time <- Sys.time()
  gsva_scores <- gsva(gsva_param, verbose = TRUE)
  end_time <- Sys.time()

  time_elapsed <- difftime(end_time, start_time, units = "mins")

  # Save checkpoint
  gsva_checkpoint <- list(
    scores = gsva_scores,
    gene_sets_filtered = all_gene_sets_filtered,
    parameters = list(
      kcdf = config$kcdf,
      minSize = config$minSize,
      maxSize = config$maxSize,
      n_pathways = length(all_gene_sets_filtered),
      n_samples = ncol(logCPM),
      time_elapsed_mins = as.numeric(time_elapsed),
      timestamp = Sys.time()
    )
  )

  saveRDS(gsva_checkpoint, gsva_checkpoint_file)
  message(sprintf("  ✓ GSVA complete in %.1f minutes", time_elapsed))
  message(sprintf("  ✓ Cached to: %s", basename(gsva_checkpoint_file)))
}

message(sprintf("✓ GSVA scores: %d pathways × %d samples",
                nrow(gsva_scores), ncol(gsva_scores)))

# ============================================================================ #
# 5. Transform to Long Format and Add Metadata                                #
# ============================================================================ #

message("\n📊 Transforming to long format...")

# Convert to long format
gsva_long <- gsva_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway_id") %>%
  tidyr::pivot_longer(
    cols = -pathway_id,
    names_to = "Sample",
    values_to = "GSVA_Score"
  )

# Add sample metadata
gsva_long <- gsva_long %>%
  left_join(annot_df, by = "Sample")

# Add pathway metadata
gsva_long <- gsva_long %>%
  left_join(gene_set_metadata %>% select(pathway_id, database, pathway_name,
                                         n_genes_total, n_genes_in_expression),
            by = "pathway_id")

message(sprintf("  ✓ Long format: %d observations", nrow(gsva_long)))

# ============================================================================ #
# 6. Calculate Group Statistics                                               #
# ============================================================================ #

message("\n📊 Calculating group statistics...")

gsva_group_stats <- gsva_long %>%
  group_by(pathway_id, database, pathway_name, Genotype, Day, Timepoint) %>%
  summarize(
    Mean_GSVA = mean(GSVA_Score, na.rm = TRUE),
    SD_GSVA = sd(GSVA_Score, na.rm = TRUE),
    SE_GSVA = SD_GSVA / sqrt(n()),
    N = n(),
    .groups = "drop"
  )

# Calculate baselines and divergence metrics
ctrl_d35_baseline <- gsva_group_stats %>%
  filter(Genotype == "Ctrl", Day == 35) %>%
  select(pathway_id, Baseline_CtrlD35 = Mean_GSVA)

ctrl_by_day <- gsva_group_stats %>%
  filter(Genotype == "Ctrl") %>%
  select(pathway_id, Day, Ctrl_Mean_Same_Day = Mean_GSVA)

gsva_group_stats <- gsva_group_stats %>%
  left_join(ctrl_d35_baseline, by = "pathway_id") %>%
  left_join(ctrl_by_day, by = c("pathway_id", "Day")) %>%
  mutate(
    Expression_vs_CtrlD35 = Mean_GSVA - Baseline_CtrlD35,
    Divergence_vs_Ctrl = ifelse(Genotype == "Ctrl", 0, Mean_GSVA - Ctrl_Mean_Same_Day)
  )

message(sprintf("  ✓ Group statistics: %d groups", nrow(gsva_group_stats)))

# ============================================================================ #
# 7. Statistical Testing                                                      #
# ============================================================================ #

message("\n📊 Performing t-tests (mutants vs controls)...")

perform_ttest <- function(pathway_id, genotype, day, gsva_data) {
  if (genotype == "Ctrl") {
    return(data.frame(
      pathway_id = pathway_id,
      Genotype = genotype,
      Day = day,
      t_statistic = NA,
      p_value = NA,
      p_adjusted = NA,
      significant = FALSE
    ))
  }

  mutant_samples <- gsva_data %>%
    filter(pathway_id == !!pathway_id, Genotype == !!genotype, Day == !!day) %>%
    pull(GSVA_Score)

  ctrl_samples <- gsva_data %>%
    filter(pathway_id == !!pathway_id, Genotype == "Ctrl", Day == !!day) %>%
    pull(GSVA_Score)

  if (length(mutant_samples) >= 2 && length(ctrl_samples) >= 2) {
    test_result <- t.test(mutant_samples, ctrl_samples)
    return(data.frame(
      pathway_id = pathway_id,
      Genotype = genotype,
      Day = day,
      t_statistic = test_result$statistic,
      p_value = test_result$p.value,
      p_adjusted = NA,
      significant = FALSE
    ))
  } else {
    return(data.frame(
      pathway_id = pathway_id,
      Genotype = genotype,
      Day = day,
      t_statistic = NA,
      p_value = NA,
      p_adjusted = NA,
      significant = FALSE
    ))
  }
}

test_combinations <- gsva_group_stats %>%
  filter(Genotype != "Ctrl") %>%
  select(pathway_id, Genotype, Day) %>%
  distinct()

message(sprintf("  Running %d comparisons...", nrow(test_combinations)))

ttest_results <- purrr::pmap_dfr(
  test_combinations,
  function(pathway_id, Genotype, Day) {
    perform_ttest(pathway_id, Genotype, Day, gsva_long)
  }
)

# Adjust p-values
ttest_results <- ttest_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = !is.na(p_adjusted) & p_adjusted < 0.05
  )

message(sprintf("  ✓ Significant comparisons (p.adj < 0.05): %d/%d",
                sum(ttest_results$significant, na.rm = TRUE),
                sum(!is.na(ttest_results$p_adjusted))))

# ============================================================================ #
# 8. Create Master GSVA Table                                                 #
# ============================================================================ #

message("\n📊 Creating master GSVA table...")

master_gsva <- gsva_group_stats %>%
  left_join(ttest_results, by = c("pathway_id", "Genotype", "Day")) %>%
  left_join(gene_set_metadata %>% select(pathway_id, n_genes_total, n_genes_in_expression),
            by = "pathway_id") %>%
  mutate(
    Trajectory_Category = case_when(
      Day == 35 & Genotype != "Ctrl" ~ "Early",
      Day == 65 & Genotype != "Ctrl" ~ "Late",
      Genotype == "Ctrl" ~ "Reference",
      TRUE ~ "Other"
    ),
    Contrast_Label = case_when(
      Day == 35 & Genotype == "G32A" ~ "Early_G32A",
      Day == 35 & Genotype == "R403C" ~ "Early_R403C",
      Day == 65 & Genotype == "G32A" ~ "Late_G32A",
      Day == 65 & Genotype == "R403C" ~ "Late_R403C",
      Day == 35 & Genotype == "Ctrl" ~ "Ctrl_D35",
      Day == 65 & Genotype == "Ctrl" ~ "Ctrl_D65",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    pathway_id, database, pathway_name,
    n_genes_total, n_genes_in_expression,
    Genotype, Day, Timepoint, Trajectory_Category, Contrast_Label,
    Mean_GSVA, SD_GSVA, SE_GSVA, N,
    Baseline_CtrlD35, Expression_vs_CtrlD35,
    Ctrl_Mean_Same_Day, Divergence_vs_Ctrl,
    t_statistic, p_value, p_adjusted, significant
  ) %>%
  arrange(database, pathway_name, Day, Genotype)

message(sprintf("  ✓ Master table: %d rows × %d columns",
                nrow(master_gsva), ncol(master_gsva)))

# ============================================================================ #
# 9. Export Master Table                                                      #
# ============================================================================ #

message("\n💾 Exporting master GSVA table...")

output_file <- file.path(config$out_root, "master_gsva_all_table.csv")
write.csv(master_gsva, output_file, row.names = FALSE)

message(sprintf("  ✓ %s", basename(output_file)))
message(sprintf("    - %d rows (pathway × genotype × timepoint)", nrow(master_gsva)))
message(sprintf("    - %d unique pathways", length(unique(master_gsva$pathway_id))))
message(sprintf("    - %d databases", length(unique(master_gsva$database))))

# ============================================================================ #
# 10. Summary Statistics                                                      #
# ============================================================================ #

message("\n📊 Summary by database:")

db_summary <- master_gsva %>%
  group_by(database) %>%
  summarize(
    n_pathways = n_distinct(pathway_id),
    n_significant = sum(significant, na.rm = TRUE),
    pct_significant = 100 * n_significant / sum(!is.na(significant))
  ) %>%
  arrange(desc(n_pathways))

print(db_summary)

message("\n📊 Top significant pathways (by p-value):")

top_sig <- master_gsva %>%
  filter(significant == TRUE) %>%
  arrange(p_adjusted) %>%
  head(10) %>%
  select(database, pathway_name, Genotype, Day, Divergence_vs_Ctrl, p_adjusted)

if (nrow(top_sig) > 0) {
  print(top_sig)
} else {
  message("  No significant pathways found (may need more replicates or larger effects)")
}

# ============================================================================ #
# 11. Complete                                                                #
# ============================================================================ #

message("\n" |> paste0(rep("=", 78) |> paste(collapse = "")))
message("✅ COMPREHENSIVE GSVA ANALYSIS COMPLETE")
message(rep("=", 78) |> paste(collapse = ""))

message("\n📁 Output file:")
message(sprintf("  %s", output_file))

message("\n📊 Summary:")
message(sprintf("  - Pathways analyzed: %d", length(unique(master_gsva$pathway_id))))
message(sprintf("  - Databases: %s", paste(unique(master_gsva$database), collapse = ", ")))
message(sprintf("  - Total comparisons: %d", sum(!is.na(master_gsva$p_value))))
message(sprintf("  - Significant comparisons: %d (%.1f%%)",
                sum(master_gsva$significant, na.rm = TRUE),
                100 * sum(master_gsva$significant, na.rm = TRUE) / sum(!is.na(master_gsva$significant))))

message("\n✨ Done!")
message("This comprehensive table now matches the scope of the GSEA master table,")
message("allowing you to explore GSVA trajectories for all pathways!\n")