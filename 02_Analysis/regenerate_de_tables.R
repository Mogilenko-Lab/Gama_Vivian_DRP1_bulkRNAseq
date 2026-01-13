#!/usr/bin/env Rscript
###############################################################################
##  Regenerate all DE results tables from checkpoint                         ##
###############################################################################
##  Purpose: Create consistent DE result CSV files from contrast_tables.rds
##  Output: 9 *_results.csv files with all genes and complete statistics
###############################################################################

library(here)

# Configuration
checkpoint_dir <- here::here("03_Results/02_Analysis/checkpoints")
output_dir <- here::here("03_Results/02_Analysis/DE_results")
helper_root <- "01_Scripts/RNAseq-toolkit"

# Source helper functions
source(here::here(helper_root, "scripts/utils_plotting.R"))

# Load contrast_tables checkpoint
message("📂 Loading contrast_tables from checkpoint...")
contrast_tables <- readRDS(file.path(checkpoint_dir, "contrast_tables.rds"))

message("📋 Found ", length(contrast_tables), " contrasts:")
for (co in names(contrast_tables)) {
  n_genes <- nrow(contrast_tables[[co]])
  message("  - ", co, ": ", n_genes, " genes")
}

# Save all contrasts as CSV files with standard naming
message("\n💾 Saving CSV files...")
save_contrast_csvs(contrast_tables, output_dir, suffix = "_results.csv")

# Verify output
message("\n✓ Regeneration complete!")
message("  Output directory: ", output_dir)
output_files <- list.files(output_dir, pattern = "_results\\.csv$")
message("  Generated ", length(output_files), " files:")
for (f in output_files) {
  message("    - ", f)
}
