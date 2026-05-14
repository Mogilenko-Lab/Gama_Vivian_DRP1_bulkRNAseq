# 02_Analysis/Supp7.verify_fdr_impact.R
#
# Objective:
#   1. Verify the specific statistical claims about NNAT from the manuscript.
#   2. Quantify the number of DEGs at FDR < 0.05 vs. FDR < 0.1 for all contrasts.
#   3. Check the significance status of key calcium-related genes at both thresholds.
#
# Input:
#   - DEG result CSVs in 03_Results/02_Analysis/DE_results/
#
# Output:
#   - A summary printed to the console.

library(dplyr)
library(readr)
library(purrr)

# --- Configuration ---

# Gene list from the user prompt
key_genes <- c(
  "NNAT", "CACNG3", "CACNA1C", "CACNA1S", "ATP2A1", "RYR1",
  "MYLK3", "CASR", "VDR", "STIM1", "STIM2", "ORAI1",
  "CALB1", "CALR", "PNPO"
)

# Base path for the results
results_dir <- here::here("03_Results/02_Analysis/DE_results")

# Get all contrast result files
deg_files <- list.files(results_dir, pattern = "\.csv$", full.names = TRUE)

# --- NNAT Verification ---
cat("--- Verifying NNAT Statistics ---

")

# Function to read a file and extract stats for a specific gene
verify_gene_stats <- function(contrast_name, gene_symbol) {
  file_path <- file.path(results_dir, paste0(contrast_name, "_results.csv"))
  if (!file.exists(file_path)) {
    cat(sprintf("SKIPPED: Could not find file for contrast '%s'
", contrast_name))
    return(NULL)
  }
  
  results <- read_csv(file_path, show_col_types = FALSE) %>%
    filter(gene == gene_symbol)
  
  if (nrow(results) > 0) {
    cat(sprintf("Stats for '%s' in contrast '%s':
", gene_symbol, contrast_name))
    cat(sprintf("  Log2FC: %.4f
", results$logFC))
    cat(sprintf("  adj.P.Val (FDR): %e
", results$adj.P.Val))
    cat("
")
  } else {
    cat(sprintf("Gene '%s' not found in contrast '%s'

", gene_symbol, contrast_name))
  }
}

verify_gene_stats("G32A_vs_Ctrl_D35", "NNAT")
verify_gene_stats("R403C_vs_Ctrl_D35", "NNAT")

# --- Impact Analysis ---
cat("
--- Analyzing Impact of FDR Thresholds ---

")

# Function to process a single DEG file
analyze_fdr_impact <- function(file_path) {
  contrast_name <- sub("_results.csv", "", basename(file_path), fixed = TRUE)
  
  results_df <- read_csv(file_path, show_col_types = FALSE)
  
  # Count DEGs at different thresholds
  degs_at_0.05 <- sum(results_df$adj.P.Val < 0.05, na.rm = TRUE)
  degs_at_0.1 <- sum(results_df$adj.P.Val < 0.1, na.rm = TRUE)
  
  # Check key genes that fall between the thresholds
  marginal_genes <- results_df %>%
    filter(gene %in% key_genes & adj.P.Val >= 0.05 & adj.P.Val < 0.1) %>%
    pull(gene)
  
  list(
    contrast = contrast_name,
    degs_fdr_0.05 = degs_at_0.05,
    degs_fdr_0.1 = degs_at_0.1,
    marginal_key_genes = if(length(marginal_genes) > 0) paste(marginal_genes, collapse=", ") else "None"
  )
}

# Apply the function to all DEG files and compile the results
impact_summary <- map_dfr(deg_files, analyze_fdr_impact)

cat("DEG Counts at Different FDR Thresholds:
")
print(as.data.frame(impact_summary))

cat("
Analysis complete.
")
