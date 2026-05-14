#!/usr/bin/env Rscript
# Supp10a.export_gsva_modules.R
# ----------------------------------------------------------------------------
# Companion to Supp10.replicate_level_gsva.py — exports the per-sample GSVA
# module scores from the binary RDS checkpoint into a tidy long-format CSV
# that the Python plotting script can read.
#
# Inputs (read-only)
#   - 03_Results/02_Analysis/checkpoints/gsva_module_scores.rds
#   - 03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv
#
# Output
#   - 03_Results/02_Analysis/Plots/Supplementary_10/replicate_level_gsva_per_sample.csv
#     Columns: sample, genotype, days, rep, module, gsva_score
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

PROJECT_ROOT <- here::here()
if (!nzchar(PROJECT_ROOT)) PROJECT_ROOT <- getwd()

rds_path  <- file.path(PROJECT_ROOT,
  "03_Results/02_Analysis/checkpoints/gsva_module_scores.rds")
meta_path <- file.path(PROJECT_ROOT,
  "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv")
out_dir   <- file.path(PROJECT_ROOT,
  "03_Results/02_Analysis/Plots/Supplementary_10")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv   <- file.path(out_dir, "replicate_level_gsva_per_sample.csv")

message("[load] ", rds_path)
obj <- readRDS(rds_path)
stopifnot(is.list(obj) && "scores" %in% names(obj))
scores <- obj$scores
message("  scores matrix: ", nrow(scores), " modules x ", ncol(scores), " samples")

# Wide -> long: rows are modules, columns are samples.
long <- as.data.frame(scores) |>
  tibble::rownames_to_column("module") |>
  pivot_longer(-module, names_to = "sample", values_to = "gsva_score")

message("[load] ", meta_path)
meta <- read.csv(meta_path, sep = ";", stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8-BOM")
meta <- meta[, c("sample", "genotype", "days", "rep")]
# Normalise the genotype label so it matches the plotting palette ("Ctrl").
meta$genotype[meta$genotype == "Control"] <- "Ctrl"

out <- long |> left_join(meta, by = "sample")
missing_meta <- sum(is.na(out$genotype))
if (missing_meta > 0L) {
  warning("Samples without metadata: ", missing_meta)
}

# Stable ordering for downstream readers.
out <- out[order(out$module, out$days, out$genotype, out$sample),
           c("sample", "genotype", "days", "rep", "module", "gsva_score")]

write.csv(out, out_csv, row.names = FALSE)
message("[write] ", out_csv)
message("  rows: ", nrow(out))
