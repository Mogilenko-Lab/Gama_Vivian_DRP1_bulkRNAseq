###############################################################################
##  Export per-sample GSVA scores in long format                              ##
###############################################################################
##                                                                            ##
##  PURPOSE: Melt the per-sample GSVA matrix in gsva_all_pathways.rds         ##
##           into long format so the interactive bump-chart dashboard         ##
##           (Phase 2b) can embed per-replicate enrichment scores under       ##
##           each pathway_id.                                                 ##
##                                                                            ##
##  INPUTS:  03_Results/02_Analysis/checkpoints/gsva_all_pathways.rds         ##
##           03_Results/02_Analysis/checkpoints/qc_variables.rds (annot)      ##
##                                                                            ##
##  OUTPUT:  03_Results/02_Analysis/replicate_level_gsva_long.csv             ##
##           Columns: pathway_id, sample_id, genotype, day, gsva_score        ##
##           (genotype recoded Control -> Ctrl to match GSEA naming)          ##
##                                                                            ##
###############################################################################

suppressMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

message("Loading GSVA checkpoint and sample annotation...")

gsva_ck <- readRDS(here("03_Results/02_Analysis/checkpoints/gsva_all_pathways.rds"))
qc <- readRDS(here("03_Results/02_Analysis/checkpoints/qc_variables.rds"))

scores <- gsva_ck$scores
annot  <- as.data.frame(qc$annot) %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(genotype = gsub("Control", "Ctrl", genotype)) %>%
  rename(day = days) %>%
  select(sample_id, genotype, day)

message(sprintf("  scores: %d pathways x %d samples", nrow(scores), ncol(scores)))
message(sprintf("  annot:  %d samples", nrow(annot)))

# Sanity: every column in scores must be in annot
missing_samples <- setdiff(colnames(scores), annot$sample_id)
if (length(missing_samples) > 0) {
  stop("Samples in GSVA scores not found in annot: ",
       paste(missing_samples, collapse = ", "))
}

# Melt to long format
long <- as.data.frame(scores) %>%
  tibble::rownames_to_column("pathway_id") %>%
  tidyr::pivot_longer(-pathway_id,
                      names_to  = "sample_id",
                      values_to = "gsva_score") %>%
  left_join(annot, by = "sample_id") %>%
  select(pathway_id, sample_id, genotype, day, gsva_score)

message(sprintf("  long-format rows: %d", nrow(long)))
message(sprintf("  unique pathways:  %d", dplyr::n_distinct(long$pathway_id)))
message(sprintf("  unique samples:   %d", dplyr::n_distinct(long$sample_id)))

# Join sanity check against GSEA universe
master_gsea_ids <- readr::read_csv(
  here("03_Results/02_Analysis/master_gsea_table.csv"),
  show_col_types = FALSE,
  col_select = "pathway_id"
)$pathway_id |> unique()

gsva_ids <- unique(long$pathway_id)
in_both    <- length(intersect(master_gsea_ids, gsva_ids))
gsea_only  <- setdiff(master_gsea_ids, gsva_ids)
gsva_only  <- setdiff(gsva_ids, master_gsea_ids)

message("\nPathway-ID join sanity check (GSEA vs GSVA):")
message(sprintf("  GSEA universe:        %d", length(master_gsea_ids)))
message(sprintf("  GSVA universe:        %d", length(gsva_ids)))
message(sprintf("  Intersection:         %d", in_both))
message(sprintf("  GSEA-only (no GSVA):  %d", length(gsea_only)))
message(sprintf("  GSVA-only (no GSEA):  %d", length(gsva_only)))

# Write
out_file <- here("03_Results/02_Analysis/replicate_level_gsva_long.csv")
readr::write_csv(long, out_file)
message(sprintf("\nWrote: %s (%.1f MB)",
                out_file,
                file.info(out_file)$size / 1024^2))
