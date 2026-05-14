#!/usr/bin/env Rscript
## 5b.extract_shared_maturation_degs.R
##
## Extracts the 38 genes that are differentially expressed (FDR < 0.05) 
## over time in BOTH G32A and R403C mutants (Time_G32A and Time_R403C).
##
## Saves the result to 03_Results/02_Analysis/Tables/shared_DEGs_G32A_R403C.csv

suppressPackageStartupMessages(library(here))

de_dir <- here::here("03_Results/02_Analysis/DE_results")
out_dir <- here::here("03_Results/02_Analysis/Tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

g32a  <- read.csv(file.path(de_dir, "Time_G32A_results.csv"),  row.names = 1)
r403c <- read.csv(file.path(de_dir, "Time_R403C_results.csv"), row.names = 1)

sig_at <- function(df, alpha) rownames(df)[df$adj.P.Val < alpha & !is.na(df$adj.P.Val)]

set_g32a  <- sig_at(g32a,  0.05)
set_r403c <- sig_at(r403c, 0.05)

shared_genes <- intersect(set_g32a, set_r403c)

# Build the summary table
res <- data.frame(
  Gene = shared_genes,
  G32A_logFC = round(g32a[shared_genes, "logFC"], 3),
  G32A_FDR = g32a[shared_genes, "adj.P.Val"],
  G32A_dir = ifelse(g32a[shared_genes, "logFC"] > 0, "Up", "Down"),
  R403C_logFC = round(r403c[shared_genes, "logFC"], 3),
  R403C_FDR = r403c[shared_genes, "adj.P.Val"],
  R403C_dir = ifelse(r403c[shared_genes, "logFC"] > 0, "Up", "Down"),
  stringsAsFactors = FALSE
)
res$Mean_logFC <- (res$G32A_logFC + res$R403C_logFC) / 2

# Sort by Mean_logFC descending (most upregulated to most downregulated)
res <- res[order(res$Mean_logFC, decreasing = TRUE), ]

out_file <- file.path(out_dir, "shared_DEGs_G32A_R403C.csv")
write.csv(res, out_file, row.names = FALSE)
message("Saved ", length(shared_genes), " shared DEGs to ", out_file)
