#!/usr/bin/env Rscript
## 5b.extract_maturation_DEGs.R
##
## Extracts the union of Time_G32A and Time_R403C maturation DEGs at
## FDR < 0.05, labels each gene as shared / G32A_only / R403C_only, and
## emits per-contrast logFC + FDR for both mutations on every row.
##
## Output: 03_Results/02_Analysis/Tables/maturation_DEGs_G32A_R403C.csv

suppressPackageStartupMessages(library(here))

de_dir  <- here::here("03_Results/02_Analysis/DE_results")
out_dir <- here::here("03_Results/02_Analysis/Tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

g32a  <- read.csv(file.path(de_dir, "Time_G32A_DE_results.csv"),  row.names = 1)
r403c <- read.csv(file.path(de_dir, "Time_R403C_DE_results.csv"), row.names = 1)

sig_at <- function(df, alpha) rownames(df)[df$adj.P.Val < alpha & !is.na(df$adj.P.Val)]

set_g32a  <- sig_at(g32a,  0.05)
set_r403c <- sig_at(r403c, 0.05)

shared     <- intersect(set_g32a, set_r403c)
g32a_only  <- setdiff(set_g32a,   set_r403c)
r403c_only <- setdiff(set_r403c,  set_g32a)
union_set  <- union(set_g32a, set_r403c)

membership_lookup <- setNames(
  c(rep("shared",     length(shared)),
    rep("G32A_only",  length(g32a_only)),
    rep("R403C_only", length(r403c_only))),
  c(shared, g32a_only, r403c_only)
)

direction <- function(x) ifelse(is.na(x), NA_character_, ifelse(x > 0, "Up", "Down"))

res <- data.frame(
  Gene        = union_set,
  Membership  = membership_lookup[union_set],
  G32A_logFC  = round(g32a[union_set,  "logFC"],     3),
  G32A_FDR    = signif(g32a[union_set,  "adj.P.Val"], 3),
  G32A_dir    = direction(g32a[union_set,  "logFC"]),
  R403C_logFC = round(r403c[union_set, "logFC"],     3),
  R403C_FDR   = signif(r403c[union_set, "adj.P.Val"], 3),
  R403C_dir   = direction(r403c[union_set, "logFC"]),
  stringsAsFactors = FALSE
)
res$Mean_logFC <- rowMeans(res[, c("G32A_logFC", "R403C_logFC")], na.rm = TRUE)

membership_order <- factor(res$Membership, levels = c("shared", "G32A_only", "R403C_only"))
res <- res[order(membership_order, -res$Mean_logFC), ]

out_file <- file.path(out_dir, "maturation_DEGs_G32A_R403C.csv")
write.csv(res, out_file, row.names = FALSE)

message(sprintf(
  "Saved %d genes (shared=%d, G32A_only=%d, R403C_only=%d) to %s",
  nrow(res), length(shared), length(g32a_only), length(r403c_only), out_file
))
