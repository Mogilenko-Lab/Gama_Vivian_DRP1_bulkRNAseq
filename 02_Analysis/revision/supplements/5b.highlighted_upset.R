#!/usr/bin/env Rscript
## 5b.highlighted_upset.R
##
## Concern: add a highlighted bar to the nine-contrast UpSet
## plot (Supplementary Fig. S6) AND a maturation-only three-contrast UpSet
## as an alternative view. Both outputs go to a NEW directory — the original
## UpSet_plot_all_contrasts.pdf is NOT modified.
##
## Reads checkpoints only; does NOT rerun the pipeline.
##
## Run:
##   Rscript 02_Analysis/5b.highlighted_upset.R

suppressPackageStartupMessages({
  library(UpSetR)
  library(grid)
})

project_root <- "/workspaces/Gama_Vivian_DRP1_bulkRNAseq"
ck_dir       <- file.path(project_root, "03_Results/02_Analysis/checkpoints")
out_dir      <- file.path(project_root,
                          "03_Results/02_Analysis/Plots/General/Supplementary_5b")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load cached fit + decideTests matrix -----------------------------------
fit         <- readRDS(file.path(ck_dir, "fit_object.rds"))
de_results  <- readRDS(file.path(ck_dir, "de_results.rds"))

# Build significant-gene lists per contrast (mirrors 1.1.main_pipeline.R:482-492)
contrast_names <- colnames(de_results)
sig_genes_list <- list()
for (co in contrast_names) {
  sig_genes_list[[co]] <- rownames(fit$coefficients)[which(de_results[, co] != 0)]
}

highlight_color <- "#E08214"  # Okabe-Ito orange, matches project palette

# --- Option 1: nine-contrast UpSet with highlighted G32A∩R403C bar ----------
upset_all <- upset(
  fromList(sig_genes_list),
  order.by = "freq",
  nsets    = length(sig_genes_list),
  queries  = list(
    list(query  = intersects,
         params = list("Time_G32A", "Time_R403C"),
         color  = highlight_color,
         active = TRUE,
         query.name = "Mutant-shared maturation DEGs")
  ),
  query.legend = "top",
  mainbar.y.label = "Intersection size (DEGs)",
  sets.x.label    = "DEGs per contrast"
)

grDevices::pdf(file.path(out_dir, "UpSet_plot_all_contrasts_highlighted.pdf"),
               width = 12, height = 8)
print(upset_all)
grDevices::dev.off()

grDevices::png(file.path(out_dir, "UpSet_plot_all_contrasts_highlighted.png"),
               width = 12, height = 8, units = "in", res = 300)
print(upset_all)
grDevices::dev.off()
message("Nine-contrast highlighted UpSet written.")

# --- Option 2: maturation-only three-contrast UpSet -------------------------
# If any of the three maturation contrasts is missing we skip the plot.
mat_contrasts <- c("Time_Ctrl", "Time_G32A", "Time_R403C")
if (all(mat_contrasts %in% names(sig_genes_list))) {
  mat_sets <- sig_genes_list[mat_contrasts]
  # Drop 0-size sets to avoid UpSetR "empty input" crash
  mat_sets <- mat_sets[sapply(mat_sets, length) > 0]

  if (length(mat_sets) >= 2) {
    upset_mat <- upset(
      fromList(mat_sets),
      order.by = "freq",
      nsets    = length(mat_sets),
      queries  = list(
        list(query  = intersects,
             params = list("Time_G32A", "Time_R403C"),
             color  = highlight_color, active = TRUE,
             query.name = "Mutant-shared maturation DEGs")
      ),
      query.legend    = "top",
      mainbar.y.label = "Intersection size (DEGs)",
      sets.x.label    = "DEGs per contrast"
    )

    grDevices::pdf(file.path(out_dir, "UpSet_plot_maturation_contrasts.pdf"),
                   width = 7, height = 5)
    print(upset_mat)
    # Annotate the missing Time_Ctrl visually (0 DEGs)
    grid.text(
      "Time_Ctrl has 0 DEGs at FDR < 0.05 (not shown as a set)",
      x = 0.5, y = 0.02,
      gp = gpar(col = "grey30", fontsize = 8, fontface = "italic")
    )
    grDevices::dev.off()

    grDevices::png(file.path(out_dir, "UpSet_plot_maturation_contrasts.png"),
                   width = 7, height = 5, units = "in", res = 300)
    print(upset_mat)
    grid.text(
      "Time_Ctrl has 0 DEGs at FDR < 0.05 (not shown as a set)",
      x = 0.5, y = 0.02,
      gp = gpar(col = "grey30", fontsize = 8, fontface = "italic")
    )
    grDevices::dev.off()
    message("Maturation-only UpSet (highlighted) written.")
  } else {
    warning("Fewer than 2 non-empty maturation contrasts; skipping three-set UpSet.")
  }
}

# --- QC: print intersection sizes to stdout ---------------------------------
cat("\n--- Intersection sizes (FDR < 0.05 decideTests matrix) ---\n")
cat("Time_G32A  DEGs: ",  length(sig_genes_list[["Time_G32A"]]),  "\n", sep = "")
cat("Time_R403C DEGs: ",  length(sig_genes_list[["Time_R403C"]]), "\n", sep = "")
cat("Time_Ctrl  DEGs: ",  length(sig_genes_list[["Time_Ctrl"]]),  "\n", sep = "")
cat("G32A ∩ R403C:   ",
    length(intersect(sig_genes_list[["Time_G32A"]],
                     sig_genes_list[["Time_R403C"]])), "\n", sep = "")

cat("\nOutput directory:\n  ", out_dir, "\n", sep = "")
