###############################################################################
## 5a.filtered_volcano_supplement.R
##
## Concern: "apply an expression filter (e.g. avg. expression > 0) to further
## highlight the significance of NNAT and other well-expressed genes in Fig. 4a."
##
## This script is ADDITIVE. It produces *new* supplementary volcano PDFs and an
## AveExpr-vs-significance scatter, without modifying any existing artifact.
##
## Inputs  (all already cached):
##   - 03_Results/02_Analysis/checkpoints/contrast_tables.rds
##   - 03_Results/02_Analysis/DE_results/*.csv           (for Task-1 diagnostics)
##
## Outputs  (all in NEW locations):
##   - 03_Results/02_Analysis/Supplementary/5a_aveexpr_diagnostics.csv
##   - 03_Results/02_Analysis/Supplementary/5a_aveexpr_calcium_genes.csv
##   - 03_Results/02_Analysis/Plots/Volcano_Supplementary_MinAveExpr0/
##        * vertical_fdr_calcium/
##        * vertical_p_calcium/
##        * vertical_fdr/
##        * vertical_p/
##        * horizontal_p/
##        * horizontal_fdr/
##        * scatter_aveexpr_significance/
##
## Visual-only filter: AveExpr > 0 log2-CPM (approximately 1 CPM linear).
## The underlying DE statistics, CSVs, GSEA, and checkpoints are untouched.
##
## Toolkit helpers (plot_standard_volcano.R, volcano_helpers.R) are SOURCED,
## not modified; the AveExpr filter is applied to the contrast tables upstream
## of the helpers (data is mutated before being passed in).
##
## Run from project root:
##     Rscript 02_Analysis/5a.filtered_volcano_supplement.R
###############################################################################

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

setwd(here::here())

# -------------------------------------------------------------------- #
# Config (mirrors config in 1.1.main_pipeline.R lines 16-39)           #
# -------------------------------------------------------------------- #
config <- list(
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",
  fdr_cutoff    = 0.05,
  p_cutoff      = 0.05,
  fc_cutoff     = 2,
  # Exactly as in 1.1.main_pipeline.R lines 27-30
  calcium_genes = c(
    "NNAT","CACNG3","CACNA1S","ATP2A1",
    "RYR1","MYLK3","VDR","STIM1","STIM2",
    "ORAI1_1","CALB1","CALR","PNPO"
  ),
  # The AveExpr visual cut-off (log2-CPM)
  min_ave_expr = 0
)

# Full curated calcium list including genes confirmed absent / outside pipeline
calcium_full <- c(
  "NNAT","CACNG3","CACNA1C","CACNA1S","ATP2A1","RYR1","MYLK3","CASR","VDR",
  "STIM1","STIM2","ORAI1","ORAI1_1","CALB1","CALR","PNPO"
)

# -------------------------------------------------------------------- #
# Source toolkit helpers (READ-ONLY; not modified)                     #
# -------------------------------------------------------------------- #
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found -> ", path)
  }
}

source_if_present(config$helper_root, "scripts/DE/plot_standard_volcano.R")
source_if_present(config$helper_root, "scripts/DE/volcano_helpers.R")
source_if_present(config$helper_root, "scripts/utils_plotting.R")

# -------------------------------------------------------------------- #
# TASK 1 - Evidence gathering (AveExpr diagnostics)                    #
# -------------------------------------------------------------------- #
message("\nTask 1: AveExpr diagnostics from existing DE CSVs ...")

de_dir <- file.path(config$out_root, "DE_results")
supp_dir <- file.path(config$out_root, "Supplementary")
ensure_dir(supp_dir)

contrast_files <- list(
  G32A_vs_Ctrl_D35          = "G32A_vs_Ctrl_D35_results.csv",
  R403C_vs_Ctrl_D35         = "R403C_vs_Ctrl_D35_results.csv",
  G32A_vs_Ctrl_D65          = "G32A_vs_Ctrl_D65_results.csv",
  R403C_vs_Ctrl_D65         = "R403C_vs_Ctrl_D65_results.csv",
  Time_Ctrl                 = "Time_Ctrl_results.csv",
  Time_G32A                 = "Time_G32A_results.csv",
  Time_R403C                = "Time_R403C_results.csv",
  Maturation_G32A_specific  = "Maturation_G32A_specific_results.csv",
  Maturation_R403C_specific = "Maturation_R403C_specific_results.csv"
)

# Read CSVs for diagnostics (doesn't replace checkpoint-based contrast_tables)
csv_tables <- list()
for (co in names(contrast_files)) {
  path <- file.path(de_dir, contrast_files[[co]])
  csv_tables[[co]] <- read.csv(path, row.names = 1, check.names = FALSE,
                               stringsAsFactors = FALSE)
}

summary_rows <- list()
for (co in names(csv_tables)) {
  df <- csv_tables[[co]]
  n_total  <- nrow(df)
  n_lt_0   <- sum(df$AveExpr < 0, na.rm = TRUE)
  n_lt_0p5 <- sum(df$AveExpr < 0.5, na.rm = TRUE)
  n_lt_1   <- sum(df$AveExpr < 1, na.rm = TRUE)
  n_fdr05  <- sum(df$adj.P.Val < 0.05, na.rm = TRUE)
  sig      <- df[df$adj.P.Val < 0.05 & !is.na(df$adj.P.Val), , drop = FALSE]
  n_fdr05_lowexpr <- sum(sig$AveExpr <= 0, na.rm = TRUE)
  ave_range <- range(df$AveExpr, na.rm = TRUE)
  summary_rows[[co]] <- data.frame(
    contrast = co,
    n_total = n_total,
    n_AveExpr_lt_0 = n_lt_0,
    pct_AveExpr_lt_0 = round(100 * n_lt_0 / n_total, 2),
    n_AveExpr_lt_0p5 = n_lt_0p5,
    n_AveExpr_lt_1 = n_lt_1,
    AveExpr_min = round(ave_range[1], 3),
    AveExpr_max = round(ave_range[2], 3),
    n_sig_FDR05 = n_fdr05,
    n_sig_FDR05_AveExpr_le_0 = n_fdr05_lowexpr,
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# Calcium-gene AveExpr diagnostics (using the full curated list)
calcium_rows <- list()
for (co in names(csv_tables)) {
  df <- csv_tables[[co]]
  hit <- intersect(calcium_full, rownames(df))
  if (length(hit) == 0) next
  for (g in hit) {
    calcium_rows[[paste(co, g, sep = "::")]] <- data.frame(
      contrast = co,
      gene = g,
      AveExpr = round(df[g, "AveExpr"], 3),
      logFC = round(df[g, "logFC"], 3),
      P.Value = signif(df[g, "P.Value"], 3),
      adj.P.Val = signif(df[g, "adj.P.Val"], 3),
      passes_AveExpr_gt_0 = df[g, "AveExpr"] > 0,
      stringsAsFactors = FALSE
    )
  }
}
calcium_df <- do.call(rbind, calcium_rows)
rownames(calcium_df) <- NULL

write.csv(summary_df,
          file.path(supp_dir, "5a_aveexpr_diagnostics.csv"),
          row.names = FALSE)
write.csv(calcium_df,
          file.path(supp_dir, "5a_aveexpr_calcium_genes.csv"),
          row.names = FALSE)

message("  Wrote 5a_aveexpr_diagnostics.csv")
message("  Wrote 5a_aveexpr_calcium_genes.csv")

cat("\nPer-contrast AveExpr distribution summary:\n")
print(summary_df, row.names = FALSE)

cat("\nCalcium genes that would be removed by AveExpr > 0 filter (if any):\n")
print(calcium_df[!calcium_df$passes_AveExpr_gt_0, ], row.names = FALSE)

# -------------------------------------------------------------------- #
# Load cached contrast_tables for volcano/scatter generation           #
# -------------------------------------------------------------------- #
message("\nLoading cached contrast_tables.rds ...")
checkpoint_file <- file.path(config$out_root, "checkpoints", "contrast_tables.rds")
if (!file.exists(checkpoint_file)) {
  stop("Missing ", checkpoint_file, "\nRun 02_Analysis/1.1.main_pipeline.R first.")
}
contrast_tables <- readRDS(checkpoint_file)
message("  Loaded ", length(contrast_tables), " contrasts.")

# -------------------------------------------------------------------- #
# PRE-CALCULATION: FDR Boundaries from UNFILTERED data                 #
# -------------------------------------------------------------------- #
# Problem: When we filter for AveExpr > 0, we often remove the genes
# that sit right at the FDR boundary. This makes the horizontal line in
# a volcano plot (which marks the FDR cutoff on a P-value axis) jump 
# to the next available significant gene, which looks inconsistent.
#
# Solution: We calculate the raw P-value threshold that corresponds to
# FDR <= 0.1 on the FULL dataset. We will then pass this "fixed" boundary
# to the plotting functions for the filtered data.
# -------------------------------------------------------------------- #
message("Calculating stable FDR boundaries from full dataset ...")
fdr_boundaries <- list()
for (co in names(contrast_tables)) {
  tbl <- contrast_tables[[co]]
  # Find the largest raw P.Value that satisfies adj.P.Val <= 0.1
  sig_pvals <- tbl$P.Value[tbl$adj.P.Val <= 0.1 & !is.na(tbl$adj.P.Val)]
  if (length(sig_pvals) > 0) {
    fdr_boundaries[[co]] <- max(sig_pvals)
  } else {
    fdr_boundaries[[co]] <- NULL
  }
}

# -------------------------------------------------------------------- #
# Upstream filter helper - apply AveExpr > min_ave_expr                #
# -------------------------------------------------------------------- #
filter_by_ave_expr <- function(tbl, cutoff = 0) {
  if (!"AveExpr" %in% colnames(tbl)) {
    warning("AveExpr column missing; returning unfiltered table.")
    return(tbl)
  }
  keep <- tbl$AveExpr > cutoff & !is.na(tbl$AveExpr)
  tbl[keep, , drop = FALSE]
}

filtered_contrast_tables <- lapply(contrast_tables,
                                   filter_by_ave_expr,
                                   cutoff = config$min_ave_expr)

# QC: report filter impact
filter_impact <- data.frame(
  contrast = names(contrast_tables),
  n_before = vapply(contrast_tables, nrow, integer(1)),
  n_after  = vapply(filtered_contrast_tables, nrow, integer(1)),
  stringsAsFactors = FALSE
)
filter_impact$n_removed <- filter_impact$n_before - filter_impact$n_after
filter_impact$pct_removed <- round(100 * filter_impact$n_removed / filter_impact$n_before, 2)
rownames(filter_impact) <- NULL
cat("\nAveExpr > 0 filter impact on contrast_tables checkpoint:\n")
print(filter_impact, row.names = FALSE)

write.csv(filter_impact,
          file.path(supp_dir, "5a_aveexpr_filter_impact.csv"),
          row.names = FALSE)

# QC: confirm retention of NNAT & PNPO everywhere
cat("\nPer-contrast retention check for NNAT & PNPO:\n")
for (co in names(filtered_contrast_tables)) {
  t1 <- filtered_contrast_tables[[co]]
  cat(sprintf("  %-28s NNAT %s  PNPO %s\n",
              co,
              if ("NNAT" %in% rownames(t1)) "YES" else "NO",
              if ("PNPO" %in% rownames(t1)) "YES" else "NO"))
}

# -------------------------------------------------------------------- #
# TASK 2 - Supplementary filtered volcano set                          #
# -------------------------------------------------------------------- #
message("\nTask 2: Generating supplementary volcano plots (AveExpr > 0) ...")

supp_volcano_root <- file.path(config$out_root,
                               "Plots/Volcano_Supplementary_MinAveExpr0")
ensure_dir(supp_volcano_root)

# Helper: save plot both as PDF and PNG @ 300 DPI
save_plot_pdf_png <- function(plot, path_no_ext, width = 7, height = 7) {
  ensure_dir(dirname(path_no_ext))
  ggsave(paste0(path_no_ext, ".pdf"), plot = plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(paste0(path_no_ext, ".png"), plot = plot,
         width = width, height = height, dpi = 300, device = "png")
  invisible(NULL)
}

# Contrast groups mirror 01_Scripts/R_scripts/generate_vertical_volcanos.R
contrast_groups <- list(
  group1 = c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35"),
  group2 = c("G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"),
  group3 = c("G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35",
             "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"),
  group4 = c("Time_Ctrl", "Time_G32A", "Time_R403C"),
  group5 = c("Maturation_G32A_specific", "Maturation_R403C_specific")
)
group_names <- c(
  "D35_comparisons",
  "D65_comparisons",
  "all_disease_vs_control",
  "time_effects",
  "maturation_effects"
)

# ---------- Vertical volcano generator (AveExpr-filtered data) ----------
make_vertical_set <- function(mode, out_dir, highlight_genes = NULL) {
  ensure_dir(out_dir)
  decision_by <- if (mode == "p") "p" else "fdr"
  p_cutoff    <- if (mode == "p") config$p_cutoff else 0.1

  # Group panels
  for (i in seq_along(contrast_groups)) {
    group      <- contrast_groups[[i]]
    group_name <- group_names[i]
    valid      <- intersect(group, names(filtered_contrast_tables))
    if (length(valid) == 0) next

    volcano_list <- list()
    for (contrast in valid) {
      volcano_list[[contrast]] <- create_vertical_volcano(
        filtered_contrast_tables[[contrast]],
        decision_by    = decision_by,
        p_cutoff       = p_cutoff,
        fc_cutoff      = config$fc_cutoff,
        label_method   = "top",
        highlight_gene = highlight_genes,
        title          = contrast
      )
    }
    combined <- combine_volcano_row(volcano_list, keep_first_caption = FALSE)
    width  <- 3 * length(volcano_list)
    height <- 6
    save_plot_pdf_png(combined,
                      file.path(out_dir, paste0(group_name, "_vertical_minAveExpr0")),
                      width = width, height = height)
    # individual contrast panels
    for (contrast in names(volcano_list)) {
      save_plot_pdf_png(volcano_list[[contrast]],
                        file.path(out_dir,
                                  paste0(contrast, "_vertical_minAveExpr0")),
                        width = 6, height = 7)
    }
  }

  # All-contrasts panel
  volcano_list <- list()
  for (contrast in names(filtered_contrast_tables)) {
    volcano_list[[contrast]] <- create_vertical_volcano(
      filtered_contrast_tables[[contrast]],
      decision_by    = decision_by,
      p_cutoff       = p_cutoff,
      fc_cutoff      = config$fc_cutoff,
      label_method   = "top",
      highlight_gene = highlight_genes,
      title          = contrast
    )
  }
  save_plot_pdf_png(
    combine_volcano_row(volcano_list, keep_first_caption = FALSE),
    file.path(out_dir, "all_contrasts_vertical_minAveExpr0"),
    width = 3 * length(volcano_list), height = 6
  )
}

# ---------- Horizontal volcano generator ----------
make_horizontal_set <- function(mode, out_dir) {
  ensure_dir(out_dir)
  decision_by <- if (mode == "p") "p" else "fdr"
  p_cutoff <- if (mode == "p") config$p_cutoff else 0.1
  subtitle <- if (mode == "p") "Highlighting by unadjusted p-value"
              else "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
  title_suffix <- if (mode == "p")
                    sprintf("(Threshold: p <= %.2f; AveExpr > 0)", config$p_cutoff)
                  else
                    "(Threshold: FDR <= 0.1; AveExpr > 0)"

  for (co in names(filtered_contrast_tables)) {
    # Use the stable boundary from the full dataset if in FDR mode
    fixed_p <- if (mode == "fdr") fdr_boundaries[[co]] else NULL

    plot <- create_standard_volcano(
      filtered_contrast_tables[[co]],
      p_cutoff       = p_cutoff,
      fc_cutoff      = config$fc_cutoff,
      decision_by    = decision_by,
      label_method   = "top",
      highlight_gene = config$calcium_genes,
      x_breaks       = 2,
      fixed_p_boundary = fixed_p,
      title          = paste(co, title_suffix),
      subtitle       = subtitle
    )
    save_plot_pdf_png(plot,
                      file.path(out_dir,
                                paste0(co, "_standard_minAveExpr0")),
                      width = 8, height = 7)
  }
}

# ---- Vertical sets (FDR + p; with/without calcium highlight) -----------
message("  vertical_fdr_calcium/")
make_vertical_set("fdr",
                  file.path(supp_volcano_root, "vertical_fdr_calcium"),
                  highlight_genes = config$calcium_genes)

message("  vertical_p_calcium/")
make_vertical_set("p",
                  file.path(supp_volcano_root, "vertical_p_calcium"),
                  highlight_genes = config$calcium_genes)

message("  vertical_fdr/ (no calcium highlight)")
make_vertical_set("fdr",
                  file.path(supp_volcano_root, "vertical_fdr"),
                  highlight_genes = NULL)

message("  vertical_p/ (no calcium highlight)")
make_vertical_set("p",
                  file.path(supp_volcano_root, "vertical_p"),
                  highlight_genes = NULL)

# ---- Horizontal sets (FDR + p) -----------------------------------------
message("  horizontal_fdr/")
make_horizontal_set("fdr", file.path(supp_volcano_root, "horizontal_fdr"))

message("  horizontal_p/")
make_horizontal_set("p", file.path(supp_volcano_root, "horizontal_p"))

message("Supplementary volcano set complete.")

# -------------------------------------------------------------------- #
# TASK 3 - AveExpr vs -log10(FDR) scatter                              #
# -------------------------------------------------------------------- #
message("\nTask 3: AveExpr-vs-significance scatter ...")

scatter_dir <- file.path(supp_volcano_root, "scatter_aveexpr_significance")
ensure_dir(scatter_dir)

scatter_contrasts <- c(
  "G32A_vs_Ctrl_D35",
  "R403C_vs_Ctrl_D35",
  "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D65"
)

make_aveexpr_scatter <- function(de, contrast_label, calcium_genes_vec, fixed_p_boundary = NULL) {
  de$gene <- rownames(de)
  de$highlight <- de$gene %in% calcium_genes_vec
  de$neglog10pval <- -log10(pmax(de$P.Value, .Machine$double.xmin))

  labeled <- de[de$highlight & de$AveExpr > 0, , drop = FALSE]

  p <- ggplot(de, aes(AveExpr, neglog10pval)) +
    geom_point(
      data = de[!de$highlight, ],
      colour = "grey70", alpha = 0.35, size = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
    
  if (!is.null(fixed_p_boundary)) {
    p <- p + geom_hline(yintercept = -log10(fixed_p_boundary), linetype = "dashed", colour = "grey40")
    cap_fdr <- "Dashed horizontal: FDR \u2264 0.1 boundary (raw p-value)."
  } else {
    cap_fdr <- ""
  }

  p <- p +
    geom_point(
      data = de[de$highlight, ],
      colour = "#D55E00", size = 2.5, alpha = 0.95) +
    ggrepel::geom_label_repel(
      data = labeled,
      aes(label = gene),
      size = 3, max.overlaps = 20,
      box.padding = 0.4, point.padding = 0.3,
      fontface = "bold", colour = "black",
      fill = scales::alpha("white", 0.85),
      segment.colour = "grey30") +
    labs(
      x = "Average expression (log2-CPM)",
      y = expression(-log[10](p-value)),
      title = paste0("AveExpr vs significance - ", contrast_label),
      caption = paste0("Dashed vertical: AveExpr = 0. ", cap_fdr, " Calcium genes in orange.")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.4),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  p
}

for (co in scatter_contrasts) {
  if (!co %in% names(contrast_tables)) {
    message("  skipping scatter for ", co, " (not in contrast_tables)")
    next
  }
  p <- make_aveexpr_scatter(contrast_tables[[co]], co, config$calcium_genes, fixed_p_boundary = fdr_boundaries[[co]])
  save_plot_pdf_png(p,
                    file.path(scatter_dir, paste0("aveexpr_scatter_", co)),
                    width = 6, height = 5.5)
  message("  scatter_aveexpr_significance/aveexpr_scatter_", co, ".pdf/.png")
}

# Also build a combined 2x2 panel for completeness
panels <- lapply(scatter_contrasts, function(co) {
  if (!co %in% names(contrast_tables)) return(NULL)
  make_aveexpr_scatter(contrast_tables[[co]], co, config$calcium_genes, fixed_p_boundary = fdr_boundaries[[co]]) +
    theme(plot.caption = element_blank())
})
panels <- panels[!vapply(panels, is.null, logical(1))]
if (length(panels) >= 2) {
  combined <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "AveExpr vs -log10(FDR) - all four mutation-vs-control contrasts",
      caption = "Dashed vertical: AveExpr = 0 (approx 1 CPM). Dotted horizontal: FDR = 0.05. Calcium genes in orange.")
  save_plot_pdf_png(combined,
                    file.path(scatter_dir, "aveexpr_scatter_all_mutation_contrasts_2x2"),
                    width = 12, height = 11)
  message("  aveexpr_scatter_all_mutation_contrasts_2x2.pdf/.png")
}

# -------------------------------------------------------------------- #
# Done                                                                  #
# -------------------------------------------------------------------- #
message("\n5a.filtered_volcano_supplement.R complete.")
message("Outputs written under:")
message("  ", supp_dir)
message("  ", supp_volcano_root)
