#!/usr/bin/env Rscript
#' Generate Rasterized FDR Volcano Plots
#'
#' Creates volcano plots with rasterized points but vector text and plot elements.
#' This allows for easy editing in vector graphics editors while handling plots
#' with many data points efficiently.
#'
#' @description
#' This script:
#' 1. Loads the contrast tables from the main pipeline
#' 2. Creates a modified version of create_standard_volcano that uses ggrastr
#' 3. Generates all FDR volcano plots with rasterized points
#' 4. Saves them to the fdr_raster subdirectory

# Setup ----
library(here)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)  # For rasterization

# Source the toolkit functions
source(here("01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R"))
source(here("01_Scripts/RNAseq-toolkit/scripts/utils_plotting.R"))

# Configuration (matching main pipeline)
config <- list(
  out_root      = "03_Results/02_Analysis",
  fc_cutoff     = 2,         # |log2FC| >= 2 (4-fold change) for volcano visualization
  calcium_genes = c(
      "NNAT","CACNG3","CACNA1S","ATP2A1",
      "RYR1","MYLK3","VDR","STIM1","STIM2",
      "ORAI1_1","CALB1","CALR","PNPO")
)

# Load contrast tables from checkpoint
message("Loading contrast tables...")
contrast_tables <- readRDS(here(config$out_root, "checkpoints/contrast_tables.rds"))
message("✓ Loaded ", length(contrast_tables), " contrasts")

# Create rasterized volcano function ----
#' @description
#' This is a modified version of create_standard_volcano that uses
#' ggrastr::geom_point_rast() instead of ggplot2::geom_point()
#' to rasterize only the points while keeping all other elements as vectors.
#'
#' @param raster_dpi Resolution for rasterized points (default 300)
create_standard_volcano_raster <- function(
    de_results,
    decision_by   = c("fdr", "p"),
    p_cutoff      = 0.05,
    fc_cutoff     = 2,
    top_n         = 5,
    highlight_gene= NULL,
    label_method  = "top",
    x_breaks      = 1,
    title         = "Volcano plot",
    subtitle      = NULL,
    caption       = NULL,
    color_palette = c(
      "NS"               = "#7F7F7F",   # grey
      "Log2FC"           = "#0173B2",   # blue
      "p-value"          = "#029E73",   # green
      "p-value & Log2FC" = "#D55E00"    # orange
    ),
    show_grid     = FALSE,
    max.overlaps  = 10,
    raster_dpi    = 300,
    ...
) {
  decision_by <- match.arg(decision_by)

  # ───────────────────────────────── helpers ──────────────────────────
  shade <- function(hex, factor = .6) {
    rgb <- grDevices::col2rgb(hex)/255 * factor
    grDevices::rgb(pmax(pmin(rgb,1),0)[1],
                   pmax(pmin(rgb,1),0)[2],
                   pmax(pmin(rgb,1),0)[3])
  }
  text_col <- function(hex) {
    lum <- sum(grDevices::col2rgb(hex) * c(0.299,0.587,0.114))/255
    ifelse(lum > .55, "black", "white")
  }
  custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
    source("scripts/custom_minimal_theme.R", local = TRUE)
    custom_minimal_theme_with_grid
  } else {
    function() ggplot2::theme_minimal()
  }

  # ──────────────────────── sanity checks ────────────────────────────
  stopifnot(all(c("logFC","P.Value","adj.P.Val") %in% colnames(de_results)))

  if (isTRUE(list(...)[["use_fdr"]])) {
    warning("`use_fdr` is deprecated: please use `decision_by = \"fdr\"` instead.")
  }

  # ────────────────── 1. annotate significance ───────────────────────
  if (decision_by == "fdr") {
    sig_stat   <- de_results$adj.P.Val
    stat_name  <- "FDR"
    sig_logic  <- sig_stat <= p_cutoff

    sig_pvals <- de_results$P.Value[sig_logic]
    if (length(sig_pvals) > 0) {
      p_thresh <- max(sig_pvals, na.rm = TRUE)
      horiz_line <- -log10(p_thresh)
      draw_horiz_line <- TRUE
    } else {
      horiz_line <- NA
      draw_horiz_line <- FALSE
    }
    legend_sig <- sprintf("FDR ≤ %.2g", p_cutoff)
  } else {  # decision_by == "p"
    sig_stat   <- de_results$P.Value
    stat_name  <- "p-value"
    sig_logic  <- sig_stat <= p_cutoff
    horiz_line <- -log10(p_cutoff)
    draw_horiz_line <- TRUE
    legend_sig <- sprintf("p ≤ %.2g", p_cutoff)
  }

  df <- dplyr::mutate(de_results,
    sig_fc  = abs(logFC) >= fc_cutoff,
    sig_dec = sig_logic,
    significance_value = sig_stat,
    cat = dplyr::case_when(sig_fc & sig_dec ~ "p-value & Log2FC",
                         sig_fc          ~ "Log2FC",
                         sig_dec         ~ "p-value",
                         TRUE            ~ "NS"))

  # ────────────────── 2. label selection ─────────────────────────────
  get_top <- function(side) {
    if (side == "up") {
      df |>
        dplyr::filter(.data$logFC > 0) |>
        dplyr::arrange(.data$significance_value) |>
        dplyr::slice_head(n = top_n)
    } else {
      df |>
        dplyr::filter(.data$logFC < 0) |>
        dplyr::arrange(.data$significance_value) |>
        dplyr::slice_head(n = top_n)
    }
  }

  if (label_method == "top") {
    lab_df <- dplyr::bind_rows(get_top("up"), get_top("down"))
  } else if (label_method == "sig") {
    lab_df <- df[df$cat == "p-value & Log2FC", ]
  } else if (label_method == "p") {
    lab_df <- df[df$sig_dec, ]
  } else if (label_method == "log2fc") {
    lab_df <- df[df$sig_fc, ]
  } else {
    lab_df <- df[0, ]
  }

  if (!is.null(highlight_gene)) {
    lab_df <- dplyr::bind_rows(lab_df,
                               df[rownames(df) %in% highlight_gene, ]) |>
              dplyr::distinct()
  }

  # ────────────────── 3. axis limits & colours ───────────────────────
  xmax <- ceiling(max(abs(df$logFC))/x_breaks)*x_breaks
  ymax <- ceiling(max(-log10(df$P.Value)))
  dark_pal <- vapply(color_palette, shade, character(1))

  # ────────────────── 4. build ggplot with RASTERIZED points ─────────
  # KEY CHANGE: Using ggrastr::geom_point_rast instead of geom_point
  g <- ggplot2::ggplot(df, ggplot2::aes(logFC, -log10(P.Value), colour = cat)) +
       ggrastr::geom_point_rast(size = 2, alpha = .65, raster.dpi = raster_dpi) +
       ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       ggplot2::scale_colour_manual(name = NULL,
         values = color_palette,
         breaks = names(color_palette),
         labels = c(
           "p-value & Log2FC" = sprintf("%s & |log2FC| ≥ %.1f", legend_sig, fc_cutoff),
           "Log2FC"           = sprintf("|log2FC| ≥ %.1f", fc_cutoff),
           "p-value"          = legend_sig,
           "NS"               = "NS")) +
       ggplot2::scale_x_continuous(breaks = seq(-xmax, xmax, by = x_breaks),
                                   limits = c(-xmax, xmax)) +
       ggplot2::coord_cartesian(ylim = c(0, ymax)) +
       ggplot2::labs(x = "log2(FC)",
                     y = expression(-log[10](p-value)),
                     title = title,
                     subtitle = subtitle,
                     caption = if (is.null(caption)) {
                       if (decision_by == "fdr") {
                         if (draw_horiz_line) {
                           sprintf("Dashed lines: horiz. – FDR ≤ %.2g (p ≤ %.2g); vert. – |log2FC| ≥ %.1f",
                                   p_cutoff, signif(10^(-horiz_line),2), fc_cutoff)
                         } else {
                           sprintf("No genes pass FDR ≤ %.2g. Dashed lines: vert. – |log2FC| ≥ %.1f",
                                   p_cutoff, fc_cutoff)
                         }
                       } else {
                         sprintf("Dashed lines: horiz. – p ≤ %.2g; vert. – |log2FC| ≥ %.1f",
                                 p_cutoff, fc_cutoff)
                       }
                     } else caption) +
       custom_minimal_theme_with_grid()

  # Add horizontal line only if there are significant genes
  if (draw_horiz_line) {
    g <- g + ggplot2::geom_hline(yintercept = horiz_line, linetype = "dashed")
  } else if (decision_by == "fdr") {
    g <- g + ggplot2::annotate("text",
                               x = xmax * 0.5,
                               y = ymax * 0.95,
                               label = sprintf("No genes pass FDR ≤ %.2g", p_cutoff),
                               size = 4,
                               color = "darkred",
                               fontface = "italic")
  }

  if (!show_grid) {
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank())
  }

  # ────────────────── 5. labels (VECTOR) ──────────────────────────────
  if (nrow(lab_df)) {
    if (!is.null(highlight_gene) && any(rownames(lab_df) %in% highlight_gene)) {
      calcium_labs <- lab_df[rownames(lab_df) %in% highlight_gene, , drop = FALSE]
      regular_labs <- lab_df[!rownames(lab_df) %in% highlight_gene, , drop = FALSE]
    } else {
      calcium_labs <- lab_df[0, , drop = FALSE]
      regular_labs <- lab_df
    }

    # Layer 1: Regular labels
    if (nrow(regular_labs) > 0) {
      g <- g + ggrepel::geom_text_repel(
        data            = regular_labs,
        ggplot2::aes(label = rownames(regular_labs)),
        colour          = dark_pal[regular_labs$cat],
        fontface        = "plain",
        size            = 3.5,
        box.padding     = .4,
        point.padding   = .3,
        max.overlaps    = max.overlaps,
        min.segment.length = 0,
        show.legend     = FALSE)
    }

    # Layer 2: Highlighted genes
    if (nrow(calcium_labs) > 0) {
      g <- g + ggrepel::geom_text_repel(
        data            = calcium_labs,
        ggplot2::aes(label = rownames(calcium_labs)),
        colour          = "black",
        fontface        = "bold",
        size            = 3.5,
        box.padding     = .5,
        point.padding   = .3,
        max.overlaps    = Inf,
        force           = 5,
        min.segment.length = 0,
        show.legend     = FALSE)
    }
  }

  return(g)
}

# Generate rasterized FDR volcano plots ----
message("\n🌋 Generating rasterized FDR volcano plots...")

# Create output directory
fdr_raster_dir <- here(config$out_root, "Plots/Volcano/fdr_raster")
dir.create(fdr_raster_dir, recursive = TRUE, showWarnings = FALSE)

# FDR parameters
fdr_params <- list(
  p_cutoff = 0.1,
  decision_by = "fdr",
  label_method = "top",
  title_suffix = "(Threshold: FDR ≤ 0.1)",
  subtitle = "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
)

# Generate plots
n_plots <- 0
for (co in names(contrast_tables)) {
  message("  Processing: ", co)

  plot <- create_standard_volcano_raster(
    contrast_tables[[co]],
    p_cutoff = fdr_params$p_cutoff,
    fc_cutoff = config$fc_cutoff,
    decision_by = fdr_params$decision_by,
    label_method = fdr_params$label_method,
    highlight_gene = config$calcium_genes,
    x_breaks = 2,
    title = paste(co, fdr_params$title_suffix),
    subtitle = fdr_params$subtitle,
    raster_dpi = 300  # High quality rasterization
  )

  # Save as PDF with rasterized points
  save_plot(plot,
            file.path(fdr_raster_dir, paste0(co, "_standard.pdf")),
            width = 8, height = 7)

  n_plots <- n_plots + 1
}

message("✓ Saved ", n_plots, " rasterized FDR volcano plots to:")
message("  ", fdr_raster_dir)
message("\nNOTE: These PDFs have:")
message("  - Rasterized points (300 DPI) for efficient rendering")
message("  - Vector text, axes, and plot elements for easy editing")
