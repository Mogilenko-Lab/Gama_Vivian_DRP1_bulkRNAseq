#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Normalised SynGO running-sum plots for Fig. 5D
#
# Builds four GSEA running-sum figures (G32A/R403C × D35/D65) that share:
#   * the same six SynGO synaptic-localized pathways (Tol/Wong colorblind-safe
#     palette, identical to the Fig. 5 chord diagram so the figure family
#     reads as one),
#   * journal-grade typography (larger fonts, thicker enrichment lines and
#     axes) sized for a small printed panel.
#
# Output: 03_Results/02_Analysis/Plots/SynGO_running_sum_normalised/
#   <contrast>_SynGO_running_sum_normalised.pdf  (one per contrast)
#   SynGO_running_sum_normalised_grid.pdf        (2 x 2 composite)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)   # rasterise heavy line / tick / area geoms
})

# Resolution for rasterised running-sum lines and rank ticks. Vector axes,
# fonts, and panel borders are preserved; only the dense per-gene geometry is
# rasterised so PDFs stay light enough to manipulate.
RASTER_DPI <- 350

proj_root <- here::here()

# Unified GSEA running-sum helper (we override line widths / theme below)
source(file.path(
  proj_root,
  "01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R"
))

# -----------------------------------------------------------------------------
# Shared SynGO pathway palette - matched 1:1 to viz_chord_diagrams.py so the
# running-sum panels carry the same color identity as the chord diagrams.
# Palette: Tol muted (synaptic compartments) + Tol vibrant (focal ribosomes).
# -----------------------------------------------------------------------------
PATHWAY_IDS <- c(
  "SYNGO:presyn_ribosome",
  "SYNGO:postsyn_ribosome",
  "GO:0045202",
  "GO:0099523",
  "GO:0099524",
  "GO:0014069",
  "GO:0045211"
)

PATHWAY_COLORS <- c(
  "SYNGO:presyn_ribosome"  = "#CC6677",  # muted rose
  "SYNGO:postsyn_ribosome" = "#AA4499",  # muted purple
  "GO:0045202"             = "#DDCC77",  # sand (SynGO "synapse" compartment)
  "GO:0099523"             = "#117733",  # forest green
  "GO:0099524"             = "#44AA99",  # teal
  "GO:0014069"             = "#999933",  # olive
  "GO:0045211"             = "#88CCEE"   # sky blue
)

PATHWAY_LABELS <- c(
  "SYNGO:presyn_ribosome"  = "Presynaptic ribosome",
  "SYNGO:postsyn_ribosome" = "Postsynaptic ribosome",
  "GO:0045202"             = "Synapse",
  "GO:0099523"             = "Presynaptic cytosol",
  "GO:0099524"             = "Postsynaptic cytosol",
  "GO:0014069"             = "Postsynaptic density",
  "GO:0045211"             = "Postsynaptic membrane"
)

CONTRASTS <- c(
  "G32A_vs_Ctrl_D35",
  "G32A_vs_Ctrl_D65",
  "R403C_vs_Ctrl_D35",
  "R403C_vs_Ctrl_D65"
)

CONTRAST_TITLES <- c(
  "G32A_vs_Ctrl_D35"  = "G32A vs Ctrl - D35",
  "G32A_vs_Ctrl_D65"  = "G32A vs Ctrl - D65",
  "R403C_vs_Ctrl_D35" = "R403C vs Ctrl - D35",
  "R403C_vs_Ctrl_D65" = "R403C vs Ctrl - D65"
)

out_dir <- file.path(
  proj_root,
  "03_Results/02_Analysis/Plots/SynGO_running_sum_normalised"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

syngo <- readRDS(file.path(
  proj_root,
  "03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds"
))

# -----------------------------------------------------------------------------
# First pass: build raw gseaplot2 outputs only so we can read off the data
# ranges and force a SHARED Y-axis scale across all four contrasts. With the
# same limits applied to every panel, the zero line lands at the same vertical
# position in every plot - the visual property the cross-mutation /
# cross-timepoint comparison hinges on.
# -----------------------------------------------------------------------------
collect_layer_vals <- function(p, candidates = c("runningScore", "y", "value")) {
  vals <- numeric()
  for (layer in p$layers) {
    d <- layer$data
    if (is.null(d) || inherits(d, "waiver")) d <- p$data
    if (is.null(d) || !is.data.frame(d)) next
    for (col in candidates) {
      if (col %in% names(d)) {
        v <- suppressWarnings(as.numeric(d[[col]]))
        vals <- c(vals, v[is.finite(v)])
      }
    }
  }
  vals
}

pad_range <- function(r, frac = 0.04) {
  span <- max(diff(r), 1e-6)
  c(r[1] - span * frac, r[2] + span * frac)
}

message("Pre-computing shared Y ranges across the four contrasts")
raw_plots <- list()
for (ct in CONTRASTS) {
  obj <- syngo[[ct]]
  ids <- PATHWAY_IDS[PATHWAY_IDS %in% obj@result$ID]
  palette_ct <- unname(PATHWAY_COLORS[ids])
  raw_plots[[ct]] <- enrichplot::gseaplot2(
    obj,
    geneSetID    = ids,
    subplots     = c(1, 2, 3),
    pvalue_table = FALSE,
    rel_heights  = c(1.6, 0.5, 0.5),
    color        = palette_ct
  )
}

vals_running <- unlist(lapply(raw_plots, function(rp)
  collect_layer_vals(rp[[1]], c("runningScore", "y"))))
vals_ranked  <- unlist(lapply(raw_plots, function(rp)
  collect_layer_vals(rp[[3]], c("y", "value"))))

# Always include 0 in the running-enrichment range so the zero line is on-plot.
YLIM_RUNNING <- pad_range(range(c(0, vals_running), na.rm = TRUE))
YLIM_RANKED  <- pad_range(range(vals_ranked,         na.rm = TRUE))

message(sprintf("  running enrichment limits: [%.2f, %.2f]",
                YLIM_RUNNING[1], YLIM_RUNNING[2]))
message(sprintf("  ranked metric    limits: [%.2f, %.2f]",
                YLIM_RANKED[1],  YLIM_RANKED[2]))

# -----------------------------------------------------------------------------
# Walk a ggplot layer list and bump the linewidth of line / segment geoms so
# they remain visible when the panel is scaled down to a small journal print.
# -----------------------------------------------------------------------------
bump_linewidths <- function(p, lw_line = 1.4, lw_seg = 0.7) {
  p$layers <- lapply(p$layers, function(layer) {
    cls <- class(layer$geom)[1]
    if (cls %in% c("GeomLine", "GeomPath", "GeomStep")) {
      layer$aes_params$linewidth <- lw_line
      layer$aes_params$size      <- lw_line  # ggplot < 3.4 fallback
    } else if (cls == "GeomSegment") {
      layer$aes_params$linewidth <- lw_seg
      layer$aes_params$size      <- lw_seg
    } else if (cls == "GeomHline" || cls == "GeomVline") {
      layer$aes_params$linewidth <- 0.5
    }
    layer
  })
  p
}

# Wrap dense gene-resolution geoms (running-sum line, rank ticks, ranked-metric
# area / column) in ggrastr::rasterise so they embed as PNG tiles rather than
# tens of thousands of vector primitives. Axes, text, and panel borders stay
# vector.
rasterise_heavy_layers <- function(p, dpi = RASTER_DPI) {
  heavy_classes <- c("GeomLine", "GeomPath", "GeomStep",
                     "GeomSegment",
                     "GeomCol", "GeomBar", "GeomArea", "GeomPolygon",
                     "GeomRect", "GeomRibbon")
  p$layers <- lapply(p$layers, function(layer) {
    if (class(layer$geom)[1] %in% heavy_classes) {
      ggrastr::rasterise(layer, dpi = dpi, dev = "ragg")
    } else {
      layer
    }
  })
  p
}

# -----------------------------------------------------------------------------
# Style helper: theme tuned for a small printed panel (base_size = 17).
# -----------------------------------------------------------------------------
stylise <- function(p,
                    show_x, show_y, show_legend,
                    base_size = 17,
                    legend_pos = c(0.99, 0.99)) {
  p +
    ggplot2::labs(color = NULL) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position      = if (show_legend) legend_pos else "none",
      legend.justification = c(1, 1),
      legend.title         = ggplot2::element_blank(),
      legend.text          = ggplot2::element_text(size = base_size * 0.72),
      legend.background    = ggplot2::element_rect(
                                fill = scales::alpha("white", 0.92),
                                colour = "grey80", linewidth = 0.35),
      legend.margin        = ggplot2::margin(4, 7, 4, 7),
      legend.key.height    = ggplot2::unit(0.75, "lines"),
      legend.key.width     = ggplot2::unit(1.6, "lines"),
      panel.background     = ggplot2::element_rect(fill = "white", color = NA),
      plot.title           = ggplot2::element_text(size = base_size * 1.05,
                                                    face = "bold",
                                                    margin = ggplot2::margin(b = 4)),
      axis.line            = ggplot2::element_line(colour = "black",
                                                    linewidth = 0.6),
      axis.ticks           = ggplot2::element_line(colour = "black",
                                                    linewidth = 0.6),
      axis.ticks.length    = ggplot2::unit(0.18, "cm"),
      axis.title.x         = if (show_x)
                                ggplot2::element_text(size = base_size,
                                                       margin = ggplot2::margin(t = 4))
                             else ggplot2::element_blank(),
      axis.title.y         = if (show_y)
                                ggplot2::element_text(size = base_size * 0.95)
                             else ggplot2::element_blank(),
      axis.text.x          = if (show_x)
                                ggplot2::element_text(size = base_size * 0.85,
                                                       colour = "black")
                             else ggplot2::element_blank(),
      axis.text.y          = if (show_y)
                                ggplot2::element_text(size = base_size * 0.85,
                                                       colour = "black")
                             else ggplot2::element_blank(),
      plot.margin          = ggplot2::margin(6, 14, 6, 6)
    )
}

# -----------------------------------------------------------------------------
# Build one running-sum patchwork for a single contrast.
#
# show_y_axis = FALSE blanks the y axis text/title on panels 1 and 3 (used for
# the right column of the grid composite so the y axis is shared visually with
# the left column).
# ylim_running / ylim_ranked apply shared y limits so the zero line aligns
# across panels.
# -----------------------------------------------------------------------------
make_running_sum <- function(gsea_obj, title, base_size = 17,
                             include_legend = TRUE,
                             show_y_axis = TRUE,
                             ylim_running = NULL,
                             ylim_ranked = NULL) {

  available_ids <- intersect(PATHWAY_IDS, gsea_obj@result$ID)
  if (length(available_ids) == 0) {
    stop("None of the focal SynGO pathways are present in this contrast.")
  }

  # Keep the global ordering so colors stay aligned across panels.
  ids     <- PATHWAY_IDS[PATHWAY_IDS %in% available_ids]
  palette <- unname(PATHWAY_COLORS[ids])
  labels  <- unname(PATHWAY_LABELS[ids])

  # CRITICAL: enrichplot::gseaplot2() maps colors via aes(color = Description)
  # and scale_color_manual(values = color) with an UNNAMED vector, so colors
  # bind to factor levels in alphabetical order of Description - NOT in the
  # order ids were provided. Without correction the rose / purple curves get
  # painted onto the wrong pathways (e.g. postsynaptic cytosol instead of
  # presynaptic ribosome). We fix this below by overriding the scale with a
  # palette named by Description.
  descriptions <- gsea_obj@result$Description[match(ids, gsea_obj@result$ID)]
  if (any(is.na(descriptions))) {
    stop("Could not look up Description for ids: ",
         paste(ids[is.na(descriptions)], collapse = ", "))
  }
  named_palette <- setNames(palette, descriptions)

  p_raw <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID    = ids,
    title        = title,
    subplots     = c(1, 2, 3),
    pvalue_table = FALSE,
    rel_heights  = c(1.6, 0.5, 0.5),
    color        = palette       # overridden by named scale below
  )

  # Override the color scale on the line panel AND the rank-tick panel so the
  # rose / purple / sand / green / teal / olive / sky-blue colors land on the
  # correct Description. `breaks` also forces the legend order to follow
  # PATHWAY_IDS (focal ribosomes first, then umbrella synapse, then the four
  # non-ribosome compartments), instead of alphabetical Description.
  p_raw[[1]] <- p_raw[[1]] +
    ggplot2::scale_color_manual(
      values = named_palette,
      breaks = descriptions,
      labels = labels,
      na.value = "grey60"
    )
  p_raw[[2]] <- p_raw[[2]] +
    ggplot2::scale_color_manual(
      values   = named_palette,
      breaks   = descriptions,
      na.value = "grey60",
      guide    = "none"
    )

  # Apply shared Y limits BEFORE rasterise so the rasterised tile maps onto
  # the new coordinate range. Use coord_cartesian (a zoom, not a data filter)
  # so curves never get clipped at panel edges.
  if (!is.null(ylim_running)) {
    p_raw[[1]] <- suppressMessages(
      p_raw[[1]] + ggplot2::coord_cartesian(ylim = ylim_running, expand = FALSE)
    )
  }
  if (!is.null(ylim_ranked)) {
    p_raw[[3]] <- suppressMessages(
      p_raw[[3]] + ggplot2::coord_cartesian(ylim = ylim_ranked, expand = FALSE)
    )
  }

  p1 <- rasterise_heavy_layers(bump_linewidths(p_raw[[1]], lw_line = 1.5))
  p2 <- rasterise_heavy_layers(bump_linewidths(p_raw[[2]], lw_seg = 0.55))
  p3 <- rasterise_heavy_layers(bump_linewidths(p_raw[[3]], lw_line = 1.1))

  # Color scale was already fixed above (named by Description). Just thicken
  # the legend line swatches for readability at small print sizes.
  p1 <- p1 +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linewidth = 2.2)))

  p1 <- stylise(p1, show_x = FALSE, show_y = show_y_axis,
                show_legend = include_legend, base_size = base_size)
  p2 <- stylise(p2, show_x = FALSE, show_y = FALSE,
                show_legend = FALSE,                base_size = base_size)
  p3 <- stylise(p3, show_x = TRUE,  show_y = show_y_axis,
                show_legend = FALSE,                base_size = base_size)

  patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(2.2, 0.45, 0.55))
}

# -----------------------------------------------------------------------------
# Render one PDF per contrast (large fonts so a ~3-inch print still reads).
# -----------------------------------------------------------------------------
individual_plots <- list()

for (ct in CONTRASTS) {
  obj <- syngo[[ct]]
  if (is.null(obj)) {
    warning("Missing contrast in syngo_gsea_results.rds: ", ct)
    next
  }
  message("Building running-sum plot for ", ct)
  pl <- make_running_sum(obj, CONTRAST_TITLES[ct], base_size = 17,
                         include_legend = TRUE,
                         show_y_axis    = TRUE,
                         ylim_running   = YLIM_RUNNING,
                         ylim_ranked    = YLIM_RANKED)
  individual_plots[[ct]] <- pl

  out <- file.path(out_dir, paste0(ct, "_SynGO_running_sum_normalised.pdf"))
  ggplot2::ggsave(out, pl, width = 5.5, height = 7.0,
                  device = cairo_pdf)
  message("  -> ", out)
}

# -----------------------------------------------------------------------------
# Build a stand-alone legend panel (used by the 2 x 2 composite). The legend
# sits outside the grid - to the right of all four contrasts - so the running-
# enrichment panels stay uncluttered.
# -----------------------------------------------------------------------------
build_legend_panel <- function(base_size = 14) {
  df <- data.frame(
    x       = 1,
    pathway = factor(PATHWAY_IDS, levels = PATHWAY_IDS)
  )
  dummy <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = pathway,
                                             color = pathway)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x + 1, yend = pathway),
                          linewidth = 2.4) +
    ggplot2::scale_color_manual(
      name   = "SynGO pathway",
      values = unname(PATHWAY_COLORS[PATHWAY_IDS]),
      labels = unname(PATHWAY_LABELS[PATHWAY_IDS])
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = base_size * 1.05,
                                                face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size * 0.95),
      legend.key.height = ggplot2::unit(1.1, "lines"),
      legend.key.width  = ggplot2::unit(2.2, "lines"),
      legend.spacing.y  = ggplot2::unit(0.3, "lines")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(linewidth = 2.6)))
  # Extract legend grob (ggplot2 >= 3.5: get_guide_data / get_plot_component)
  if (utils::packageVersion("ggplot2") >= "3.5.0" &&
      exists("get_plot_component", envir = asNamespace("ggplot2"))) {
    leg <- ggplot2::get_plot_component(dummy, "guide-box-right",
                                        return_all = FALSE)
  } else {
    g   <- ggplot2::ggplotGrob(dummy)
    leg <- g$grobs[[which(vapply(g$grobs,
                                 function(x) x$name, character(1)) == "guide-box")]]
  }
  patchwork::wrap_elements(full = leg)
}

# -----------------------------------------------------------------------------
# 2 x 2 composite: columns = mutation (G32A left, R403C right), rows = age
# (D35 top, D65 bottom). A vertical spacer column sits between the two
# mutations so the visual grouping reads cleanly. Shared legend grob placed
# to the right of the grid, spanning both rows.
# -----------------------------------------------------------------------------
message("Building 2 x 2 composite for Fig. 5D")

# Build one running-sum panel per contrast. Left column keeps Y axis
# (show_y_axis = TRUE); right column drops it (text + title) so the shared
# Y scale reads from the left column only.
composite_panels <- list(
  G32A_vs_Ctrl_D35  = make_running_sum(syngo[["G32A_vs_Ctrl_D35"]],
                                       CONTRAST_TITLES["G32A_vs_Ctrl_D35"],
                                       base_size      = 14,
                                       include_legend = FALSE,
                                       show_y_axis    = TRUE,
                                       ylim_running   = YLIM_RUNNING,
                                       ylim_ranked    = YLIM_RANKED),
  R403C_vs_Ctrl_D35 = make_running_sum(syngo[["R403C_vs_Ctrl_D35"]],
                                       CONTRAST_TITLES["R403C_vs_Ctrl_D35"],
                                       base_size      = 14,
                                       include_legend = FALSE,
                                       show_y_axis    = FALSE,
                                       ylim_running   = YLIM_RUNNING,
                                       ylim_ranked    = YLIM_RANKED),
  G32A_vs_Ctrl_D65  = make_running_sum(syngo[["G32A_vs_Ctrl_D65"]],
                                       CONTRAST_TITLES["G32A_vs_Ctrl_D65"],
                                       base_size      = 14,
                                       include_legend = FALSE,
                                       show_y_axis    = TRUE,
                                       ylim_running   = YLIM_RUNNING,
                                       ylim_ranked    = YLIM_RANKED),
  R403C_vs_Ctrl_D65 = make_running_sum(syngo[["R403C_vs_Ctrl_D65"]],
                                       CONTRAST_TITLES["R403C_vs_Ctrl_D65"],
                                       base_size      = 14,
                                       include_legend = FALSE,
                                       show_y_axis    = FALSE,
                                       ylim_running   = YLIM_RUNNING,
                                       ylim_ranked    = YLIM_RANKED)
)

legend_pn <- build_legend_panel(base_size = 14)

# Letter slots:
#   A = top-left     : G32A  D35     B = top-right    : R403C D35
#   C = bottom-left  : G32A  D65     D = bottom-right : R403C D65
#   S = thin vertical strip between mutation columns
#   L = shared legend (spans both rows)
# Cell counts per row: A(6) + S(1) + B(6) + L(4) = 17  ->  strip = 1/17 ~ 6%
design <- "
AAAAAASBBBBBBLLLL
AAAAAASBBBBBBLLLL
AAAAAASBBBBBBLLLL
CCCCCCSDDDDDDLLLL
CCCCCCSDDDDDDLLLL
CCCCCCSDDDDDDLLLL
"

composite <- patchwork::wrap_plots(
  A = composite_panels[["G32A_vs_Ctrl_D35"]],
  B = composite_panels[["R403C_vs_Ctrl_D35"]],
  C = composite_panels[["G32A_vs_Ctrl_D65"]],
  D = composite_panels[["R403C_vs_Ctrl_D65"]],
  S = patchwork::plot_spacer(),
  L = legend_pn,
  design = design
) +
  patchwork::plot_annotation(
    title = "SynGO running-sum enrichment - shared synaptic compartments",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold",
                                          margin = ggplot2::margin(b = 6))
    )
  )

composite_out <- file.path(out_dir,
                           "SynGO_running_sum_normalised_grid.pdf")
ggplot2::ggsave(composite_out, composite,
                width = 14, height = 12.5,
                device = cairo_pdf)
message("  -> ", composite_out)

message("\nDone. ", length(individual_plots),
        " per-contrast PDFs + composite written to:\n  ", out_dir)
