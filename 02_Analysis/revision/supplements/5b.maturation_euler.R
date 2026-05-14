#!/usr/bin/env Rscript
## 5b.maturation_euler.R
##
## Concern: replace Fig. 4c UpSet with a three-set Euler diagram
## of the maturation contrasts (Time_Ctrl, Time_G32A, Time_R403C).
##
## Produces Euler diagrams at FDR<0.05 (primary Fig. 4c replacement) and
## FDR<0.10 (relaxed-threshold robustness panel for Supplementary Fig. S6).
##
## Strictly additive: reads DE CSVs only, writes to NEW directory:
##   03_Results/02_Analysis/Plots/General/Supplementary_5b/
## Does NOT modify or overwrite any existing artifact.
##
## Run:
##   Rscript 02_Analysis/5b.maturation_euler.R

suppressPackageStartupMessages({
  library(here)
})

# Package availability: prefer eulerr, fall back to VennDiagram
has_eulerr <- requireNamespace("eulerr", quietly = TRUE)
has_venn   <- requireNamespace("VennDiagram", quietly = TRUE)

if (!has_eulerr) {
  message("eulerr not installed; attempting install...")
  install.packages("eulerr", repos = "https://cloud.r-project.org", quiet = TRUE)
  has_eulerr <- requireNamespace("eulerr", quietly = TRUE)
}

project_root <- "/workspaces/Gama_Vivian_DRP1_bulkRNAseq"
de_dir  <- file.path(project_root, "03_Results/02_Analysis/DE_results")
out_dir <- file.path(project_root, "03_Results/02_Analysis/Plots/General/Supplementary_5b")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load DE results for 3 maturation contrasts -----------------------------
g32a  <- read.csv(file.path(de_dir, "Time_G32A_results.csv"),  row.names = 1)
r403c <- read.csv(file.path(de_dir, "Time_R403C_results.csv"), row.names = 1)
ctrl  <- read.csv(file.path(de_dir, "Time_Ctrl_results.csv"),  row.names = 1)

sig_at <- function(df, alpha) rownames(df)[df$adj.P.Val < alpha]

# Okabe-Ito palette (colorblind-safe); matches project color_config scheme
col_ctrl  <- "#999999"   # grey — Time_Ctrl
col_g32a  <- "#D55E00"   # vermilion — G32A
col_r403c <- "#0072B2"   # blue — R403C

# --- Plotting helper --------------------------------------------------------
plot_euler <- function(sets, out_pdf, out_png, title_text) {
  ns <- sapply(sets, length)
  labels <- c(
    sprintf("Time_Ctrl\n(n=%d)",  ns[["Time_Ctrl"]]),
    sprintf("Time_G32A\n(n=%d)",  ns[["Time_G32A"]]),
    sprintf("Time_R403C\n(n=%d)", ns[["Time_R403C"]])
  )
  if (has_eulerr) {
    fit <- eulerr::euler(sets, shape = "ellipse")
    p <- plot(fit,
              fills      = list(fill  = c(col_ctrl, col_g32a, col_r403c),
                                alpha = 0.5),
              labels     = list(labels = labels, fontfamily = "sans", cex = 1.0),
              quantities = list(fontfamily = "sans", cex = 0.9),
              main       = title_text)
    grDevices::pdf(out_pdf, width = 5.5, height = 4.2)
    print(p); grDevices::dev.off()
    grDevices::png(out_png, width = 5.5, height = 4.2, units = "in", res = 300)
    print(p); grDevices::dev.off()
    return("eulerr")
  }
  if (has_venn) {
    # Fallback: static Venn (circles, not area-proportional)
    VennDiagram::venn.diagram(
      x              = sets,
      filename       = out_pdf,
      imagetype      = "pdf",
      height         = 5, width = 5,
      category.names = labels,
      fill           = c(col_ctrl, col_g32a, col_r403c),
      alpha          = 0.5, cex = 1.1, fontfamily = "sans",
      cat.cex = 1.0, cat.fontfamily = "sans",
      main           = title_text
    )
    return("VennDiagram")
  }
  stop("Neither eulerr nor VennDiagram available.")
}

# --- FDR < 0.05 primary panel (Fig. 4c replacement) -------------------------
sets_05 <- list(
  Time_Ctrl  = sig_at(ctrl,  0.05),
  Time_G32A  = sig_at(g32a,  0.05),
  Time_R403C = sig_at(r403c, 0.05)
)
engine_used <- plot_euler(
  sets    = sets_05,
  out_pdf = file.path(out_dir, "Maturation_DEG_Euler_FDR05.pdf"),
  out_png = file.path(out_dir, "Maturation_DEG_Euler_FDR05.png"),
  title_text = "Maturation-responsive DEGs (FDR < 0.05)"
)
message("Primary Euler (FDR<0.05) written; engine: ", engine_used)

# --- FDR < 0.10 relaxed-threshold panel (Supp Fig. S6 robustness) -----------
sets_10 <- list(
  Time_Ctrl  = sig_at(ctrl,  0.10),
  Time_G32A  = sig_at(g32a,  0.10),
  Time_R403C = sig_at(r403c, 0.10)
)
plot_euler(
  sets    = sets_10,
  out_pdf = file.path(out_dir, "Maturation_DEG_Euler_FDR10.pdf"),
  out_png = file.path(out_dir, "Maturation_DEG_Euler_FDR10.png"),
  title_text = "Maturation-responsive DEGs (FDR < 0.10, robustness)"
)
message("Relaxed Euler (FDR<0.10) written.")

# --- Always produce a VennDiagram fallback for both thresholds --------------
# Useful because eulerr renders 0-member sets as degenerate points which some
# readers may find unusual; the Venn makes the empty-set explicit. VennDiagram
# does not support imagetype="pdf" directly, so we build the grob and route it
# through grDevices::pdf manually.
make_venn_grob <- function(sets, title_text) {
  ns <- sapply(sets, length)
  grid::grid.newpage()
  VennDiagram::venn.diagram(
    x              = sets,
    filename       = NULL,        # returns a grid grob instead of rendering
    category.names = c(
      sprintf("Time_Ctrl (n=%d)",  ns[["Time_Ctrl"]]),
      sprintf("Time_G32A (n=%d)",  ns[["Time_G32A"]]),
      sprintf("Time_R403C (n=%d)", ns[["Time_R403C"]])),
    fill     = c(col_ctrl, col_g32a, col_r403c), alpha = 0.5,
    cex      = 1.2, fontfamily = "sans",
    cat.cex  = 1.0, cat.fontfamily = "sans",
    main     = title_text
  )
}

if (has_venn) {
  venn_pdf_05 <- file.path(out_dir, "Maturation_DEG_Venn_FDR05.pdf")
  venn_pdf_10 <- file.path(out_dir, "Maturation_DEG_Venn_FDR10.pdf")
  tryCatch({
    g05 <- make_venn_grob(sets_05, "Maturation-responsive DEGs (FDR < 0.05) — Venn fallback")
    grDevices::pdf(venn_pdf_05, width = 5, height = 5); grid::grid.draw(g05); grDevices::dev.off()

    g10 <- make_venn_grob(sets_10, "Maturation-responsive DEGs (FDR < 0.10) — Venn fallback")
    grDevices::pdf(venn_pdf_10, width = 5, height = 5); grid::grid.draw(g10); grDevices::dev.off()

    # Clean up any VennDiagram log files it may drop in the working dir
    file.remove(list.files(out_dir, pattern = "VennDiagram.*\\.log$",
                           full.names = TRUE))
    file.remove(list.files(".", pattern = "VennDiagram.*\\.log$",
                           full.names = TRUE))
    message("VennDiagram fallbacks written.")
  }, error = function(e) {
    message("VennDiagram fallback failed (non-fatal): ", conditionMessage(e))
  })
}

cat("\nDone.\n")
cat("Output directory:\n  ", out_dir, "\n", sep = "")
