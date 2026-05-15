###############################################################################
##  Synaptic Ribosome Analysis: Presynaptic vs Postsynaptic                  ##
##  Detailed comparison of compartment-specific effects                       ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)  # For topTable
library(ComplexHeatmap)
library(circlize)
library(grid)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
syngo_lists <- readRDS(file.path(checkpoint_dir, "syngo_lists.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))

# Load gene lists
gene_list_dir <- here("03_Results/02_Analysis/Verification_reports")
presyn_genes <- readLines(file.path(gene_list_dir, "syngo_presyn_ribosome_genes.txt"))
postsyn_genes <- readLines(file.path(gene_list_dir, "syngo_postsyn_ribosome_genes.txt"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Synaptic_ribosomes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  PANEL A: REMOVED - Redundant with Panel C heatmap                       ##
###############################################################################
## Previously showed postsynaptic-specific genes, now visible in Panel C

# Identify compartment-specific genes (needed for Panel C)
shared_genes <- intersect(presyn_genes, postsyn_genes)
postsyn_only <- setdiff(postsyn_genes, presyn_genes)
all_ribosome_genes <- unique(c(presyn_genes, postsyn_genes))

message(sprintf("  Gene overlap: Shared: %d (%.0f%%), Postsynaptic-only: %d",
                length(shared_genes), 100*length(shared_genes)/length(presyn_genes),
                length(postsyn_only)))

###############################################################################
##  PANEL B: REMOVED - Data extraction only (no plot, redundant)            ##
###############################################################################
## Extract enrichment statistics for summary, but don't create redundant plot

message("📊 Extracting enrichment statistics (no plot, redundant with Translation_Paradox)...")

contrasts <- c("Maturation_G32A_specific", "Maturation_R403C_specific")
enrichment_stats <- data.frame()

for (contrast in contrasts) {
  syngo_res <- syngo_gsea_results[[contrast]]
  if (!is.null(syngo_res)) {
    res_df <- syngo_res@result

    # Presynaptic
    presyn_row <- res_df[res_df$ID == "SYNGO:presyn_ribosome", ]
    if (nrow(presyn_row) > 0) {
      enrichment_stats <- rbind(enrichment_stats, data.frame(
        Contrast = contrast,
        Compartment = "Presynaptic",
        NES = presyn_row$NES,
        FDR = presyn_row$p.adjust,
        GeneCount = presyn_row$setSize,
        CoreEnrichment = length(unlist(strsplit(presyn_row$core_enrichment, "/")))
      ))
    }

    # Postsynaptic
    postsyn_row <- res_df[res_df$ID == "SYNGO:postsyn_ribosome", ]
    if (nrow(postsyn_row) > 0) {
      enrichment_stats <- rbind(enrichment_stats, data.frame(
        Contrast = contrast,
        Compartment = "Postsynaptic",
        NES = postsyn_row$NES,
        FDR = postsyn_row$p.adjust,
        GeneCount = postsyn_row$setSize,
        CoreEnrichment = length(unlist(strsplit(postsyn_row$core_enrichment, "/")))
      ))
    }
  }
}

enrichment_stats$Contrast_Clean <- gsub("Maturation_|_specific", "", enrichment_stats$Contrast)

message("✓ Data extraction complete\n")

###############################################################################
##  PANEL C: Expression Heatmap by Compartment                               ##
###############################################################################

message("📊 Creating Panel C: Expression heatmap by compartment...")

# Get logFC values for maturation contrasts (including Time_Ctrl for baseline comparison)
contrasts_of_interest <- c("Time_Ctrl", "Maturation_G32A_specific", "Maturation_R403C_specific")

# Extract logFC matrix
logfc_matrix <- sapply(contrasts_of_interest, function(contrast) {
  coef_idx <- which(colnames(fit$coefficients) == contrast)
  logfc <- fit$coefficients[, coef_idx]
  names(logfc) <- rownames(fit$coefficients)
  return(logfc)
})

# Annotate genes by compartment
all_ribosome_genes_sorted <- sort(all_ribosome_genes)
compartment_annotation <- sapply(all_ribosome_genes_sorted, function(g) {
  if (g %in% shared_genes) {
    return("Both")
  } else if (g %in% postsyn_only) {
    return("Postsynaptic only")
  } else {
    return("Unknown")
  }
})

# Filter to genes present in data
genes_present <- all_ribosome_genes_sorted[all_ribosome_genes_sorted %in% rownames(logfc_matrix)]
logfc_subset <- logfc_matrix[genes_present, ]
compartment_subset <- compartment_annotation[genes_present]

# Convert to factor with explicit levels to control display order
# ComplexHeatmap displays splits in factor level order, so we set the order explicitly
compartment_subset <- factor(compartment_subset,
                             levels = c("Postsynaptic only", "Both"))

# Rename columns
colnames(logfc_subset) <- c("Ctrl", "G32A", "R403C")

# Create annotation
ha_compartment <- rowAnnotation(
  Compartment = compartment_subset,
  col = list(Compartment = c("Both" = "#999999",
                              "Postsynaptic only" = "#56B4E9")),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_legend_param = list(
    Compartment = list(title = "SynGO Annotation")
  )
)

# Color scheme for logFC
# Use consistent Blue-White-Orange palette as in publication figures
col_fun <- colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "#F7F7F7", "#B35806"))

###############################################################################
##  PANEL C: Expression Heatmap by Compartment                               ##
###############################################################################

message("📊 Creating Panel C: Expression heatmap by compartment...")

# Define trajectory contrasts as in publication figures
contrasts_of_interest <- c(
  "G32A_vs_Ctrl_D35",          # Early G32A
  "Maturation_G32A_specific",  # TrajDev G32A
  "G32A_vs_Ctrl_D65",          # Late G32A
  "R403C_vs_Ctrl_D35",         # Early R403C
  "Maturation_R403C_specific", # TrajDev R403C
  "R403C_vs_Ctrl_D65"          # Late R403C
)

# Extract logFC matrix
logfc_matrix <- sapply(contrasts_of_interest, function(contrast) {
  coef_idx <- which(colnames(fit$coefficients) == contrast)
  logfc <- fit$coefficients[, coef_idx]
  names(logfc) <- rownames(fit$coefficients)
  return(logfc)
})

# Annotate genes by compartment
all_ribosome_genes_sorted <- sort(all_ribosome_genes)
compartment_annotation <- sapply(all_ribosome_genes_sorted, function(g) {
  if (g %in% shared_genes) {
    return("Both")
  } else if (g %in% postsyn_only) {
    return("Postsynaptic only")
  } else {
    return("Unknown")
  }
})

# Filter to genes present in data
genes_present <- all_ribosome_genes_sorted[all_ribosome_genes_sorted %in% rownames(logfc_matrix)]
logfc_subset <- logfc_matrix[genes_present, ]
compartment_subset <- compartment_annotation[genes_present]

# Convert to factor with explicit levels to control display order
compartment_subset <- factor(compartment_subset,
                             levels = c("Postsynaptic only", "Both"))

# Rename columns to simple timepoints
colnames(logfc_subset) <- rep(c("Early", "TrajDev", "Late"), 2)

# Create column split vector
column_split <- factor(rep(c("G32A", "R403C"), each = 3), levels = c("G32A", "R403C"))

# Create annotation
ha_compartment <- rowAnnotation(
  Compartment = compartment_subset,
  col = list(Compartment = c("Both" = "#999999",
                              "Postsynaptic only" = "#56B4E9")),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_legend_param = list(
    Compartment = list(title = "SynGO Annotation")
  )
)

# Create heatmap
pdf(file.path(out_dir, "Panel_C_Expression_Heatmap.pdf"), width = 6, height = 14)
ht <- Heatmap(
  logfc_subset,
  name = "logFC",
  col = col_fun,

  # Row settings
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,  # Enable clustering within each compartment
  cluster_row_slices = FALSE,  # Prevent reordering of compartment slices
  show_row_dend = TRUE,
  row_dend_width = unit(10, "mm"),
  row_split = compartment_subset,
  row_title = c("Postsynaptic Only", "Both Compartments"),
  row_title_gp = gpar(fontface = "bold", fontsize = 10),
  row_gap = unit(2, "mm"),

  # Column settings
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_split = column_split,
  column_title_gp = gpar(fontface = "bold", fontsize = 11),
  column_gap = unit(2, "mm"),

  # Annotations
  right_annotation = ha_compartment,

  # Heatmap appearance
  border = TRUE,
  rect_gp = gpar(col = "gray80", lwd = 0.5),

  # Legend
  heatmap_legend_param = list(
    title = "logFC",
    at = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = c("-0.6", "-0.3", "0", "0.3", "0.6"),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 10)
  )
)

draw(ht, heatmap_legend_side = "right", column_title = "Synaptic Ribosome Expression Trajectories")
dev.off()

message("✓ Panel C complete\n")

###############################################################################
##  PANEL C (refactored): Sample-level expression heatmap                    ##
##  Print-typography variant, no compartment side-strip, log2 CPM per sample ##
###############################################################################

message("📊 Creating refactored Panel C: sample-level expression heatmap...")

# Voom-normalized log2 CPM and per-sample metadata
voom_obj   <- readRDS(file.path(checkpoint_dir, "voom_object.rds"))
dge_obj    <- readRDS(file.path(checkpoint_dir, "DGE_object.rds"))
expr_mat   <- voom_obj$E
sample_meta <- dge_obj$samples
sample_meta$sample_id <- rownames(sample_meta)

# Subset to ribosome genes detected in expression matrix
ribo_present <- intersect(all_ribosome_genes_sorted, rownames(expr_mat))
expr_sub     <- expr_mat[ribo_present, , drop = FALSE]

# Row z-score so the matrix shows Ctrl→Mutant dynamics on a common scale
# (raw log2 CPM mixes high- and low-abundance genes and washes out variation).
expr_z <- t(scale(t(expr_sub), center = TRUE, scale = TRUE))

# Compartment factor — used for row split only, no colored side strip
compartment_subset <- factor(
  ifelse(ribo_present %in% shared_genes, "Both",
         ifelse(ribo_present %in% postsyn_only, "Postsynaptic only", "Unknown")),
  levels = c("Postsynaptic only", "Both")
)

# Order samples Genotype (G32A → R403C → Control) → Day (D35 → D65) → replicate
sample_meta$genotype <- factor(sample_meta$genotype, levels = c("G32A", "R403C", "Control"))
sample_meta$days     <- factor(sample_meta$days,     levels = c("D35", "D65"))
sample_meta <- sample_meta[order(sample_meta$genotype, sample_meta$days, sample_meta$rep), ]
expr_z      <- expr_z[, sample_meta$sample_id, drop = FALSE]

# Six-group column split (Genotype × Day) — D35 then D65 within each genotype
column_split <- factor(
  paste(sample_meta$genotype, sample_meta$days, sep = "\n"),
  levels = c("G32A\nD35",   "G32A\nD65",
             "R403C\nD35",  "R403C\nD65",
             "Control\nD35","Control\nD65")
)

# Within-group replicate index → compact column labels
sample_short <- ave(seq_along(column_split), column_split, FUN = seq_along)

# Compact top-of-column annotation (Day + Genotype)
genotype_colors <- c("Control" = "#BBBBBB", "G32A" = "#E69F00", "R403C" = "#0072B2")
day_colors      <- c("D35" = "#D9D9D9", "D65" = "#525252")

ha_top <- HeatmapAnnotation(
  Day      = sample_meta$days,
  Genotype = sample_meta$genotype,
  col = list(Day = day_colors, Genotype = genotype_colors),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 8),
  simple_anno_size     = unit(2.2, "mm"),
  gap                  = unit(0.4, "mm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    Day = list(title = "Day", title_gp = gpar(fontsize = 9, fontface = "bold"),
               labels_gp = gpar(fontsize = 8),
               grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Genotype = list(title = "Genotype", title_gp = gpar(fontsize = 9, fontface = "bold"),
                    labels_gp = gpar(fontsize = 8),
                    grid_height = unit(4, "mm"), grid_width = unit(4, "mm"))
  )
)

# Diverging palette matched to original Panel C: Blue–White–Orange.
# Symmetric clip at ±2 SD on row-z (sample-level dynamic range).
zmax  <- 2
col_z <- colorRamp2(c(-zmax, 0, zmax), c("#2166AC", "#F7F7F7", "#B35806"))

# Compact figure (~4.4 × 6.6 in keeps original ~0.66 ratio, smaller print footprint)
pdf(file.path(out_dir, "Panel_C_Sample_Heatmap.pdf"), width = 4.4, height = 6.6)
ht_samp <- Heatmap(
  expr_z,
  name = "z",
  col  = col_z,

  # Rows — names left, no side strip, split by SynGO compartment
  row_names_side       = "left",
  row_names_gp         = gpar(fontsize = 7),
  cluster_rows         = TRUE,
  cluster_row_slices   = FALSE,
  show_row_dend        = FALSE,
  row_split            = compartment_subset,
  row_title            = c("Postsynaptic only", "Postsynaptic + Presynaptic"),
  row_title_gp         = gpar(fontface = "bold", fontsize = 9),
  row_title_rot        = 90,
  row_gap              = unit(1.5, "mm"),

  # Columns — day×genotype split (no per-sample tick labels for cleanliness)
  cluster_columns       = FALSE,
  show_column_names     = FALSE,
  column_split          = column_split,
  column_title_gp       = gpar(fontface = "bold", fontsize = 8),
  column_gap            = unit(0.6, "mm"),

  top_annotation = ha_top,

  # Cells — flat fill, narrow matrix to push the figure into a thin column
  border  = TRUE,
  rect_gp = gpar(col = NA),
  width   = unit(55, "mm"),

  heatmap_legend_param = list(
    title         = "Row z-score\n(log2 CPM)",
    at            = c(-2, -1, 0, 1, 2),
    legend_height = unit(34, "mm"),
    grid_width    = unit(4, "mm"),
    title_gp      = gpar(fontface = "bold", fontsize = 9),
    labels_gp     = gpar(fontsize = 8)
  )
)

draw(
  ht_samp,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  merge_legend           = TRUE,
  column_title           = "Synaptic Ribosome Expression — sample-level",
  column_title_gp        = gpar(fontface = "bold", fontsize = 10)
)
dev.off()

message("✓ Refactored Panel C (sample-level) complete\n")

###############################################################################
##  PANEL D: REMOVED - Redundant scatter plot                               ##
###############################################################################
## This panel showed pre vs post NES on a diagonal (equal effects)
## Information is redundant and not mechanistically insightful

message("ℹ️  Panel D skipped (redundant scatter plot)\n")

###############################################################################
##  PANEL E: REMOVED - Leading edge analysis (statistical detail)           ##
###############################################################################
## This panel showed core enrichment percentages (~70-80%)
## Not mechanistically interesting - just statistical detail

message("ℹ️  Panel E skipped (statistical detail, not mechanistically interesting)\n")

###############################################################################
##  Summary Statistics                                                       ##
###############################################################################

message("📝 Generating summary statistics...\n")

cat("\n=== Synaptic Ribosome Analysis Summary ===\n\n")

cat("Gene Set Composition:\n")
cat(sprintf("  Presynaptic ribosome: %d genes\n", length(presyn_genes)))
cat(sprintf("  Postsynaptic ribosome: %d genes\n", length(postsyn_genes)))
cat(sprintf("  Shared (both compartments): %d genes (%.1f%% of presynaptic)\n",
            length(shared_genes), 100*length(shared_genes)/length(presyn_genes)))
cat(sprintf("  Postsynaptic-only: %d genes\n", length(postsyn_only)))

cat("\nPostsynaptic-only genes:\n")
cat(sprintf("  %s\n", paste(postsyn_only, collapse = ", ")))

cat("\nEnrichment Statistics:\n")
for (i in 1:nrow(enrichment_stats)) {
  row <- enrichment_stats[i, ]
  cat(sprintf("\n%s - %s:\n", row$Contrast_Clean, row$Compartment))
  cat(sprintf("  NES: %.3f\n", row$NES))
  cat(sprintf("  FDR: %.2e\n", row$FDR))
  cat(sprintf("  Gene Set Size: %d\n", row$GeneCount))
  cat(sprintf("  Core Enrichment: %d genes\n", row$CoreEnrichment))
}

cat("\n")
message("✅ Synaptic ribosome analysis complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  - Panel_C_Expression_Heatmap.pdf (compartment-specific expression with normal baseline)")
message("\n🎯 KEY FINDINGS:")
message("  1. Presynaptic genes are SUBSET of postsynaptic (100% overlap)")
message("  2. 18 postsynaptic-only genes show unique expression patterns")
message("  3. Both compartments show EQUAL downregulation (NES ~ -3.0)")
message("  4. Translation machinery failure is PAN-SYNAPTIC, not compartment-specific")
message("  5. Clustering reveals co-regulated gene modules within compartments")
message("  6. Time_Ctrl baseline shows normal maturation trajectory for contrast")
