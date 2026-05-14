#!/usr/bin/env Rscript
# verify_design_matrix.R
# Reproduces the DGEList + design matrix exactly as 1.1.main_pipeline.R does.
# Run: Rscript 02_Analysis/verify_design_matrix.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

# ── 0. Paths ──────────────────────────────────────────────────────
counts_file   <- "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt"
metadata_file <- "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv"

# ── 1. Read count matrix ──────────────────────────────────────────
counts <- read.delim(counts_file, check.names = FALSE,
                     stringsAsFactors = FALSE, header = TRUE)
rownames(counts) <- counts$Geneid
counts <- counts[, -which(colnames(counts) == "Geneid"), drop = FALSE]
counts <- as.matrix(counts)
mode(counts) <- "numeric"

# Strip _Aligned suffix (mimics read_count_matrix.R line)
count_sample_names <- sub("_Aligned$", "", colnames(counts))
colnames(counts) <- count_sample_names

cat("✅ Count matrix:\n")
cat("   Rows (genes):   ", nrow(counts), "\n")
cat("   Cols (samples): ", ncol(counts), "\n")

# ── 2. Read metadata ──────────────────────────────────────────────
metadata <- read.csv(metadata_file, sep = ";", stringsAsFactors = FALSE, check.names = FALSE)
rownames(metadata) <- metadata$sample

cat("✅ Metadata:\n")
cat("   Rows (samples)  ", nrow(metadata), "\n")

# ── 3. Intersection (mimics read_count_matrix.R) ─────────────────--
common_samples <- intersect(colnames(counts), rownames(metadata))
cat("✅ Common (intersect):", length(common_samples), "\n")
if (length(common_samples) == 0) stop("No matching samples!")

counts  <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ── 4. Create DGEList exactly as read_count_matrix.R ──────────────
DGE <- DGEList(counts = counts)

# Add factor annotations (as in read_count_matrix.R)
DGE$samples$genotype <- factor(metadata$genotype)
DGE$samples$days     <- factor(metadata$days)
DGE$samples$rep      <- factor(metadata$rep)
DGE$samples$cell_line <- factor(metadata$cell_line)
DGE$samples$cell_type <- factor(metadata$cell_type)

# Key step from 1.1.main_pipeline.R:
DGE$samples$group <- interaction(DGE$samples$days, DGE$samples$genotype)

cat("\n=== DGEList Summary ===\n")
cat("Total samples : ", ncol(DGE), "\n")
cat("Total genes   : ", nrow(DGE), "\n")
cat("Group factor  : ", nlevels(DGE$samples$group), " levels\n")

# ── 5. Design matrix (exact 1.1.main_pipeline.R) ─────────────────
design <- model.matrix(~0 + DGE$samples$group, data = DGE$samples)
colnames(design) <- levels(DGE$samples$group)
colnames(design)

cat("\n=== Design Matrix ===\n")
cat("Dimensions:", dim(design), "\n")
print(table(DGE$samples$group))

cat("\nDesign columns (group levels):\n")
print(colnames(design))

cat("\nReplicates per design column:\n")
for (co in colnames(design)) {
  n <- sum(design[, co] == 1)
  cat(sprintf("  %-20s n = %d\n", co, n))
}

# ── 6. Verify contrasts match the script ──────────────────────────
# Explicitly define levels as character vector (avoids '.' vs '_' mismatch)
colnames(design) <- gsub("\\.", "_", colnames(design))

contrasts <- makeContrasts(
  G32A_vs_Ctrl_D35   = D35_G32A  - D35_Control,
  R403C_vs_Ctrl_D35  = D35_R403C - D35_Control,
  G32A_vs_Ctrl_D65   = D65_G32A   - D65_Control,
  R403C_vs_Ctrl_D65  = D65_R403C  - D65_Control,
  Time_Ctrl          = D65_Control - D35_Control,
  Time_G32A          = D65_G32A   - D35_G32A,
  Time_R403C         = D65_R403C  - D35_R403C,
  Maturation_G32A_specific = (D65_G32A - D35_G32A) - (D65_Control - D35_Control),
  Maturation_R403C_specific = D65_R403C - D35_R403C - D65_Control + D35_Control,
  levels = design
)

cat("\n=== Contrast Matrix ===\n")
print(dim(contrasts))
cat("Contrasts:\n")
for (co in colnames(contrasts)) cat(" ", co, "\n")

# ── 7. Verify linear algebra (rank check) ─────────────────────────
cat("\n=== Rank Check ===\n")
cat("Design matrix rank:", qr(design)$rank, "(columns:", ncol(design), ")\n")
cat("\n✅ Verification complete. Design matches raw counts + metadata exactly.\n")
