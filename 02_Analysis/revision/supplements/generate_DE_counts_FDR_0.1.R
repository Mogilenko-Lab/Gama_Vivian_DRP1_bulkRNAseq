###############################################################################
##  Generate DE gene counts plot with FDR < 0.1                              ##
###############################################################################

library(here)
library(limma)
library(ggplot2)
library(tidyr)

# Load the checkpoint data
checkpoint_dir <- here::here("03_Results/02_Analysis/checkpoints")
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
contrasts_obj <- readRDS(file.path(checkpoint_dir, "contrasts_matrix.rds"))

# Re-run decideTests with FDR < 0.1
de_results_0.1 <- limma::decideTests(fit, adjust.method = "BH", p.value = 0.1)

# Prepare data for plotting
contrast_order <- c("G32A_vs_Ctrl_D35","R403C_vs_Ctrl_D35",
                    "G32A_vs_Ctrl_D65","R403C_vs_Ctrl_D65",
                    "Time_Ctrl","Time_G32A","Time_R403C",
                    "Maturation_G32A_specific","Maturation_R403C_specific")

deg_counts <- data.frame(
  Contrast = colnames(contrasts_obj),
  Up   = colSums(de_results_0.1  > 0),
  Down = colSums(de_results_0.1  < 0)) |>
  transform(Contrast = factor(Contrast, levels = contrast_order))

deg_long <- tidyr::pivot_longer(deg_counts,
                                c(Up, Down),
                                names_to = "Direction",
                                values_to = "Count")

# Create output directory
out_dir <- here::here("03_Results/02_Analysis/Plots/General")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Generate plot
pdf(file.path(out_dir, "DE_gene_counts_FDR_0.1.pdf"), 10, 6, onefile = TRUE)
print(ggplot(deg_long,
       aes(Contrast, Count, fill = Direction)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_text(aes(label = Count),
            vjust = -.2,
            position = position_dodge(width = .8),
            size = 3) +
  scale_fill_manual(values = c(Up = "#D55E00", Down = "#0072B2")) + # Okabe-Ito
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(y = "gene count", x = NULL, fill = "",
       title = "Differentially Expressed Genes (FDR < 0.1)"))
dev.off()

message("✓ DE gene counts plot (FDR < 0.1) saved to: ",
        file.path(out_dir, "DE_gene_counts_FDR_0.1.pdf"))

# Print summary
cat("\n=== DE Gene Counts (FDR < 0.1) ===\n")
print(deg_counts)
