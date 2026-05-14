
library(dplyr)

df <- read.csv("03_Results/02_Analysis/master_gsea_table.csv")

# Check unique values of ever_significant
cat("Unique values of ever_significant:", paste(unique(df$ever_significant), collapse=", "), "\n")

# Create wide format
pathways <- df %>%
  group_by(pathway_id, database, ID, Description) %>%
  summarise(
    Pattern_G32A = first(Pattern_G32A),
    Pattern_R403C = first(Pattern_R403C),
    NES_Early_G32A = first(NES_Early_G32A),
    NES_Early_R403C = first(NES_Early_R403C),
    NES_TrajDev_G32A = first(NES_TrajDev_G32A),
    NES_TrajDev_R403C = first(NES_TrajDev_R403C),
    NES_Late_G32A = first(NES_Late_G32A),
    NES_Late_R403C = first(NES_Late_R403C),
    p_Early_G32A = first(p.adjust_Early_G32A),
    p_Early_R403C = first(p.adjust_Early_R403C),
    p_TrajDev_G32A = first(p.adjust_TrajDev_G32A),
    p_TrajDev_R403C = first(p.adjust_TrajDev_R403C),
    p_Late_G32A = first(p.adjust_Late_G32A),
    p_Late_R403C = first(p.adjust_Late_R403C),
    ever_significant = first(ever_significant),
    .groups = 'drop'
  )

# Test different universe definitions
u_both <- pathways %>% filter(
    (p_Early_G32A < 0.05 | p_TrajDev_G32A < 0.05 | p_Late_G32A < 0.05) &
    (p_Early_R403C < 0.05 | p_TrajDev_R403C < 0.05 | p_Late_R403C < 0.05)
)
u_either <- pathways %>% filter(
    (p_Early_G32A < 0.05 | p_TrajDev_G32A < 0.05 | p_Late_G32A < 0.05) |
    (p_Early_R403C < 0.05 | p_TrajDev_R403C < 0.05 | p_Late_R403C < 0.05)
)

cat("Universe size (Intersection - BOTH):", nrow(u_both), "\n")
cat("Universe size (Union - EITHER):", nrow(u_either), "\n")
sum(pathways$ever_significant == "True")
cat(
  "Universe size (ever_significant == TRUE):", 
  sum(pathways$ever_significant == "True"), 
  "\n"
  )
