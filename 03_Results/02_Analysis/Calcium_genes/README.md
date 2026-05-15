# 03_Results/02_Analysis/Calcium_genes — focused calcium signaling gene analysis

## Overview

This directory holds the focused analysis of 13 calcium signaling genes relevant to DRP1 mutation phenotypes. DRP1 mutations impair mitochondrial calcium buffering (via reduced mitochondrial membrane potential), ATP-dependent SERCA pump activity, and store-operated calcium entry (STIM/ORAI), potentially linking the energetic crisis to neuronal hyperexcitability. The outputs here serve as extended data for the manuscript's calcium phenotype discussion. **Generating script:** `02_Analysis/2.6.viz_calcium_genes.R` (sourced from `02_Analysis/1.1.main_pipeline.R`).

---

## Biological Rationale

**Why calcium genes?**

DRP1 mutations affect mitochondrial dynamics, which impacts:
1. **Mitochondrial calcium buffering** - Mitochondria buffer cytosolic Ca2+
2. **ATP-dependent calcium pumps** - SERCA (ATP2A1) requires ATP
3. **Store-operated calcium entry** - STIM/ORAI pathway
4. **Neuronal excitability** - Calcium signaling underlies action potentials

---

## Target Gene Panel

| Gene | Function | Category |
|------|----------|----------|
| **NNAT** | Neuronatin - ER calcium regulation | ER/calcium |
| **CACNG3** | Calcium channel auxiliary subunit | Channel |
| **CACNA1S** | L-type calcium channel | Channel |
| **ATP2A1** | SERCA calcium ATPase | Pump |
| **RYR1** | Ryanodine receptor (calcium release) | Channel |
| **MYLK3** | Myosin light chain kinase | Kinase |
| **VDR** | Vitamin D receptor | Transcription |
| **STIM1** | Stromal interaction molecule 1 | SOCE |
| **STIM2** | Stromal interaction molecule 2 | SOCE |
| **ORAI1_1** | Calcium release-activated channel | SOCE |
| **CALB1** | Calbindin - calcium buffer | Buffer |
| **CALR** | Calreticulin - ER calcium buffer | Buffer |
| **PNPO** | Pyridoxine phosphate oxidase | Metabolism |

**Note:** CACNA1C, CASR, and ORAI1 were not detected in the dataset. ORAI1_1 is an alternative symbol.

---

## File Inventory

| File | Description |
|------|-------------|
| `calcium_genes_DE_results.csv` | DE statistics for all 13 calcium genes across all 9 contrasts |
| `calcium_genes_boxplots.pdf` | Log2-CPM expression boxplots per gene, grouped by genotype and timepoint |
| `calcium_genes_expression_heatmap.pdf` | Z-score expression heatmap with hierarchical clustering (empty if insufficient variance) |
| `G32A_vs_Ctrl_D35_calcium_volcano.pdf` | Volcano plot (G32A vs Ctrl at D35) with calcium genes highlighted in red |
| `G32A_vs_Ctrl_D65_calcium_volcano.pdf` | Volcano plot (G32A vs Ctrl at D65) with calcium genes highlighted |
| `R403C_vs_Ctrl_D35_calcium_volcano.pdf` | Volcano plot (R403C vs Ctrl at D35) with calcium genes highlighted |
| `R403C_vs_Ctrl_D65_calcium_volcano.pdf` | Volcano plot (R403C vs Ctrl at D65) with calcium genes highlighted |
| `Time_Ctrl_calcium_volcano.pdf` | Volcano plot (D65 − D35 in Controls) with calcium genes highlighted |
| `Time_G32A_calcium_volcano.pdf` | Volcano plot (D65 − D35 in G32A) with calcium genes highlighted |
| `Time_R403C_calcium_volcano.pdf` | Volcano plot (D65 − D35 in R403C) with calcium genes highlighted |
| `Maturation_G32A_specific_calcium_volcano.pdf` | Volcano plot (G32A-specific TrajDev interaction) with calcium genes highlighted |
| `Maturation_R403C_specific_calcium_volcano.pdf` | Volcano plot (R403C-specific TrajDev interaction) with calcium genes highlighted |

---

## Data File Structure

### calcium_genes_DE_results.csv

| Column | Description |
|--------|-------------|
| gene | Gene symbol |
| contrast | Comparison name |
| logFC | Log2 fold-change |
| AveExpr | Average expression |
| t | Moderated t-statistic |
| P.Value | Raw p-value |
| adj.P.Val | FDR-adjusted p-value |
| B | B-statistic |
| significant | TRUE if adj.P.Val < 0.05 |

---

## Column Dictionary: calcium_genes_DE_results.csv

| Column | Description |
|--------|-------------|
| `gene` | HGNC gene symbol |
| `contrast` | Limma contrast name |
| `logFC` | Log2 fold-change |
| `AveExpr` | Average log2-CPM expression |
| `t` | Moderated t-statistic |
| `P.Value` | Raw p-value |
| `adj.P.Val` | BH FDR-adjusted p-value |
| `B` | B-statistic (log-odds of DE) |
| `significant` | TRUE if adj.P.Val < 0.05 |

## Reading Guide

Volcano PDFs show log2FC vs −log10(P.Value) for all 21,122 genes; the 13 calcium panel genes are overlaid in red with labels regardless of significance. Key observation: **NNAT** (neuronatin) is the most strongly downregulated calcium gene (logFC ≈ −4 in multiple contrasts), consistent with ER calcium dysregulation. Three candidate genes (CACNA1C, CASR, ORAI1) were absent from the expression dataset. Calcium pathway-level GSEA statistics are in `../Verification_reports/calcium_pathways_*.csv`.
