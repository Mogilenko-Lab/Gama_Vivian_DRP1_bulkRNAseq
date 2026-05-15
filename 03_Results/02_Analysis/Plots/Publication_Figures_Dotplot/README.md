# Publication Figures Dotplot — Primary Submitted Figures (Fig 6A)

This directory contains dotplot versions of the main publication figures. Dotplots encode statistical significance through dot size and NES through color, making it easier to identify significant pathways than the heatmap versions. **Fig4_Semantic_Pathway_Overview_dotplot.pdf** backs **Fig 6A** of the manuscript.

## File inventory

| File | Description |
|---|---|
| `Fig1_Ribosome_Paradox_dotplot.pdf` | 3-panel dotplot: cytoplasmic biogenesis (Compensation), presynaptic ribosome (Sign_reversal), postsynaptic ribosome (Sign_reversal) across trajectory stages |
| `Fig1_Ribosome_Paradox_dotplot.png` | Raster version (300 DPI) |
| `Fig2_MitoCarta_Trajectory_Patterns_dotplot.pdf` | Dotplot of all MitoCarta pathways across Early / TrajDev / Late for G32A and R403C |
| `Fig2_MitoCarta_Trajectory_Patterns_dotplot.png` | Raster version |
| `Fig3b_SynGO_Trajectory_Patterns_dotplot.pdf` | Dotplot of all SynGO pathways across trajectory stages |
| `Fig3b_SynGO_Trajectory_Patterns_dotplot.png` | Raster version |
| `Fig3_Pattern_Classification_Summary.pdf` | Bar chart of pattern distribution counts across databases (identical to heatmap version — it is a bar chart, not a heatmap) |
| `Fig3_Pattern_Classification_Summary.png` | Raster version |
| `Fig4_Semantic_Pathway_Overview_dotplot.pdf` | Fig 6A: dotplot of semantic pathway categories across all databases and contrasts |
| `Fig4_Semantic_Pathway_Overview_dotplot.png` | Raster version |

## Generating script

`02_Analysis/3.2.publication_figures_dotplot.py`

```bash
python3 02_Analysis/3.2.publication_figures_dotplot.py
```

Input: `03_Results/02_Analysis/master_gsea_table.csv`; semantic categories from `01_Scripts/Python/semantic_categories.py`; pattern definitions from `01_Scripts/Python/pattern_definitions.py`; rendering via `01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting_python/dotplot_renderer.py`.

## Reading guide

**Dot encoding**: Each dot carries three pieces of information simultaneously:
1. **Color** — NES (blue = downregulated, orange = upregulated, white = no enrichment; Blue-White-Orange colorblind-safe diverging scale)
2. **Size** — statistical significance (−log10(FDR)); larger = more significant. Size legend included in all figures showing example dots at FDR = 0.001, 0.01, 0.05, 0.1
3. **Border** — black outline = FDR < 0.05 (significant); gray outline = FDR ≥ 0.05

**Trajectory framework**: Columns ordered Early (D35 mutation effect) → TrajDev (mutation-specific maturation interaction) → Late (D65 mutation effect), separately for G32A and R403C.

**Fig 6A (Fig4_Semantic_Pathway_Overview_dotplot.pdf)**: Each row is a semantic category aggregating related pathways across multiple databases. Columns are trajectory stages for each mutation. Large, black-outlined blue dots in the Synaptic/Translation row at TrajDev identify the amplitude extreme of the ribosome Sign_reversal; large orange dots in Mitochondrial rows show the compensatory Compensation pattern.

**Filtering applied**: CGP database excluded (cancer-focused). Pathways must have at least 3 trajectory data points and at least one significant result (FDR < 0.05).

## Manuscript figure caption (Fig 6A, `Fig4_Semantic_Pathway_Overview_dotplot.pdf`)

**Semantic-category pathway enrichment overview across the developmental trajectory in DRP1 mutant cortical neurons.** Each row is a semantic category that aggregates related pathways across multiple enrichment databases (e.g. mitochondrial translation, cytoplasmic translation, synaptic translation, OXPHOS, ATP synthesis, mtDNA, cell-cycle/DNA-replication, cilium); category-to-pathway mapping is curated in `01_Scripts/Python/semantic_categories.py`. Each column is one trajectory cell of one mutation (Early = D35 mutant-vs-control, TrajDev = mutation-specific maturation interaction, Late = D65 mutant-vs-control), drawn separately for G32A and R403C. Each dot encodes three values simultaneously: colour = aggregated NES on a blue–white–orange diverging scale (blue = downregulated, orange = upregulated; values clipped at NES = ±4); size = −log10(BH-FDR), with a fixed-position legend showing reference dots at FDR = 0.001, 0.01, 0.05, 0.1; border = black outline if BH-FDR < 0.05, grey outline otherwise. Statistical thresholds inherited from the upstream fgsea pipeline (10,000 permutations); the CGP database is excluded (cancer-focused) and categories must have at least three trajectory data points and at least one significant aggregate. The panel resolves the translation paradox at semantic-category scale: the synaptic-translation row carries large, black-outlined blue dots in the TrajDev and Late cells of both mutations (the amplitude extreme of the Sign_reversal program documented in Fig 5F and Supp Fig S9), while the mitochondrial-translation and OXPHOS rows carry large orange dots in the same cells (the Compensation arm documented in Fig 6B and the GSVA panels). The cell-cycle / DNA-replication row carries the second-largest TrajDev signal, consistent with the mutation-independent genome-maintenance response reported in the main text.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/3.2.publication_figures_dotplot.py`
