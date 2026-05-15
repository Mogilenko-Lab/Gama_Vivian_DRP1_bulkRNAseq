# Complex V Analysis — ATP Synthase Deep-Dive

This directory contains a focused analysis of ATP synthase (OXPHOS Complex V) gene expression and pathway enrichment in DRP1 mutant cortical neurons. Complex V is the final step of oxidative phosphorylation and the direct source of the ATP whose synaptic depletion underlies the translation paradox. These figures are supplementary to the main manuscript.

## File inventory

| File | Description |
|---|---|
| `Complex_V_Gene_Expression_Heatmap.pdf` | Log2FC heatmap for Complex V genes across trajectory stages (Early / TrajDev / Late) for G32A and R403C |
| `Complex_V_Pathway_Enrichment.pdf` | GSEA enrichment dotplot for all MitoCarta Complex V pathways across all 9 contrasts |
| `Complex_V_Genes.csv` | Gene list extracted from significant GSEA core enrichment |
| `Complex_V_Pathway_Data.csv` | Full GSEA statistics (NES, FDR, core enrichment) for Complex V pathways |

## Generating script

`02_Analysis/2.5.viz_complex_v_analysis.R`

```bash
Rscript 02_Analysis/2.5.viz_complex_v_analysis.R
```

Requires checkpoints `checkpoints/mitocarta_gsea_results.rds` and `checkpoints/fit_object.rds`.

## Reading guide

### Complex_V_Gene_Expression_Heatmap.pdf

Side-by-side heatmaps (G32A left, R403C right). Rows are individual Complex V subunits and assembly factors; columns are trajectory stages:

- **Early**: D35 mutation vs control (G32A_vs_Ctrl_D35 or R403C_vs_Ctrl_D35)
- **TrajDev**: Mutation-specific maturation interaction term (Maturation_G32A/R403C_specific)
- **Late**: D65 mutation vs control

Color: blue = downregulated, white = no change, red/orange = upregulated (scale ±0.6 log2FC). Both F0 and F1 sector subunits typically show early suppression followed by a TrajDev compensatory upregulation, illustrating the cell's adaptive response to ATP shortage. The assembly factors (ATPAF1, ATPAF2, TMEM70) cluster with their target subunits.

### Complex_V_Pathway_Enrichment.pdf

Dotplot with pathways on the y-axis and contrasts on the x-axis. Dot color encodes NES (red = upregulated, blue = downregulated); dot size encodes −log10(FDR). Pathways were retrieved from MitoCarta 3.0 using keywords: `complex.*v`, `cv_`, `atp.*synth`, `atp5`.

Gene set used: curated KEGG hsa00190 Complex V subunits (F0 sector: ATP5PB, ATP5MC1–3, ATP5ME/F/G, ATP5PD/F/O; F1 sector: ATP5F1A–E; peripheral stalk: ATP5MJ/K; assembly factors: ATPAF1, ATPAF2, TMEM70, DMAC2L, ATPSCKMT). This is distinct from GO:CC "ATPase complex" which includes SWI/SNF chromatin remodelers.

## Manuscript supplementary caption (`Complex_V_Gene_Expression_Heatmap.pdf` + `Complex_V_Pathway_Enrichment.pdf`)

**ATP synthase (Complex V) gene-level and pathway-level enrichment in DRP1 mutant cortical neurons.** (Left) Gene-level log2FC heatmap of curated KEGG hsa00190 Complex V subunits and assembly factors across three trajectory stages (Early = D35 mutant-vs-control, TrajDev = mutation-specific maturation interaction, Late = D65 mutant-vs-control) for G32A and R403C. Values are limma-voom coefficients; the blue–white–orange diverging scale spans ±0.6 log2FC. Genes are arranged by sub-complex (F1 catalytic core, FO membrane sector, peripheral stalk, assembly factors). (Right) GSEA dotplot of all MitoCarta Complex V pathways across the nine contrasts. Dot colour = NES (blue = downregulated, orange = upregulated); dot size = −log10(FDR); pathways are retained when retrieved by keywords `complex.*v`, `cv_`, `atp.*synth`, `atp5` in MitoCarta 3.0. Statistical significance is from fgsea (10,000 permutations, BH-corrected FDR; significant at FDR < 0.05). Both panels show the Compensation signature for Complex V in G32A — early downregulation (blue Early column) followed by positive TrajDev and near-baseline Late — consistent with the cell mounting a compensatory upregulation of the ATP synthesis machinery in response to mitochondrial dysfunction. The R403C signal is weaker but directionally identical. This panel supplies the gene-level evidence for the OXPHOS compensation arm of the translation paradox documented at the GSVA level in `../Critical_period_trajectories/Panel_G_OXPHOS_*.pdf`.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/2.5.viz_complex_v_analysis.R`
