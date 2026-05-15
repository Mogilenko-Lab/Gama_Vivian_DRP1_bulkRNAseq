# Synaptic Ribosomes — Expression Heatmap (Fig 5F)

This directory contains a gene-level log2FC heatmap for synaptic ribosome genes from SynGO annotations, backing **Fig 5F** of the manuscript. The figure documents the pan-synaptic collapse of ribosome programs during the D35 → D65 critical period in both DRP1 mutations.

## File inventory

| File | Description |
|---|---|
| `Panel_C_Expression_Heatmap.pdf` | Fig 5F: log2FC heatmap for 67 synaptic ribosome genes across trajectory stages (Early / TrajDev / Late) for G32A and R403C |
| `Panel_C_Sample_Heatmap.pdf` | Companion per-sample expression heatmap: rows = 67 genes, columns = individual biological replicates ordered by genotype × timepoint; values = log2-CPM z-scores. Shows replicate-level concordance behind the contrast-level Panel_C signal. |

## Generating script

`02_Analysis/2.3.viz_synaptic_ribosomes.R`

```bash
Rscript 02_Analysis/2.3.viz_synaptic_ribosomes.R
```

Requires checkpoints: `checkpoints/syngo_gsea_results.rds`, `checkpoints/syngo_lists.rds`, `checkpoints/fit_object.rds`, and gene-list text files in `Verification_reports/`.

## Gene sets

- **Presynaptic ribosome (SynGO)**: 52 genes; GSEA TrajDev NES = −2.90 (G32A), −2.71 (R403C); FDR ≈ 10⁻¹²–10⁻¹⁰
- **Postsynaptic ribosome (SynGO)**: 70 genes; GSEA TrajDev NES = −3.02 (G32A), −2.89 (R403C); FDR ≈ 10⁻¹⁵–10⁻¹²
- All 52 presynaptic genes are a subset of the 70 postsynaptic genes (100% overlap). 18 postsynaptic-only genes are annotated exclusively in the dendritic compartment.
- Heatmap shows 67 genes after filtering for detected expression (14 postsynaptic-only + 53 shared).
- All 70 postsynaptic-ribosome genes are members of the curated cytoplasmic ribosomal gene set (Jaccard = 0.574), establishing this as a cytoplasmic structural-ribosome collapse rather than a synapse-specific phenomenon.

## Reading guide

**Structure**: Rows = 67 ribosomal genes split into two sections — "Postsynaptic Only" (14 genes, top) and "Both Compartments" (53 genes, bottom). Columns = trajectory stages: Early (D35 mutation vs control), TrajDev (mutation-specific maturation interaction), Late (D65 mutation vs control), separately for G32A (left) and R403C (right).

**Values displayed**: Log2FC coefficients directly from the limma-voom fitted model (`fit$coefficients`). NOT normalized expression or z-scores.

**Color scale**: Blue-White-Orange diverging palette; range ±0.6 log2FC. Blue = downregulated vs control; white = no change; orange = upregulated vs control.

**Key patterns**: Early columns (D35) are mostly white or light-colored — minimal mutation effect before the critical period. TrajDev and Late columns are deep blue — strong and persistent downregulation emerges during D35 → D65 maturation. Both mutations and both compartment groups show nearly identical profiles, indicating a mutation-independent, pan-synaptic ribosome failure.

**Clustering**: Hierarchical clustering (Euclidean distance, complete linkage) applied within each compartment slice independently. Compartment order is fixed (Postsynaptic Only on top).

## Manuscript figure caption (Fig 5F)

**(F) Pan-synaptic ribosomal-gene collapse emerges during the D35 → D65 critical period in both DRP1 mutations.** Heatmap of 67 SynGO-annotated ribosomal protein genes (rows; 14 postsynaptic-only, top split; 53 shared presynaptic+postsynaptic, bottom split) across six trajectory cells (columns): Early (D35 mutant vs control), TrajDev (mutation-specific maturation interaction term), and Late (D65 mutant vs control), shown separately for G32A and R403C. Values are limma-voom coefficients (log2FC) drawn directly from `fit$coefficients`. The blue–white–orange diverging colour scale spans ±0.6 log2FC. Within-slice hierarchical clustering uses Euclidean distance and complete linkage. Significance for the parent pathways (from fgsea, 10,000 permutations, BH-FDR) is documented elsewhere: presynaptic ribosome NES_TrajDev = −2.90 (G32A) and −2.71 (R403C); postsynaptic ribosome NES_TrajDev = −3.02 (G32A) and −2.89 (R403C); both at FDR < 10⁻¹⁰ in each contrast. The figure supports the Sign_reversal trajectory classification for the synaptic-ribosome program: minimal effect at Early, deep blue at TrajDev and Late in both compartment slices and both mutations, indicating a mutation-independent collapse that emerges only during maturation. The full SynGO presynaptic-ribosome gene set (52 genes) is a strict subset of the postsynaptic-ribosome gene set (70 genes), and all 70 are members of the curated cytoplasmic ribosomal proteome (Jaccard = 0.574 vs GO:CC cytosolic ribosome) — see `../Supplementary_9/` for the cross-compartment quantification.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/2.3.viz_synaptic_ribosomes.R`
