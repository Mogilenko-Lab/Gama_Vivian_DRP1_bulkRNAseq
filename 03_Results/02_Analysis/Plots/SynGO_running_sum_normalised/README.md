# SynGO Running-Sum (Normalised) — Fig 5D Companion

GSEA running-sum enrichment panels for the seven SynGO synaptic-localized pathways across the four mutation-vs-control contrasts (G32A and R403C × D35 and D65). These panels use the same Tol/Wong colorblind-safe palette as the Fig 5E chord diagrams so the synaptic-translation figure family reads as a single visual identity.

## File inventory

| File | Description |
|---|---|
| `G32A_vs_Ctrl_D35_SynGO_running_sum_normalised.pdf` | Running sum panel: G32A vs Control at D35 |
| `G32A_vs_Ctrl_D65_SynGO_running_sum_normalised.pdf` | Running sum panel: G32A vs Control at D65 |
| `R403C_vs_Ctrl_D35_SynGO_running_sum_normalised.pdf` | Running sum panel: R403C vs Control at D35 |
| `R403C_vs_Ctrl_D65_SynGO_running_sum_normalised.pdf` | Running sum panel: R403C vs Control at D65 |
| `SynGO_running_sum_normalised_grid.pdf` | 2×2 composite of the four per-contrast panels (primary submission file) |

## Generating script

`02_Analysis/3.10.viz_syngo_running_sum_normalised.R`

```bash
Rscript 02_Analysis/3.10.viz_syngo_running_sum_normalised.R
```

Requires checkpoint `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds`. Uses the unified running-sum helper `01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R` with rasterised line geometry (`RASTER_DPI = 350`) so the dense per-gene tick marks stay light while axes, fonts, and panel borders remain vector.

## Pathways shown

Seven SynGO cellular-component sets (palette matched 1:1 to `02_Analysis/3.7.viz_chord_diagrams.py`):

| Pathway ID | Display name | Tol/Wong colour |
|---|---|---|
| `SYNGO:presyn_ribosome` | Presynaptic ribosome | muted rose `#CC6677` |
| `SYNGO:postsyn_ribosome` | Postsynaptic ribosome | muted purple `#AA4499` |
| `GO:0045202` | Synapse | sand `#DDCC77` |
| `GO:0099523` | Presynaptic cytosol | forest green `#117733` |
| `GO:0099524` | Postsynaptic cytosol | teal `#44AA99` |
| `GO:0014069` | Postsynaptic density | olive `#999933` |
| `GO:0045211` | Postsynaptic membrane | sky blue `#88CCEE` |

## How to read this folder

Each PDF carries one contrast (e.g. `G32A_vs_Ctrl_D35`); the composite `SynGO_running_sum_normalised_grid.pdf` is the main figure file for submission. The x-axis is the rank position in the limma t-statistic-ordered gene list (most-upregulated in the contrast at left, most-downregulated at right). The y-axis is the running enrichment score; the peak (or trough, for negative-NES pathways) is the enrichment score (ES). Vertical tick marks along each line show positions of gene-set members in the ranked list. Statistical thresholds are inherited from the upstream fgsea run in `02_Analysis/1.1.main_pipeline.R` (10,000 permutations, BH-corrected FDR; significance at FDR < 0.05).

## Manuscript figure caption (Fig 5D)

**(D) GSEA running-sum enrichment plots for seven SynGO synaptic-localized cellular-component pathways across the four mutation-vs-control contrasts (G32A and R403C, D35 and D65; 2×2 composite).** Each panel plots the cumulative running enrichment score (y-axis) along the limma t-statistic-ranked gene list (x-axis; left = up-regulated in mutant, right = down-regulated). Vertical tick marks along the bottom of each curve denote rank positions of gene-set members. Pathway lines share the Tol/Wong colorblind-safe palette used in the Fig 5E chord diagrams: focal presynaptic and postsynaptic ribosome (muted rose, muted purple) on top, with synaptic-compartment sets (sand, forest, teal, olive, sky) beneath. Statistical significance is from fgsea (10,000 permutations, BH-corrected FDR < 0.05) as reported in `master_gsea_table.csv`. The figure documents the biphasic synaptic-ribosome trajectory: at D35 the presynaptic and postsynaptic ribosome pathways enrich with strong positive NES (early peak in the running sum), while at D65 the same pathways flip to large-magnitude negative NES (late trough) — the Sign_reversal signature that anchors the synaptic translation paradox.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/3.10.viz_syngo_running_sum_normalised.R`
