# Mito Translation Cascade — Fig 6B Heatmaps

Three-module gene-level log2FC heatmap covering ATP synthase (Complex V), curated calcium-signaling genes, and the mitochondrial central-dogma machinery (mtDNA maintenance, transcription, mitoribosome assembly). The main heatmap backs **Fig 6B** of the manuscript and complements the Fig 5F synaptic-ribosome heatmap by showing the mitochondrial compensation arm of the translation paradox.

## File inventory

| File | Description |
|---|---|
| `Mechanistic_Cascade_Heatmap.pdf` | Fig 6B: 57-gene log2FC heatmap across Early / TrajDev / Late for G32A and R403C; rows split into three modules (ATP Synthase, Calcium Signaling, Mt Central Dogma) with within-module hierarchical clustering |
| `Module_Summary_Heatmap.pdf` | Supplementary: per-module mean log2FC summary with numeric annotations |

## Generating script

`02_Analysis/2.2.viz_mito_translation_cascade.R`

```bash
Rscript 02_Analysis/2.2.viz_mito_translation_cascade.R
```

Requires checkpoints `03_Results/02_Analysis/checkpoints/fit_object.rds` (limma model fit) and `03_Results/02_Analysis/checkpoints/mitocarta_gsea_results.rds` (MitoCarta GSEA).

## Modules and gene sets

| Module | n genes | Source | Notes |
|---|---|---|---|
| **ATP Synthase (Complex V)** | 20 | Curated KEGG hsa00190 Complex V subunits | F1 (ATP5F1A–E), FO (ATP5MC1–3, ATP5ME/F/G, ATP5MJ/K, ATP5PB/D/F/O), assembly (ATPAF1, ATPAF2, TMEM70). Distinct from GO:CC ATPase complex (which contains chromatin-remodeling ATPases). |
| **Calcium Signaling** | 12 | Hypothesis-driven panel including NNAT and PNPO | NNAT, PNPO, CACNG3, CACNA1S, ATP2A1, RYR1, MYLK3, VDR, STIM1, STIM2, CALB1, CALR. ORAI1 was not detected in this dataset. |
| **Mt Central Dogma** | 25 | Top-25 |logFC| genes from MitoCarta Mitochondrial_central_dogma core enrichment for Maturation_G32A_specific (parent pathway: NES = 2.41, p.adj = 2.14×10⁻¹³, set size 216) | mtDNA maintenance (POLQ, DNA2, RECQL4, PIF1, UNG); replication (PRIMPOL, EXOG, ENDOG); transcription (MTERF3, MTG1); mt-tRNA modification (PUS1, MRM2, TRMT61B, DARS2); mitoribosome (MRPL1/11/15, MRPS11/17/18C); other (DDX28, GUF1, MPV17L2, PDF, RBFA). |

## How to read this folder

Open `Mechanistic_Cascade_Heatmap.pdf` for the panel that goes into Fig 6B. `Module_Summary_Heatmap.pdf` is a numeric-summary supplement that aggregates the per-module pattern (one row per module, three columns per mutation) for at-a-glance comparison and is referenced in the supplementary text. Gene-level statistical detail is not on the heatmaps — refer to `../../DE_results/` for per-gene FDR.

## Manuscript figure caption (Fig 6B)

**(B) Gene-level log2 fold-change trajectories for ATP synthase, calcium-signaling, and mitochondrial central-dogma genes in DRP1 mutant cortical neurons.** Heatmap of 57 genes (rows) across six trajectory cells (columns): Early (D35 mutant vs control), TrajDev (mutation-specific maturation interaction term, `(MutD65 − MutD35) − (CtrlD65 − CtrlD35)`), and Late (D65 mutant vs control), shown separately for G32A and R403C. Rows are split into three biologically curated modules (ATP synthase / Complex V, n = 20, curated KEGG hsa00190 subunits and assembly factors; calcium signaling, n = 12, hypothesis-driven panel; mitochondrial central dogma, n = 25, top-effect genes from MitoCarta `Mitochondrial_central_dogma` core enrichment) with hierarchical clustering (Euclidean distance, complete linkage) applied within each module. Color encodes log2FC from the limma-voom fit (`fit$coefficients`); the blue–white–orange diverging scale spans ±0.6 (values outside this range are clipped to the endpoints). Significance thresholds inherited from the upstream limma-voom analysis (BH-FDR < 0.05); module-level enrichment for the central-dogma module reached p.adj = 2.14×10⁻¹³ in the Maturation_G32A_specific contrast (n = 25 top genes selected by |logFC|, original set size 216). The panel supports the conclusion that DRP1 mutants mount a compensatory upregulation of mitochondrial biogenesis (TrajDev = orange in ATP synthase and central-dogma modules) while the calcium-signaling module shows progressive downregulation deepening from Early to Late (NNAT and PNPO drive the signal). Together with Fig 5F (synaptic-ribosome collapse), this panel establishes the dissociation between mitochondrial translation compensation and synaptic translation failure.

### Module Summary supplementary caption

**Module-level summary of expression trajectories.** Mean log2FC values averaged across all genes within each module are displayed in cells. Positive TrajDev values for ATP Synthase (G32A: +0.14) and Mt Central Dogma (G32A: +0.69) indicate compensatory upregulation during maturation; Calcium Signaling is uniformly negative (Early: −0.39, TrajDev: −0.52, Late: −0.91 for G32A) indicating progressive deficit without compensation.

## Related visualisations

| Sibling folder | Relationship |
|---|---|
| `../Synaptic_ribosomes/Panel_C_Expression_Heatmap.pdf` | Sister Fig 5F panel, same style, synaptic translation genes |
| `../Complex_V_analysis/` | Deep-dive on the Complex V subunit GSEA |
| `../Critical_period_trajectories/` | GSVA-level view of the same modules |

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/2.2.viz_mito_translation_cascade.R`
