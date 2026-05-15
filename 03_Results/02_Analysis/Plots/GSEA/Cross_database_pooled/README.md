# GSEA/Cross_database_pooled — Cross-Database Pooled Dotplots

Per-contrast GSEA dotplots that pool significant pathways across multiple enrichment databases into a single unified view. Used in supplementary discussion of the cross-database convergence of the translation-paradox finding.

## File inventory

| Path | Description |
|---|---|
| `focused/{contrast}_pooled_focused.pdf` | Manuscript-ready pooled dotplot per contrast (4 databases: KEGG, Reactome, SynGO, MitoCarta) |
| `comprehensive/{contrast}_pooled_comprehensive.pdf` | Exploratory pooled dotplot per contrast (11 databases: focused + Hallmark, GO:BP/CC/MF, WikiPathways, Canonical, TF) |
| `csv_data/{contrast}_raw_{scope}.csv`, `{contrast}_filtered_{scope}.csv` | Intermediate raw and post-filter tables used to build each dotplot |
| `csv_data/FILTERING_STRATEGY.txt` | Plain-text record of the neuronal-pathway include / exclude rules |
| `test_output/TEST_*.pdf` | Diagnostic dotplots used to verify the colour-gradient and significance-outline encoding (kept for reviewer reproducibility, not cited) |
| `pooled_dotplots_summary.txt` | Legacy summary text |

Each `focused/` and `comprehensive/` subdir contains nine PDFs (one per contrast: G32A_vs_Ctrl_{D35,D65}, R403C_vs_Ctrl_{D35,D65}, Time_{Ctrl,G32A,R403C}, Maturation_{G32A,R403C}_specific).

## Generating script

`02_Analysis/3.9.viz_pooled_dotplots.R`

```bash
Rscript 02_Analysis/3.9.viz_pooled_dotplots.R
```

Requires checkpoints `checkpoints/all_gsea_results.rds`, `checkpoints/syngo_gsea_results.rds`, `checkpoints/mitocarta_gsea_results.rds`. Helper functions in `01_Scripts/R_scripts/gsea_dotplot_helpers.R`; colour palette from `01_Scripts/R_scripts/color_config.R`.

## Encoding

| Channel | Variable |
|---|---|
| Y-axis | Pathway description, sorted by gene ratio |
| X-axis | Gene ratio (leading-edge genes / set size) |
| Dot colour | NES (blue–white–orange diverging; `#2166AC` → `#F7F7F7` → `#B35806`) |
| Dot size | −log10(FDR) |
| Dot border | Black outline = BH-FDR < 0.05; no border = 0.05 ≤ FDR < 0.10 |

## Filter and pooling rules

- FDR cutoff: < 0.05 for inclusion in the pool (also drives black-outline flag)
- Per-database top N: 10 (focused) / 5 (comprehensive), ranked by FDR
- Neuronal-pathway filter: excludes pathways with non-neuronal tissue keywords (pancreatic, cardiac, kidney, intestinal, lung, liver, breast/prostate, adipogenesis, xenobiotic, coagulation); always retains SynGO pathways and pathways with neural/synaptic/mitochondrial/cell-cycle keywords.
- NES colour-scale cap: auto-calculated up to a max of |NES| = 4 to keep the scale comparable across panels.

## How to read this folder

Open the `focused/` PDFs for the manuscript-grade per-contrast snapshot. Use `comprehensive/` PDFs for unbiased exploration including the GO and TF databases. The `csv_data/` tables document the raw vs filtered pathway lists so a reviewer can audit any filter decision. Statistical significance is from the upstream fgsea pipeline (10,000 permutations, BH-FDR; see `../../README.md` for the canonical method description).

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/3.9.viz_pooled_dotplots.R`
