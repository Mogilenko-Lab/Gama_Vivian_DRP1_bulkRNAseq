# Supplementary_6b — Per-Database Pattern Heatmap (Supp Fig S7)

This directory contains the per-database trajectory-pattern heatmap that backs **Supplementary Fig S7** of the manuscript. It summarizes how the eight trajectory patterns are distributed across the 12 enrichment databases for both mutations.

## File inventory

| File | Description |
|---|---|
| `per_database_pattern_heatmap.pdf` | Supp Fig S7: heatmap of pattern counts per database (G32A and R403C panels), vector format |
| `per_database_pattern_heatmap.png` | Raster version (300 DPI) |

## Generating script

`02_Analysis/revision/supplements/6b.per_database_pattern_summary.py`

```bash
python3 02_Analysis/revision/supplements/6b.per_database_pattern_summary.py
```

Input: `03_Results/02_Analysis/master_gsea_table.csv`

## Reading guide

Rows = enrichment databases (Hallmark, KEGG, Reactome, GO:BP, GO:CC, GO:MF, WikiPathways, Canonical, CGP, TF, SynGO, MitoCarta). Columns = trajectory patterns (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient, Complex). Color encodes count per cell; two side-by-side panels show G32A and R403C independently.

Key patterns to read:
- **Compensation** is the dominant non-Complex pattern across all databases; larger absolute counts in databases with more pathways (GO:BP, Reactome, CGP).
- **Sign_reversal** is enriched in SynGO and MitoCarta relative to their database size — reflecting the structural-ribosome collapse and synaptic failure narrative.
- **Complex** typically dominates in large databases (GO:BP, CGP) because the strict trajectory-classification thresholds leave most pathways unclassified.

This figure is complementary to `../Pattern_Summary_Normalized/pattern_summary_normalized.pdf` (Fig 5A), which shows proportions rather than absolute counts and excludes Complex.

## Manuscript supplementary caption (Supp Fig S7)

**Trajectory-pattern composition by enrichment database in DRP1 mutant cortical neurons.** Two-panel count heatmap: rows are the 12 enrichment databases (Hallmark, KEGG, Reactome, GO:BP, GO:CC, GO:MF, WikiPathways, Canonical, CGP, TF, SynGO, MitoCarta 3.0); columns are the eight trajectory pattern classes (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient, Complex) as defined in `01_Scripts/Python/pattern_definitions.py`. Each cell encodes the absolute count of pathways from that database assigned to that pattern; G32A and R403C are shown side-by-side. Statistical thresholds inherited from the upstream fgsea pipeline (10,000 permutations, BH-FDR < 0.05 for the Active_* super-categories' TrajDev requirement). The panel supports three observations cited in the supplementary text: (i) Compensation is the dominant non-Complex pattern across every database for both mutations, with absolute counts scaling with database size (GO:BP, Reactome, CGP carry the highest counts); (ii) Sign_reversal is enriched in SynGO and MitoCarta relative to their database size, reflecting the synaptic-ribosome collapse and the parallel structural-ribosome dynamics; (iii) Complex predominates in large databases (GO:BP, CGP) because the strict trajectory-classifier thresholds leave most pathways unclassified. The proportions counterpart of this panel (excluding Complex) is `../Pattern_Summary_Normalized/pattern_summary_normalized.pdf` (Fig 5A).

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/revision/supplements/6b.per_database_pattern_summary.py`
