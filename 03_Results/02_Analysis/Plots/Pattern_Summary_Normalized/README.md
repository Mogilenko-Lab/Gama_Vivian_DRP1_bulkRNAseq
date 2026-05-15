# Pattern Summary Normalized — Trajectory Pattern Distribution (Fig 5A)

This directory contains normalized pattern-distribution bar charts summarizing temporal trajectory classifications across all 12,221 GSEA pathways tested. The main figure backs **Fig 5A** of the manuscript.

**Important scope note**: These figures show proportions across **all 12,221 tested pathways** (the all-contrasts universe), regardless of statistical significance. The "Complex" pattern category is excluded from the bars to focus on pathways with interpretable trajectories; the `n=` annotation on each bar row gives the count of non-Complex pathways per database. This is distinct from the 5,267-pathway ever-significantly-enriched universe used in `Supplementary_6a/`.

## File inventory

| File | Description |
|---|---|
| `pattern_summary_normalized.pdf` | Main figure: 100% stacked bars for both mutations (G32A and R403C), one row per database, Complex excluded |
| `pattern_summary_normalized.png` | Raster version (300 DPI) |
| `pattern_comparison_dual_G32A.pdf` | Dual-panel: normalized proportions (left) and absolute counts (right) for G32A |
| `pattern_comparison_dual_G32A.png` | Raster version |
| `pattern_comparison_dual_R403C.pdf` | Dual-panel: normalized proportions (left) and absolute counts (right) for R403C |
| `pattern_comparison_dual_R403C.png` | Raster version |

### Subdirectory: Supplementary_6a/

Contains geometric NES_Early × NES_Late scatter plots using the **5,267-pathway ever-significantly-enriched universe**. See `Supplementary_6a/README.md`.

## Generating script

`02_Analysis/3.4.pattern_summary_normalized.py`

```bash
python3 02_Analysis/3.4.pattern_summary_normalized.py
```

Input: `03_Results/02_Analysis/master_gsea_table.csv`

## Reading guide

**pattern_summary_normalized.pdf**: Each horizontal bar represents one database. Bar segments show proportions of the eight trajectory patterns (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient — Complex excluded). Proportions are computed within the non-Complex subset for that database; the `n=` label on the right gives the absolute count of non-Complex pathways.

Pattern taxonomy:
- **Compensation**: TrajDev significantly opposes the Early defect (most common interpretable pattern, ~12–13% of all tested)
- **Sign_reversal**: Sign flips from Early to Late (~2–3%)
- **Progressive**: TrajDev amplifies the Early defect (rare, < 1%)
- **Natural_improvement / Natural_worsening**: Improvement or worsening without significant TrajDev
- **Late_onset**: No Early effect; defect emerges only at Late (~< 1%)
- **Transient**: Early defect resolved by Late (very rare)

**Dual-panel figures**: Left panel = normalized proportions (same as main figure but single-mutation); right panel = absolute counts, showing the raw frequency context for each pattern.

**Interpretation caveat**: ~78% of all tested pathways are classified as Complex (do not meet the criteria for any of the eight defined patterns). Bars represent only the remaining ~22% (G32A) or ~28% (R403C). Do not interpret "X% of pathways show Compensation" without specifying "of non-Complex pathways."

## Manuscript figure caption (Fig 5A, `pattern_summary_normalized.pdf`)

**Trajectory-pattern distribution across 12 enrichment databases in DRP1 mutant cortical neurons.** Horizontal 100% stacked bar chart of pathway trajectory patterns, one row per enrichment database (Hallmark, KEGG, Reactome, GO:BP, GO:CC, GO:MF, WikiPathways, Canonical, CGP, TF, SynGO, MitoCarta), separately for G32A and R403C panels. Each bar shows the proportion of non-Complex pathways assigned to each of seven pattern classes (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient; Complex is excluded from the bars but its count is reported in the `n=` annotation alongside each row). Pattern assignments follow the criteria in `01_Scripts/Python/pattern_definitions.py`. Statistical thresholds are inherited from the upstream fgsea pipeline (10,000 permutations, BH-corrected FDR; significance at FDR < 0.05 for the Active_Compensation / Active_Reversal / Active_Progression super-categories; non-significant TrajDev defines the Passive Natural_* patterns). The universe (n = 12,221 pathways × 2 mutations) is all pathways tested for GSEA regardless of significance; ~78% are unclassifiable as Complex and are excluded from the proportions. The panel supports two manuscript claims: (i) Compensation is the dominant interpretable pattern across every database in both mutations, and (ii) Sign_reversal is preferentially enriched in SynGO and MitoCarta relative to their database size — the database-level fingerprint of the translation paradox. The geometric companion panel in `Supplementary_6a/` plots the same pathways in (NES_Early, NES_Late) coordinates against the 5,267-pathway ever-significantly-enriched subset.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/3.4.pattern_summary_normalized.py`
