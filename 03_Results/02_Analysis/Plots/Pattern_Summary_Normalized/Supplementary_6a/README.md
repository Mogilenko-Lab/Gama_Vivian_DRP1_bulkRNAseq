# Plots/Pattern_Summary_Normalized/Supplementary_6a — geometric NES_Early × NES_Late scatter

## Overview

Geometric interpretation of the trajectory-pattern classifier: each ever-significantly enriched pathway is plotted as a single point in (NES_Early, NES_Late) space, colored by super-category (Active_Compensation, Active_Reversal, Active_Progression, Passive, Late_onset, Other). The diagonal dotted line marks `Early = Late` (no trajectory change); the four quadrants show where each major pattern class lives. This panel turns the reviewer's "thresholds are arbitrary" objection into a verifiable visual partition: pattern boundaries map to interpretable regions of NES space rather than to opaque rules.

## Reviewer concern addressed

**concern 6a** — interaction-GSEA vs trajectories. See:
- `Manuscript/current_submission/compbio-REVIEW.md` lines 15–17 (reviewer comment)
- `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/concept.md` (motivation)
- `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/plan.md` (Choice 3 — geometric scatter)

## Universe / denominator — IMPORTANT

These figures use the **5,267 ever-significantly enriched pathway universe** (FDR<0.05 in any of the 9 GSEA contrasts: G32A_vs_Ctrl_{D35,D65}, R403C_vs_Ctrl_{D35,D65}, Maturation_{G32A,R403C}_specific, Time_{Ctrl,G32A,R403C}). Pathways with missing NES_Early or NES_Late are dropped (per-mutation n shown in panel titles); the actual point count therefore differs slightly between G32A and R403C panels.

**This subfolder's universe is different from its parent folder.** The parent
[`../README.md`](../README.md) explicitly scopes the parent folder to the **12,221-pathway universe** ("ALL tested pathways regardless of statistical significance"). The 5,267-restricted subset used here was added by the 2026-04-24 external audit so the geometric scatter aligns with RESULTS_combio.md L11/13 percentages and with the new dual-denominator table. See [`../../../Tables/README_pattern_summary_denominators.md`](../../../Tables/README_pattern_summary_denominators.md) for the full universe vocabulary.

## Files

| File | Purpose | Format | Size |
|---|---|---|---|
| `geometric_scatter_G32A.pdf` | Single-mutation panel (G32A); top 10 pathways by \|NES_TrajDev\| labelled | PDF | ~115 KB |
| `geometric_scatter_G32A.png` | Raster version of the same | PNG | ~610 KB |
| `geometric_scatter_R403C.pdf` | Single-mutation panel (R403C); same labelling rule | PDF | ~115 KB |
| `geometric_scatter_R403C.png` | Raster version | PNG | ~720 KB |
| `geometric_scatter_both_mutations.pdf` | Side-by-side combined panel (G32A left, R403C right); intended as Fig 5A right-bottom panel A4 and as Supplementary Fig. S9 | PDF | ~200 KB |
| `geometric_scatter_both_mutations.png` | Raster version | PNG | ~1.3 MB |

## Cross-references

- **Manuscript anchors that cite these figures:**
  - `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/for-the-paper.md` Edit 4 (Fig 5A caption) — names `geometric_scatter_both_mutations.pdf` as panel A4 and as Supp Fig S9
  - `for-the-paper.md` Edit 6 (Supp Fig S9 legend) — full legend for the supplementary version
  - `for-the-paper.md` reviewer-response letter — cites this scatter as the "geometric reframe" answer to the "arbitrary thresholds" critique
- **Composite-figure assembly:** [`../../Publication_Figures/Supplementary_6a/Fig5A_composite_plan.md`](../../Publication_Figures/Supplementary_6a/Fig5A_composite_plan.md)
- **Generator script:** `02_Analysis/6a.geometric_scatter.py`
- **Color / super-category source:** `01_Scripts/Python/pattern_definitions.py` (`SUPER_CATEGORY_COLORS`)

## How to regenerate

```bash
python3 02_Analysis/6a.geometric_scatter.py
```

Reads `master_gsea_table.csv`. Deterministic, ~5 seconds. Overwrites the six files above.

## Read-only constraints

- This script reads `master_gsea_table.csv` and `pattern_definitions.py` only — neither is modified.
- The parent `Pattern_Summary_Normalized/` folder's existing 12,221-universe figures are not touched; this subfolder is strictly additive.
