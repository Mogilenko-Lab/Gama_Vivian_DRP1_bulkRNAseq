# Pattern_Summary_Normalized/Supplementary_6a — Geometric NES_Early × NES_Late Scatter

## Overview

Geometric reinterpretation of the trajectory-pattern classifier: each ever-significantly enriched pathway is plotted as a single point in (NES_Early, NES_Late) space, coloured by super-category (Active_Compensation, Active_Reversal, Active_Progression, Passive, Late_onset, Other). The diagonal `Early = Late` reference line marks "no trajectory change"; the four quadrants partition pattern-class boundaries into interpretable regions of NES space. This panel addresses reviewer concern **6a** ("classifier thresholds are arbitrary") by mapping each pattern label to a geometrically interpretable region, and is the panel A4 source for the manuscript's revised Fig 5A composite. The manuscript revision letter cites it as the geometric companion to the classifier-summary panels in the parent `../` folder.

## Universe / denominator — IMPORTANT

These figures use the **5,267 ever-significantly-enriched pathway universe** (FDR < 0.05 in any of the 9 GSEA contrasts: G32A_vs_Ctrl_{D35,D65}, R403C_vs_Ctrl_{D35,D65}, Maturation_{G32A,R403C}_specific, Time_{Ctrl,G32A,R403C}). Pathways with missing NES_Early or NES_Late are dropped (per-mutation n is shown in panel titles); the actual point count therefore differs slightly between G32A and R403C panels.

**This subfolder's universe differs from its parent folder.** The parent `../README.md` scopes its pattern-summary bars to the **12,221-pathway "all tested" universe**. The 5,267-restricted subset used here was added by the 2026-04-24 external audit so the geometric scatter aligns with the RESULTS percentages cited in the manuscript and with the dual-denominator table. See `../../../Tables/README_pattern_summary_denominators.md` for the universe vocabulary.

## File inventory

| File | Purpose |
|---|---|
| `geometric_scatter_G32A.pdf` / `.png` | Single-mutation panel, G32A; top-10 pathways by |NES_TrajDev| labelled |
| `geometric_scatter_R403C.pdf` / `.png` | Single-mutation panel, R403C; same labelling rule |
| `geometric_scatter_both_mutations.pdf` / `.png` | Side-by-side combined panel (G32A left, R403C right); manuscript Fig 5A panel A4 source |

PDFs are vector; PNGs are 300 dpi rasters.

## Generating script

`02_Analysis/revision/supplements/6a.geometric_scatter.py`

```bash
python3 02_Analysis/revision/supplements/6a.geometric_scatter.py
```

Reads `03_Results/02_Analysis/master_gsea_table.csv` and pattern super-category colours from `01_Scripts/Python/pattern_definitions.py` (`SUPER_CATEGORY_COLORS`). Deterministic, ~5 seconds. Overwrites the six files above.

## How to read this folder

Open `geometric_scatter_both_mutations.pdf` for the manuscript-ready combined panel. The per-mutation PDFs are used in supplementary breakouts. Pattern boundaries map to NES-space regions as follows: Active_Compensation lies on the side of the diagonal where the Late NES has moved toward zero from the Early NES; Active_Reversal lies across-quadrant from Early; Passive lies along the diagonal (no significant TrajDev); Late_onset lies along the NES_Late axis with near-zero NES_Early. The dotted `Early = Late` line is the geometric locus of "no trajectory change". Member-pathway counts and threshold definitions are in the parent `../README.md` and in `01_Scripts/Python/pattern_definitions.py`.

## Manuscript figure caption (Fig 5A panel A4 / supplementary geometric companion)

**Geometric reinterpretation of the trajectory-pattern classifier in NES space.** Each point represents one ever-significantly-enriched GSEA pathway (FDR < 0.05 in at least one of the nine experimental contrasts; total universe 5,267; per-panel n shown in title), plotted at its (NES_Early, NES_Late) coordinates and coloured by trajectory super-category from `pattern_definitions.py` (Active_Compensation, Active_Reversal, Active_Progression, Passive, Late_onset, Other). The dotted diagonal marks `NES_Late = NES_Early` — the geometric locus of "no trajectory deviation" — so that vertical distance from the diagonal encodes the magnitude of the trajectory-deviation (TrajDev) contrast directly. Active_Compensation populates the region where Late NES has moved toward zero from a non-zero Early NES; Active_Reversal sits across the origin from its Early position (sign-flipped); Passive entries cluster along the diagonal (Early ≈ Late, no significant TrajDev); Late_onset entries lie along the NES_Late axis with near-zero NES_Early. The top ten pathways by |NES_TrajDev| are labelled per panel. Side-by-side G32A (left) and R403C (right) panels show that the ribosome-translation Sign_reversal pathways (postsynaptic and presynaptic ribosome) occupy the extreme Active_Reversal region in both mutations, while the dominant Active_Compensation cluster is populated by mitochondrial-translation and ATP-synthesis pathways. The geometric mapping demonstrates that the classifier's discrete pattern labels correspond to interpretable regions of a continuous two-dimensional trajectory space rather than to opaque rule firings.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/revision/supplements/6a.geometric_scatter.py`
