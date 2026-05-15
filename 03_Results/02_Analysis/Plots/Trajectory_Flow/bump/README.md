# Trajectory_Flow/bump — Bump Chart Variants (Non-Key Figures)

This subdirectory contains alternative bump chart visualizations used during analysis to explore different visual encodings. The primary paper figures (`bump_curved_nes_significant.pdf` for Fig 5B and `bump_focused_FINAL_paper_combined.pdf` for Fig 5C) remain in the parent `../` directory.

## File naming convention

```
bump_{scope}_{style}_{metric}.{ext}
```

- **scope**: `focused` (MitoCarta + SynGO, ~100–200 pathways) or `significant` (all non-Complex patterns, ~2,000–4,000 pathways)
- **style**: `uniform`, `weighted`, `highlight`, `curved`, `FINAL_paper_combined`
- **metric**: `nes` (Normalized Enrichment Score on y-axis) or `rank` (relative rank position on y-axis)
- **ext**: `pdf` or `png`

## File inventory

| File | Description |
|---|---|
| `bump_curved_nes_focused.pdf/.png` | Focused scope, curved lines (TrajDev Bezier), NES y-axis |
| `bump_curved_nes_significant.pdf/.png` | Significant scope, curved lines, NES y-axis (note: version with paper labels is in parent dir) |
| `bump_focused_highlight_nes.pdf/.png` | Focused scope, weighted lines with key pathway labels, NES y-axis |
| `bump_focused_uniform_nes.pdf/.png` | Focused scope, uniform line width, NES y-axis |
| `bump_focused_uniform_rank.pdf/.png` | Focused scope, uniform line width, rank y-axis |
| `bump_focused_weighted_nes.pdf/.png` | Focused scope, weighted line width, NES y-axis |
| `bump_focused_weighted_rank.pdf/.png` | Focused scope, weighted line width, rank y-axis |
| `bump_significant_curved_nes.pdf/.png` | Significant scope, curved lines, NES y-axis |
| `bump_significant_curved_rank.pdf/.png` | Significant scope, curved lines, rank y-axis |
| `bump_significant_FINAL_paper_combined.pdf/.png` | Significant scope version of combined weighted+curved+labels figure |
| `bump_significant_highlight_nes.pdf/.png` | Significant scope with pathway labels, NES y-axis |
| `bump_significant_uniform_nes.pdf/.png` | Significant scope, uniform line width, NES y-axis |
| `bump_significant_uniform_rank.pdf/.png` | Significant scope, uniform line width, rank y-axis |
| `bump_significant_weighted_nes.pdf/.png` | Significant scope, weighted line width, NES y-axis |
| `bump_significant_weighted_rank.pdf/.png` | Significant scope, weighted line width, rank y-axis |
| `bump_uniform_nes_focused.pdf/.png` | Alias variant (uniform/NES/focused) |
| `bump_uniform_nes_significant.pdf/.png` | Alias variant (uniform/NES/significant) |
| `bump_uniform_rank_focused.pdf/.png` | Alias variant (uniform/rank/focused) |
| `bump_uniform_rank_significant.pdf/.png` | Alias variant (uniform/rank/significant) |
| `bump_weighted_nes_focused.pdf/.png` | Alias variant (weighted/NES/focused) |
| `bump_weighted_nes_significant.pdf/.png` | Alias variant (weighted/NES/significant) |
| `bump_weighted_rank_focused.pdf/.png` | Alias variant (weighted/rank/focused) |
| `bump_weighted_rank_significant.pdf/.png` | Alias variant (weighted/rank/significant) |

## Generating scripts

- `02_Analysis/3.7.viz_bump_chart.py` — static variants
- `02_Analysis/3.8.viz_interactive_bump.py` — interactive HTML variants (if present)

## Visual encoding notes

**uniform**: All lines same width and opacity — baseline with no pattern-frequency weighting.

**weighted**: Line width inversely proportional to pattern frequency. Dominant patterns (> 30%, e.g., Compensation) are thin/low-opacity; rare patterns (< 1%, e.g., Sign_reversal, Progressive) are thick/full-opacity. This ensures rare trajectories are visible against the mass.

**highlight**: Weighted lines with key pathway labels added using adjustText for non-overlapping placement. Labels prioritize semantic categories (Mitochondrial > Synaptic > Ribosomal > Other).

**curved**: Lines rendered as quadratic Bezier curves. The control point displacement from the straight line encodes the TrajDev NES magnitude and direction. Curves only shown for pathways with significant TrajDev (FDR < 0.05); others are straight lines.

**NES vs rank**: NES y-axis preserves absolute effect magnitudes and allows reading enrichment strength directly. Rank y-axis normalizes for cross-database NES differences and emphasizes relative position changes.

## How to read this folder

For the manuscript figures jump up one level (`../bump_curved_nes_significant.pdf` for Fig 5B and `../bump_focused_FINAL_paper_combined.pdf` for Fig 5C). The variants in this folder are diagnostic / alternative encodings used during figure selection — they are kept in the repository so a reviewer can audit the visual choices that led to the manuscript versions, and so the `weighted_rank_*` and `uniform_rank_*` alternatives can be substituted if the editor prefers a rank-y-axis view. None of the files here are independently cited in the main text or supplement.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/3.7.viz_bump_chart.py`
