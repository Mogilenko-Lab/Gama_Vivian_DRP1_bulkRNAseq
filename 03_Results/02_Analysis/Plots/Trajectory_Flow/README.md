# Trajectory Flow — Bump Charts and Alluvial Diagrams (Fig 5B, Fig 5C)

This directory holds visualizations of pathway trajectory dynamics across the Early (D35) → TrajDev → Late (D65) developmental axis. Two primary figure types are present: bump charts (static and interactive) and alluvial/Sankey diagrams. The key paper figures are:

- **Fig 5B** → `bump_curved_nes_significant.pdf`
- **Fig 5C** → `bump_focused_FINAL_paper_combined.pdf`

## File inventory (root directory)

| File | Description |
|---|---|
| `interactive_bump_dashboard.html` | Comprehensive interactive explorer: filter by pattern/database/significance, toggle NES vs rank, color by pattern, hover for per-pathway stats |
| `bump_curved_nes_significant.pdf` | **Fig 5B**: all significant pathways (non-Complex) with Bezier TrajDev curves, NES y-axis |
| `bump_curved_nes_significant.png` | Raster version |
| `bump_focused_FINAL_paper_combined.pdf` | **Fig 5C**: MitoCarta + SynGO focused pathways, weighted lines + curves + key pathway labels combined |
| `bump_focused_FINAL_paper_combined.png` | Raster version |
| `bump_focused_curved_nes.pdf` | Focused modules with Bezier curves, NES y-axis (unfocused variant of Fig 5C) |
| `bump_focused_curved_nes.png` | Raster version |

### bump/ subdirectory

Contains bump chart variants used during analysis but not selected as primary figures. See `bump/README.md` for the naming convention and variant inventory.

### alluvial/ subdirectory

Contains alluvial/Sankey diagrams (Plotly HTML and static PDF) showing how early defects branch into different trajectory outcomes. See `alluvial/README.md`.

## Generating scripts

| Script | Output |
|---|---|
| `02_Analysis/3.7.viz_bump_chart.py` | Static bump charts (all variants in root + bump/) |
| `02_Analysis/3.8.viz_interactive_bump_dashboard.py` | `interactive_bump_dashboard.html` |
| `02_Analysis/3.8.viz_interactive_bump.py` | Interactive variants in bump/ |
| `02_Analysis/3.5.viz_trajectory_flow.py` | Alluvial HTML diagrams in alluvial/ |
| `02_Analysis/3.6.viz_alluvial_ggalluvial.R` | Classical ggalluvial PDFs in alluvial/ |

## Reading guide

**Bump charts (Fig 5B, Fig 5C)**: X-axis spans Early (D35) to Late (D65). Y-axis = NES (or rank in rank-mode figures). Each line is one pathway. Line curvature encodes the TrajDev magnitude and direction (upward bulge = positive TrajDev = pathway was upregulated more in mutants during maturation; downward bulge = negative TrajDev). Curves only appear for pathways with significant TrajDev (FDR < 0.05). Straight lines indicate no significant developmental trajectory deviation.

**Line weight encoding (weighted/FINAL variants)**: Line width is inversely proportional to pattern frequency — dominant patterns (Compensation, > 30% of pathways) are thin and low-opacity (background); rare patterns (Sign_reversal, Progressive, < 1%) are thick and fully opaque (foreground), ensuring rare trajectories are visible against the mass.

**Color by pattern**: Compensation = green; Sign_reversal = orange/red; Progressive = purple; Natural_improvement = light blue; Late_onset = yellow. See `01_Scripts/Python/pattern_definitions.py` for canonical color mapping.

**Key finding**: The Sign_reversal lines for SynGO synaptic-ribosome pathways (postsynaptic and presynaptic ribosome) arch visibly downward through TrajDev — the largest-amplitude downward curves in the focused panel — while MitoCarta mitochondrial pathways arch upward (Compensation). This divergence is the visual summary of the translation paradox.

**Interactive dashboard**: Click any pathway line to see per-pathway stats (NES at each stage, FDR, pattern, database). Use the filter panel to isolate specific databases or patterns.

## Manuscript figure captions

### Fig 5B (`bump_curved_nes_significant.pdf`)

**Pathway-level trajectory bump chart for all significantly enriched non-Complex pathways across the Early → TrajDev → Late developmental axis in DRP1 mutant cortical neurons.** Each line represents one GSEA pathway (n shown in panel header), with x-axis = trajectory stage (Early = D35 mutant-vs-control, TrajDev = mutation-specific maturation interaction, Late = D65 mutant-vs-control) and y-axis = GSEA NES. Lines are drawn as quadratic Bezier curves whose control-point displacement encodes the TrajDev NES magnitude and direction; pathways without significant TrajDev (FDR ≥ 0.05) are drawn as straight Early-to-Late segments. Line colour encodes trajectory pattern (Compensation = green; Sign_reversal = orange/red; Progressive = purple; Natural_improvement = light blue; Late_onset = yellow; canonical palette in `01_Scripts/Python/pattern_definitions.py`). Line width is inversely proportional to pattern frequency in the full universe so that rare pattern classes (Sign_reversal, Progressive; < 1% of pathways) remain visible against the Compensation mass (~30% of pathways). G32A and R403C panels are drawn side-by-side. The figure provides the system-level visual summary of the trajectory-pattern landscape that the rest of Fig 5 dissects at the per-module level.

### Fig 5C (`bump_focused_FINAL_paper_combined.pdf`)

**Trajectory bump chart restricted to the focused MitoCarta + SynGO pathway panel.** Same axis encoding and line-style conventions as Fig 5B but with pathway scope restricted to the MitoCarta mitochondrial pathway set (n = 149) and the SynGO synaptic-localisation set (~300 pathways); the combined "weighted + curved + key-pathway-label" rendering is the manuscript-version composite. Key pathway labels (postsynaptic ribosome, presynaptic ribosome, mitochondrial ribosome, OXPHOS, ribosome biogenesis) are placed using adjustText with semantic priority (Mitochondrial > Synaptic > Ribosomal > Other). The panel makes visible the bidirectional split that defines the translation paradox: SynGO synaptic-ribosome pathways arch sharply downward through TrajDev (the largest-amplitude downward curves in the panel, NES_TrajDev ≈ −3.0), while MitoCarta mitochondrial-translation and OXPHOS pathways arch upward (Compensation pattern, NES_TrajDev ≈ +1.5 to +2.5). Statistical thresholds inherited from the upstream fgsea pipeline (10,000 permutations, BH-FDR < 0.05) and pattern assignments from `master_gsea_table.csv` per the criteria in `01_Scripts/Python/pattern_definitions.py`.

---

**Last Updated**: 2026-05-15
**Generating scripts**: `02_Analysis/3.7.viz_bump_chart.py`, `02_Analysis/3.8.viz_interactive_bump_dashboard.py`, `02_Analysis/3.5.viz_trajectory_flow.py`, `02_Analysis/3.6.viz_alluvial_ggalluvial.R`
