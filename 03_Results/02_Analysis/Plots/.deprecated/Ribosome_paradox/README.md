# Ribosome Paradox — Core Finding Plots

This directory contains the core ribosome paradox visualization: during neuronal maturation (D35 → D65), DRP1 mutations drive cytoplasmic ribosome biogenesis upward while synaptic ribosome programs collapse. These figures underpin the central narrative and are referenced throughout the manuscript.

## File inventory

| File | Description |
|---|---|
| `Ribosome_Paradox_Three_Pools.pdf` | NES trajectory plot for three ribosome pools (cytoplasmic biogenesis, synaptic pre/post, mitochondrial) across Early / TrajDev / Late for G32A and R403C |
| `Ribosome_Temporal_Trajectory.pdf` | Line plot of ribosome biogenesis NES over the developmental trajectory, showing the V-shape compensation arc |
| `Ribosome_Paradox_Data.csv` | Per-pathway NES, FDR, and pool assignments for all three pools across all contrasts |

## Generating script

`02_Analysis/2.1.viz_ribosome_paradox.R`

```bash
Rscript 02_Analysis/2.1.viz_ribosome_paradox.R
```

Requires checkpoints `checkpoints/all_gsea_results.rds` and `checkpoints/syngo_gsea_results.rds`.

## The three ribosome pools

| Pool | Source | Trajectory pattern | TrajDev NES (G32A) | TrajDev FDR |
|---|---|---|---|---|
| 1. Cytoplasmic biogenesis | GO:BP RIBOSOME_BIOGENESIS | Compensation | +2.25 | 5.5 × 10⁻¹² |
| 2. Presynaptic ribosome | SynGO presyn_ribosome | Sign_reversal | −2.90 | 3.2 × 10⁻¹² |
| 2. Postsynaptic ribosome | SynGO postsyn_ribosome | Sign_reversal | −3.02 | 1.9 × 10⁻¹⁵ |
| 3. Mitochondrial ribosome | MitoCarta Mitochondrial_ribosome | Compensation | +1.89 | 7.9 × 10⁻⁴ |

The synaptic ribosome Sign_reversal is the amplitude extreme of a broader cytoplasmic structural-ribosome collapse — not a synapse-specific phenomenon. The SynGO 70-gene postsynaptic ribosome set is 98.6% contained in the cytoplasmic ribosomal proteome (see `../Publication_Figures/Fig1b_Ribosome_Gene_Overlap_UpSet.pdf` for the UpSet illustration).

## Reading guide

**Ribosome_Paradox_Three_Pools.pdf**: X-axis = trajectory stage (Early / TrajDev / Late); Y-axis = NES. Each line connects one pathway across the three stages. Color and shape encode pool membership. Solid points = FDR < 0.05; hollow = not significant. Separate facets for G32A and R403C. The divergence at TrajDev — Pools 1 and 3 going up while Pool 2 goes down — is the core visual representation of the paradox.

**Ribosome_Temporal_Trajectory.pdf**: Focused view of the cytoplasmic ribosome biogenesis NES, showing the V-shape (DOWN at Early, UP at TrajDev, near-normal at Late). Shaded region indicates significance. This panel emphasises that compensation is strongest during the D35–D65 critical period and subsides by Late stage.

**Ribosome_Paradox_Data.csv**: Columns include Pool (1/2/3), Pathway_Short, Contrast, Mutation, Timepoint (Early/TrajDev/Late), NES, p.adjust, pvalue, setSize.

---

**Last Updated**: 2026-05-14
**Generating script**: `02_Analysis/2.1.viz_ribosome_paradox.R`
