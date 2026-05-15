# Plots — Visualization Overview

This directory collects all figure outputs from the DRP1 mutation bulk RNA-seq analysis pipeline. Plots are organized by analysis type; each subdirectory has its own README. The overarching finding documented here is a **synaptic-ribosome translation paradox**: DRP1 mutations drive cytoplasmic ribosome biogenesis upward (Compensation, `ctrl_driven`) while simultaneously collapsing the assembled synaptic ribosome program (Sign_reversal, genuine sample-level crossover). All 70 genes of the SynGO synaptic-ribosome set are members of the curated cytoplasmic ribosomal gene set (Jaccard = 0.574), establishing this as the amplitude extreme of a broader structural-ribosome collapse rather than a synapse-specific phenomenon.

---

## Directory Inventory

| Subdirectory | Manuscript figure(s) | Generating script(s) |
|---|---|---|
| `Chord_Diagrams/` | Fig 5E | `02_Analysis/3.7.viz_chord_diagrams.py` |
| `Complex_V_analysis/` | Supplementary | `02_Analysis/2.5.viz_complex_v_analysis.R` |
| `Critical_period_trajectories/` | Supplementary | `02_Analysis/2.4.viz_critical_period_trajectories_gsva.R` |
| `General/` | QC; Supp 5b | `02_Analysis/1.1.main_pipeline.R`; `02_Analysis/revision/supplements/5b.*.R` |
| `GSEA/` | Fig 5D; multiple supp. | `02_Analysis/1.1.main_pipeline.R`, `02_Analysis/1.3.add_mitocarta.R`, `02_Analysis/3.9.viz_pooled_dotplots.R` |
| `Mito_translation_cascade/` | Fig 6B | `02_Analysis/2.2.viz_mito_translation_cascade.R` |
| `Pattern_Summary_Normalized/` | Fig 5A | `02_Analysis/3.4.pattern_summary_normalized.py`; `02_Analysis/revision/supplements/6a.geometric_scatter.py` |
| `Publication_Figures/` | Internal reference heatmap versions; Supp Fig S6 (UpSet) | `02_Analysis/3.1.publication_figures.py`, `02_Analysis/3.3.ribosome_upset_plot.py` |
| `Publication_Figures_Dotplot/` | Fig 6A | `02_Analysis/3.2.publication_figures_dotplot.py` |
| `Supplementary_6b/` | Supp Fig S7 | `02_Analysis/revision/supplements/6b.per_database_pattern_summary.py` |
| `Supplementary_8/` | Supp Fig S8 | `02_Analysis/revision/supplements/Supp8.focused_panel_classifications.py` |
| `Supplementary_9/` | Supp Fig S9 (cross-compartment) | `02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py` |
| `Supplementary_10/` | Supp Fig S10 | `02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py` |
| `Synaptic_ribosomes/` | Fig 5F | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| `SynGO_running_sum_normalised/` | Fig 5D companion (normalised running-sum panels) | `02_Analysis/3.10.viz_syngo_running_sum_normalised.R` |
| `Trajectory_Flow/` | Fig 5B, Fig 5C | `02_Analysis/3.7.viz_bump_chart.py`, `02_Analysis/3.8.viz_interactive_bump_dashboard.py`, `02_Analysis/3.5.viz_trajectory_flow.py`, `02_Analysis/3.6.viz_alluvial_ggalluvial.R` |
| `Volcano/` | Supp Fig S2E (standard set) | `02_Analysis/1.1.main_pipeline.R`, `01_Scripts/R_scripts/generate_vertical_volcanos.R` |
| `Volcano_Supplementary_MinAveExpr0/` | Supp Fig S2E, S2F | `01_Scripts/R_scripts/generate_vertical_volcanos.R`, `01_Scripts/R_scripts/generate_fdr_raster_volcanos.R` |

---

## Experimental design summary

- **Cell type**: iPSC-derived cortical neurons
- **Mutations**: DRP1 G32A (GTPase domain) and R403C (stalk domain)
- **Timepoints**: Day 35 (D35, early maturation) and Day 65 (D65, late maturation)
- **Samples**: 25 total (n = 3–6 per group)
- **Contrasts (9)**: G32A/R403C_vs_Ctrl at D35 and D65; Time_Ctrl/G32A/R403C; Maturation_G32A/R403C_specific (interaction term)
- **Databases (12)**: Hallmark, KEGG, Reactome, GO:BP, GO:CC, GO:MF, WikiPathways, Canonical, CGP, TF, SynGO, MitoCarta 3.0

## Trajectory pattern taxonomy

Eight mutually exclusive labels classify how a pathway's enrichment evolves across the Early → TrajDev → Late axis:

| Pattern | Definition |
|---|---|
| Compensation | Significant TrajDev opposing the Early defect |
| Sign_reversal | Sign flips between Early and Late; TrajDev crosses zero |
| Progressive | TrajDev amplifies the Early defect |
| Natural_improvement | Improvement without significant TrajDev |
| Natural_worsening | Worsening without significant TrajDev |
| Late_onset | No Early effect; new defect at Late |
| Transient | Early defect resolved by Late |
| Complex | Multiphasic or inconsistent |

Pattern labels describe the interaction-contrast view. The GSVA driver classification (`mutant_driven` / `ctrl_driven` / `both_moving` / `neither_moving`) is the replicate-level companion and can assign different biological meaning to the same contrast-level label (see `Supplementary_10/`).

## Key narrative anchor

The synaptic-ribosome Sign_reversal (NES TrajDev ≈ −2.9 to −3.0) is the amplitude extreme of a broader cytoplasmic structural-ribosome collapse, not a synapse-specific phenomenon. All 70 SynGO synaptic-ribosome genes are members of the curated cytoplasmic ribosomal gene set (Jaccard = 0.574). Cytoplasmic ribosome biogenesis runs opposite at the contrast level (Compensation), but replicate-level GSVA shows this label is `ctrl_driven` — controls descend toward mutant baseline rather than mutants rebounding.

---

**Last Updated**: 2026-05-15
