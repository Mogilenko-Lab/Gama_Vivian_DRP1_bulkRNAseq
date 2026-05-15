# Supplementary_8 — Focused Panel Classifications (Supp Fig S8)

This directory contains trajectory-pattern classification visualizations for the focused pathway panel (MitoCarta and SynGO databases), backing **Supplementary Fig S8** of the manuscript.

## File inventory

| File | Description |
|---|---|
| `focused_panel_classifications.pdf` | Supp Fig S8: pattern assignment summary for focused MitoCarta + SynGO pathways, both mutations, vector format |
| `focused_panel_classifications.png` | Raster version (300 DPI) |

## Generating script

`02_Analysis/revision/supplements/Supp8.focused_panel_classifications.py`

```bash
python3 02_Analysis/revision/supplements/Supp8.focused_panel_classifications.py
```

Input: `03_Results/02_Analysis/master_gsea_table.csv`

## Reading guide

The figure shows how the high-confidence focused pathway panel (MitoCarta mitochondrial pathways + SynGO synaptic pathways) distributes across the eight trajectory patterns for G32A and R403C mutations. Because this panel contains the pathways most central to the DRP1 phenotype narrative, the pattern distribution here is more interpretable than the cross-database summary in `../Supplementary_6b/`.

Key observations:
- **Mitochondrial pathways (MitoCarta)**: predominantly Compensation — the cell's adaptive response to DRP1-driven ATP shortage
- **Synaptic pathways (SynGO)**: include high-amplitude Sign_reversal entries for ribosome pathways (presynaptic and postsynaptic ribosome programs that collapse during maturation), alongside Compensation entries for synaptic-vesicle and scaffolding pathways

The contrast between the SynGO ribosome Sign_reversals and the MitoCarta Compensation pattern is the primary narrative tension captured in Fig 5B/C (bump charts) and the GSVA panels in `../Critical_period_trajectories/`.

## Manuscript supplementary caption (Supp Fig S8)

**Trajectory-pattern composition of the focused MitoCarta + SynGO pathway panel.** Two-panel summary of pattern assignments for the high-confidence focused pathway panel: MitoCarta 3.0 mitochondrial pathways (n = 149) and SynGO synaptic-localisation pathways (~300), drawn separately for G32A (left) and R403C (right). Each bar shows the count of pathways in each of the eight trajectory pattern classes (Compensation, Sign_reversal, Progressive, Natural_improvement, Natural_worsening, Late_onset, Transient, Complex) as classified by the criteria in `01_Scripts/Python/pattern_definitions.py` against the master GSEA table. Statistical thresholds are inherited from the upstream fgsea pipeline (10,000 permutations, BH-FDR; Active_* super-categories require BH-FDR < 0.05 for the TrajDev contrast). Because this panel contains the pathways most central to the DRP1 phenotype narrative, the pattern distribution here is more interpretable than the cross-database summary in `../Supplementary_6b/`. Two observations support the manuscript narrative: (i) Compensation exceeds passive-recovery patterns in both mutations and Progressive is rare-to-absent in both — consistent with an active adaptive response rather than progressive maladaptation; (ii) SynGO ribosomes (presynaptic and postsynaptic ribosome) contribute the Sign_reversal entries that distinguish the synaptic-translation collapse from the broader MitoCarta Compensation arm. R403C consistently shows more Compensation than G32A, with Compensation constituting a strict majority (> 50%) of classifiable G32A pathways.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/revision/supplements/Supp8.focused_panel_classifications.py`
