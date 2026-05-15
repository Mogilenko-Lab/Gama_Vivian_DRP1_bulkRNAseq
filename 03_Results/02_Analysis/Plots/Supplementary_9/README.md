# Supplementary_9 — Cross-Compartment Ribosome Trajectories (Supp Fig S9)

Cross-compartment ribosome trajectory analysis backing **Supplementary Fig S9**. The figure reveals that the translation paradox is a **structural-vs-biogenesis split** rather than a mitochondrial-vs-cytoplasmic or synaptic-vs-bulk distinction.

## File inventory

| File | Description |
|---|---|
| `cross_compartment_ribosome_trajectory.pdf` | Supp Fig S9: mean NES trajectories for five ribosome compartment-level categories at Early (D35) and Late (D65), G32A and R403C panels |
| `cross_compartment_ribosome_trajectory.png` | Raster version (300 dpi) |
| `cross_compartment_ribosome_trajectory_aggregate.csv` | Aggregated mean NES values per category × mutation × stage (used to build the trajectory lines) |
| `cross_compartment_ribosome_trajectory_per_pathway.csv` | Per-pathway NES values, category membership, and source database (full membership audit) |

## Generating script

`02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py`

```bash
python3 02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py
```

Input: `03_Results/02_Analysis/master_gsea_table.csv`.

## The five categories

| Category | Source pathways | n |
|---|---|---|
| Mitochondrial Ribosome — structural | MitoCarta `Mitochondrial_ribosome`; GO:CC mitochondrial large / small ribosomal subunit; GO:CC organellar ribosome | 4 |
| Mitochondrial Ribosome — biogenesis | MitoCarta `Mitochondrial_ribosome_assembly`; GO:BP mitochondrial ribosome assembly | 2 |
| Cytoplasmic Ribosome — biogenesis | GO:BP ribosome biogenesis / assembly / large- and small-subunit biogenesis / subunit export from nucleus; GO:CC pre-ribosome compartments | 9 |
| Cytoplasmic Ribosome — structural | GO:CC cytosolic ribosome and its large / small subunits, GO:CC ribosome, GO:MF structural constituent of ribosome | 7 |
| Synaptic Ribosome | SynGO `SYNGO:presyn_ribosome`, `SYNGO:postsyn_ribosome` | 2 |

## How to read this folder

`cross_compartment_ribosome_trajectory.pdf` is the panel that goes into Supp Fig S9 — open this first. The two CSVs back the figure: `_aggregate.csv` gives the category-level means (the trajectory lines) and `_per_pathway.csv` gives the underlying per-pathway NES values (the small dots in the figure) plus the category membership audit a reviewer would want to inspect. To regenerate, rerun the Python script above; no checkpoint dependencies.

## Manuscript figure caption (Supp Fig S9)

**Cross-compartment ribosome trajectories: a structural-vs-biogenesis split that maps unevenly onto compartments.** Each panel plots the mean GSEA NES across all member pathways of five compartment-level categories at Early (D35; 35 DIV) and Late (D65; 65 DIV) stages for the G32A (A) and R403C (B) mutations, respectively. Shaded bands are ±1 SE around the category mean; small semi-transparent dots are the underlying per-pathway NES values (membership listed in `cross_compartment_ribosome_trajectory_per_pathway.csv`); the white trajectory lines connect the Early and Late category means. The five categories are: *Mitochondrial Ribosome — structural* (MitoCarta `Mitochondrial_ribosome`, GO:CC mitochondrial large/small ribosomal subunit, GO:CC organellar ribosome; n = 4); *Mitochondrial Ribosome — biogenesis* (MitoCarta `Mitochondrial_ribosome_assembly`, GO:BP mitochondrial ribosome assembly; n = 2); *Cytoplasmic Ribosome — biogenesis* (GO:BP ribosome biogenesis / assembly / large- and small-subunit biogenesis / subunit export from nucleus; GO:CC pre-ribosome and pre-ribosome precursor compartments; n = 9); *Cytoplasmic Ribosome — structural* (GO:CC cytosolic ribosome and its large/small subunits, GO:CC large ribosomal subunit and ribosomal subunit, GO:CC ribosome, GO:MF structural constituent of ribosome; n = 7); and *Synaptic Ribosome* (SynGO `SYNGO:presyn_ribosome`, `SYNGO:postsyn_ribosome`; n = 2). The three categories shown in cool/green tones (mitochondrial ribosome — structural, mitochondrial ribosome — biogenesis, cytoplasmic ribosome — biogenesis) all trace a Compensation-like arc — Early negative enrichment recovering toward zero or positive enrichment at Late — together with the broader mitochondrial compensatory programme. The two categories shown in warm tones (cytoplasmic ribosome — structural, synaptic ribosome) trace the opposite Sign_reversal arc, with the SynGO synaptic ribosome category at the largest amplitude. The split is therefore structural-vs-biogenesis rather than mitochondrial-vs-cytoplasmic or synaptic-vs-bulk: the assembled cytoplasmic ribosome (and its synaptic-localised subset) collapses, while its biogenesis machinery, the assembled mitochondrial ribosome, and the mitochondrial biogenesis programme all recover. Statistical thresholds inherited from the upstream fgsea pipeline (10,000 permutations, BH-FDR; per-pathway and aggregate values in the companion CSVs). Colours follow the Wong (2011) colorblind-safe palette.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py`
