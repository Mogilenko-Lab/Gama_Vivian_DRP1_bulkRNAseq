# Supplementary Figure S9 — caption

**Cross-compartment ribosome trajectories: a structural-vs-biogenesis split
that maps unevenly onto compartments.** Each panel shows the mean NES across
all member pathways of five compartment-level categories at Early (35 DIV)
and Late (65 DIV) stages, for the G32A (A) and R403C (B) mutations
respectively. Shaded bands are ±1 SE around the mean; small semi-transparent
dots are the underlying per-pathway NES values. The five categories are
defined as follows (full member lists in
`cross_compartment_ribosome_trajectory_per_pathway.csv`):
*Mitochondrial Ribosome — structural* (MitoCarta `Mitochondrial_ribosome`,
GO:CC mitochondrial large / small ribosomal subunit, GO:CC organellar
ribosome; n=4); *Mitochondrial Ribosome — biogenesis* (MitoCarta
`Mitochondrial_ribosome_assembly`, GO:BP mitochondrial ribosome assembly;
n=2); *Cytoplasmic Ribosome — biogenesis* (GO:BP ribosome biogenesis /
assembly / large- and small-subunit biogenesis / subunit export from
nucleus; GO:CC pre-ribosome and pre-ribosome precursor compartments; n=9);
*Cytoplasmic Ribosome — structural* (GO:CC cytosolic ribosome and its large
/ small subunits, GO:CC large ribosomal subunit and ribosomal subunit,
GO:CC ribosome, GO:MF structural constituent of ribosome; n=7); and
*Synaptic Ribosome* (SynGO `SYNGO:presyn_ribosome`, `SYNGO:postsyn_ribosome`;
n=2). The three categories shown in cool / green tones all trace a
Compensation-like arc — early negative enrichment recovering toward zero or
positive late enrichment — together with the broader mitochondrial
compensatory programme. The two categories shown in warm tones trace the
opposite Sign-reversal arc, with the SynGO synaptic-ribosome category at the
largest amplitude. The split is therefore structural-vs-biogenesis rather
than mitochondrial-vs-cytoplasmic or synaptic-vs-bulk: the assembled
cytoplasmic ribosome (and its synaptic-localised subset) collapses, while
its biogenesis machinery, the assembled mitochondrial ribosome, and the
mitochondrial biogenesis programme all recover. Aggregated and per-pathway
numerical values, including the category × mutation × stage means used here,
are provided in `cross_compartment_ribosome_trajectory_aggregate.csv` and
`cross_compartment_ribosome_trajectory_per_pathway.csv`. Source script:
`02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py`.
Colors are from the Wong (2011) colorblind-safe palette.
