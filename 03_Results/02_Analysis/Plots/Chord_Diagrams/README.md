# Chord Diagrams — Synaptic Ribosome Gene-Pathway Membership

This directory contains GOChord-style chord diagrams linking RPL/RPS ribosomal protein
genes to their GSEA pathway memberships for the two DRP1 mutations (G32A, R403C) in
cortical neurons at D35 and D65. Three diagram variants are generated per mutation, each
expanding the pathway scope to contextualize the synaptic ribosome signal within the
broader cytoplasmic translation program.

---

## Diagram variants

### 1. focused

**Files:** `chord_diagram_G32A.{pdf,png}`, `chord_diagram_R403C.{pdf,png}`

Shows the six SynGO synaptic-localized pathways only:

| Pathway ID | Display name |
|---|---|
| SYNGO:presyn_ribosome | Presyn. Ribosome |
| SYNGO:postsyn_ribosome | Postsyn. Ribosome |
| GO:0099523 | Presyn. Cytosol |
| GO:0099524 | Postsyn. Cytosol |
| GO:0014069 | Postsyn. Density |
| GO:0045211 | Postsyn. Membrane |

No outer grouping brackets are drawn (single group only).

---

### 2. expanded_clean

**Files:** `chord_diagram_G32A_expanded_clean.{pdf,png}`, `chord_diagram_R403C_expanded_clean.{pdf,png}`

Adds **GO:CC Cytosolic Ribosome** to the six SynGO pathways (seven pathways total):

| Pathway ID | Display name | Group |
|---|---|---|
| SYNGO:presyn_ribosome | Presyn. Ribosome | SynGO synaptic-localized |
| SYNGO:postsyn_ribosome | Postsyn. Ribosome | SynGO synaptic-localized |
| GO:0099523 | Presyn. Cytosol | SynGO synaptic-localized |
| GO:0099524 | Postsyn. Cytosol | SynGO synaptic-localized |
| GO:0014069 | Postsyn. Density | SynGO synaptic-localized |
| GO:0045211 | Postsyn. Membrane | SynGO synaptic-localized |
| GOCC_CYTOSOLIC_RIBOSOME | GO:CC Cytosolic Ribosome | Broader cytoplasmic ribosome / translation |

The warm-vs-cool palette conveys group nesting at a glance; no outer-arc brackets are
drawn (the palette encodes nesting sufficiently). The two-layer ribbon structure (see
**Ribbon semantics** below) makes the 100% gene-set containment of the SynGO synaptic
ribosome sets within GO:CC cytosolic ribosome directly visible.

**Recommended for main figure panel.**

---

### 3. expanded_full

**Files:** `chord_diagram_G32A_expanded_full.{pdf,png}`, `chord_diagram_R403C_expanded_full.{pdf,png}`

Adds cytosolic ribosome plus three further pathway sets (ten pathways total):

| Pathway ID | Display name | Group |
|---|---|---|
| SYNGO:presyn_ribosome | Presyn. Ribosome | SynGO synaptic-localized |
| SYNGO:postsyn_ribosome | Postsyn. Ribosome | SynGO synaptic-localized |
| GO:0099523 | Presyn. Cytosol | SynGO synaptic-localized |
| GO:0099524 | Postsyn. Cytosol | SynGO synaptic-localized |
| GO:0014069 | Postsyn. Density | SynGO synaptic-localized |
| GO:0045211 | Postsyn. Membrane | SynGO synaptic-localized |
| GOCC_CYTOSOLIC_RIBOSOME | GO:CC Cytosolic Ribosome | Broader cytoplasmic ribosome / translation |
| KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION | KEGG Translation Initiation | Broader cytoplasmic ribosome / translation |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | REACTOME Translation Elong. | Broader cytoplasmic ribosome / translation |
| GOBP_RIBOSOME_BIOGENESIS | GO:BP Ribosome Biogenesis | Biogenesis (counter-direction) |

Three outer-ring brackets are drawn: "SynGO", "Broader cyto.", and "Biogenesis". The
biogenesis arc runs in the opposite NES direction from the structural translation
pathways, representing the compensation arm.

**Recommended for supplemental figure.**

---

## Layout and visual encoding

### Symmetric layout

The pathway arc is centered on the right (angular midpoint = 0°) and the gene arc is
centered on the left (angular midpoint = 180°), with equal top and bottom gaps. This
eliminates label collisions and produces a balanced, publication-ready layout.

### Gene ring order (all variants)

Each gene occupies a stacked pair of concentric colored rectangles on the left arc.
Reading outward from the center of the diagram:

- **Inner ring** = D65 log2FC
- **Outer ring** = D35 log2FC

Concentric reading outer to inner therefore corresponds to temporal progression from
early (D35) to late (D65). This makes the biphasic sign inversion — strong positive
enrichment at D35 flipping to strong negative at D65 — visible as a color reversal
between the two rings. The color scale is a blue-white-orange diverging map
(range: −3 to +3 log2FC); gray indicates a missing or near-zero value.

### Pathway arcs (right side)

Each pathway occupies a colored arc on the right semicircle. Pathway labels appear
outside the arc, rotated radially. NES values are not printed on the figure; they are
available in `chord_pathways_*.csv` and in the manuscript text.

### Ribbon semantics

**This is the most important encoding in the figure.** Each gene-pathway pair can be
represented by one of two ribbon types, drawn in two layers:

| Ribbon type | Alpha | Z-order | Meaning |
|---|---|---|---|
| Bold (per-tier alpha, foreground) | 0.70 (focal) / 0.30 (recessive) / 0.25 (envelope/biogenesis) | Top | Gene is in the GSEA **leading edge** for that pathway in this contrast (D35 contrast, falling back to D65 if D35 had no enrichment) |
| Faint (background) | 0.13 | Bottom | Gene is **annotated** to the pathway (gene-set membership in MSigDB / SynGO) but was NOT captured in that pathway's leading edge for this contrast |

The faint background layer is essential for biological honesty. Without it, a gene
appearing with a cool ribbon but no warm GO:CC ribbon could be misread as evidence that
the gene is absent from GO:CC cytosolic ribosome. In reality, SynGO synaptic-ribosome
gene sets are 100% contained within GO:CC cytosolic ribosome (70/70 genes). The absence
of a bold warm ribbon for some genes reflects GSEA leading-edge statistics in that
contrast — not gene-set biology.

**Concrete example:** In R403C, RPL14 is in the SynGO synaptic-ribosome leading edge
but only in the gene-set annotation (not the leading edge) of GO:CC cytosolic ribosome.
Its connections therefore appear as bold cool ribbons and one faint warm ribbon. This is
a GSEA statistical artifact, not biology: RPL14 is part of the cytoplasmic ribosomal
proteome. The same pattern holds for RPLP1 and RPS5 in G32A.

The legend includes a dedicated "Ribbons" section with two swatches illustrating the
distinction (saturated = leading edge; faint = annotation only).

### Three-tier color palette

A colorblind-safe palette encodes the manuscript's visual hierarchy. Tiers are drawn
in z-order from bottom to top:

| Tier | Pathways | Hex codes | Ribbon α | Visual role |
|---|---|---|---|---|
| **Focal** (drawn on top) | SYNGO:presyn_ribosome, SYNGO:postsyn_ribosome | `#882255` (wine), `#332288` (indigo) | 0.70 | Primary message — synaptic ribosome pathways |
| **Recessive** | Presyn. Cytosol, Postsyn. Cytosol, Postsyn. Density, Postsyn. Membrane | `#117733`, `#44AA99`, `#999933`, `#88CCEE` | 0.30 | Supporting synaptic compartment pathways |
| **Envelope** (drawn underneath) | GO:CC Cytosolic Ribosome, KEGG Translation Initiation, REACTOME Translation Elong. | `#DDCC77`, `#E69F00`, `#D55E00` | 0.25 | Broader cytoplasmic context; warm tones contrast with cool SynGO tones |
| **Counter-direction** | GO:BP Ribosome Biogenesis | `#6699CC` (desaturated steel) | 0.25 | Runs opposite to structural-ribosome cluster; biogenesis compensation arm |

### Group brackets

Outer-arc brackets with inward tick marks and short text labels are drawn only for
`expanded_full` (three groups: "SynGO", "Broader cyto.", "Biogenesis"). They are omitted
from `expanded_clean` (palette alone encodes nesting) and from `focused` (single group).

### Legend elements

- **Top-left box**: Gene ring explanation (outer = D35, inner = D65; temporal reading direction).
- **Bottom-right box**: Pathway legend grouped by category, with colored swatches.
- **Left-side colorbar**: log2 fold change scale for the gene ring color map.
- **Ribbons section**: Two swatches — saturated (leading edge) and faint (annotation only).

---

## Output file naming convention

```
chord_diagram_{MUT}.{pdf,png}                         # focused variant
chord_diagram_{MUT}_expanded_clean.{pdf,png}          # expanded_clean variant
chord_diagram_{MUT}_expanded_full.{pdf,png}           # expanded_full variant

chord_genes_{MUT}.csv                                 # focused gene table
chord_genes_{MUT}_expanded_clean.csv
chord_genes_{MUT}_expanded_full.csv

chord_pathways_{MUT}.csv                              # focused pathway table (includes NES)
chord_pathways_{MUT}_expanded_clean.csv
chord_pathways_{MUT}_expanded_full.csv

chord_connections_{MUT}.csv                           # focused leading-edge gene-pathway links
chord_connections_{MUT}_expanded_clean.csv
chord_connections_{MUT}_expanded_full.csv

chord_connections_geneset_only_{MUT}.csv              # gene-set member but NOT in leading edge
chord_connections_geneset_only_{MUT}_expanded_clean.csv
chord_connections_geneset_only_{MUT}_expanded_full.csv

gsea_pathway_data.csv                                 # all-contrast GSEA statistics
                                                      # (shared input, single file)
```

Where `{MUT}` is `G32A` or `R403C`. PDF files are saved at 600 dpi for print
publication; PNG files at 300 dpi.

---

## Data CSV contents

| File pattern | Key columns |
|---|---|
| `chord_genes_*.csv` | gene name, logfc_d35, logfc_d65, padj_d35, padj_d65, n_pathways, is_ribosomal, mutation |
| `chord_pathways_*.csv` | pathway id, display_name, group, nes_d35, nes_d65, padj_d35, padj_d65, core_genes, n_genes, mutation |
| `chord_connections_*.csv` | gene, pathway, mutation (leading-edge membership only) |
| `chord_connections_geneset_only_*.csv` | gene, pathway, mutation (gene-set annotation, NOT leading edge) |
| `gsea_pathway_data.csv` | ID, NES, pvalue, p.adjust, core_enrichment, contrast |
| `pathway_gene_set_memberships.csv` | pathway × gene long-format membership table used by the gene-set-only ribbon layer |

---

## Generating script

**Script:** `02_Analysis/3.7.viz_chord_diagrams.py`

**Command (run from project root):**

```bash
python3 02_Analysis/3.7.viz_chord_diagrams.py
```

The script calls R internally via `subprocess` to extract NES and leading-edge
memberships from the checkpoint RDS files:

- `03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds` — SynGO pathway results
- `03_Results/02_Analysis/checkpoints/all_gsea_results.rds` — MSigDB / KEGG / REACTOME

No command-line arguments are required. Both mutations and all three variants are
generated in a single run. Total output: 12 PDF/PNG pairs + 15 CSV files.

The CONFIGURATION block at the top of the script (`PATHWAY_DEFS`, `VARIANTS`,
`PATHWAY_COLOR_MAP`, `GROUP_COLORS`) is the single point of control for pathway
membership, variant scope, and color assignments.

---

## Scientific reading

The focused variant establishes that SynGO-annotated synaptic-localized ribosomal
pathways carry the largest-amplitude enrichment signal in the dataset, with a biphasic
trajectory: strong positive NES at D35 reversing to strong negative NES at D65.

The expanded variants place this signal in its genomic context. The two SynGO
synaptic-ribosome pathways are 100% contained within GO:CC cytosolic ribosome at the
gene-set level (70/70 genes). The same D35-to-D65 sign inversion propagates across the
full cytoplasmic structural-ribosome cluster: GO:CC cytosolic ribosome, KEGG translation
initiation, and REACTOME eukaryotic translation elongation all mirror the synaptic
pattern. The synaptic fraction sits at the amplitude extreme of this broader program
rather than as a distinct phenomenon.

The two-layer ribbon structure makes this containment visible: bold cool ribbons (leading
edge) connect most genes to the SynGO arcs, while faint warm ribbons (annotation only)
document that the same genes are GO:CC members regardless of whether the GSEA leading
edge captured them in that contrast. The `chord_connections_geneset_only_*.csv` files
enumerate these discordant pairs for manuscript enumeration.

The biogenesis group (GO:BP ribosome biogenesis, `expanded_full` only) runs in the
opposite direction — positive NES at D65 — representing a transcriptional upregulation
of the ribosome-assembly program as a compensatory response. This counter-direction arc
is spatially separated from the structural-translation cluster and labeled "Biogenesis"
in the outer bracket.

---

**Last updated:** 2026-05-15
**Generating script:** `02_Analysis/3.7.viz_chord_diagrams.py`

**Manuscript figure mapping:** the `expanded_clean` variant backs **Fig 5E** of the main text; the `expanded_full` variant is the matched supplementary panel.

---

## Suggested manuscript captions

### Caption A — panel uses `expanded_clean` (recommended for main figure)

**(E) Chord diagrams linking ribosomal protein genes to synaptic and cytoplasmic
ribosomal pathways for G32A (left) and R403C (right) mutations.** Gene arcs (left
semicircle) show each RPL/RPS ribosomal protein as two concentric colored rectangles:
the outer rectangle encodes the D35 log2 fold change and the inner encodes the D65
log2 fold change (blue-white-orange diverging scale, range ±3 log2FC), so reading
outward to inward follows temporal progression from early (D35) to late (D65). Pathway
arcs (right semicircle) represent six SynGO synaptic-localized pathways (presynaptic
ribosome, postsynaptic ribosome, and four synaptic-compartment sets) rendered in cool
jewel tones, and GO:CC cytosolic ribosome rendered in warm amber — a colorblind-safe
palette in which the focal presynaptic and postsynaptic ribosome arcs (wine `#882255`
and indigo `#332288`) are drawn on top at full saturation. Bezier ribbons encode
gene-pathway relationships in two layers: bold ribbons (higher opacity) indicate that
the gene was in the GSEA leading edge for that pathway in the D35 contrast; faint
ribbons (α = 0.13, drawn underneath) indicate that the gene is annotated to the pathway
by gene-set membership but was not captured in the leading edge for this contrast. The
100% gene-set containment of SynGO synaptic-ribosome sets within GO:CC cytosolic
ribosome (70/70 genes) is visible as a dense overlap between the cool SynGO ribbon
bundles and faint warm GO:CC ribbons: where a cool ribbon has no accompanying bold warm
ribbon, the faint warm ribbon restores biological truth — the gene is part of the
cytoplasmic ribosomal proteome regardless of leading-edge statistics. The same biphasic
enrichment pattern — strong positive NES at D35 reversing to strong negative NES at D65
— runs uniformly across the structural-ribosome cluster in both mutations.

---

### Caption B — panel uses `expanded_full` (recommended for supplemental figure)

**(E) Chord diagrams linking ribosomal protein genes to synaptic, cytoplasmic-structural,
and biogenesis pathways for G32A (left) and R403C (right) mutations.** Gene arcs (left
semicircle) show each RPL/RPS ribosomal protein as two concentric colored rectangles;
the outer rectangle encodes D35 log2 fold change and the inner encodes D65 log2 fold
change (blue-white-orange diverging scale, ±3 log2FC), so reading outward to inward
corresponds to temporal progression from early (D35) to late (D65). Pathway arcs (right
semicircle) span three groups demarcated by outer-arc brackets: "SynGO" (six
synaptic-localized pathways, cool tones — focal presynaptic and postsynaptic ribosome
arcs in wine and indigo drawn on top), "Broader cyto." (GO:CC cytosolic ribosome, KEGG
translation initiation, REACTOME eukaryotic translation elongation — warm earth tones,
drawn underneath), and "Biogenesis" (GO:BP ribosome biogenesis — desaturated steel,
counter-direction). Bezier ribbons are drawn in two layers: bold ribbons (per-tier
opacity) connect a gene to pathways in whose GSEA leading edge it was captured for the
D35 contrast; faint ribbons (α = 0.13) connect a gene to pathways to which it is
annotated by gene-set membership but whose leading edge did not include it in this
contrast. This two-layer encoding exposes the 100% gene-set containment of SynGO
synaptic-ribosome sets within GO:CC cytosolic ribosome (70/70 genes): genes present with
bold cool ribbons but only faint warm ribbons are leading-edge members of SynGO pathways
whose GO:CC membership is a biological fact independent of GSEA statistics. The
structural-translation cluster shares a nearly identical ribbon bundle across all three
groups, and the same D35-positive to D65-negative sign inversion runs uniformly across
both mutations. The biogenesis arc runs in the opposite direction — elevated at D65 —
representing transcriptional upregulation of ribosome assembly machinery as a
compensatory response to structural-ribosome depletion.
