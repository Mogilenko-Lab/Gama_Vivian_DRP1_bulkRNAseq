# GSEA — Gene Set Enrichment Analysis Results

This directory holds all GSEA visualizations organized by contrast and pathway database. The SynGO results back **Fig 5D** (running sum plots) of the manuscript. Results span 9 experimental contrasts and 12 databases (~1,000 PDF figures total).

## Directory structure

```
GSEA/
├── Cross_database_pooled/        # Combined visualizations across databases
│   ├── comprehensive/            # All 12 databases pooled
│   ├── focused/                  # Hallmark + KEGG + Reactome only
│   ├── csv_data/                 # Filtered and raw CSV exports
│   └── test_output/              # Verification plots
│
├── G32A_vs_Ctrl_D35/            # Per-contrast folder (× 9)
│   ├── hallmark/ kegg/ reactome/
│   ├── gobp/ gocc/ gomf/
│   ├── wiki/ canon/ cgp/ tf/
│   ├── SynGO/                   # SynGO cellular-component pathways
│   └── MitoCarta/               # MitoCarta 3.0 mitochondrial pathways
│
├── G32A_vs_Ctrl_D65/
├── R403C_vs_Ctrl_D35/
├── R403C_vs_Ctrl_D65/
├── Time_Ctrl/
├── Time_G32A/
├── Time_R403C/
├── Maturation_G32A_specific/
└── Maturation_R403C_specific/
```

## Databases

| Database | Folder | Approximate gene sets | Notes |
|---|---|---|---|
| Hallmark | hallmark/ | 50 | MSigDB curated hallmark states |
| KEGG | kegg/ | 186 | Metabolic and signaling pathways |
| Reactome | reactome/ | 1,615 | Curated reaction pathways |
| GO:BP | gobp/ | 7,658 | Biological Process |
| GO:CC | gocc/ | 1,006 | Cellular Component |
| GO:MF | gomf/ | 1,738 | Molecular Function |
| WikiPathways | wiki/ | 664 | Community-curated |
| Canonical | canon/ | 2,922 | MSigDB C2 canonical |
| CGP | cgp/ | 3,358 | Chemical & Genetic Perturbations |
| TF | tf/ | 1,137 | Transcription Factor targets |
| SynGO | SynGO/ | ~300 | Synaptic Gene Ontology (CC namespace) |
| MitoCarta | MitoCarta/ | 149 | Mitochondrial pathways (MitoCarta 3.0) |

## File types per database folder

| File pattern | Description |
|---|---|
| `[contrast]_[db]_running_sum.pdf` | GSEA running sum enrichment plot (Fig 5D for SynGO) |
| `[contrast]_[db]_combined.pdf` | Multi-panel: dotplot + bar + running sum |
| `[contrast]_[db]_facet.pdf` | Dotplot faceted by up/down regulation |
| `[contrast]_[db]_nes_bar.pdf` | NES bar chart (top pathways) |
| `[contrast]_[db]_up_dot.pdf` | Upregulated pathways dotplot (MSigDB) |
| `[contrast]_[db]_down_dot.pdf` | Downregulated pathways dotplot (MSigDB) |
| `[contrast]_[db]_dot.pdf` | Combined dotplot (SynGO and MitoCarta) |
| `[contrast]_[db]_bar.pdf` | Bar plot (SynGO and MitoCarta) |
| `[contrast]_[db]_results.txt` | Text summary of enrichment statistics |
| `GSEA_[db]_result.rds` | R object with full fgsea output (SynGO, MitoCarta) |

## Generating scripts

- `02_Analysis/1.1.main_pipeline.R` — MSigDB + SynGO GSEA
- `02_Analysis/1.3.add_mitocarta.R` — MitoCarta GSEA (separate pass)
- `02_Analysis/3.9.viz_pooled_dotplots.R` — Cross_database_pooled/ visualizations

## Reading guide

**Dotplots** (`*_up_dot.pdf`, `*_down_dot.pdf`, `*_dot.pdf`): Y-axis = pathway name; X-axis = gene ratio (fraction of set in leading edge); dot color = NES (orange/red = upregulated, blue = downregulated); dot size = −log10(FDR). Larger, more saturated dots indicate stronger and more significant enrichments.

**Running sum plots** (`*_running_sum.pdf`, Fig 5D): X-axis = rank position in the gene list (ordered by t-statistic); Y-axis = running enrichment score. The peak is the enrichment score (ES). Vertical tick marks show positions of gene-set members. A sharp, early peak for a downregulated pathway indicates leading-edge genes cluster at the top of the ranked list (most downregulated in the contrast). For Fig 5D (SynGO synaptic-ribosome pathways in G32A/R403C vs Ctrl at D35 and D65), the running sum illustrates the amplitude and leading-edge composition of the synaptic-ribosome collapse.

**NES interpretation**: |NES| > 2.0 = strong effect; 1.5–2.0 = moderate; < 1.5 = mild. Significance threshold FDR < 0.05 throughout.

**Contrast interpretation**: A positive NES means the pathway is enriched toward the first condition listed (e.g., G32A for G32A_vs_Ctrl). For interaction contrasts (Maturation_G32A_specific), a positive NES means G32A maturation exceeds control maturation for that pathway.

## Manuscript figure mapping

| File path | Manuscript role |
|---|---|
| `*/SynGO/*_SynGO_running_sum.pdf` | Per-contrast SynGO running-sum panels — back **Fig 5D** (composite version in `../SynGO_running_sum_normalised/`) |
| `Cross_database_pooled/focused/` | Pooled cross-database dotplots used in supplementary text |
| All other contrast × database folders | Extended-data per-contrast enrichment outputs cited in Supp Tables S2 and the methods section |

---

**Last Updated**: 2026-05-15
**Generating scripts**: `02_Analysis/1.1.main_pipeline.R`, `02_Analysis/1.3.add_mitocarta.R`, `02_Analysis/3.9.viz_pooled_dotplots.R`
