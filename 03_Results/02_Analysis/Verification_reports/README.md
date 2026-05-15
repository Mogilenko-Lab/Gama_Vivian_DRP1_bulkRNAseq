# 03_Results/02_Analysis/Verification_reports — QC reports, ribosome and calcium gene lists

## Overview

This directory holds the quality-control and verification outputs that support the manuscript's core ribosome and calcium findings: (a) the three-compartment SynGO gene lists that confirm the synaptic ribosome scope, (b) GSEA statistics for ribosome and calcium pathways to demonstrate cross-database consistency, and (c) core-enrichment gene lists from the TrajDev GSEA leading edge. The `verification_summary_report.md` provides a concise human-readable summary of all validation checks.

## File Inventory

| File | Description | Generating Script |
|------|-------------|-------------------|
| `verification_summary_report.md` | Narrative QC summary: gene detection rates, calcium gene coverage, ribosome pathway direction, cross-database consistency | Manual/viz scripts |
| `calcium_genes_DE_summary.csv` | DE statistics (logFC, adj.P.Val, Significant flag) for the 13 calcium genes × 9 contrasts | `02_Analysis/2.6.viz_calcium_genes.R` |
| `calcium_pathways_all.csv` | All GSEA results for calcium-annotated pathways across all databases and contrasts | `02_Analysis/2.6.viz_calcium_genes.R` |
| `calcium_pathways_significant.csv` | Subset of `calcium_pathways_all.csv` filtered to FDR < 0.05 | `02_Analysis/2.6.viz_calcium_genes.R` |
| `calcium_pathways_summary_by_contrast.csv` | Count of significant calcium pathways per contrast | `02_Analysis/2.6.viz_calcium_genes.R` |
| `ribosome_pathway_statistics.csv` | GSEA statistics for key ribosome pathways across databases and contrasts | `02_Analysis/2.1.viz_ribosome_paradox.R` |
| `ribosome_pathway_significant.csv` | Subset filtered to FDR < 0.05 ribosome pathways | `02_Analysis/2.1.viz_ribosome_paradox.R` |
| `syngo_presyn_ribosome_genes.txt` | 52 presynaptic ribosome genes from SynGO CC (one symbol per line) | `02_Analysis/1.8.extract_syngo_ribosome_genes.R` |
| `syngo_postsyn_ribosome_genes.txt` | 70 postsynaptic ribosome genes from SynGO CC (one symbol per line) | `02_Analysis/1.8.extract_syngo_ribosome_genes.R` |
| `syngo_all_synaptic_ribosome_genes.txt` | Union of pre- and postsynaptic SynGO ribosome genes | `02_Analysis/1.8.extract_syngo_ribosome_genes.R` |
| `syngo_ribosome_gene_membership.csv` | Gene-level compartment flags: `Gene`, `Presynaptic`, `Postsynaptic`, `Both` | `02_Analysis/1.8.extract_syngo_ribosome_genes.R` |
| `Maturation_G32A_specific_presyn_core_genes.txt` | GSEA leading-edge genes for G32A TrajDev × presynaptic ribosome pathway | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| `Maturation_G32A_specific_postsyn_core_genes.txt` | GSEA leading-edge genes for G32A TrajDev × postsynaptic ribosome pathway | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| `Maturation_R403C_specific_presyn_core_genes.txt` | GSEA leading-edge genes for R403C TrajDev × presynaptic ribosome pathway | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |
| `Maturation_R403C_specific_postsyn_core_genes.txt` | GSEA leading-edge genes for R403C TrajDev × postsynaptic ribosome pathway | `02_Analysis/2.3.viz_synaptic_ribosomes.R` |

---

## Reading Guide

### calcium_pathways_all.csv / ribosome_pathway_statistics.csv

Both share the same GSEA column schema: `ID`, `Description`, `setSize`, `enrichmentScore`, `NES`, `pvalue`, `p.adjust`, `qvalue`, `rank`, `leading_edge`, `core_enrichment`, followed by either a `SearchTerm` (calcium files) or a `Contrast` column. `NES` > 0 = upregulated gene set in numerator condition; `p.adjust` < 0.05 = significant after BH correction.

### syngo_ribosome_gene_membership.csv

Columns: `Gene` (HGNC symbol), `Presynaptic` (TRUE/FALSE), `Postsynaptic` (TRUE/FALSE), `Both` (TRUE/FALSE). Gene counts: 52 presynaptic, 70 postsynaptic, 52 shared (all presynaptic genes are also postsynaptic), 18 postsynaptic-only. This asymmetry reflects that postsynaptic densities accumulate a superset of the presynaptic ribosome complement and is key context for interpreting Supplementary Table S3 (Jaccard 0.561 synaptic↔cytoplasmic).

### Core enrichment gene lists (`*_core_genes.txt`)

One HGNC symbol per line. These are the GSEA leading-edge genes — the subset that drives the enrichment score for the TrajDev contrast in the SynGO synaptic ribosome pathway. They represent genes most downregulated during the G32A/R403C maturation deviation relative to controls.
