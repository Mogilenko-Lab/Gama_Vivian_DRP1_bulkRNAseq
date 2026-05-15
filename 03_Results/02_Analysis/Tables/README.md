# 03_Results/02_Analysis/Tables — aggregate pathway-summary and ribosome-compartment tables

## Overview

Compact CSV summaries that aggregate pathway-level GSEA results into paper-citable tables. This directory backs **Supplementary Table S3** (`ribosome_compartment_summary.csv`, three-compartment Jaccard comparison) and **Supplementary Table S4** (`per_database_pattern_summary.csv` and the preferred dual-denominator companion `per_database_pattern_dual_denominator.csv`, per-collection trajectory-pattern frequencies for all 12 databases). A two-mutation shared-DEG list is also stored here. All tables are derived artifacts: regenerate from `master_gsea_table.csv` and the DE results in `../DE_results/` rather than editing CSVs by hand.

## File inventory

| File | Rows | Purpose | Manuscript anchor | Generator |
|---|---|---|---|---|
| `per_database_pattern_dual_denominator.csv` | 208 | Per-`(database, mutation, pattern)` counts with all three denominators side-by-side; preferred citation source | Supp Table S4 | `02_Analysis/revision/supplements/6a.per_database_pattern_dual_denominator.py` |
| `per_database_pattern_summary.csv` | 24 | Legacy per-`(database, mutation)` pattern counts; single denominator (5,267 universe) | Supp Table S4 (legacy view) | `02_Analysis/revision/supplements/6b.per_database_pattern_summary.py` |
| `ribosome_compartment_summary.csv` | 95 | Synaptic / mitochondrial / cytoplasmic ribosome pathways with Jaccard, per-stage NES, and pattern labels | Supp Table S3 | `02_Analysis/revision/supplements/6c.ribosome_compartment_summary.py` |
| `shared_DEGs_G32A_R403C.csv` | 38 | Genes DE in both G32A and R403C at Day 35 | Reviewer response 5b | `02_Analysis/revision/supplements/5b.extract_shared_maturation_degs.R` |

---

## Denominator vocabulary

The trajectory-pattern files use the manuscript-aligned **5,267 ever-significant pathway universe** by default (FDR < 0.05 in any of 9 GSEA contrasts, drawn from the 12,221-pathway full universe). The dual-denominator CSV carries all three denominators in explicitly named columns so any percentage in the paper can be audited against the right base.

| Denominator key | Meaning | Global analog |
|---|---|---|
| `n_db_total` / `pct_of_db_total` | All pathways in this database, full GSEA universe | 12,221 |
| `n_db_eversig` / `pct_of_db_eversig` | Pathways in this database with FDR < 0.05 in any of 9 GSEA contrasts | 5,267 |
| `n_db_classifiable` / `pct_of_db_classifiable` | Ever-significant pathways minus `Complex` (per mutation); used for Methods majority/plurality claims | varies |

Universe-level totals (`database == 'ALL'` rows in the dual-denominator CSV): 12,221 total, 5,267 ever-significant, 2,533 classifiable (G32A) / 3,086 classifiable (R403C).

### Manuscript anchor → file → denominator

| Manuscript anchor | Denominator | Example wording |
|---|---|---|
| `RESULTS_combio.md` L11 | 5,267 universe | "Complex (2,734/5,267 pathways in G32A [52%], 2,181/5,267 in R403C [42%])" |
| `RESULTS_combio.md` L13 | 5,267 universe | "Compensation … most common (1,462/5,267 [28%] in G32A, 1,612/5,267 [31%] in R403C)" |
| `RESULTS_combio.md` L17 | active-pattern subset (sig TrajDev within 5,267) | "Compensation was the dominant pattern (1,462 pathways [72%] in G32A; 1,612 [72%] in R403C)" |
| `for-the-paper.md` Methods (Edit 1) | classifiable subset | "Compensation … 54.5–58.5% of classifiable pathways in G32A … 46.9–54.5% in R403C" |
| `for-the-paper.md` Supp Fig S8 legend (Edit 5) | classifiable subset | same as Methods |
| Reviewer-response letter | 81-combination range, 5,267 universe | cross-references `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv` |

When adding a new percentage to the paper: compute it from one of the files in this folder (or from `../master_gsea_table.csv`), cite the file path AND the denominator (e.g. "% of `db_eversig`"), and extend this README rather than introducing a new ambiguous column.

---

## Table captions

### `per_database_pattern_dual_denominator.csv` — Supp Table S4 (preferred)

Per-database, per-mutation, per-pattern counts of GSEA pathways classified into the 8-pattern trajectory taxonomy (Compensation, Sign_reversal, Progressive, Late_onset, Transient, Natural_improvement, Natural_worsening, Complex). One row per `(database, mutation, pattern)` triple, plus 16 universe-level rows where `database == 'ALL'`. Each row carries the pattern count `n_pattern` along with all three denominators (`n_db_total`, `n_db_eversig`, `n_db_classifiable`) and the matching `pct_of_db_*` percentages, so a reader can independently verify any percentage cited in the manuscript without recomputing. Classifier thresholds: NES_EFFECT = 0.5, NES_STRONG = 1.0, IMPROVEMENT_RATIO = 0.7, WORSENING_RATIO = 1.3. Universe-level summary: Compensation is the dominant active pattern (1,462/2,533 = 57.7% of classifiable G32A pathways; 1,612/3,086 = 52.2% R403C). For threshold-sensitivity ranges across 81 combinations see `../Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`.

**Key columns:**

| Column | Description |
|---|---|
| `database` | One of 12 GSEA collections (or `ALL` for the universe-level rollup) |
| `mutation` | `G32A` or `R403C` |
| `pattern` | One of the 8 trajectory-pattern labels |
| `n_pattern` | Count of pathways in this database × mutation × pattern cell |
| `n_db_total` | Pathways in this database, full universe |
| `n_db_eversig` | Pathways in this database with FDR < 0.05 in any of 9 contrasts |
| `n_db_classifiable` | Ever-significant pathways minus `Complex` for this mutation |
| `pct_of_db_total`, `pct_of_db_eversig`, `pct_of_db_classifiable` | Percentages against each denominator (`pct_of_db_classifiable` is `NA` for the `Complex` row) |

### `per_database_pattern_summary.csv` — Supp Table S4 (legacy single-denominator view)

Wide-format per-`(database, mutation)` table with one row per pair. Reports `n_total_gene_sets` (per-db analog of 12,221), `n_pathways_ever_sig` (per-db analog of 5,267), and `n_<Pattern>` / `pct_<Pattern>` columns for each of the 8 trajectory patterns. **All `pct_<Pattern>` columns are computed against `n_pathways_ever_sig` (5,267-aligned), not against `n_total_gene_sets`** — the column name does not say so explicitly, which is the reason the dual-denominator file above exists. Preserved verbatim for backward compatibility; the dual-denominator CSV is strictly additive, not a replacement.

### `ribosome_compartment_summary.csv` — Supp Table S3

Per-pathway summary of ribosome / translation-machinery gene sets across three cellular compartments (synaptic via SynGO, mitochondrial via MitoCarta, cytoplasmic via GO:BP/CC/MF). Demonstrates that the three compartments were analyzed simultaneously and have quantitatively distinct trajectories. Reference-set Jaccard overlaps: synaptic ↔ cytoplasmic = 0.574, synaptic ↔ mitochondrial = 0.000, cytoplasmic ↔ mitochondrial = 0.000 — confirming non-redundancy. Compartment-level pattern signatures: synaptic ribosomes show Sign_reversal in both mutations (NES_TrajDev ≈ −2.9 to −3.0); mitochondrial ribosomes show Compensation; cytoplasmic structural ribosome subunits show Late_onset; cytoplasmic ribosome biogenesis shows Compensation.

**Key columns:**

| Column | Description |
|---|---|
| `compartment` | `Synaptic`, `Mitochondrial`, or `Cytoplasmic` |
| `source_database` | SynGO, MitoCarta, GO:BP/CC/MF, Reactome, or KEGG |
| `pathway_id`, `pathway_name`, `n_genes` | Pathway metadata |
| `jaccard_with_synaptic` / `_with_mito` / `_with_cytoplasmic` | Gene-set overlap with the canonical reference set for each compartment |
| `NES_D35_<mut>`, `NES_TrajDev_<mut>`, `NES_D65_<mut>` | Stage-specific NES values for G32A and R403C |
| `pattern_G32A`, `pattern_R403C` | Trajectory-pattern label per mutation |

### `shared_DEGs_G32A_R403C.csv`

38 genes called DE in **both** G32A and R403C at Day 35 (FDR < 0.05 in `G32A_vs_Ctrl_D35` and `R403C_vs_Ctrl_D35`). Supports the reviewer-response argument that the two stalk/GTPase mutations converge on a shared early transcriptional signature.

**Columns:** `Gene`, `G32A_logFC`, `G32A_FDR`, `G32A_dir`, `R403C_logFC`, `R403C_FDR`, `R403C_dir`, `Mean_logFC`.

---

## How to regenerate

```bash
# Trajectory-pattern tables (read-only on master_gsea_table.csv, run in seconds)
python3 02_Analysis/revision/supplements/6a.per_database_pattern_dual_denominator.py
python3 02_Analysis/revision/supplements/6b.per_database_pattern_summary.py

# Ribosome compartment summary (Supp Table S3)
python3 02_Analysis/revision/supplements/6c.ribosome_compartment_summary.py

# Shared DEGs (reads DE_results CSVs)
Rscript 02_Analysis/revision/supplements/5b.extract_shared_maturation_degs.R
```

All four generators are deterministic.

---

## Cross-references

- **Upstream data:** `../master_gsea_table.csv` (canonical pathway × contrast results, ~110K rows); `../DE_results/*.csv` (per-contrast DE statistics); `01_Scripts/Python/pattern_definitions.py` (classifier thresholds and `classify_pattern` function).
- **Sibling sensitivity digest:** `../Supplementary/6a_sensitivity_stability_digest_5267universe.tsv` (81-combination threshold robustness on the 5,267 universe).
- **Pattern taxonomy reference:** `../README.md` and `docs/PATTERN_CLASSIFICATION.md`.

## Read-only constraint

These tables are derived artifacts. Do not hand-edit any `.csv` file in this folder — regenerate from the upstream master table instead.
