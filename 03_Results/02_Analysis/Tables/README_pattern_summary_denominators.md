# Trajectory-pattern summary tables — denominator reference

Last updated: 2026-04-24 (external audit of reviewer-concern 6a).

This README disambiguates which **denominator** each pattern-summary table in this
folder uses, since the same `pct_<Pattern>` column name can mean different things in
different files. Always check this README before citing a percentage in the paper.

## Files in this folder

### `per_database_pattern_summary.csv` (legacy, kept for backward compatibility)

- Generator: `02_Analysis/6b.per_database_pattern_summary.py`
- Default thresholds: NES_EFFECT=0.5, NES_STRONG=1.0, IMPROVEMENT_RATIO=0.7, WORSENING_RATIO=1.3.
- Reports one row per `(database, mutation)` pair.
- Columns:
  - `n_total_gene_sets`    — per-database analog of the 12,221 universe
  - `n_pathways_ever_sig`  — per-database analog of the 5,267 universe (pathways with
                              `ever_significant=True`, i.e. FDR<0.05 in any of 9 GSEA contrasts)
  - `n_<Pattern>`          — count of pathways in that pattern
  - `pct_<Pattern>`        — **percentage of `n_pathways_ever_sig`** (i.e. 5,267-aligned, NOT
                              of `n_total_gene_sets`). This is the only denominator this file
                              exposes; the column name does not say so explicitly, hence the
                              dual-denominator file below.

### `per_database_pattern_dual_denominator.csv` (preferred for paper citations)

- Generator: `02_Analysis/6a.per_database_pattern_dual_denominator.py`
- Default thresholds: same as above.
- Reports one row per `(database, mutation, pattern)` triple, plus a `database='ALL'`
  universe-level set of rows aggregating across all 12 databases.
- Every row carries all three denominators side-by-side, with column names that say
  exactly which denominator each percentage is computed against:
  - `n_db_total` and `pct_of_db_total`            — 12,221-universe analog at per-db level
  - `n_db_eversig` and `pct_of_db_eversig`        — 5,267-universe analog at per-db level
  - `n_db_classifiable` and `pct_of_db_classifiable` — classifiable subset (used for
                                                       Methods majority/plurality claims)
- Markdown companion: `per_database_pattern_dual_denominator.md`

### `ribosome_compartment_summary.csv` (separate concern, listed for completeness)

- Generator: `02_Analysis/6c.ribosome_compartment_summary.py`
- This file pertains to reviewer-concern 6c (ribosome compartment Jaccard analysis); it
  does NOT use the trajectory-pattern denominators described above. See its own MD
  companion for its provenance.

## Global denominators (cross-reference)

| Universe | n | Used by | Generating script |
|---|---|---|---|
| 12,221 (all GSEA pathways) | 12,221 | `Sensitivity_Analysis/sensitivity_results.csv`; `Supplementary/6a_sensitivity_stability_digest.tsv` | `02_Analysis/Supp4.sensitivity_analysis.py` |
| 5,267 (ever-significant) | 5,267 | `RESULTS_combio.md` L11/13/17; `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`; `Tables/per_database_pattern_summary.csv`; this dual-denominator table | `02_Analysis/6a.sensitivity_5267universe.py`; `02_Analysis/6b.per_database_pattern_summary.py`; `02_Analysis/6a.per_database_pattern_dual_denominator.py` |
| classifiable subset (per mutation) | varies | `for-the-paper.md` Methods Edit 1 (majority/plurality claims) | downstream of the above |

## Where the manuscript cites each universe

| Manuscript anchor | Denominator | Example wording (from RESULTS_combio.md / for-the-paper.md) |
|---|---|---|
| RESULTS_combio.md L11 | 5,267 universe | "Complex (2,734/5,267 pathways in G32A [52%], 2,181/5,267 in R403C [42%])" |
| RESULTS_combio.md L13 | 5,267 universe | "Compensation … most common (1,462/5,267 pathways [28%] in G32A, 1,612/5,267 [31%] in R403C)" |
| RESULTS_combio.md L17 | active-pattern subset (sig TrajDev within 5,267) | "Compensation was the dominant pattern (1,462 pathways [72%] in G32A; 1,612 [72%] in R403C)" |
| for-the-paper.md Edit 1 (Methods) | 5,267 classifiable subset | "Compensation … 54.5–58.5% of classifiable pathways in G32A … 46.9–54.5% in R403C" |
| for-the-paper.md Edit 5 (Supp Fig S8 legend) | 5,267 classifiable subset | same as Methods |

If you add a new percentage to the paper, please:
  1. Compute it from one of the files in this folder (or from the master table).
  2. Cite the file path AND the denominator (e.g. "% of db_eversig").
  3. If the denominator is not already enumerated above, extend this README rather
     than introducing a new ambiguous column.
