# 03_Results/02_Analysis/Tables — pattern-summary and ribosome-compartment tables

## Overview

CSV / Markdown tables that aggregate pathway-level GSEA results into compact, paper-citable
summaries. This pass (2026-04-24) documents the **concern-6a** entries fully; concern-6c
entries are listed at filename level only (their detailed provenance is captured in their
own MD companions and in `Manuscript/current_submission/docs/6c_ribosome_scope/`).

**Sibling reference:** [`README_pattern_summary_denominators.md`](README_pattern_summary_denominators.md) — disambiguates which denominator each `pct_*` column refers to in the trajectory-pattern files. Always consult this file before citing a percentage in the manuscript.

## Universe / denominator

The trajectory-pattern files in this folder use the manuscript-aligned **5,267 ever-significantly enriched pathway universe** (FDR<0.05 in any of 9 GSEA contrasts) by default. The newer `per_database_pattern_dual_denominator.csv` carries all three denominators (`db_total` / `db_eversig` / `db_classifiable`) side-by-side. See `README_pattern_summary_denominators.md` for the global vocabulary.

## Files (concern 6a — fully documented)

| File | Purpose | Generator | Universe / denominator |
|---|---|---|---|
| `per_database_pattern_summary.csv` | Legacy per-database × per-mutation pattern counts (one row per `(database, mutation)`); `pct_<Pattern>` is `% of n_pathways_ever_sig` (i.e. 5,267-aligned) — the column name does not say so explicitly, hence the dual-denominator file below | `02_Analysis/6b.per_database_pattern_summary.py` | 5,267 universe (per-db ever-sig denominator) |
| `per_database_pattern_summary.md` | Markdown companion to the legacy CSV; same content, friendlier for paper-supplement copy-paste | same generator | same |
| `per_database_pattern_dual_denominator.csv` | **Preferred for paper citations.** One row per `(database, mutation, pattern)` triple plus 16 universe-level `database='ALL'` rows. Every row carries `n_pattern` along with three explicitly named denominators (`n_db_total`, `n_db_eversig`, `n_db_classifiable`) and the matching `pct_of_db_*` columns | `02_Analysis/6a.per_database_pattern_dual_denominator.py` | All three side-by-side |
| `per_database_pattern_dual_denominator.md` | Markdown view of the same; broken out per (database, mutation) section with all three percentages in each pattern row | same | same |
| `README_pattern_summary_denominators.md` | Topic-specific reference defining the three denominators, mapping them to manuscript anchors, listing the generators for each universe | hand-written, 2026-04-24 audit | — |

## Files (concern 6c — listed for completeness)

| File | Concern | Generator |
|---|---|---|
| `ribosome_compartment_summary.csv` | Ribosome compartment Jaccard (concern 6c) | `02_Analysis/6c.ribosome_compartment_summary.py` |
| `ribosome_compartment_summary.md` | Markdown companion | same |

For 6c provenance, see `Manuscript/current_submission/docs/6c_ribosome_scope/`.

## Cross-references

- **Manuscript anchors that cite files in this folder:**
  - `RESULTS_combio.md` L11/13/17 — trajectory-pattern percentages (5,267 universe)
  - `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/for-the-paper.md` Edit 1 (Methods), Edit 5 (Supp Fig S8 legend), reviewer-response letter
  - `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/status.md` Q6 / Q7 audit-log entries
- **Sibling READMEs:** [`README_pattern_summary_denominators.md`](README_pattern_summary_denominators.md)
- **Upstream data:** `03_Results/02_Analysis/master_gsea_table.csv` (canonical pathway × contrast results); `01_Scripts/Python/pattern_definitions.py` (classifier thresholds and `classify_pattern` function).

## How to regenerate

```bash
# Legacy single-denominator per-database table
python3 02_Analysis/6b.per_database_pattern_summary.py

# Dual-denominator paired table (preferred for new paper citations)
python3 02_Analysis/6a.per_database_pattern_dual_denominator.py
```

Both scripts are deterministic, read-only on `master_gsea_table.csv`, and run in seconds.

## Read-only constraints

- These tables are **derived artifacts**. Do not hand-edit any `.csv` or `.md` file in this folder; regenerate from the upstream master table instead.
- The legacy `per_database_pattern_summary.csv` (Apr 23) is preserved verbatim — the dual-denominator file is **strictly additive**, not a replacement.
