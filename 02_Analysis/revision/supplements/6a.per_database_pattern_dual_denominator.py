#!/usr/bin/env python3
"""
6a.per_database_pattern_dual_denominator.py
============================================

Per-database trajectory-pattern frequencies with all three denominators side-by-side
-------------------------------------------------------------------------------------

PURPOSE
-------
Produce a single per-database × per-mutation × per-pattern table where every count is
paired with **all three denominators** the manuscript cites elsewhere, so a reader can
go straight from any sentence in Methods, Results, or a figure legend to the exact cell
in this table without having to redo the arithmetic.

WHY THIS SCRIPT EXISTS (motivation, 2026-04-24)
-----------------------------------------------
The existing per-database pattern summary `03_Results/02_Analysis/Tables/per_database_pattern_summary.csv`
(produced by `02_Analysis/6b.per_database_pattern_summary.py`) is correct, but it reports
percentages under a single denominator — the per-database `ever_significant` subset (a
per-database analog of the 5,267 ever-significantly-enriched universe used by RESULTS).

Two issues surfaced during the 2026-04-24 external audit of reviewer-concern 6a:

  1. The column `pct_<Pattern>` does not say *which* denominator it is computed against.
     A reader looking at hallmark G32A `pct_Compensation = 37.5` does not know whether
     the denominator is `n_total_gene_sets` (=50, giving 30%), `n_pathways_ever_sig`
     (=40, giving 37.5%, the actual choice), or `classifiable` (=25, giving 60%). All
     three percentages are quoted in different parts of the paper at the universe level;
     a per-database table that only carries one of them invites cross-reference confusion.

  2. The for-the-paper.md Methods edit cites Compensation-as-fraction-of-classifiable
     ranges from the sensitivity grid (54.5–58.5% G32A; 46.9–54.5% R403C across 81 threshold
     combinations on the 5,267 universe). A reviewer who tries to reproduce that ratio
     per-database from the existing CSV would have to compute it themselves; the dual-
     denominator table produced here exposes it as a column.

This script does NOT modify the existing `per_database_pattern_summary.csv`. It writes
strictly additive new artifacts that supersede it for paper-citation purposes.

WHAT THIS SCRIPT DOES
---------------------
For each (gene-set database, mutation, trajectory pattern):
  - Counts the pathways assigned to that pattern.
  - Reports three denominators explicitly:
      * n_db_total       = pathways in this database in the full GSEA universe (per-database
                           analog of the 12,221 universe that `sensitivity_results.csv` uses).
      * n_db_eversig     = pathways in this database with FDR<0.05 in any of the 9 GSEA
                           contrasts (per-database analog of the 5,267 ever-significantly-
                           enriched universe used by RESULTS_combio.md).
      * n_db_classifiable = pathways in this database with `ever_significant=True` AND
                           Pattern_<mut> != "Complex" (per-database analog of the
                           "classifiable" denominator the Methods edit cites for the
                           majority/plurality claims; computed PER MUTATION because Complex
                           membership depends on the mutation).
  - Reports each pattern count as percentages under all three denominators in a single row.

Plus a "universe row" per mutation aggregating across all 12 databases, so the global
12,221 / 5,267 / classifiable numbers appear in the same table as their per-database
breakdowns.

WHAT THIS SCRIPT DOES *NOT* DO
------------------------------
- Does not regenerate `master_gsea_table.csv`, `pattern_definitions.py`, or any GSEA
  output. Pattern labels are read directly from the `Pattern_G32A` / `Pattern_R403C`
  columns already present in the master table at default classifier thresholds.
- Does not re-run the threshold sensitivity grid. For sensitivity-aware ranges see
  `02_Analysis/6a.sensitivity_5267universe.py` and its digest.
- Does not modify or replace `Tables/per_database_pattern_summary.csv`; the existing
  file remains valid and is referenced for backward compatibility in the new README.

REPRODUCIBILITY
---------------
Inputs (read-only):
  - `03_Results/02_Analysis/master_gsea_table.csv`
Outputs:
  - `03_Results/02_Analysis/Tables/per_database_pattern_dual_denominator.csv`
  - `03_Results/02_Analysis/Tables/per_database_pattern_dual_denominator.md`
  - `03_Results/02_Analysis/Tables/README_pattern_summary_denominators.md`
Determinism: fully deterministic (no RNG, single pass over master table).
Runtime: a few seconds.

REFERENCED BY (paper-facing)
----------------------------
- `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/for-the-paper.md`
  (Methods Edit 1; Supp Fig S8 legend Edit 5; per-database citations)
- `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/status.md`
  (Q9 audit-log entry; documents this artifact's denominator semantics)

DENOMINATOR VOCABULARY (used throughout the repo from 2026-04-24 onward)
------------------------------------------------------------------------
At the GLOBAL level (all 12 databases combined):
  * "12,221 universe"   = all pathways in the GSEA universe; used by
                          `Sensitivity_Analysis/sensitivity_results.csv` and the original
                          `Supplementary/6a_sensitivity_stability_digest.tsv`.
  * "5,267 universe"    = pathways with FDR<0.05 in any of the 9 GSEA contrasts; used by
                          `RESULTS_combio.md` lines 11/13/17 and by
                          `Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`.
  * "classifiable subset" = within either universe, the pathways receiving any non-Complex
                          label; this is the denominator the Methods edit cites for
                          Compensation-as-strict-majority and Natural-improvement-as-
                          alternative claims. Computed per mutation.

At the PER-DATABASE level (this table):
  * "n_db_total"        = analog of 12,221 universe at per-database level
  * "n_db_eversig"      = analog of  5,267 universe at per-database level (matches the
                          denominator used in the legacy `per_database_pattern_summary.csv`)
  * "n_db_classifiable" = analog of "classifiable subset" at per-database level

Author: External audit pass, 2026-04-24
"""

from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

# -----------------------------------------------------------------------------
# Paths anchored to repo root (one level up from 02_Analysis/).
# -----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
MASTER_GSEA_CSV = REPO_ROOT / "03_Results" / "02_Analysis" / "master_gsea_table.csv"
TABLES_DIR = REPO_ROOT / "03_Results" / "02_Analysis" / "Tables"
OUT_CSV = TABLES_DIR / "per_database_pattern_dual_denominator.csv"
OUT_MD = TABLES_DIR / "per_database_pattern_dual_denominator.md"
OUT_README = TABLES_DIR / "README_pattern_summary_denominators.md"

PATTERNS = [
    "Compensation",
    "Sign_reversal",
    "Progressive",
    "Late_onset",
    "Transient",
    "Natural_improvement",
    "Natural_worsening",
    "Complex",
]
MUTATIONS = [("G32A", "Pattern_G32A"), ("R403C", "Pattern_R403C")]

# Database display order — same as the legacy per_database_pattern_summary.csv to ease
# side-by-side comparison.
DB_ORDER = [
    "hallmark", "kegg", "reactome", "wiki", "canon", "cgp", "tf",
    "gobp", "gocc", "gomf", "MitoCarta", "SynGO",
]


def fmt_pct(num: int, denom: int) -> str:
    """Format a percentage as 'X.X%' with NaN handling for zero denominators."""
    if denom == 0:
        return "NA"
    return f"{100 * num / denom:.2f}"


def main() -> int:
    print(f"[6a.per_database_pattern_dual_denominator] Loading {MASTER_GSEA_CSV} …")
    df = pd.read_csv(MASTER_GSEA_CSV, low_memory=False)
    uniq = df.drop_duplicates(subset="pathway_id").copy()
    n_master_total = len(uniq)
    n_master_eversig = int(uniq["ever_significant"].sum())
    print(f"   Master table: {len(df):,} long-form rows; {n_master_total:,} unique pathways; "
          f"{n_master_eversig:,} ever-significant.")
    if n_master_total != 12221 or n_master_eversig != 5267:
        # Surface drift loudly — both numbers are cited verbatim in the paper.
        print(f"WARNING: universe sizes drifted ({n_master_total}, {n_master_eversig}); "
              f"manuscript expects (12221, 5267). Investigate before citing.",
              file=sys.stderr)

    rows: list[dict] = []

    # -- Per-database rows --
    for db in DB_ORDER:
        sub_total = uniq[uniq["database"] == db]
        sub_eversig = sub_total[sub_total["ever_significant"] == True]
        n_db_total = len(sub_total)
        n_db_eversig = len(sub_eversig)

        for mut, pat_col in MUTATIONS:
            # 'classifiable' = ever-sig AND not labelled Complex for THIS mutation.
            n_db_classifiable = int((sub_eversig[pat_col] != "Complex").sum())
            for pat in PATTERNS:
                n_pat = int((sub_eversig[pat_col] == pat).sum())
                rows.append({
                    "scope":                  "per_database",
                    "database":               db,
                    "mutation":               mut,
                    "pattern":                pat,
                    "n_pattern":              n_pat,
                    "n_db_total":             n_db_total,
                    "n_db_eversig":           n_db_eversig,
                    "n_db_classifiable":      n_db_classifiable,
                    # Triple-denominator percentages, explicit names so the column
                    # answers the question "% of WHAT".
                    "pct_of_db_total":        fmt_pct(n_pat, n_db_total),
                    "pct_of_db_eversig":      fmt_pct(n_pat, n_db_eversig),
                    "pct_of_db_classifiable": (fmt_pct(n_pat, n_db_classifiable)
                                               if pat != "Complex" else "NA"),
                })

    # -- Universe rows (all 12 databases combined; ALL.universe scope) --
    # These are the global numbers cited at the top of for-the-paper.md and in
    # RESULTS_combio.md L11/13/17. Including them in the same file as the per-db
    # breakdown lets a reader audit cross-table consistency in one place.
    for mut, pat_col in MUTATIONS:
        n_master_classif = int((uniq.loc[uniq["ever_significant"] == True, pat_col] != "Complex").sum())
        for pat in PATTERNS:
            n_pat = int((uniq.loc[uniq["ever_significant"] == True, pat_col] == pat).sum())
            rows.append({
                "scope":                  "ALL_universe",
                "database":               "ALL",
                "mutation":               mut,
                "pattern":                pat,
                "n_pattern":              n_pat,
                "n_db_total":             n_master_total,    # 12,221
                "n_db_eversig":           n_master_eversig,  # 5,267
                "n_db_classifiable":      n_master_classif,
                "pct_of_db_total":        fmt_pct(n_pat, n_master_total),
                "pct_of_db_eversig":      fmt_pct(n_pat, n_master_eversig),
                "pct_of_db_classifiable": (fmt_pct(n_pat, n_master_classif)
                                           if pat != "Complex" else "NA"),
            })

    out = pd.DataFrame(rows)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"   Wrote: {OUT_CSV}")

    # ---- Markdown view ------------------------------------------------------
    # Shows one section per database, plus the universe-level section at the end.
    # Each section pivots so that the three denominator columns appear together
    # for every pattern; this is the format the user asked to see "alongside".
    md_lines: list[str] = [
        "# Per-database trajectory-pattern frequencies — dual-denominator view",
        "",
        f"Source: `{OUT_CSV.relative_to(REPO_ROOT)}` (generated by "
        f"`02_Analysis/6a.per_database_pattern_dual_denominator.py`)",
        "",
        f"All counts and percentages are computed at the **default** classifier thresholds "
        f"(NES_EFFECT=0.5, NES_STRONG=1.0, IMPROVEMENT_RATIO=0.7, WORSENING_RATIO=1.3). For "
        f"sensitivity ranges across 81 threshold combinations on the 5,267 universe see "
        f"`Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`.",
        "",
        "## Denominator vocabulary",
        "",
        "| Denominator key | Meaning | Global analog |",
        "|---|---|---|",
        "| `n_db_total`        | All pathways in this database, full GSEA universe | 12,221 |",
        "| `n_db_eversig`      | Pathways in this database with FDR<0.05 in any of 9 GSEA contrasts | 5,267 |",
        "| `n_db_classifiable` | Pathways in this database that are ever-significant AND not labelled Complex (per-mutation) | varies |",
        "",
        f"_(Universe-level numbers: {n_master_total:,} total; {n_master_eversig:,} ever-significant.)_",
        "",
    ]
    cols = ["pattern", "n_pattern", "pct_of_db_total", "pct_of_db_eversig", "pct_of_db_classifiable"]
    for db in DB_ORDER + ["ALL"]:
        sub = out[out["database"] == db]
        if sub.empty:
            continue
        for mut in ("G32A", "R403C"):
            sub_mut = sub[sub["mutation"] == mut]
            if sub_mut.empty:
                continue
            anchor = sub_mut.iloc[0]
            md_lines.append(f"## {db} — {mut} "
                            f"(n_db_total={anchor['n_db_total']}, "
                            f"n_db_eversig={anchor['n_db_eversig']}, "
                            f"n_db_classifiable={anchor['n_db_classifiable']})")
            md_lines.append("")
            md_lines.append("| pattern | n | % of db_total | % of db_eversig | % of db_classifiable |")
            md_lines.append("|---|---|---|---|---|")
            for _, r in sub_mut.iterrows():
                md_lines.append(
                    f"| {r['pattern']} | {r['n_pattern']} | "
                    f"{r['pct_of_db_total']} | {r['pct_of_db_eversig']} | {r['pct_of_db_classifiable']} |"
                )
            md_lines.append("")

    OUT_MD.write_text("\n".join(md_lines))
    print(f"   Wrote: {OUT_MD}")

    # ---- Tables/-level README documenting all related artifacts and their
    #      denominators. Strictly additive — does not modify the legacy MD.
    readme = f"""# Trajectory-pattern summary tables — denominator reference

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
"""
    OUT_README.write_text(readme)
    print(f"   Wrote: {OUT_README}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
