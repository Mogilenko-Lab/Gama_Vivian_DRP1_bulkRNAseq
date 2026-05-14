#!/usr/bin/env python3
"""
6a.sensitivity_5267universe.py — Universe-aligned sensitivity grid for trajectory-pattern classifier
=====================================================================================================

PURPOSE
-------
Re-run the combination threshold sensitivity analysis on the **5,267 ever-significantly enriched pathway subset** that the main-text Results section uses, instead of the full 12,221-pathway GSEA universe used by `02_Analysis/Supp4.sensitivity_analysis.py`.

WHY THIS SCRIPT EXISTS (motivation, 2026-04-24)
-----------------------------------------------
During the Review process one of the reviewers expressed a concern
"Would presenting and plotting the GSEA analysis on the interaction-term contrast not directly accomplish this? In that case there would not be the requirement to rely on a classification system that the authors term 'descriptive and based on heuristic rules with thresholds, that we deemed, strict, but arbitrary in their nature" to the phrase I used in the Manuscript about the "arbitrary" thresholds used for pathway pattern classification. 

Compensation-as-fraction-of-classifiable ranges of 52.3–55.3% (G32A) and 43.6–49.4% (R403C) as based on `03_Results/02_Analysis/Sensitivity_Analysis/sensitivity_results.csv` and from the digest TSV `03_Results/02_Analysis/Supplementary/6a_sensitivity_stability_digest.tsv`.

Both of those upstream artifacts compute pattern frequencies on **all 12,221 pathways** in the GSEA universe, including the ~6,954 pathways that are not significantly enriched in any contrast (those are auto-classified as Complex). The main-text Results section cites percentages on a different denominator: the **5,267 pathway subset** that was significantly enriched (FDR < 0.05) in at least one of the 9 GSEA contrasts.

Because the two universes differ in the Complex count (and hence in the "classifiable" = non-Complex denominator), the percentage ranges differ too:

  5,267 universe (matches Results)         12,221 universe (sensitivity_results.csv)
  --------------------------------         -----------------------------------------
  G32A Comp/classifiable: 54.5–58.5%       G32A Comp/classifiable: 52.3–55.3%
  R403C Comp/classifiable: 46.9–54.5%      R403C Comp/classifiable: 43.6–49.4%
  R403C strict-majority: 54/81 combos       R403C strict-majority: 0/81 combos

A reviewer who reads RESULTS line 11 (Complex G32A = 2,734 / 5,267) and computes
classifiable = 2,533 → Comp/classifiable = 1,462/2,533 = 57.7%, will find that 57.7% is outside the cited [52.3, 55.3] range. The cited range is correct for the 12,221 universe but incorrect for the universe the reviewer is computing in.

The fix is to re-run the sensitivity grid on the 5,267 universe and cite the resulting ranges (Methods Edit, Supp Fig S8 legend Edit 5, and Reviewer Response letter)

WHAT THIS SCRIPT DOES
---------------------
1. Reads the canonical pathway-level GSEA results: `03_Results/02_Analysis/master_gsea_table.csv`.
2. Identifies the 5,267 ever-significantly-enriched pathway IDs as the union of pathwayswith `p.adjust < 0.05` in any of the 9 GSEA contrasts present in the master table:
      - G32A_vs_Ctrl_D35, G32A_vs_Ctrl_D65
      - R403C_vs_Ctrl_D35, R403C_vs_Ctrl_D65
      - Maturation_G32A_specific, Maturation_R403C_specific
      - Time_Ctrl, Time_G32A, Time_R403C
3. Imports the canonical parameterized classifier `classify_pattern_parameterized` from `02_Analysis/Supp4.sensitivity_analysis.py` to ensure byte-for-byte agreement with the published sensitivity analysis. 
4. Re-classifies each of the 5,267 pathways at every one of the 81 threshold combinations
   in the sensitivity grid, for each mutation (G32A and R403C):
      - NES_EFFECT       ∈ {0.4, 0.5 (default), 0.6}
      - NES_STRONG       ∈ {0.8, 1.0 (default), 1.2}
      - IMPROVEMENT_RATIO ∈ {0.6, 0.7 (default), 0.8}
      - WORSENING_RATIO   ∈ {1.25, 1.3 (default), 1.4}
5. Aggregates pattern counts per (combination × mutation), computes the Comp- and Natural_improvement-as-fraction-of-classifiable percentages, and writes:
      a. Full grid: `03_Results/02_Analysis/Supplementary/sensitivity_5267universe.csv`
         (one row per combination × mutation; intended for replication / audit only, not directly cited in the paper).
      b. Compact digest: `03_Results/02_Analysis/Supplementary/6a_sensitivity_stability_digest_5267universe.tsv`
         (reader-facing; cited in for-the-paper.md Methods Edit 1, Supp Fig S8 legend Edit 5, status.md correction block, and the reviewer-response letter).

WHAT THIS SCRIPT DOES *NOT* DO
------------------------------
- Does not modify `master_gsea_table.csv`, `pattern_definitions.py`, `Supp4.sensitivity_analysis.py`, `sensitivity_results.csv`, or any existing PDF/PNG. The original 12,221-universe sensitivity outputs remain valid (they answer a different but related question).
- Does not re-render the published Supp Fig S8 heatmaps (those continue to reflect the 12,221 universe; the legend Edit 5 cites the 5,267-universe ranges in text and points to this script's TSV digest for the underlying numbers).
- Does not regenerate any GSEA result.

REPRODUCIBILITY
---------------
Inputs (read-only):
  - 03_Results/02_Analysis/master_gsea_table.csv
  - 02_Analysis/Supp4.sensitivity_analysis.py (for the canonical classifier)
Outputs:
  - 03_Results/02_Analysis/Supplementary/sensitivity_5267universe.csv
  - 03_Results/02_Analysis/Supplementary/6a_sensitivity_stability_digest_5267universe.tsv
Runtime: ~3 minutes (81 combinations × 2 mutations × 5,267 row-wise classifier calls).
Determinism: fully deterministic; same inputs → identical outputs.
"""

from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Configuration — paths anchored to repo root (one level up from 02_Analysis/).
# Resolving via __file__ keeps the script runnable from any cwd.
# -----------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
MASTER_GSEA_CSV = REPO_ROOT / "03_Results" / "02_Analysis" / "master_gsea_table.csv"
SENSITIVITY_PY = REPO_ROOT / "02_Analysis" / "Supp4.sensitivity_analysis.py"
OUT_DIR = REPO_ROOT / "03_Results" / "02_Analysis" / "Supplementary"
OUT_GRID_CSV = OUT_DIR / "sensitivity_5267universe.csv"
OUT_DIGEST_TSV = OUT_DIR / "6a_sensitivity_stability_digest_5267universe.tsv"

# Sensitivity grid — must mirror Supp4.sensitivity_analysis.py defaults exactly.
NES_EFFECT_VALUES = [0.4, 0.5, 0.6]
NES_STRONG_VALUES = [0.8, 1.0, 1.2]
IMPROVEMENT_RATIO_VALUES = [0.6, 0.7, 0.8]
WORSENING_RATIO_VALUES = [1.25, 1.3, 1.4]
DEFAULT_THRESHOLDS = (0.5, 1.0, 0.7, 1.3)
MUTATIONS = ["G32A", "R403C"]


# -----------------------------------------------------------------------------
# Import the canonical parameterized classifier from Supp4.sensitivity_analysis.py.
# This is the SAME function used to produce sensitivity_results.csv on the 12,221
# universe. By importing it (rather than re-implementing), we guarantee that the
# only difference between this script's output and the published sensitivity is
# the universe being classified.
#
# The Supp4 module has top-level imports of seaborn/matplotlib that are heavy
# but harmless; we tolerate them rather than refactor the upstream module.
# -----------------------------------------------------------------------------
def _load_classifier():
    spec = importlib.util.spec_from_file_location("supp4_sensitivity", SENSITIVITY_PY)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {SENSITIVITY_PY}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["supp4_sensitivity"] = mod
    spec.loader.exec_module(mod)
    return mod.classify_pattern_parameterized


# -----------------------------------------------------------------------------
# Identify the 5,267 ever-significantly-enriched pathway IDs.
#
# The master GSEA table is in long format: each row is one (pathway, contrast)
# pair with its NES and p.adjust. A pathway is "ever-significant" if ANY of its
# nine rows (one per contrast) has p.adjust < 0.05. This matches the criterion
# implicit in RESULTS_combio.md line 9 ("union of pathways significantly
# enriched at one of any key contrasts").
# -----------------------------------------------------------------------------
def identify_ever_significant_pathways(master_df: pd.DataFrame) -> set:
    sig_paths = set(master_df.loc[master_df["p.adjust"] < 0.05, "pathway_id"].unique())
    if len(sig_paths) != 5267:
        # Sanity check: this number is the manuscript-cited 5,267. If it drifts,
        # either master_gsea_table.csv was regenerated or the contrast list changed.
        # Surface it loudly rather than silently embedding an incorrect denominator.
        print(
            f"WARNING: ever-significant pathway count = {len(sig_paths)} "
            f"(manuscript states 5,267). Investigate before citing.",
            file=sys.stderr,
        )
    return sig_paths


# -----------------------------------------------------------------------------
# Load the wide-format pathway-level table (one row per pathway, NES/padj for
# each of Early/TrajDev/Late × G32A/R403C as columns). The wide columns are
# duplicated across long-format rows (one per contrast), so dedup by pathway_id.
# -----------------------------------------------------------------------------
WIDE_COLS = [
    "pathway_id",
    "NES_Early_G32A", "NES_TrajDev_G32A", "NES_Late_G32A",
    "NES_Early_R403C", "NES_TrajDev_R403C", "NES_Late_R403C",
    "p.adjust_Early_G32A", "p.adjust_TrajDev_G32A", "p.adjust_Late_G32A",
    "p.adjust_Early_R403C", "p.adjust_TrajDev_R403C", "p.adjust_Late_R403C",
]


def load_pathway_table(master_df: pd.DataFrame, ever_sig: set) -> pd.DataFrame:
    uniq = master_df[WIDE_COLS].drop_duplicates(subset="pathway_id")
    es = uniq[uniq["pathway_id"].isin(ever_sig)].copy().reset_index(drop=True)
    return es


# -----------------------------------------------------------------------------
# Per-combination classification. Returns a dict of pattern counts plus derived
# percentages (Comp/classifiable, NI/classifiable). 'Classifiable' = all pathways
# receiving any non-Complex label, restricted to the 5,267-pathway universe.
# -----------------------------------------------------------------------------
def classify_grid_cell(
    es: pd.DataFrame, mutation: str, classifier,
    ne: float, ns: float, ir: float, wr: float,
) -> dict:
    s = f"_{mutation}"
    pattern_series = es.apply(
        lambda row: classifier(
            row[f"NES_Early{s}"],   row[f"p.adjust_Early{s}"],
            row[f"NES_TrajDev{s}"], row[f"p.adjust_TrajDev{s}"],
            row[f"NES_Late{s}"],    row[f"p.adjust_Late{s}"],
            nes_effect=ne, nes_strong=ns,
            improvement_ratio=ir, worsening_ratio=wr,
        )[0],
        axis=1,
    )
    counts = pattern_series.value_counts().to_dict()
    n_total = len(es)
    n_complex = counts.get("Complex", 0)
    n_classif = n_total - n_complex
    n_comp = counts.get("Compensation", 0)
    n_ni = counts.get("Natural_improvement", 0)
    return {
        "ne": ne, "ns": ns, "ir": ir, "wr": wr, "mut": mutation,
        "n_total": n_total,
        "Comp": n_comp,
        "SR": counts.get("Sign_reversal", 0),
        "Prog": counts.get("Progressive", 0),
        "NI": n_ni,
        "NW": counts.get("Natural_worsening", 0),
        "LO": counts.get("Late_onset", 0),
        "Tr": counts.get("Transient", 0),
        "Complex": n_complex,
        "Classif": n_classif,
        "Comp_class_pct": 100 * n_comp / max(n_classif, 1),
        "NI_class_pct": 100 * n_ni / max(n_classif, 1),
    }


# -----------------------------------------------------------------------------
# Aggregate the full 81 × 2 grid into a digest TSV that's safe to cite directly
# in for-the-paper.md. Keys mirror the original 6a_sensitivity_stability_digest.tsv
# so that downstream consumers can swap the digest without changing key names.
# -----------------------------------------------------------------------------
def build_digest(grid: pd.DataFrame) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = [
        ("metric", "value"),
        ("universe", "5,267 ever-significantly enriched pathways (matches RESULTS section denominator)"),
        ("definition_ever_significant",
         "FDR<0.05 in any of 9 GSEA contrasts: G32A_vs_Ctrl_{D35,D65}, "
         "R403C_vs_Ctrl_{D35,D65}, Maturation_{G32A,R403C}_specific, Time_{Ctrl,G32A,R403C}"),
        ("n_combinations", "81"),
    ]
    g32 = grid[grid["mut"] == "G32A"].reset_index(drop=True)
    r403 = grid[grid["mut"] == "R403C"].reset_index(drop=True)

    for mut, sub in [("G32A", g32), ("R403C", r403)]:
        out.append((f"comp_dominates_classifiable_{mut}_true",
                    f"{(sub['Comp_class_pct'] > 50).sum()}/81"))
        out.append((f"comp_of_classifiable_{mut}_range",
                    f"{sub['Comp_class_pct'].min():.1f} - {sub['Comp_class_pct'].max():.1f}%"))
        out.append((f"NI_of_classifiable_{mut}_range",
                    f"{sub['NI_class_pct'].min():.1f} - {sub['NI_class_pct'].max():.1f}%"))
        out.append((f"{mut}_Compensation_range_across_grid",
                    f"{sub['Comp'].min()} - {sub['Comp'].max()}"))
        out.append((f"{mut}_Sign_reversal_range_across_grid",
                    f"{sub['SR'].min()} - {sub['SR'].max()}"))
        ne_def, ns_def, ir_def, wr_def = DEFAULT_THRESHOLDS
        d = sub[(sub["ne"] == ne_def) & (sub["ns"] == ns_def)
                & (sub["ir"] == ir_def) & (sub["wr"] == wr_def)].iloc[0]
        out.append((f"default_Compensation_pct_{mut}",
                    f"{100 * d['Comp'] / d['n_total']:.2f}% (n={int(d['Comp'])})"))
        out.append((f"default_Sign_reversal_pct_{mut}",
                    f"{100 * d['SR'] / d['n_total']:.2f}% (n={int(d['SR'])})"))
        out.append((f"default_Complex_pct_{mut}",
                    f"{100 * d['Complex'] / d['n_total']:.2f}% (n={int(d['Complex'])})"))
        out.append((f"default_Comp_pct_classifiable_{mut}",
                    f"{d['Comp_class_pct']:.1f}% (n_classif={int(d['Classif'])})"))

    # Stability of the load-bearing claims (universe-aligned recomputation).
    for mut, sub in [("G32A", g32), ("R403C", r403)]:
        out.append((f"comp_exceeds_passive_{mut}_true",
                    f"{((sub['Comp'] > (sub['NI'] + sub['NW']))).sum()}/81"))
        out.append((f"progressive_rare_{mut}_true",
                    f"{(sub['Prog'] == 0).sum()}/81"))
    out.append(("R403C_more_compensation_true",
                f"{(r403['Comp'].values > g32['Comp'].values).sum()}/81"))
    out.append(("source_script", "02_Analysis/6a.sensitivity_5267universe.py"))
    out.append(("classifier_source",
                "02_Analysis/Supp4.sensitivity_analysis.py:classify_pattern_parameterized"))
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[6a.sensitivity_5267universe] Loading {MASTER_GSEA_CSV} …")
    master = pd.read_csv(MASTER_GSEA_CSV, low_memory=False)

    classifier = _load_classifier()
    ever_sig = identify_ever_significant_pathways(master)
    es = load_pathway_table(master, ever_sig)
    print(f"[6a.sensitivity_5267universe] Classifying {len(es)} pathways "
          f"× {len(MUTATIONS)} mutations × {3*3*3*3} combos …")

    rows: list[dict] = []
    for ne in NES_EFFECT_VALUES:
        for ns in NES_STRONG_VALUES:
            for ir in IMPROVEMENT_RATIO_VALUES:
                for wr in WORSENING_RATIO_VALUES:
                    for mut in MUTATIONS:
                        rows.append(classify_grid_cell(es, mut, classifier, ne, ns, ir, wr))
    grid = pd.DataFrame(rows)
    grid.to_csv(OUT_GRID_CSV, index=False)
    print(f"[6a.sensitivity_5267universe] Full grid → {OUT_GRID_CSV}")

    digest = build_digest(grid)
    with open(OUT_DIGEST_TSV, "w") as f:
        for k, v in digest:
            f.write(f"{k}\t{v}\n")
    print(f"[6a.sensitivity_5267universe] Digest → {OUT_DIGEST_TSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
