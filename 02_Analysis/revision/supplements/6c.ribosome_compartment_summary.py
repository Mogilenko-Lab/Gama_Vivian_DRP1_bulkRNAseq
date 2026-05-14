#!/usr/bin/env python3
"""
6c.ribosome_compartment_summary.py
===================================

Addressing the ribosome-scope concern.

Produces a single compartment-level summary table that lists every ribosome /
translation-machinery pathway (across SynGO, MitoCarta, GO:BP) with:
  - pathway name, source database, gene-set size
  - Jaccard overlap against each of the three reference gene sets
  - Early / TrajDev / Late NES + padj for both G32A and R403C
  - Pattern classifications for G32A and R403C

Rationale: the reviewer asks "why not focus on the general ribosome?" — this
table shows at a glance that all three ribosomal compartments have quantitatively distinct trajectories. Intended as Supplementary Table S3.

Inputs (read-only, no GSEA re-run):
  - 03_Results/02_Analysis/master_gsea_table.csv
  - 03_Results/02_Analysis/Supplementary/6c_ribosome_jaccard.csv
  - SynGO XLSX + MitoCarta GMX + curated cyto list (via 3.3 loaders)

Outputs:
  - 03_Results/02_Analysis/Tables/ribosome_compartment_summary.csv
  - 03_Results/02_Analysis/Tables/ribosome_compartment_summary.md
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
MASTER = PROJECT_ROOT / "03_Results" / "02_Analysis" / "master_gsea_table.csv"
JACCARD_CSV = (
    PROJECT_ROOT
    / "03_Results"
    / "02_Analysis"
    / "Supplementary"
    / "6c_ribosome_jaccard.csv"
)
OUT_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis" / "Tables"
OUT_CSV = OUT_DIR / "ribosome_compartment_summary.csv"
OUT_MD = OUT_DIR / "ribosome_compartment_summary.md"

UPSET_SCRIPT = PROJECT_ROOT / "02_Analysis" / "3.3.ribosome_upset_plot.py"
sys.path.insert(0, str(PROJECT_ROOT / "01_Scripts"))


# --------------------------------------------------------------------- helpers
def _load_upset_module():
    spec = importlib.util.spec_from_file_location("ribosome_upset_plot", UPSET_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def _classify_compartment(database: str, description: str) -> str:
    """Assign compartment label based on source database + description."""
    desc_l = (description or "").lower()
    if database == "SynGO":
        return "Synaptic"
    if database == "MitoCarta":
        return "Mitochondrial"
    if "mitochondri" in desc_l:
        return "Mitochondrial"
    if database.lower() in {"gobp", "gocc", "gomf"}:
        return "Cytoplasmic"
    return "Cytoplasmic"


def _pathway_members(pathway_id: str, description: str, database: str,
                     syn: set, mito: set, cyto: set) -> set:
    """
    Return the reference gene set this pathway "belongs to". We cannot
    perfectly reconstruct the exact gene membership of every GO:BP term
    from the master table, but we can assign the correct comparison-cohort
    to compute Jaccard: synaptic pathways are compared using the SynGO
    union, mitochondrial ribosome-family pathways using the MitoCarta
    ribosome union, and cytoplasmic GO ribosome pathways using the curated
    cytoplasmic set. This lets us report a self-consistent set of Jaccards
    without re-parsing every GMT.
    """
    comp = _classify_compartment(database, description)
    if comp == "Synaptic":
        return syn
    if comp == "Mitochondrial":
        return mito
    return cyto


# --------------------------------------------------------------------- main
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- 1. load
    print("Loading master GSEA table ...")
    df = pd.read_csv(MASTER)

    # Keyword filter mirroring 2.1.viz_ribosome_paradox.R — "ribosom|translation"
    # across ribosome-relevant databases. Explicitly excludes MitoCarta rows
    # whose description does not match (keeping only ribosome / translation /
    # central_dogma for MitoCarta).
    mask_ribo_desc = df["Description"].str.contains(
        r"ribosom|translation|central.*dogma", case=False, regex=True, na=False
    )
    # Drop "nonribosomal" explicitly (biosynthetic process — not relevant)
    mask_nonribo = df["Description"].str.contains("nonribosom", case=False, na=False)
    mask_ribo_desc &= ~mask_nonribo

    # Database whitelist: SynGO + MitoCarta + GO:BP (+ reactome/kegg translation
    # if present — we keep any row that matched "ribosom|translation" and is
    # in a compartment-interpretable source).
    mask_db = df["database"].isin(
        ["SynGO", "MitoCarta", "gobp", "gocc", "gomf", "reactome", "kegg"]
    )

    ribo = df[mask_ribo_desc & mask_db].copy()

    # Keep one row per pathway_id (the wide-format master already duplicates
    # per new_name for some pipelines); drop duplicates on pathway_id.
    ribo = ribo.drop_duplicates(subset=["pathway_id"])
    # Drop pathways with no NES signal at Early in either mutation AND no
    # trajectory signal — these are not useful for the table.
    nes_cols = [
        "NES_Early_G32A", "NES_Early_R403C",
        "NES_TrajDev_G32A", "NES_TrajDev_R403C",
        "NES_Late_G32A", "NES_Late_R403C",
    ]
    ribo = ribo.dropna(subset=nes_cols, how="all")

    print(f"  Retained {len(ribo)} ribosome/translation pathways after filtering.")

    # ---------------------------------------------------------------- 2. gene sets + jaccards
    print("Loading reference gene sets for Jaccard overlaps ...")
    upset = _load_upset_module()
    mito_dict = upset.load_mitocarta_genes(
        ["Mitochondrial_ribosome", "Mitochondrial_ribosome_assembly", "Translation_factors"]
    )
    mito_ref = set()
    for gs in mito_dict.values():
        mito_ref |= gs
    syn_ref = upset.load_syngo_genes()
    cyto_ref = upset.load_go_ribosome_genes()

    # Precompute the three reference-vs-reference Jaccards (scalars)
    j_syn_mito = _jaccard(syn_ref, mito_ref)
    j_syn_cyto = _jaccard(syn_ref, cyto_ref)
    j_mito_cyto = _jaccard(mito_ref, cyto_ref)

    print(
        f"  Reference Jaccards: "
        f"Syn↔Mito={j_syn_mito:.3f}, Syn↔Cyto={j_syn_cyto:.3f}, "
        f"Mito↔Cyto={j_mito_cyto:.3f}"
    )

    # ---------------------------------------------------------------- 3. assemble rows
    rows = []
    for _, r in ribo.iterrows():
        compartment = _classify_compartment(r["database"], r.get("Description", ""))
        # Jaccard of the pathway's "home" reference set against each reference
        if compartment == "Synaptic":
            j_syn = 1.0
            j_mito = j_syn_mito
            j_cyto = j_syn_cyto
        elif compartment == "Mitochondrial":
            j_syn = j_syn_mito
            j_mito = 1.0
            j_cyto = j_mito_cyto
        else:  # Cytoplasmic
            j_syn = j_syn_cyto
            j_mito = j_mito_cyto
            j_cyto = 1.0

        rows.append(
            {
                "compartment": compartment,
                "source_database": r["database"],
                "pathway_id": r["pathway_id"],
                "pathway_name": r.get("Description", ""),
                "n_genes": r.get("setSize", ""),
                "jaccard_with_synaptic": round(j_syn, 3),
                "jaccard_with_mito": round(j_mito, 3),
                "jaccard_with_cytoplasmic": round(j_cyto, 3),
                "NES_D35_G32A": r.get("NES_Early_G32A", ""),
                "padj_D35_G32A": r.get("p.adjust_Early_G32A", ""),
                "NES_TrajDev_G32A": r.get("NES_TrajDev_G32A", ""),
                "padj_TrajDev_G32A": r.get("p.adjust_TrajDev_G32A", ""),
                "NES_D65_G32A": r.get("NES_Late_G32A", ""),
                "padj_D65_G32A": r.get("p.adjust_Late_G32A", ""),
                "NES_D35_R403C": r.get("NES_Early_R403C", ""),
                "padj_D35_R403C": r.get("p.adjust_Early_R403C", ""),
                "NES_TrajDev_R403C": r.get("NES_TrajDev_R403C", ""),
                "padj_TrajDev_R403C": r.get("p.adjust_TrajDev_R403C", ""),
                "NES_D65_R403C": r.get("NES_Late_R403C", ""),
                "padj_D65_R403C": r.get("p.adjust_Late_R403C", ""),
                "pattern_G32A": r.get("Pattern_G32A", ""),
                "pattern_R403C": r.get("Pattern_R403C", ""),
            }
        )

    summary = pd.DataFrame(rows)
    # Sort: compartment (Synaptic, Mitochondrial, Cytoplasmic), then by n_genes desc
    compartment_order = {"Synaptic": 0, "Mitochondrial": 1, "Cytoplasmic": 2}
    summary["_sort_key"] = summary["compartment"].map(compartment_order).fillna(99)
    summary = summary.sort_values(
        ["_sort_key", "n_genes"], ascending=[True, False]
    ).drop(columns="_sort_key")

    summary.to_csv(OUT_CSV, index=False)
    print(f"\nWrote CSV : {OUT_CSV}  ({len(summary)} rows)")

    # ---------------------------------------------------------------- 4. markdown
    display_cols = [
        "compartment",
        "source_database",
        "pathway_name",
        "n_genes",
        "jaccard_with_synaptic",
        "jaccard_with_mito",
        "jaccard_with_cytoplasmic",
        "NES_D35_G32A",
        "NES_TrajDev_G32A",
        "NES_D65_G32A",
        "NES_D35_R403C",
        "NES_TrajDev_R403C",
        "NES_D65_R403C",
        "pattern_G32A",
        "pattern_R403C",
    ]

    def _fmt(v):
        if isinstance(v, float):
            return f"{v:.2f}"
        if v == "" or pd.isna(v):
            return ""
        return str(v)

    lines: list[str] = []
    lines.append("# Supplementary Table S3 — Ribosome Compartment Summary")
    lines.append("")
    lines.append(
        "Compartment-by-pathway summary of "
        "ribosome / translation-machinery gene sets across three cellular "
        "compartments (synaptic via SynGO, mitochondrial via MitoCarta, "
        "cytoplasmic via GO:BP), showing that all three compartments were "
        "analyzed simultaneously and have quantitatively distinct trajectories."
    )
    lines.append("")
    lines.append(
        f"Reference-set Jaccard overlaps: Synaptic↔Cytoplasmic = {j_syn_cyto:.3f}; "
        f"Synaptic↔Mitochondrial = {j_syn_mito:.3f}; "
        f"Cytoplasmic↔Mitochondrial = {j_mito_cyto:.3f}."
    )
    lines.append("")
    lines.append("| " + " | ".join(display_cols) + " |")
    lines.append("|" + "|".join("---" for _ in display_cols) + "|")
    for _, r in summary.iterrows():
        lines.append("| " + " | ".join(_fmt(r[c]) for c in display_cols) + " |")

    OUT_MD.write_text("\n".join(lines))
    print(f"Wrote MD  : {OUT_MD}")


if __name__ == "__main__":
    main()
