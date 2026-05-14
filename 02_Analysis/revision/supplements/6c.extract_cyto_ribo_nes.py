#!/usr/bin/env python3
"""
6c.extract_cyto_ribo_nes.py
============================

Pulls exact NES / p.adjust values for cytoplasmic (non-mitochondrial) ribosome
pathways from GO:BP out of master_gsea_table.csv.

Output:
  03_Results/02_Analysis/Supplementary/6c_cytoplasmic_ribo_nes.csv
"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
MASTER = PROJECT_ROOT / "03_Results" / "02_Analysis" / "master_gsea_table.csv"
OUT_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis" / "Supplementary"
OUT_CSV = OUT_DIR / "6c_cytoplasmic_ribo_nes.csv"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MASTER)

    cyto = df[
        (df["database"].str.lower().isin(["gobp", "go:bp"]))
        & (df["Description"].str.contains("ribosom", case=False, na=False))
        & (~df["Description"].str.contains("mitochondri", case=False, na=False))
    ][
        [
            "pathway_id",
            "Description",
            "database",
            "setSize",
            "NES_Early_G32A",
            "NES_Late_G32A",
            "NES_TrajDev_G32A",
            "NES_Early_R403C",
            "NES_Late_R403C",
            "NES_TrajDev_R403C",
            "p.adjust_Early_G32A",
            "p.adjust_Late_G32A",
            "p.adjust_TrajDev_G32A",
            "p.adjust_Early_R403C",
            "p.adjust_Late_R403C",
            "p.adjust_TrajDev_R403C",
            "Pattern_G32A",
            "Pattern_R403C",
        ]
    ].drop_duplicates(subset="Description")

    cyto = cyto.sort_values(["Pattern_G32A", "Description"])
    cyto.to_csv(OUT_CSV, index=False)

    print(f"Found {len(cyto)} cytoplasmic GO:BP ribosome-related pathways.")
    print("-" * 70)
    # Keep the stdout compact but informative
    display_cols = [
        "Description",
        "setSize",
        "NES_Early_G32A",
        "NES_Late_G32A",
        "NES_TrajDev_G32A",
        "NES_Early_R403C",
        "NES_Late_R403C",
        "NES_TrajDev_R403C",
        "Pattern_G32A",
        "Pattern_R403C",
    ]
    with pd.option_context(
        "display.max_rows", None,
        "display.max_colwidth", 70,
        "display.width", 220,
    ):
        print(cyto[display_cols].to_string(index=False))

    print("-" * 70)
    print(f"Wrote: {OUT_CSV}")


if __name__ == "__main__":
    main()
