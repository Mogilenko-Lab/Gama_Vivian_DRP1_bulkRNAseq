#!/usr/bin/env python3
"""
Supp8.focused_panel_classifications.py
======================================

Two-panel supplementary figure supporting the focused-panel paragraph that
opens the MitoCarta + SynGO subsections in RESULTS:

    Panel A: MitoCarta (n=72) classification tile — G32A vs R403C
    Panel B: SynGO     (n=32) classification tile — G32A vs R403C

Sorting in each panel groups Compensation / Sign_reversal rows at the top so
that the cross-mutation themes referenced in the paragraph are visually
explicit. The two synaptic ribosome pathways (presyn, postsyn) are outlined
in black in Panel B.

Inputs (read-only):
    03_Results/02_Analysis/master_gsea_table.csv

Outputs:
    03_Results/02_Analysis/Plots/Supplementary_8/
        focused_panel_classifications.pdf
        focused_panel_classifications.png

Run:
    python3 02_Analysis/revision/supplements/Supp8.focused_panel_classifications.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


# =============================================================================
# Paths and constants
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[3]
IN_CSV       = PROJECT_ROOT / "03_Results/02_Analysis/master_gsea_table.csv"
OUT_DIR      = PROJECT_ROOT / "03_Results/02_Analysis/Plots/Supplementary_8"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF      = OUT_DIR / "focused_panel_classifications.pdf"
OUT_PNG      = OUT_DIR / "focused_panel_classifications.png"

# Pattern colors — from 01_Scripts/Python/color_config.py
PATTERN_COLORS = {
    "Compensation":        "#009E73",
    "Sign_reversal":       "#9467BD",
    "Natural_improvement": "#56B4E9",
    "Natural_worsening":   "#E69F00",
    "Late_onset":          "#CC79A7",
    "Transient":           "#0072B2",
    "Complex":             "#F0E442",
}
PATTERN_ORDER = [
    "Compensation",
    "Sign_reversal",
    "Late_onset",
    "Transient",
    "Natural_improvement",
    "Natural_worsening",
    "Complex",
]

RIBOSOME_PATHWAYS = {
    "SynGO::SYNGO:postsyn_ribosome": "postsynaptic ribosome",
    "SynGO::SYNGO:presyn_ribosome":  "presynaptic ribosome",
}

# =============================================================================
# Font sizes — single source of truth so the figure scales legibly
# =============================================================================
FS_SUPTITLE   = 16
FS_TITLE      = 14
FS_AXIS_TICK  = 12       # G32A / R403C column labels
FS_ROW_LABEL  = 11       # pathway names on the y-axis
FS_LEGEND     = 12
FS_LEGEND_TIT = 13

# PDF settings — embed fonts (AGENTS.md guidance)
mpl.rcParams["pdf.fonttype"]  = 42
mpl.rcParams["ps.fonttype"]   = 42
mpl.rcParams["font.family"]   = "DejaVu Sans"


# =============================================================================
# Data loading
# =============================================================================
def load_master() -> pd.DataFrame:
    df = pd.read_csv(IN_CSV)
    keep = [
        "pathway_id", "database", "Description",
        "Pattern_G32A", "Pattern_R403C",
        "ever_significant",
    ]
    df = df[keep].drop_duplicates(subset=["pathway_id"]).reset_index(drop=True)
    return df


def short_label(row: pd.Series) -> str:
    """Compact display label per pathway."""
    if row["database"] == "MitoCarta":
        leaf = row["pathway_id"].split("::", 1)[-1].split(".")[-1]
        return leaf.replace("_", " ")
    if row["database"] == "SynGO":
        return str(row["Description"])
    return str(row["Description"])


# =============================================================================
# Classification tile builder
# =============================================================================
def _pattern_to_idx(series: pd.Series) -> np.ndarray:
    idx_map = {p: i for i, p in enumerate(PATTERN_ORDER)}
    return series.map(idx_map).to_numpy()


def _sort_for_tile(sub: pd.DataFrame) -> pd.DataFrame:
    """Sort pathways so Compensation/Sign reversal blocks group at the top."""
    rank_map = {p: i for i, p in enumerate(PATTERN_ORDER)}
    sub = sub.copy()
    sub["_g_rank"] = sub["Pattern_G32A"].map(rank_map).fillna(99)
    sub["_r_rank"] = sub["Pattern_R403C"].map(rank_map).fillna(99)
    sub["_best"] = sub[["_g_rank", "_r_rank"]].min(axis=1)
    sub = sub.sort_values(
        by=["_best", "_g_rank", "_r_rank", "label"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)
    return sub.drop(columns=["_g_rank", "_r_rank", "_best"])


def draw_classification_tile(ax: plt.Axes, sub: pd.DataFrame,
                             title: str, highlight_ribosomes: bool) -> None:
    """Tile of pathway classifications, rows=pathways, cols=[G32A, R403C]."""
    sub = sub.copy()
    sub["label"] = sub.apply(short_label, axis=1)
    sub = _sort_for_tile(sub)

    mat = np.column_stack([
        _pattern_to_idx(sub["Pattern_G32A"]),
        _pattern_to_idx(sub["Pattern_R403C"]),
    ])

    cmap = ListedColormap([PATTERN_COLORS[p] for p in PATTERN_ORDER])
    ax.imshow(mat, aspect="auto", cmap=cmap,
              vmin=0, vmax=len(PATTERN_ORDER) - 1, interpolation="nearest")

    n_rows = len(sub)
    ax.set_xticks([-0.5, 0.5, 1.5], minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["G32A", "R403C"], fontsize=FS_AXIS_TICK,
                       fontweight="bold")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(sub["label"].tolist(), fontsize=FS_ROW_LABEL)
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", pad=10,
                 loc="left")

    if highlight_ribosomes:
        for i, pid in enumerate(sub["pathway_id"]):
            if pid in RIBOSOME_PATHWAYS:
                ax.add_patch(
                    plt.Rectangle((-0.5, i - 0.5), 2, 1,
                                  fill=False, edgecolor="black", linewidth=2.2)
                )

    for spine in ax.spines.values():
        spine.set_visible(False)


# =============================================================================
# Legend builder
# =============================================================================
def draw_pattern_legend(fig: plt.Figure, present: list[str]) -> None:
    handles = [
        Patch(facecolor=PATTERN_COLORS[p], edgecolor="white",
              label=p.replace("_", " "))
        for p in PATTERN_ORDER if p in present
    ]
    fig.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=len(handles), frameon=False,
        fontsize=FS_LEGEND,
        title="Trajectory pattern (per mutation)",
        title_fontsize=FS_LEGEND_TIT,
    )


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    print(f"[load] {IN_CSV}")
    df = load_master()

    mt = df[df["database"] == "MitoCarta"].copy()
    sy = df[df["database"] == "SynGO"].copy()
    print(f"  MitoCarta unique pathways: {len(mt)}")
    print(f"  SynGO     unique pathways: {len(sy)}")

    for db_name, sub in [("MitoCarta", mt), ("SynGO", sy)]:
        print(f"\n[{db_name}]  Pattern_G32A:")
        print(sub["Pattern_G32A"].value_counts().to_string())
        print(f"[{db_name}]  Pattern_R403C:")
        print(sub["Pattern_R403C"].value_counts().to_string())

    # Two side-by-side panels with generous horizontal gutter (wspace=1.3) so
    # the long MitoCarta and SynGO pathway names don't collide.
    fig = plt.figure(figsize=(18, 16))
    gs  = fig.add_gridspec(
        nrows=1, ncols=2,
        width_ratios=[1.0, 1.0],
        wspace=1.3,
        left=0.18, right=0.98, top=0.90, bottom=0.04,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    draw_classification_tile(
        ax_a, mt,
        f"A. MitoCarta pathways (n={len(mt)})",
        highlight_ribosomes=False,
    )
    draw_classification_tile(
        ax_b, sy,
        f"B. SynGO pathways (n={len(sy)})",
        highlight_ribosomes=True,
    )

    present = sorted(
        set(mt["Pattern_G32A"]) | set(mt["Pattern_R403C"])
        | set(sy["Pattern_G32A"]) | set(sy["Pattern_R403C"]),
        key=PATTERN_ORDER.index,
    )
    draw_pattern_legend(fig, present)

    fig.suptitle(
        "Focused panel: MitoCarta + SynGO trajectory-pattern classifications",
        fontsize=FS_SUPTITLE, fontweight="bold", y=0.995,
    )

    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[write] {OUT_PDF}")
    print(f"[write] {OUT_PNG}")


if __name__ == "__main__":
    main()
