#!/usr/bin/env python3
"""
6a.geometric_scatter.py
-----------------------

Geometric visualization of trajectory patterns as (NES_Late, NES_Early)
scatters, colored by super-category.

The three dimensions of the classifier —
    X  = NES_Late    (D65 mutant vs. control — developmental outcome)
    Y  = NES_Early   (D35 mutant vs. control — initial state)
    color = super-category, which encodes the TrajDev gate
            (Active_* = TrajDev significant; Passive = TrajDev not significant)
— are all present in the 2D image.  Color is the third axis: Active and Passive
patterns occupy the same (NES_Late, NES_Early) quadrant but are distinguished
by whether TrajDev was statistically significant.

Input:
    03_Results/02_Analysis/master_gsea_table.csv

Output:
    03_Results/02_Analysis/Plots/Pattern_Summary_Normalized/Supplementary_6a/
        geometric_scatter_G32A.{pdf,png}
        geometric_scatter_R403C.{pdf,png}
        geometric_scatter_both_mutations.{pdf,png}

Run:
    python3 02_Analysis/revision/supplements/6a.geometric_scatter.py
"""

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "01_Scripts/Python"))

from pattern_definitions import SUPER_CATEGORY_COLORS, SUPER_CATEGORY_MAP  # noqa: E402

IN_CSV  = PROJECT_ROOT / "03_Results/02_Analysis/master_gsea_table.csv"
OUT_DIR = (
    PROJECT_ROOT
    / "03_Results/02_Analysis/Plots/Pattern_Summary_Normalized/Supplementary_6a"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Canonical plotting order — same as pattern_definitions.SUPER_CATEGORY_ORDER
CAT_ORDER = [
    "Active_Compensation",
    "Active_Reversal",
    "Active_Progression",
    "Passive",
    "Late_onset",
    "Other",
    "Insufficient_data",
]


# ──────────────────────────────────────────────────────────────────────────────
# Low-level helpers
# ──────────────────────────────────────────────────────────────────────────────

def _prepare(df: pd.DataFrame, mut: str) -> pd.DataFrame:
    """Extract and annotate per-mutation ever-significant pathways."""
    col_e = f"NES_Early_{mut}"
    col_l = f"NES_Late_{mut}"
    col_p = f"Pattern_{mut}"
    col_t = f"NES_TrajDev_{mut}"

    sub = df.drop_duplicates(subset="pathway_id").copy()
    sub = sub[sub["ever_significant"] == True]  # noqa: E712
    sub = sub.dropna(subset=[col_e, col_l, col_p])
    sub["_super_cat"] = sub[col_p].map(SUPER_CATEGORY_MAP).fillna("Other")
    sub["_color"]     = sub["_super_cat"].map(SUPER_CATEGORY_COLORS).fillna("#999999")
    sub["_abs_traj"]  = sub[col_t].abs().fillna(0)
    sub["_early"]     = sub[col_e]
    sub["_late"]      = sub[col_l]
    sub["_pathway"]   = sub["Description"].fillna(sub["pathway_id"])
    sub["_database"]  = sub["database"].fillna("unknown")

    # Label-priority score: pathways from curated biological DBs with large
    # effect on BOTH axes get a 100-pt bonus so they rank above generic
    # high-TrajDev pathways from CGP/tf/cancer signatures.
    _PREF_DBS  = {"MitoCarta", "SynGO", "gobp", "gocc"}
    _PREF_NES  = 1.0
    sub["_label_priority"] = (
        sub["_database"].isin(_PREF_DBS) &
        (sub["_early"].abs() > _PREF_NES) &
        (sub["_late"].abs()  > _PREF_NES)
    ).astype(float) * 100.0 + sub["_abs_traj"]

    return sub


# Preferred databases for label selection (curated biological scope)
_PREF_DBS = {"MitoCarta", "SynGO", "gobp", "gocc"}


def _place_labels(ax, sub: pd.DataFrame, top_n: int, lim: float) -> None:
    """
    Label the top-N pathways by label-priority score.

    Priority score = 100 (if preferred DB + |Early|>1 + |Late|>1) + |TrajDev|.
    This ensures MitoCarta/SynGO/gobp/gocc pathways with strong effects on both
    axes are labelled in preference to generic high-TrajDev CGP/cancer entries.

    All labels for a quadrant are anchored to a fixed vertical column
    (volcano-plot style): no overlap, no y-axis collision.
    Annotations use clip_on=True so text boxes never extend outside the axes.
    """
    score_col = "_label_priority" if "_label_priority" in sub.columns else "_abs_traj"
    top = sub.nlargest(top_n, score_col)

    quad_rows: dict = defaultdict(list)
    for _, r in top.iterrows():
        qx = int(float(r["_late"])  >= 0)   # Late on X
        qy = int(float(r["_early"]) >= 0)   # Early on Y
        quad_rows[(qx, qy)].append(r)

    # Labels always use ha="left" and are anchored near the respective
    # axis edge so text extends INWARD (toward plot centre), never outward.
    # Right quadrant: anchor at +30 % of lim (text grows rightward, safe).
    # Left  quadrant: anchor at -lim + margin (text grows rightward, safe).
    MARGIN   = 0.12   # gap from axis edge in data units
    MIN_DY   = 0.38
    MAX_CHARS = 24    # hard cap — at fontsize 5.5 this stays inside either half

    for (qx, qy), rows in quad_rows.items():
        rows = sorted(rows, key=lambda r: float(r["_early"]), reverse=(qy == 1))

        # Both sides: ha="left", text grows toward centre
        ha = "left"
        if qx == 1:   # right quadrant: put column at 30 % of lim
            x_col = lim * 0.30
        else:         # left quadrant: put column just inside left edge
            x_col = -lim + MARGIN

        y_start = float(np.clip(float(rows[0]["_early"]), -lim + 0.30, lim - 0.30))
        placed_y: list = [y_start]
        direction = -1 if qy == 1 else 1

        for i, r in enumerate(rows):
            raw = r["_pathway"] if isinstance(r["_pathway"], str) else str(r["_pathway"])
            db  = str(r.get("_database", ""))
            # Truncate first, then prepend DB tag for preferred databases
            short = (raw[:MAX_CHARS] + "…") if len(raw) > MAX_CHARS else raw
            label = f"[{db}] {short}" if db in _PREF_DBS else short

            ly = y_start if i == 0 else float(
                np.clip(placed_y[-1] + direction * MIN_DY, -lim + 0.15, lim - 0.15)
            )
            placed_y.append(ly)

            ax.annotate(
                label,
                xy=(float(r["_late"]), float(r["_early"])),
                xytext=(x_col, ly),
                fontsize=5.5,
                ha=ha,
                va="center",
                annotation_clip=False,   # never suppress or clip
                arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.5,
                                shrinkA=3, shrinkB=2,
                                relpos=(0.0, 0.5)),
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
            )


def _build_legend_handles() -> list:
    """Patch handles for all non-Insufficient_data super-categories."""
    return [
        mpatches.Patch(color=v, label=k)
        for k, v in SUPER_CATEGORY_COLORS.items()
        if k != "Insufficient_data"
    ]


# ──────────────────────────────────────────────────────────────────────────────
# High-level orchestration
# ──────────────────────────────────────────────────────────────────────────────

def _plot_one(
    ax,
    sub: pd.DataFrame,
    mut_label: str,
    top_n_label: int = 4,
    legend: bool = True,
    title: str | None = None,
) -> None:
    """Draw one geometric scatter panel onto *ax*."""
    lim = float(min(max(sub["_early"].abs().max(), sub["_late"].abs().max(), 3.0), 4.5))

    ax.axhline(0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.axvline(0, color="k", lw=0.7, ls="--", alpha=0.6)
    # Diagonal: NES_Late = NES_Early  →  x = y  (no trajectory change)
    ax.plot([-lim, lim], [-lim, lim], color="grey", lw=0.5, ls=":", alpha=0.5,
            label="_nolegend_")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    for cat in CAT_ORDER:
        s = sub[sub["_super_cat"] == cat]
        if s.empty:
            continue
        ax.scatter(s["_late"], s["_early"],   # x=Late, y=Early
                   c=SUPER_CATEGORY_COLORS.get(cat, "#999999"),
                   alpha=0.42, s=10, linewidths=0, label=cat)

    _place_labels(ax, sub, top_n_label, lim)

    ax.set_xlabel("NES Late (D65 mutant vs. control)",  fontsize=10)   # x
    ax.set_ylabel("NES Early (D35 mutant vs. control)", fontsize=10, labelpad=8)  # y
    ax.set_title(
        title or f"{mut_label}: trajectory patterns in (NES_Late, NES_Early) space",
        fontsize=11, fontweight="bold",
    )
    ax.set_aspect("equal", "box")

    # Quadrant labels: x=Late, y=Early
    corner_specs = [
        ( lim * 0.95,  lim * 0.95, "Sustained\nup",                          "right", "top"),
        (-lim * 0.95, -lim * 0.95, "Sustained\ndown",                        "left",  "bottom"),
        (-lim * 0.95,  lim * 0.95, "Early-up /\nLate-down\n(Sign_reversal)", "left",  "top"),
        ( lim * 0.95, -lim * 0.95, "Early-down /\nLate-up\n(Sign_reversal)", "right", "bottom"),
    ]
    for tx, ty, txt, ha, va in corner_specs:
        ax.text(tx, ty, txt, ha=ha, va=va, fontsize=7, color="grey", alpha=0.70)

    if legend:
        handles = _build_legend_handles()
        ax.legend(handles=handles, fontsize=7, loc="upper left",
                  bbox_to_anchor=(1.02, 1.0), frameon=True,
                  edgecolor="#dddddd", facecolor="white",
                  title="Super-category", title_fontsize=7.5)


def main() -> None:
    df = pd.read_csv(IN_CSV)

    # ── Single-mutation figures ──────────────────────────────────────────────
    for mut in ("G32A", "R403C"):
        sub = _prepare(df, mut)
        print(f"{mut}: {len(sub)} ever-significant pathways")
        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        _plot_one(ax, sub, mut_label=mut, legend=True)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            out = OUT_DIR / f"geometric_scatter_{mut}.{ext}"
            fig.savefig(out, bbox_inches="tight", **({"dpi": 300} if ext == "png" else {}))
            print(f"  Saved: {out}")
        plt.close(fig)

    # ── Side-by-side figure ──────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))
    for ax, mut in zip(axes, ("G32A", "R403C")):
        sub = _prepare(df, mut)
        _plot_one(ax, sub, mut_label=mut, top_n_label=3, legend=False)

    handles = _build_legend_handles()
    fig.legend(handles=handles, fontsize=7, loc="center left",
               bbox_to_anchor=(1.00, 0.50), frameon=True,
               edgecolor="#dddddd", facecolor="white",
               title="Super-category", title_fontsize=7.5)
    fig.suptitle(
        "Trajectory patterns in (NES_Late, NES_Early) space \u2014 "
        "geometric interpretation of the classifier boundaries\n"
        "(x = Late outcome; y = initial state; color = TrajDev significance gate: "
        "Active vs. Passive)",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.30)
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"geometric_scatter_both_mutations.{ext}"
        fig.savefig(out, bbox_inches="tight", **({"dpi": 300} if ext == "png" else {}))
        print(f"  Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
