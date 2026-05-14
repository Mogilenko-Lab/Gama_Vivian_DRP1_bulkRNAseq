"""
Supplementary Figure: NES Trajectory Curves for the Focused Panel
(MitoCarta 72 + SynGO 32 = 104 unique pathways)

Purpose: Direct visual support for the "Focused Panel" paragraph describing the
104 co-enriched MitoCarta + SynGO pathways, their pattern classifications,
and the two synaptic-ribosome Sign reversal pathways.

Design:
  - Curved (bezier) NES trajectories (Early → TrajDev → Late)
  - Color encodes the 8-class pattern taxonomy
  - SynGO synaptic ribosome pathways highlighted with annotations
  - MitoCarta and SynGO sections separated by a dashed divider
  - Ranks shown beside ribosome pathways

Input:  /03_Results/02_Analysis/master_gsea_table.csv
Output: 03_Results/…/Plots/revision/Supp_Fig_bump_focused_MitoCarta_SynGO.pdf
"""

# ============================================================================
# 1. NAMESPACE
# ============================================================================
import os
import warnings

import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import FancyBboxPatch
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT_ROOT = "/workspaces/Gama_Vivian_DRP1_bulkRNAseq"
MASTER_TABLE = os.path.join(
    PROJECT_ROOT, "03_Results", "02_Analysis", "master_gsea_table.csv"
)
PLOTS_DIR = os.path.join(
    PROJECT_ROOT, "03_Results", "02_Analysis", "Plots", "revision"
)
os.makedirs(PLOTS_DIR, exist_ok=True)

# ============================================================================
# 2. CONFIGURATION
# ===============================================================================

PATTERN_COLORS = {
    "Compensation":      "#FF6B35",   # warm orange-red
    "Sign reversal":     "#6B8E23",   # green
    "Late onset":        "#8B5CF6",   # purple
    "Transient":         "#E8A94A",   # gold
    "Natural improvement": "#417690",  # blue
    "Natural worsening": "#CD5C5C",   # indianred
    "Progressive":       "#00A86B",   # dark green
    "Complex":           "#9E9E9E",   # gray
}

RIBOSOME_ID = {
    "SYNGO:presyn_ribosome": "presynaptic ribosome",
    "SYNGO:postsyn_ribosome": "postsynaptic ribosome",
}


# ============================================================================
# 3. DATA
# ===============================================================================

def load_master() -> pd.DataFrame:
    """Load master table, deduplicate on pathway_id, filter to focused panel."""
    df = pd.read_csv(MASTER_TABLE)
    # Keep first occurrence per unique pathway_id (unique across all databases)
    df = df.groupby("pathway_id").first().reset_index()
    # Only the 5,267 ever-significant pathways
    df = df[df["ever_significant"] == "True"].reset_index(drop=True)
    # Focused panel
    df = df[df["database"].isin(["MitoCarta", "SynGO"])].reset_index(drop=True)
    # Ensure key columns are numeric
    for col in [
        "NES_TrajDev_G32A", "NES_TrajDev_R403C",
        "NES_Early_G32A",   "NES_Early_R403C",
        "NES_Late_G32A",    "NES_Late_R403C",
        "p.adjust_TrajDev_G32A", "p.adjust_TrajDev_R403C",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def compute_ranks(df: pd.DataFrame) -> pd.DataFrame:
    """|NES_TrajDev| rank within the 5,267-pathway universe (method='min')."""
    for mut, rank_col in [
        ("G32A", "rank_G32A"), ("R403C", "rank_R403C"),
    ]:
        n_e = df[f"NES_Early_{mut}"].fillna(0)
        n_t = df[f"NES_TrajDev_{mut}"].fillna(0)
        order = -(n_t.abs() + n_e.abs())
        df[rank_col] = order.rank(method="min").astype(int)
    return df


# ============================================================================
# 4. PLOTTING
# ===============================================================================

def bezier3(p0, p1, p2, n=50):
    """Quadratic Bézier through (p0, p1, p2)."""
    t = np.linspace(0, 1, n)
    x = (1 - t)**2 * p0[0] + 2 * (1 - t) * t * p1[0] + t**2 * p2[0]
    y = (1 - t)**2 * p0[1] + 2 * (1 - t) * t * p1[1] + t**2 * p2[1]
    return x, y


def pat_color(pat: str) -> str:
    return PATTERN_COLORS.get(pat, "#CCCCCC")


def plot_focused_panel():
    # -- data ---------------------------------
    df = load_master()
    df = compute_ranks(df)

    # Split into G32A rows
    g32a = df[df["mutation"] == "G32A"].reset_index(drop=True)
    n_mito = len(g32a[g32a["database"] == "MitoCarta"])
    n_syngo = len(g32a[g32a["database"] == "SynGO"])

    # Sort pathways within each section alphabetically
    mito_data = g32a[g32a["database"] == "MitoCarta"].sort_values("Description").reset_index(drop=True)
    syngo_data = g32a[g32a["database"] == "SynGO"].sort_values("Description").reset_index(drop=True)

    # -- figure ------------------------------
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 16), dpi=300, facecolor="white")

    y_gap = 0.50
    LINE = 1.2

    # Y starting positions
    y_top = 15.0
    mito_rows = list(range(n_mito))
    syngo_start = y_top - n_mito * y_gap - 1.1
    syngo_rows = [syngo_start - i * y_gap for i in range(n_syngo)]

    # -- draw pathway trajectories ------
    for i in range(n_mito):
        row = mito_rows.loc[i]
        y = y_top - i * y_gap
        pat = row.get("Pattern_G32A", "")
        color = pat_color(pat)
        sig = row["p.adjust_TrajDev_G32A"] < 0.05 and abs(row["NES_TrajDev_G32A"]) > 0.5
        alpha = 0.85 if sig else 0.35

        p0 = (0, y)
        p1 = (1, y + row["NES_TrajDev_G32A"])
        p2 = (2, y + row["NES_Late_G32A"])

        xc, yc = bezier3(p0, p1, p2)
        ax.plot(xc, yc, color=color, lw=LINE, alpha=alpha)
        ax.plot(p0[0], p0[1], "o", color=color, ms=2.5, alpha=alpha)
        ax.plot(p2[0], p2[1], "s", color=color, ms=2.5, alpha=alpha)

        # Label (Right side)
        desc = row["Description"]
        ax.text(2.02, y + row["NES_Late_G32A"], desc,
                fontsize=7.5, fontweight="semibold" if sig else "normal",
                va="center", ha="left", color=color,
                family="sans-serif", alpha=min(1, 0.4 + 0.6 * sig))

    # -- divider ------
    div_y = (y_top - n_mito * y_gap - syngo_start) / 2 + syngo_start - y_gap / 2
    ax.plot([-0.3, 2.5], [div_y, div_y], "gray", ls="--", lw=0.8, alpha=0.5)
    ax.text(-0.22, div_y, "SynGO", fontsize=6.5, color="gray", va="center",
            ha="center", rotation=90, family="sans-serif")

    # -- SynGO trajectories ------
    for i in range(n_syngo):
        row = syngo_data.loc[i]
        y = syngo_start - i * y_gap
        pat = row.get("Pattern_G32A", "")
        color = pat_color(pat)
        sig = row["p.adjust_TrajDev_G32A"] < 0.05 and abs(row["NES_TrajDev_G32A"]) > 0.5
        alpha = 0.85 if sig else 0.35

        p0 = (0, y)
        p1 = (1, y + row["NES_TrajDev_G32A"])
        p2 = (2, y + row["NES_Late_G32A"])

        xc, yc = bezier3(p0, p1, p2)
        ax.plot(xc, yc, color=color, lw=LINE, alpha=alpha, clip_on=False)
        ax.plot(p0[0], p0[1], "o", color=color, ms=2.5, alpha=alpha)
        ax.plot(p2[0], p2[1], "s", color=color, ms=2.5, alpha=alpha)

        # Ribosome row: annotate with name + ranks
        if row["ID"] in RIBOSOME_ID:
            label = RIBOSOME_ID[row["ID"]]
            ax.text(2.02, y + row["NES_Late_G32A"], label,
                    fontsize=8.5, fontweight="bold",
                    va="center", ha="left", color="#C62828",
                    family="sans-serif")
            r_g32 = row.get("rank_G32A", "")
            r_r403 = row.get("rank_R403C", "")
            ax.text(1.4,  y + row["NES_TrajDev_G32A"] + 0.15,
                   f"#{r_g32} G32A", fontsize=7, fontweight="bold",
                   color="#C62828", va="bottom", ha="center", family="sans-serif")
            ax.text(1.0,  y + row["NES_TrajDev_G32A"] - 0.12,
                   f"#{r_r403} R403C", fontsize=7, fontweight="bold",
                   color="#C62828", va="top", ha="center", family="sans-serif")

    # -- Ribosome highlight box ------
    if n_syngo > 0:
        ribo_y_syngo = [syngo_rows[i] for i, row in enumerate(syngo_data.itertuples())
                        if row.ID in RIBOSOME_ID]
        if ribo_y_syngo:
            lo = min(ribo_y_syngo) - y_gap / 2
            hi = max(ribo_y_syngo) + y_gap / 2
            bx = FancyBboxPatch(
                (-0.3, lo), 2.8, hi - lo,
                linewidth=1.2, edgecolor="#C62828", facecolor="#FDEDEC",
                alpha=0.12, clip_on=False,
            )
            ax.add_patch(bx)

    # -- formatting ------
    ax.set_xlim(-0.35, 2.4)
    ax.set_ylim(syngo_start - n_syngo * y_gap - 1.0, 15.3)

    # X labels
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Early (35 DIV)", "TrajDev (Δ35→65)", "Late (65 DIV)"],
                        fontsize=10, fontweight="bold", family="sans-serif")

    # Shared NES axis on the right
    ax.set_ylabel("NES (standardized)", fontsize=9, fontweight="bold", family="sans-serif")
    ax.tick_params(left=True, labelleft=True, left_labelsize=7)
    # Set y-ticks by NES range
    ax.set_yticks(range(-6, 7))
    for spine in ax.spines.values():
        spine.set_visible(False)

    # -- title ------
    fig.suptitle(
        "Focused Panel: NES Trajectories — MitoCarta (72) + SynGO (32) Pathways",
        fontsize=14, fontweight="bold", color="#1B3A5C",
        family="sans-serif", y=1.01,
    )
    ax.text(0.5, -0.02,
            "Color = trajectory pattern  |  Line opacity/width = significance  |  O = Early  □ = Late  |  "
            "Syn-gos ribosome ranks within 5,267 ever-significant pathways:\n"
            "  #7(G32A) / #23(R403C) postsynaptic   |   #17(G32A) / #49(R403C) presynaptic",
            transform=ax.transAxes, fontsize=7, color="#555555",
            ha="center", va="top", family="sans-serif", style="italic")

    # -- legend ------
    legend_lines = [
        mlines.Line2D([0], [0], color="#FF6B35", lw=LINE, label="Compensation"),
        mlines.Line2D([0], [0], color="#6B8E23", lw=LINE, label="Sign reversal"),
        mlines.Line2D([0], [0], color="#8B5CF6", lw=LINE, label="Late onset"),
        mlines.Line2D([0], [0], color="#417690", lw=LINE, label="Natural improvement"),
        mlines.Line2D([0], [0], color="#E8A94A", lw=LINE, label="Transient"),
        mlines.Line2D([0], [0], color="#CD5C5C", lw=LINE, label="Worsening"),
        mlines.Line2D([0], [0], color="#9E9E9E", lw=LINE, label="Complex"),
    ]
    lgd = ax.legend(
        handles=legend_lines, loc="upper left", ncol=3,
        framealpha=0.92, fontsize=7.5,
        title="Trajectory Pattern",
        frameon=True, edgecolor="#AAAAAA", facecolor="white",
    )

    # -- save ------
    fig.tight_layout()
    out = os.path.join(PLOTS_DIR, "Supp_Fig_bump_focused_MitoCarta_SynGO.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # -- summary ------
    print(f"✅ Written: {out}")
    print(f"   MitoCarta: {n_mito}  SynGO: {n_syngo}  Total: {n_mito + n_syngo}")
    print(f"   MitoCarta patterns (G32A): {dict(mito_rows['Pattern_G32A'].value_counts())}")
    print(f"   SynGO      patterns (G32A): {dict(syngo_data['Pattern_G32A'].value_counts())}")


# ============================================================================
if __name__ == "__main__":
    plot_focused_panel()
