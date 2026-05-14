#!/usr/bin/env python3
"""
Supp10.replicate_level_gsva.py
==============================

Per-sample GSVA module scores for the ribosome / translation compartments
that anchor the focused-panel narrative. Each panel is one curated module;
each dot is one biological-replicate sample, coloured by genotype; a thin
boxplot per (genotype × day) cell shows median + IQR. This figure provides
the replicate-level evidence underlying Supp. Fig. S9 (which plots
pathway-level NES means). For each compartment module, the
3–6 biological replicates per (genotype × day) cell are shown as
individual dots so as to show cross-sample reproducibility

Module-to-compartment mapping (matches Supp. Fig. S9 categories)
---    ------------------    ----------    ----------------    ---------
    Synaptic_Ribosomes        - Synaptic Ribosome
    Cytoplasmic_Translation   - Cytoplasmic Ribosome / translation (structural)
    Ribosome_Biogenesis       - Cytoplasmic Ribosome - biogenesis
    Mitochondrial_Ribosome    - Mitochondrial Ribosome - structural
    Mito_Ribosome_Assembly    - Mitochondrial Ribosome - biogenesis

The remaining two modules in `gsva_module_scores.rds` (mtDNA_Maintenance,
OXPHOS) are included as a sixth panel to surface the broader mitochondrial
compensatory programme.

Inputs (read-only)
--    -------    ---
    03_Results/02_Analysis/Plots/Supplementary_10/replicate_level_gsva_per_sample.csv
        - produced by Supp10a.export_gsva_modules.R

Outputs (written to the same folder)
---    --------------    ---    --------
    replicate_level_gsva.pdf
    replicate_level_gsva.png
    replicate_level_gsva_group_summary.csv
    FIGURE_CAPTION.md

Run
---
    Rscript 02_Analysis/revision/supplements/Supp10a.export_gsva_modules.R   # one-time
    python3 02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ===========================================================================
# Paths
# ===========================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = PROJECT_ROOT / "03_Results/02_Analysis/Plots/Supplementary_10"
IN_CSV = OUT_DIR / "replicate_level_gsva_per_sample.csv"
OUT_PDF = OUT_DIR / "replicate_level_gsva.pdf"
OUT_PNG = OUT_DIR / "replicate_level_gsva.png"
OUT_GROUP = OUT_DIR / "replicate_level_gsva_group_summary.csv"


# ===========================================================================
# Panels — module order + display labels
# ===========================================================================
# (1) the two collapse panels first (synaptic, cytoplasmic translation);
# (2) the cytoplasmic biogenesis panel (failed compensation);
# (3) the two mitochondrial panels (clean compensation);
# (4) a final broader-mitochondrial panel that joins mtDNA + OXPHOS.
PANELS: list[tuple[str, str]] = [
    ("Synaptic_Ribosomes", "Synaptic Ribosome"),
    ("Cytoplasmic_Translation", "Cytoplasmic Translation"),
    ("Ribosome_Biogenesis", "Cytoplasmic Ribosome — biogenesis"),
    ("Mitochondrial_Ribosome", "Mitochondrial Ribosome — structural"),
    ("Mito_Ribosome_Assembly", "Mitochondrial Ribosome — biogenesis"),
]

# second-row mitochondrial context panels (ATP Hydrolysis + OXPHOS combined).
EXTRA_PANELS: list[tuple[str, str]] = [
    ("ATP_Hydrolysis", "ATP Hydrolysis"),
    ("OXPHOS", "OXPHOS"),
]


# ===========================================================================
# Caption
# ===========================================================================
FIGURE_CAPTION = (
    "Supp. Fig. S10. Replicate-level GSVA scores confirm the structural-vs-biogenesis "
    "ribosome trajectory split. Each dot is one biological-replicate sample, coloured by "
    "genotype; the thick grey dashed line at y = 0 marks the reference level. Within each "
    "(genotype × day) cell, individual dots show the distribution of per-sample GSVA module "
    "scores, while the white ring marks the within-cell median and the shaded band spans the "
    "IQR. Panels show the five focused-compartment modules (top to bottom, left to right) "
    "plus a combined mitochondrial context panel (ATP Hydrolysis ○ + OXPHOS ◇). "
    "White trajectory lines connect the early (D35) and late (D65) medians for each genotype "
    "so the crossover between structural and biogenesis programs is immediately visible."
)

OUT_CAPTION = OUT_DIR / "FIGURE_CAPTION.md"


# ===========================================================================
# Colours — colourblind-safe
# ===========================================================================
GENOTYPE_COLORS: dict[str, str] = {
    "Ctrl": "#999999",  # gray
    "G32A": "#0072B2",  # Okabe-Ito blue
    "R403C": "#D55E00",  # Okabe-Ito vermilion
}
GENOTYPE_ORDER = ["Ctrl", "G32A", "R403C"]
DAY_ORDER = ["D35", "D65"]


# ===========================================================================
# Typography
# ===========================================================================
FS_SUPTITLE = 19
FS_TITLE = 16
FS_AXIS = 14
FS_TICK = 13
FS_LEGEND = 14
FS_LEGEND_T = 15

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 1.2


# ===========================================================================
# Plotting helpers (low-level geometry)
# ===========================================================================
def _x_position(day: str, genotype: str) -> float:
    """Project (day, genotype) into a single x coordinate.

    Layout: D35 cluster centred at x=1, D65 cluster at x=2; within each
    cluster the three genotypes are offset by ±0.28 so the Ctrl / G32A /
    R403C clusters are clearly separated and each genotype's
    Early -> Late trajectory line is visually parallel across days.
    """
    day_x = {"D35": 1.0, "D65": 2.0}[day]
    geno_off = {"Ctrl": -0.28, "G32A": 0.0, "R403C": +0.28}[genotype]
    return day_x + geno_off


def _deterministic_spread(n: int, half_width: float = 0.07) -> np.ndarray:
    """Symmetric, deterministic horizontal offsets for n sample dots.

    Using a fixed evenly-spaced pattern (instead of random jitter) guarantees
    that the visual centre of the per-sample dot cluster always coincides
    with the genotype's deterministic x_centre. This in turn ensures the
    aggregate (median) marker drawn at x_centre always sits at the visible
    centre of the cluster — across panels and across reruns.
    """
    if n <= 1:
        return np.array([0.0])
    return np.linspace(-half_width, +half_width, n)


# ===========================================================================
# Panel rendering (per-module)
# ===========================================================================
def draw_panel(ax: plt.Axes, df: pd.DataFrame, module: str, title: str) -> None:
    """Per-genotype median trajectory line across D35 -> D65 + sample dots.

    Each genotype gets its own median line connecting Early -> Late so the
    crossover that defines the structural-vs-biogenesis trajectory is
    immediately visible. Per-sample dots use a deterministic symmetric
    spread (NOT random jitter) so the cluster centre always coincides with
    the aggregate ring's x position.
    """
    sub = df[df["module"] == module]

    # Track median values per (genotype, day) for trajectory lines.
    medians: dict[str, list[float]] = {g: [] for g in GENOTYPE_ORDER}
    iqr_low: dict[str, list[float]] = {g: [] for g in GENOTYPE_ORDER}
    iqr_high: dict[str, list[float]] = {g: [] for g in GENOTYPE_ORDER}
    xs: dict[str, list[float]] = {g: [] for g in GENOTYPE_ORDER}

    for day in DAY_ORDER:
        for geno in GENOTYPE_ORDER:
            vals = sub[(sub["days"] == day) & (sub["genotype"] == geno)][
                "gsva_score"
            ].dropna().to_numpy()
            if vals.size == 0:
                continue
            x_centre = _x_position(day, geno)
            color = GENOTYPE_COLORS[geno]

            medians[geno].append(float(np.median(vals)))
            iqr_low[geno].append(float(np.quantile(vals, 0.25)))
            iqr_high[geno].append(float(np.quantile(vals, 0.75)))
            day_x = 1.0 if day == "D35" else 2.0
            xs[geno].append(day_x)

            # Deterministic symmetric spread — visual cluster centre always
            # equals x_centre regardless of N or which samples are present.
            offsets = _deterministic_spread(vals.size)
            ax.scatter(
                np.full_like(vals, x_centre) + offsets,
                vals,
                color=color,
                s=90,
                alpha=0.85,
                edgecolor="white",
                linewidth=1.0,
                zorder=4,
            )

    # Trajectory lines per genotype (median) with shaded IQR band.
    for geno in GENOTYPE_ORDER:
        if len(xs[geno]) < 2:
            continue
        color = GENOTYPE_COLORS[geno]
        ax.fill_between(
            xs[geno], iqr_low[geno], iqr_high[geno],
            color=color, alpha=0.15, linewidth=0, zorder=2
        )
        ax.plot(
            xs[geno], medians[geno], color=color, lw=2.8, alpha=0.95, zorder=3
        )
        ax.scatter(
            xs[geno], medians[geno],
            facecolors="white", edgecolors=color, linewidths=2.6, s=170, zorder=5,
        )

    ax.axhline(0, color="gray", lw=0.9, linestyle="--", alpha=0.7, zorder=1)
    ax.set_xticks([1.0, 2.0])
    ax.set_xticklabels(["Early (D35)", "Late (D65)"],
                       fontsize=FS_TICK, fontweight="bold")
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.set_ylabel("GSVA score", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(
        title, fontsize=FS_TITLE, fontweight="bold", pad=8, loc="left"
    )
    ax.set_xlim(0.4, 2.6)
    ax.grid(axis="y", alpha=0.30, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _draw_combined_extra_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    """Overlay ATP Hydrolysis and OXPHOS as a single panel."""
    # Plot both modules with slight x-offset and different marker shapes so
    # they remain distinguishable; deterministic spread within each sub-group.
    marker_per_module = {"ATP_Hydrolysis": "o", "OXPHOS": "D"}
    x_off_per_module = {"ATP_Hydrolysis": -0.05, "OXPHOS": +0.05}

    for module, _ in EXTRA_PANELS:
        sub = df[df["module"] == module]
        marker = marker_per_module[module]
        x_off = x_off_per_module[module]
        for day in DAY_ORDER:
            for geno in GENOTYPE_ORDER:
                vals = sub[
                    (sub["days"] == day) & (sub["genotype"] == geno)
                ]["gsva_score"].dropna().to_numpy()
                if vals.size == 0:
                    continue
                x_centre = _x_position(day, geno) + x_off
                color = GENOTYPE_COLORS[geno]
                offsets = _deterministic_spread(vals.size, half_width=0.035)
                ax.scatter(
                    np.full_like(vals, x_centre) + offsets,
                    vals,
                    color=color,
                    marker=marker,
                    s=80,
                    alpha=0.85,
                    edgecolor="white",
                    linewidth=0.9,
                    zorder=4,
                )
    ax.axhline(0, color="gray", lw=0.9, linestyle="--", alpha=0.7, zorder=1)
    ax.set_xticks([1.0, 2.0])
    ax.set_xticklabels(["Early (D35)", "Late (D65)"],
                       fontsize=FS_TICK, fontweight="bold")
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.set_ylabel("GSVA score", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(
        "Mitochondrial programme context\n(ATP Hydrolysis ○, OXPHOS ◇)",
        fontsize=FS_TITLE - 1, fontweight="bold", pad=8, loc="left",
    )
    ax.set_xlim(0.4, 2.6)
    ax.grid(axis="y", alpha=0.30, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ===========================================================================
# High-level orchestration helpers
# ===========================================================================
def _write_group_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-(module, genotype, day) stats and write the summary CSV.

    Returns the (in-memory) grouped dataframe so the caller can inspect it.
    """
    group = (
        df.groupby(["module", "genotype", "days"], as_index=False)
        .agg(
            n=("gsva_score", "size"),
            mean=("gsva_score", "mean"),
            sd=("gsva_score", "std"),
            median=("gsva_score", "median"),
        )
    )
    panel_modules = [m for m, _ in PANELS] + [m for m, _ in EXTRA_PANELS]
    group["module"] = pd.Categorical(
        group["module"], categories=panel_modules, ordered=True
    )
    group["genotype"] = pd.Categorical(
        group["genotype"], categories=GENOTYPE_ORDER, ordered=True
    )
    group["days"] = pd.Categorical(
        group["days"], categories=DAY_ORDER, ordered=True
    )
    group = group.sort_values(["module", "genotype", "days"])
    group.to_csv(OUT_GROUP, index=False, float_format="%.4f")
    return group


def _add_legend_and_suptitle(fig: plt.Figure) -> None:
    """Place a shared genotype legend at the top and the suptitle."""
    legend_handles = [
        plt.Line2D(
            [0], [0],
            marker="o", linestyle="None",
            markerfacecolor=GENOTYPE_COLORS[g],
            markeredgecolor="white",
            markersize=12,
            label=g,
        )
        for g in GENOTYPE_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center", bbox_to_anchor=(0.5, 0.965),
        ncol=3, frameon=False, fontsize=FS_LEGEND,
        title="Genotype",
        title_fontsize=FS_LEGEND_T,
        columnspacing=2.4, handletextpad=0.6,
    )
    fig.suptitle(
        "structural-vs-biogenesis ribosome trajectory split in Replicate-level GSVA scores",
        fontsize=FS_SUPTITLE, fontweight="bold", y=0.998,
    )


def _write_figure(fig: plt.Figure) -> None:
    """Save figure to PDF + PNG, close, and print status."""
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_caption(path: Path, text: str) -> None:
    """Write the figure-caption markdown file."""
    path.write_text(text)


# ===========================================================================
# Main
# ===========================================================================
def main() -> None:
    # --- load ---
    print(f"[load] {IN_CSV}")
    df = pd.read_csv(IN_CSV)
    print(
        f"  rows: {len(df)}  modules: {df['module'].nunique()}  "
        f"samples: {df['sample'].nunique()}"
    )

    # --- group summary ---
    group = _write_group_summary(df)
    print(f"[write] {OUT_GROUP}")

    print("\n=== Group means by (module, genotype, day) ===")
    print(
        group[["module", "genotype", "days", "n", "mean"]]
        .assign(mean=lambda d: d["mean"].round(2))
        .to_string(index=False)
    )

    # --- figure ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), sharey=False)
    axes_flat = axes.ravel()

    # Top row: the three "cytoplasmic / synaptic" panels.
    for ax, (module, title) in zip(axes_flat[:3], PANELS[:3]):
        draw_panel(ax, df, module, title)

    # Bottom row: the two mitochondrial ribosome panels + combined panel.
    for ax, (module, title) in zip(axes_flat[3:5], PANELS[3:5]):
        draw_panel(ax, df, module, title)

    _draw_combined_extra_panel(axes_flat[5], df)

    # Shared legend + suptitle.
    _add_legend_and_suptitle(fig)

    # Save figure + caption.
    _write_figure(fig)
    print(f"[write] {OUT_PDF}")
    print(f"[write] {OUT_PNG}")

    _write_caption(OUT_CAPTION, FIGURE_CAPTION)
    print(f"[write] {OUT_CAPTION}")


if __name__ == "__main__":
    main()
