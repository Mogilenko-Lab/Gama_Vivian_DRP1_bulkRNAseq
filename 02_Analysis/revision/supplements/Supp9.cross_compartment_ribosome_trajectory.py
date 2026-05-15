#!/usr/bin/env python3
"""
Supp9.cross_compartment_ribosome_trajectory.py
==============================================

Aggregated NES trajectory across five ribosome-related compartments, two
mutations. Each compartment line is the mean NES over its member pathways at
Early (D35) and Late (D65) stages, with a shaded ±1 SE band. Individual
pathway-level values are overlaid as semi-transparent dots so within-category
dispersion is visible.

Why this figure exists
----------------------
The reviewer-driven revision (`docs/6c_ribosome_scope/`) reframes the
"Ribosome Paradox" along a structural-vs-biogenesis axis rather than a
synaptic-vs-cytoplasmic axis. The five-category taxonomy used here makes that
reframe visually unambiguous in a single panel pair:

    Compensation arc  (DOWN -> UP)              Sign-reversal arc (UP -> DOWN)
    ----------------------------------          ----------------------------------
    Mitochondrial Ribosome - structural         Cytoplasmic Ribosome - structural
    Mitochondrial Ribosome - biogenesis         Synaptic Ribosome  (amplitude max)
    Cytoplasmic Ribosome - biogenesis

Mitochondrial ribosomes (structural AND biogenesis) travel with the
cytoplasmic ribosome biogenesis program towards compensation; 
cytoplasmic structural ribosomes do the opposite, and the 
synaptic-localised RPL/RPS proteins sit at the
amplitude extreme of that opposite arc.

Semantic aggregation rationale (one paragraph per category)
-----------------------------------------------------------
Each category groups GSEA gene sets that annotate the same molecular
machinery class. The same RPL/RPS proteins appear in several of these gene
sets (e.g. the SynGO synaptic-ribosome term is 98.6 %-contained in the
curated cytoplasmic ribosomal proteome) - the divergent trajectories
therefore reflect differential ranked coherence of gene subsets within the
GSEA statistic, not biologically distinct gene families. See `docs/6c_*` and
the README in this folder for the full curation argument.

1. Mitochondrial Ribosome - structural
   - MitoCarta "Mitochondrial_ribosome" 
   - GO:CC mitochondrial large/small ribosomal subunit 
   - GO:CC organellar ribosome (effectively means mitochondrial - no plastid ribosomes).
   Members are the assembled MRPL/MRPS subunits.

2. Mitochondrial Ribosome - biogenesis
   - MitoCarta "Mitochondrial_ribosome_assembly" 
   - GO:BP mitochondrial ribosome assembly term.
   Members are accessory factors (chaperones,
   modifying enzymes) that build the mito ribosome but are not subunits of
   the assembled particle.

3. Cytoplasmic Ribosome - biogenesis
   - GO:BP ribosome biogenesis and assembly, 
   - the GO:BP small/large subunit biogenesis terms, ribosomal subunit export from nucleus, 
   - the GO:CC pre-ribosome compartments. 
   Members are nucleolar / nuclear-export biogenesis factors (NOP, FBL, PES1, BOP1, DKC1, etc.), not assembled cytosolic ribosomes.

4. Cytoplasmic Ribosome - structural
   The assembled 80S compartment: GO:CC cytosolic ribosome, cytosolic
   large/small ribosomal subunit, the generic large ribosomal subunit and
   ribosomal subunit, GO:CC ribosome, and the GO:MF structural constituent
   of ribosome term. Members are RPL and RPS proteins, mostly cytoplasmic.

5. Synaptic Ribosome
   The two SynGO Cellular-Component ribosome terms (presyn_ribosome,
   postsyn_ribosome). These annotate RPL/RPS proteins by experimental
   evidence of synaptic localisation, NOT by molecular distinctness; the
   gene-set membership is essentially (98-100%) a subset of the curated
   cytoplasmic ribosomal proteome. The divergent trajectory of the SynGO
   subset is therefore a localisation-coherence signal within the same
   structural-ribosome family that category 4 captures in its non-localised
   form.

Inputs (read-only)
------------------
    03_Results/02_Analysis/master_gsea_table.csv

Outputs (this script writes only into the output folder below)
--------------------------------------------------------------
    03_Results/02_Analysis/Plots/Supplementary_9/
        cross_compartment_ribosome_trajectory.pdf
        cross_compartment_ribosome_trajectory.png
        cross_compartment_ribosome_trajectory_per_pathway.csv
        cross_compartment_ribosome_trajectory_aggregate.csv
        FIGURE_CAPTION.md

Run
---
    python3 02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# =============================================================================
# Paths
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[3]
IN_CSV       = PROJECT_ROOT / "03_Results/02_Analysis/master_gsea_table.csv"
OUT_DIR      = PROJECT_ROOT / "03_Results/02_Analysis/Plots/Supplementary_9"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF      = OUT_DIR / "cross_compartment_ribosome_trajectory.pdf"
OUT_PNG      = OUT_DIR / "cross_compartment_ribosome_trajectory.png"
OUT_LONG_CSV = OUT_DIR / "cross_compartment_ribosome_trajectory_per_pathway.csv"
OUT_AGG_CSV  = OUT_DIR / "cross_compartment_ribosome_trajectory_aggregate.csv"


# =============================================================================
# Category membership - single declarative source of truth
# =============================================================================
# Keep this dictionary in sync with the docstring above and with the project
# README in 03_Results/02_Analysis/Plots/Supplementary_9/. Changes here are
# the only edit required to add or remove a member pathway.
CATEGORIES: dict[str, list[str]] = {
    # Assembled mitochondrial ribosome: MRPL/MRPS subunits.
    "Mitochondrial Ribosome — structural": [
        "MitoCarta::Mitochondrial_central_dogma.Translation.Mitochondrial_ribosome",
        "gocc::GOCC_MITOCHONDRIAL_LARGE_RIBOSOMAL_SUBUNIT",
        "gocc::GOCC_MITOCHONDRIAL_SMALL_RIBOSOMAL_SUBUNIT",
        # "Organellar ribosome" in mammals is effectively the mito ribosome.
        "gocc::GOCC_ORGANELLAR_RIBOSOME",
    ],
    # Accessory factors that build the mitochondrial ribosome.
    "Mitochondrial Ribosome — biogenesis": [
        "MitoCarta::Mitochondrial_central_dogma.Translation.Mitochondrial_ribosome_assembly",
        "gobp::GOBP_MITOCHONDRIAL_RIBOSOME_ASSEMBLY",
    ],
    # Nucleolar / nuclear-export factors that build the cytoplasmic ribosome
    # (BOP1, PES1, NOP56/58, FBL, DKC1, etc.).
    "Cytoplasmic Ribosome — biogenesis": [
        "gobp::GOBP_RIBOSOME_BIOGENESIS",
        "gobp::GOBP_RIBOSOME_ASSEMBLY",
        "gobp::GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS",
        "gobp::GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
        "gobp::GOBP_RIBOSOMAL_SUBUNIT_EXPORT_FROM_NUCLEUS",
        "gocc::GOCC_PRERIBOSOME",
        "gocc::GOCC_90S_PRERIBOSOME",
        "gocc::GOCC_PRERIBOSOME_LARGE_SUBUNIT_PRECURSOR",
        "gocc::GOCC_PRERIBOSOME_SMALL_SUBUNIT_PRECURSOR",
    ],
    # Assembled cytoplasmic 80S ribosome: RPL/RPS proteins.
    "Cytoplasmic Ribosome — structural": [
        "gocc::GOCC_CYTOSOLIC_RIBOSOME",
        "gocc::GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT",
        "gocc::GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT",
        "gocc::GOCC_LARGE_RIBOSOMAL_SUBUNIT",
        "gocc::GOCC_RIBOSOMAL_SUBUNIT",
        "gocc::GOCC_RIBOSOME",
        "gomf::GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
    ],
    # SynGO localisation-annotated RPL/RPS proteins (98.6%-contained subset
    # of the cytoplasmic-structural category above).
    "Synaptic Ribosome": [
        "SynGO::SYNGO:presyn_ribosome",
        "SynGO::SYNGO:postsyn_ribosome",
    ],
}

# Plot order (top-to-bottom in legend). Compensation arc first, sign-reversal
# arc last, synaptic at the very end so the eye lands on the extreme amplitude.
CATEGORY_ORDER = [
    "Mitochondrial Ribosome — structural",
    "Mitochondrial Ribosome — biogenesis",
    "Cytoplasmic Ribosome — biogenesis",
    "Cytoplasmic Ribosome — structural",
    "Synaptic Ribosome",
]

# Colorblind-safe palette (Wong 2011, Nature Methods 8:441). All five colors
# are distinguishable under deuteranopia and protanopia. Compensation
# (DOWN->UP) categories use the cool/green end; Sign-reversal categories use
# the warm end; synaptic gets vermilion as the strongest signal.
CATEGORY_COLORS: dict[str, str] = {
    "Mitochondrial Ribosome — structural":  "#0072B2",  # Wong blue
    "Mitochondrial Ribosome — biogenesis":  "#56B4E9",  # Wong sky blue
    "Cytoplasmic Ribosome — biogenesis":    "#009E73",  # Wong bluish green
    "Cytoplasmic Ribosome — structural":    "#E69F00",  # Wong orange
    "Synaptic Ribosome":                    "#D55E00",  # Wong vermilion
}


# =============================================================================
# Typography - tuned for legibility as a small journal-paper inset
# =============================================================================
FS_SUPTITLE  = 19
FS_TITLE     = 17
FS_AXIS      = 16
FS_TICK      = 15
FS_LEGEND    = 14
FS_LEGEND_T  = 15

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 1.2


# =============================================================================
# Data prep
# =============================================================================
def build_long_table(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (pathway_id, mutation, stage) with NES."""
    keep = [
        "pathway_id",
        "NES_Early_G32A",  "NES_Late_G32A",
        "NES_Early_R403C", "NES_Late_R403C",
    ]
    sub = df[keep].drop_duplicates(subset="pathway_id").copy()

    long = sub.melt(
        id_vars=["pathway_id"],
        value_vars=keep[1:],
        var_name="col",
        value_name="NES",
    )
    long[["_nes", "stage", "mutation"]] = long["col"].str.split("_", n=2, expand=True)
    long = long.drop(columns=["col", "_nes"])
    return long


def categorise(long: pd.DataFrame) -> pd.DataFrame:
    """Attach the category column based on CATEGORIES membership."""
    pid_to_cat: dict[str, str] = {}
    for cat, members in CATEGORIES.items():
        for pid in members:
            pid_to_cat[pid] = cat
    long = long.copy()
    long["category"] = long["pathway_id"].map(pid_to_cat)
    long = long.dropna(subset=["category"])
    return long


def summarise(long: pd.DataFrame) -> pd.DataFrame:
    """Mean + SE of NES per category × mutation × stage."""
    agg = (
        long.groupby(["category", "mutation", "stage"], as_index=False)
            .agg(mean_NES=("NES", "mean"),
                 sd_NES=("NES", "std"),
                 n=("NES", "size"))
    )
    agg["se_NES"] = agg.apply(
        lambda r: (r["sd_NES"] / np.sqrt(r["n"])) if r["n"] > 1 else 0.0,
        axis=1,
    )
    return agg


# =============================================================================
# Plotting
# =============================================================================
STAGE_ORDER = ["Early", "Late"]
STAGE_LABEL = {"Early": "Early (D35)", "Late": "Late (D65)"}


def draw_panel(ax: plt.Axes, agg: pd.DataFrame, long: pd.DataFrame,
               mutation: str, title: str) -> None:
    x_positions = np.arange(len(STAGE_ORDER))

    for cat in CATEGORY_ORDER:
        color = CATEGORY_COLORS[cat]
        means, ses = [], []
        for stage in STAGE_ORDER:
            row = agg[(agg["category"] == cat)
                      & (agg["mutation"] == mutation)
                      & (agg["stage"] == stage)]
            if row.empty:
                means.append(np.nan); ses.append(0.0)
            else:
                means.append(float(row["mean_NES"].iloc[0]))
                ses.append(float(row["se_NES"].iloc[0]))

        means_a = np.asarray(means)
        ses_a   = np.asarray(ses)

        # SE band
        ax.fill_between(x_positions, means_a - ses_a, means_a + ses_a,
                        color=color, alpha=0.20, linewidth=0)

        # Per-pathway dots. Two adjustments to keep low-n / tightly-clustered
        # categories visible:
        #   (1) jitter width adapts to category size — small n gets wider
        #       jitter so the dots don't all land under the mean marker;
        #   (2) per-pathway dots are drawn ABOVE the mean line so a dot is
        #       never fully occluded.
        sub = long[(long["category"] == cat) & (long["mutation"] == mutation)]
        # Seed per (category, mutation) so jitter is stable across runs but
        # different categories don't get identical jitter patterns.
        rng = np.random.default_rng(abs(hash((cat, mutation))) % (2**32))
        for stage, x in zip(STAGE_ORDER, x_positions):
            vals = sub[sub["stage"] == stage]["NES"].dropna().to_numpy()
            if vals.size == 0:
                continue
            jitter_half = 0.18 if vals.size <= 3 else 0.09
            jitter = rng.uniform(-jitter_half, jitter_half, size=vals.size)
            ax.scatter(np.full_like(vals, x) + jitter, vals,
                       color=color, s=96, alpha=0.85,
                       edgecolor="white", linewidth=1.0, zorder=5)

        # Mean line — drawn as an OPEN ring (no fill) so the per-pathway
        # dots in front of it remain visible even when they coincide with
        # the mean position.
        ax.plot(x_positions, means_a, color=color, lw=3.2, alpha=0.95,
                zorder=4, label=cat)
        ax.scatter(x_positions, means_a,
                   facecolors="white", edgecolors=color, linewidths=3.0,
                   s=240, zorder=6)

    ax.axhline(0, color="gray", lw=1.0, linestyle="--", alpha=0.7)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([STAGE_LABEL[s] for s in STAGE_ORDER],
                       fontsize=FS_TICK, fontweight="bold")
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.set_ylabel("Mean NES across category", fontsize=FS_AXIS, fontweight="bold")
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", pad=12, loc="left")
    ax.grid(axis="y", alpha=0.30, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # x-axis with a bit of breathing room
    ax.set_xlim(-0.35, len(STAGE_ORDER) - 0.65)


# =============================================================================
# Caption (written as a side-output for manuscript reuse)
# =============================================================================
FIGURE_CAPTION = """\
# Supplementary Figure S9 — caption

**Cross-compartment ribosome trajectories: a structural-vs-biogenesis split
that maps unevenly onto compartments.** Each panel shows the mean NES across
all member pathways of five compartment-level categories at Early (35 DIV)
and Late (65 DIV) stages, for the G32A (A) and R403C (B) mutations
respectively. Shaded bands are ±1 SE around the mean; small semi-transparent
dots are the underlying per-pathway NES values. The five categories are
defined as follows (full member lists in
`cross_compartment_ribosome_trajectory_per_pathway.csv`):
*Mitochondrial Ribosome — structural* (MitoCarta `Mitochondrial_ribosome`,
GO:CC mitochondrial large / small ribosomal subunit, GO:CC organellar
ribosome; n=4); *Mitochondrial Ribosome — biogenesis* (MitoCarta
`Mitochondrial_ribosome_assembly`, GO:BP mitochondrial ribosome assembly;
n=2); *Cytoplasmic Ribosome — biogenesis* (GO:BP ribosome biogenesis /
assembly / large- and small-subunit biogenesis / subunit export from
nucleus; GO:CC pre-ribosome and pre-ribosome precursor compartments; n=9);
*Cytoplasmic Ribosome — structural* (GO:CC cytosolic ribosome and its large
/ small subunits, GO:CC large ribosomal subunit and ribosomal subunit,
GO:CC ribosome, GO:MF structural constituent of ribosome; n=7); and
*Synaptic Ribosome* (SynGO `SYNGO:presyn_ribosome`, `SYNGO:postsyn_ribosome`;
n=2). The three categories shown in cool / green tones all trace a
Compensation-like arc — early negative enrichment recovering toward zero or
positive late enrichment — together with the broader mitochondrial
compensatory programme. The two categories shown in warm tones trace the
opposite Sign-reversal arc, with the SynGO synaptic-ribosome category at the
largest amplitude. The split is therefore structural-vs-biogenesis rather
than mitochondrial-vs-cytoplasmic or synaptic-vs-bulk: the assembled
cytoplasmic ribosome (and its synaptic-localised subset) collapses, while
its biogenesis machinery, the assembled mitochondrial ribosome, and the
mitochondrial biogenesis programme all recover. Aggregated and per-pathway
numerical values, including the category × mutation × stage means used here,
are provided in `cross_compartment_ribosome_trajectory_aggregate.csv` and
`cross_compartment_ribosome_trajectory_per_pathway.csv`. Source script:
`02_Analysis/revision/supplements/Supp9.cross_compartment_ribosome_trajectory.py`.
Colors are from the Wong (2011) colorblind-safe palette.
"""


# =============================================================================
# Main
# =============================================================================
def main() -> None:
    print(f"[load] {IN_CSV}")
    df = pd.read_csv(IN_CSV)

    long = build_long_table(df)
    long = categorise(long)

    # Sanity: which declared pathways were not found in the master table?
    found = set(long["pathway_id"].drop_duplicates().tolist())
    declared = {p for ms in CATEGORIES.values() for p in ms}
    missing = sorted(declared - found)
    if missing:
        print(f"[warn] {len(missing)} declared pathways not in master table:")
        for m in missing:
            print(f"   - {m}")

    # Per-category, per-pathway count
    print("\n=== Per-category n_pathways ===")
    cnt = (long.drop_duplicates(["pathway_id", "category"])
                .groupby("category")["pathway_id"].nunique()
                .reindex(CATEGORY_ORDER))
    print(cnt.to_string())

    agg = summarise(long)
    print("\n=== Aggregate means (NES) ===")
    pretty = (agg.assign(value=lambda d: d["mean_NES"].round(2))
                 .pivot_table(index=["category", "mutation"], columns="stage",
                              values="value")
                 .reindex(CATEGORY_ORDER, level=0))
    print(pretty.to_string())

    # ----- Write CSV side-outputs ----------------------------------------
    long_out = long[["pathway_id", "category", "mutation", "stage", "NES"]]\
        .sort_values(["category", "mutation", "stage", "pathway_id"])
    long_out.to_csv(OUT_LONG_CSV, index=False, float_format="%.4f")

    agg_out = (agg.assign(category=lambda d: pd.Categorical(
                    d["category"], categories=CATEGORY_ORDER, ordered=True))
                  .sort_values(["category", "mutation", "stage"]))
    agg_out[["category", "mutation", "stage",
             "n", "mean_NES", "sd_NES", "se_NES"]]\
        .to_csv(OUT_AGG_CSV, index=False, float_format="%.4f")

    print(f"[write] {OUT_LONG_CSV}")
    print(f"[write] {OUT_AGG_CSV}")

    # ----- Plot ----------------------------------------------------------
    # Wider canvas + tall enough for legend space; tight_layout will scale
    # axes proportionally.
    fig, axes = plt.subplots(1, 2, figsize=(18, 10), sharey=True)
    draw_panel(axes[0], agg, long, "G32A",  "A. G32A")
    draw_panel(axes[1], agg, long, "R403C", "B. R403C")

    # Shared legend at the top - two rows to keep label font large
    handles, labels = axes[0].get_legend_handles_labels()
    seen = set()
    ordered = [(h, l) for h, l in zip(handles, labels)
               if not (l in seen or seen.add(l))]
    fig.legend(
        [h for h, _ in ordered], [l for _, l in ordered],
        loc="upper center", bbox_to_anchor=(0.5, 0.955),
        ncol=3, frameon=False, fontsize=FS_LEGEND,
        title="Ribosome / translation compartment",
        title_fontsize=FS_LEGEND_T,
        columnspacing=2.4, handletextpad=0.9, labelspacing=0.6,
    )

    fig.suptitle(
        "Cross-compartment ribosome trajectories: structural-vs-biogenesis split maps unevenly onto compartments",
        fontsize=FS_SUPTITLE, fontweight="bold", y=0.998,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.84))
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {OUT_PDF}")
    print(f"[write] {OUT_PNG}")

    # ----- Caption -------------------------------------------------------
    caption_path = OUT_DIR / "FIGURE_CAPTION.md"
    caption_path.write_text(FIGURE_CAPTION)
    print(f"[write] {caption_path}")


if __name__ == "__main__":
    main()
