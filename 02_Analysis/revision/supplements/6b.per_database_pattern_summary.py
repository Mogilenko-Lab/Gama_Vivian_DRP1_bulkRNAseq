#!/usr/bin/env python3
"""
6b.per_database_pattern_summary.py

Concern: quantify per-database trajectory-pattern frequencies to address the
SynGO "why specifically synaptic?" concern.

Input:  03_Results/02_Analysis/master_gsea_table.csv
Output: 03_Results/02_Analysis/Tables/per_database_pattern_summary.csv
        03_Results/02_Analysis/Tables/per_database_pattern_summary.md
        03_Results/02_Analysis/Plots/Supplementary_6b/per_database_pattern_heatmap.{pdf,png}

Run:  python3 02_Analysis/6b.per_database_pattern_summary.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
IN_CSV       = PROJECT_ROOT / "03_Results/02_Analysis/master_gsea_table.csv"

TABLES_DIR   = PROJECT_ROOT / "03_Results/02_Analysis/Tables"
PLOTS_DIR    = PROJECT_ROOT / "03_Results/02_Analysis/Plots/Supplementary_6b"
for d in (TABLES_DIR, PLOTS_DIR):
    d.mkdir(parents=True, exist_ok=True)

OUT_CSV  = TABLES_DIR / "per_database_pattern_summary.csv"
OUT_MD   = TABLES_DIR / "per_database_pattern_summary.md"
OUT_PDF  = PLOTS_DIR / "per_database_pattern_heatmap.pdf"
OUT_PNG  = PLOTS_DIR / "per_database_pattern_heatmap.png"

print(f"Project root identified as: {PROJECT_ROOT}")
print(f"Inputs: {IN_CSV}")
print(f"Outputs will be saved to: {TABLES_DIR} and {PLOTS_DIR}")


PATTERNS = [
    "Sign_reversal",
    "Compensation",
    "Progressive",
    "Late_onset",
    "Transient",
    "Natural_improvement",
    "Natural_worsening",
    "Complex",
]
MUTATIONS = [("G32A", "Pattern_G32A"), ("R403C", "Pattern_R403C")]

DB_ORDER = [
    "hallmark", "kegg", "reactome", "wiki", "canon", "cgp", "tf",
    "gobp", "gocc", "gomf", "MitoCarta", "SynGO",
]

md_cols = [
    "database", "mutation", "n_total_gene_sets", "n_pathways_ever_sig",
    "n_Sign_reversal", "pct_Sign_reversal",
    "n_Compensation", "pct_Compensation",
    "n_Complex", "pct_Complex",
]


def df_to_md(df_):
    """Convert a DataFrame to a Markdown table string."""
    cols = list(df_.columns)
    out_lines = ["| " + " | ".join(cols) + " |",
                 "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df_.iterrows():
        out_lines.append("| " + " | ".join(
            f"{r[c]:.2f}" if isinstance(r[c], float) else str(r[c])
            for c in cols
        ) + " |")
    return "\n".join(out_lines)


def load_and_filter(in_csv, patterns):
    """Load CSV, deduplicate, and return the ever-significant subset."""
    df = pd.read_csv(in_csv)
    uniq = df.drop_duplicates(subset="pathway_id").copy()
    print(f"Loaded {len(df):,} long-form rows; {len(uniq):,} unique pathways.")
    print(f"ever_significant universe: {int(uniq['ever_significant'].sum()):,}")
    n_total_per_db = uniq.groupby("database").size().rename("n_total_gene_sets")
    sig_df = uniq[uniq["ever_significant"] == True].copy()
    return sig_df, n_total_per_db


def build_table_df(sig_df, n_total_per_db, db_order, mutations, patterns,
                   out_path=None):
    """Build table rows and write CSV if out_path is provided."""
    rows = []
    for db in db_order:
        grp = sig_df[sig_df["database"] == db]
        n_sig = len(grp)
        for mut, pat_col in mutations:
            vc = grp[pat_col].value_counts() if n_sig else pd.Series(dtype=int)
            row = {
                "database":              db,
                "mutation":              mut,
                "n_total_gene_sets":     int(n_total_per_db.get(db, 0)),
                "n_pathways_ever_sig":   int(n_sig),
            }
            for p in patterns:
                row[f"n_{p}"]   = int(vc.get(p, 0))
                row[f"pct_{p}"] = round(100 * vc.get(p, 0) / n_sig, 2) if n_sig else 0.0
            rows.append(row)

    out = pd.DataFrame(rows)
    if out_path:
        out.to_csv(out_path, index=False)
        print(f"Saved: {out_path}")
    return out


md_df_cols = [
    "database", "mutation", "n_total_gene_sets", "n_pathways_ever_sig",
    "n_Sign_reversal", "pct_Sign_reversal",
    "n_Compensation", "pct_Compensation",
    "n_Complex", "pct_Complex",
]


def export_markdown(out_df, out_path):
    """Write per-database pattern summary markdown table."""
    md = out_df[md_df_cols].copy()
    with open(out_path, "w") as f:
        f.write("# Per-database Trajectory Pattern Summary\n\n")
        f.write("Source: `03_Results/02_Analysis/master_gsea_table.csv` "
                "(ever_significant subset, FDR<0.05 in at least one stage "
                "of either mutation).\n\n")
        f.write("## Key columns\n\n")
        f.write(df_to_md(md))
        f.write("\n\n## Full column set (extended table)\n\n")
        f.write(df_to_md(out_df))
        f.write("\n")
    print(f"Saved: {out_path}")


PAT_CONFIGS = [
    ("Sign_reversal",       "Reds"),
    ("Compensation",        "Blues"),
    ("Progressive",         "Reds"),
    ("Late_onset",          "Blues"),
    ("Transient",           "Reds"),
    ("Natural_improvement", "Blues"),
    ("Natural_worsening",   "Reds"),
    ("Complex",             "Blues"),
]


def build_heatmap(out_df, sig_df, pdf_path, png_path, db_order, patterns):
    """Build 2x4 heatmap panels and save as PDF + PNG."""
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

    g32_baselines = {}
    r40_baselines = {}
    for pat in patterns:
        g32_baselines[pat] = round(
            100 * (sig_df["Pattern_G32A"] == pat).sum() / len(sig_df), 2)
        r40_baselines[pat] = round(
            100 * (sig_df["Pattern_R403C"] == pat).sum() / len(sig_df), 2)

    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    axes = axes.flatten()

    for ax, (pat, cmap) in zip(axes, PAT_CONFIGS):
        pivot = (
            out_df.pivot(index="database", columns="mutation",
                         values=f"pct_{pat}")
               .reindex(db_order)
        )
        baseline = (g32_baselines[pat], r40_baselines[pat])

        sns.heatmap(
            pivot, ax=ax, annot=True, fmt=".1f", cmap=cmap,
            linewidths=0.5, linecolor="white", rasterized=True,
            cbar_kws={"label": "% of ever-sig pathways"},
        )
        ax.set_title(
            f"% {pat}\n"
            f"baseline: G32A {baseline[0]:.1f}% / "
            f"R403C {baseline[1]:.1f}%",
            fontsize=11, fontweight="bold",
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(["G32A", "R403C"])



    plt.suptitle(
        "Per-database trajectory-pattern frequencies (ever-significant "
        "pathways, FDR<0.05 at any stage in either mutation)",
        fontsize=12, y=1.02)
    plt.tight_layout()

    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")


def print_manuscript_key_numbers(out_df, sig_df):
    """Print manuscript-key pathway counts and baselines to stdout."""
    print("\n=== Key numbers for manuscript ===")

    def _pct(db, pat, mut):
        sub = out_df[(out_df["database"] == db) & (out_df["mutation"] == mut)]
        if sub.empty:
            return None, None, None
        return (
            float(sub[f"pct_{pat}"].iloc[0]),
            int(sub[f"n_{pat}"].iloc[0]),
            int(sub["n_pathways_ever_sig"].iloc[0]),
        )

    for db in ("SynGO", "MitoCarta", "gocc"):
        for mut in ("G32A", "R403C"):
            parts = []
            for pat in PATTERNS:
                sr = _pct(db, pat, mut)
                if sr and sr[0] > 0:
                    parts.append(
                        f"{pat}: {sr[0]:5.2f}% ({sr[1]}/{sr[2]})")
            print(f"  {db:10s} {mut:6s}  {' | '.join(parts)}")

    print("\n  Universe-wide baselines:")
    for pat in PATTERNS:
        g32 = round(100 * (sig_df["Pattern_G32A"] == pat).sum() / len(sig_df), 2)
        r40 = round(100 * (sig_df["Pattern_R403C"] == pat).sum() / len(sig_df), 2)
        print(f"     {pat:25s}: "
              f"G32A {g32}% / R403C {r40}%")


def main():
    """Entry point: orchestrate pipeline steps."""
    sig_df, n_total_per_db = load_and_filter(IN_CSV, PATTERNS)
    out_df = build_table_df(sig_df, n_total_per_db, DB_ORDER,
                            MUTATIONS, PATTERNS, OUT_CSV)
    export_markdown(out_df, OUT_MD)
    build_heatmap(out_df, sig_df, OUT_PDF, OUT_PNG, DB_ORDER, PATTERNS)
    print_manuscript_key_numbers(out_df, sig_df)


if __name__ == "__main__":
    main()
