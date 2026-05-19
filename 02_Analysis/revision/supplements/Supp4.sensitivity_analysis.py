#!/usr/bin/env python3
"""
Sensitivity Analysis for Pattern Classification Thresholds

Tests robustness of pattern classification across reasonable threshold variations. Generates supplementary figure/table showing qualitative conclusions remain stable.

Thresholds tested:
- NES_EFFECT: 0.4, 0.5 (default), 0.6
- NES_STRONG: 0.8, 1.0 (default), 1.2
- IMPROVEMENT_RATIO: 0.6, 0.7 (default), 0.8
- WORSENING_RATIO: 1.25, 1.3 (default), 1.4
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict, List
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns

# Add module path for project imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / '01_Scripts'))
from Python.pattern_definitions import MEANINGFUL_PATTERNS

# Import unified color configuration
from Python.color_config import HEATMAP_ANNOTATION_COLORS, MUTATION_COLORS, create_diverging_cmap

# =============================================================================
# PATHS
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
DATA_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis"
OUTPUT_DIR = DATA_DIR / "Sensitivity_Analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# PARAMETERIZED CLASSIFICATION FUNCTION
# =============================================================================

def classify_pattern_parameterized(
    early_nes: float,
    early_padj: float,
    trajdev_nes: float,
    trajdev_padj: float,
    late_nes: float,
    late_padj: float,
    # Thresholds as parameters
    padj_significant: float = 0.05,
    padj_trending: float = 0.10,
    nes_effect: float = 0.5,
    nes_strong: float = 1.0,
    improvement_ratio: float = 0.7,
    worsening_ratio: float = 1.3
) -> Tuple[str, Optional[str]]:
    """
    Classify trajectory pattern with parameterized thresholds.

    This allows testing different threshold combinations for sensitivity analysis.
    """
    # Handle missing values
    if pd.isna([early_nes, trajdev_nes, late_nes, early_padj]).any():
        return ('Insufficient_data', None)

    early_abs = abs(early_nes)
    late_abs = abs(late_nes)
    trajdev_abs = abs(trajdev_nes)

    # Step 1: Early defect assessment
    early_sig_defect = (early_padj < padj_significant) and (early_abs > nes_effect)
    early_trending = (early_padj < padj_trending) and (early_abs > nes_effect)
    early_strong = (early_padj < padj_significant) and (early_abs > nes_strong)
    early_no_defect = (early_padj >= padj_trending) or (early_abs <= nes_effect)

    # Step 2: Late outcome assessment
    late_sig_defect = (late_padj < padj_significant) and (late_abs > nes_strong)
    late_resolved = late_abs < nes_effect

    # Improvement/worsening ratios
    if early_abs > 0.1:
        ratio = late_abs / early_abs
        improved = (ratio < improvement_ratio) or late_resolved
        worsened = ratio > worsening_ratio
    else:
        improved = late_resolved
        worsened = late_abs > nes_strong

    # Step 3: TrajDev assessment
    trajdev_sig = (trajdev_padj < padj_significant) and (trajdev_abs > nes_effect)

    # Direction assessment
    if early_abs > 0.1:
        trajdev_opposes = np.sign(trajdev_nes) != np.sign(early_nes)
        trajdev_amplifies = np.sign(trajdev_nes) == np.sign(early_nes)
    else:
        trajdev_opposes = False
        trajdev_amplifies = False

    # Step 4: Pattern classification
    # Late_onset: no early defect, significant late defect
    if early_no_defect and late_sig_defect:
        return ('Late_onset', 'High')

    # Patterns requiring early defect
    if early_sig_defect or early_trending:
        confidence = 'High' if early_sig_defect else 'Medium'

        # Active patterns - check BEFORE Transient
        if trajdev_sig and trajdev_opposes and improved:
            return ('Compensation', confidence)

        # Sign_reversal: TrajDev opposes AND sign flipped between Early and Late
        if trajdev_sig and trajdev_opposes:
            sign_flip = np.sign(early_nes) != np.sign(late_nes)
            late_substantial = late_abs > nes_effect
            if sign_flip and late_substantial:
                return ('Sign_reversal', confidence)

        if trajdev_sig and trajdev_amplifies and worsened:
            return ('Progressive', confidence)

        # Transient: strong early, fully resolved
        if early_strong and late_resolved:
            return ('Transient', 'High')

        # Passive patterns
        if not trajdev_sig:
            if improved:
                return ('Natural_improvement', 'High' if early_sig_defect else 'Medium')
            if worsened:
                return ('Natural_worsening', 'High' if early_sig_defect else 'Medium')

    # Edge case transient
    if early_strong and late_resolved:
        return ('Transient', 'High')

    return ('Complex', None)


# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================

def run_sensitivity_analysis(df: pd.DataFrame, mutations: List[str] = ['G32A', 'R403C']) -> pd.DataFrame:
    """
    Run classification with all threshold combinations and compare results.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format GSEA results with NES and p.adjust columns
    mutations : list
        Mutations to analyze

    Returns
    -------
    pd.DataFrame
        Summary of pattern distributions across threshold combinations
    """
    # Threshold combinations to test
    nes_effect_values = [0.4, 0.5, 0.6]
    nes_strong_values = [0.8, 1.0, 1.2]
    improvement_ratio_values = [0.6, 0.7, 0.8]
    worsening_ratio_values = [1.25, 1.3, 1.4]

    results = []

    # Column mappings
    col_map = {
        'G32A': {
            'early_nes': 'NES_G32A_vs_Ctrl_D35',
            'early_padj': 'p.adjust_G32A_vs_Ctrl_D35',
            'trajdev_nes': 'NES_Maturation_G32A_specific',
            'trajdev_padj': 'p.adjust_Maturation_G32A_specific',
            'late_nes': 'NES_G32A_vs_Ctrl_D65',
            'late_padj': 'p.adjust_G32A_vs_Ctrl_D65'
        },
        'R403C': {
            'early_nes': 'NES_R403C_vs_Ctrl_D35',
            'early_padj': 'p.adjust_R403C_vs_Ctrl_D35',
            'trajdev_nes': 'NES_Maturation_R403C_specific',
            'trajdev_padj': 'p.adjust_Maturation_R403C_specific',
            'late_nes': 'NES_R403C_vs_Ctrl_D65',
            'late_padj': 'p.adjust_R403C_vs_Ctrl_D65'
        }
    }

    # Generate all combinations
    combinations = list(product(
        nes_effect_values,
        nes_strong_values,
        improvement_ratio_values,
        worsening_ratio_values
    ))

    print(f"Testing {len(combinations)} threshold combinations...")

    for nes_effect, nes_strong, imp_ratio, wors_ratio in combinations:
        # Skip invalid combinations (nes_effect should be < nes_strong)
        if nes_effect >= nes_strong:
            continue

        combo_id = f"NES_eff={nes_effect}_strong={nes_strong}_imp={imp_ratio}_wors={wors_ratio}"

        for mutation in mutations:
            cols = col_map[mutation]

            # Classify each pathway
            patterns = []
            for _, row in df.iterrows():
                pattern, conf = classify_pattern_parameterized(
                    early_nes=row[cols['early_nes']],
                    early_padj=row[cols['early_padj']],
                    trajdev_nes=row[cols['trajdev_nes']],
                    trajdev_padj=row[cols['trajdev_padj']],
                    late_nes=row[cols['late_nes']],
                    late_padj=row[cols['late_padj']],
                    nes_effect=nes_effect,
                    nes_strong=nes_strong,
                    improvement_ratio=imp_ratio,
                    worsening_ratio=wors_ratio
                )
                patterns.append(pattern)

            # Count patterns
            pattern_counts = pd.Series(patterns).value_counts()
            total = len(patterns)

            # Calculate super-category percentages
            active_comp = pattern_counts.get('Compensation', 0)
            active_rev = pattern_counts.get('Sign_reversal', 0)
            active_prog = pattern_counts.get('Progressive', 0)
            passive = (pattern_counts.get('Natural_improvement', 0) +
                      pattern_counts.get('Natural_worsening', 0))
            late_onset = pattern_counts.get('Late_onset', 0)
            transient = pattern_counts.get('Transient', 0)
            complex_pat = pattern_counts.get('Complex', 0)
            insufficient = pattern_counts.get('Insufficient_data', 0)

            results.append({
                'combination_id': combo_id,
                'NES_EFFECT': nes_effect,
                'NES_STRONG': nes_strong,
                'IMPROVEMENT_RATIO': imp_ratio,
                'WORSENING_RATIO': wors_ratio,
                'mutation': mutation,
                'n_pathways': total,
                'Compensation': active_comp,
                'Compensation_pct': active_comp / total * 100,
                'Sign_reversal': active_rev,
                'Sign_reversal_pct': active_rev / total * 100,
                'Progressive': active_prog,
                'Progressive_pct': active_prog / total * 100,
                'Natural_improvement': pattern_counts.get('Natural_improvement', 0),
                'Natural_worsening': pattern_counts.get('Natural_worsening', 0),
                'Passive_total': passive,
                'Passive_pct': passive / total * 100,
                'Late_onset': late_onset,
                'Late_onset_pct': late_onset / total * 100,
                'Transient': transient,
                'Transient_pct': transient / total * 100,
                'Complex': complex_pat,
                'Complex_pct': complex_pat / total * 100,
                'Insufficient_data': insufficient,
                'is_default': (nes_effect == 0.5 and nes_strong == 1.0 and
                              imp_ratio == 0.7 and wors_ratio == 1.3)
            })

    return pd.DataFrame(results)


def analyze_key_claims(sensitivity_df: pd.DataFrame) -> pd.DataFrame:
    """
    Check whether key biological claims remain stable across thresholds.

    Key claims to verify:
    1. Compensation is dominant classifiable pattern (among non-Complex)
    2. R403C shows more compensation than G32A
    3. Progressive patterns are rare (<5%)
    4. Passive patterns exist but are secondary to Active compensation
    """
    claims = []

    for combo_id in sensitivity_df['combination_id'].unique():
        combo_data = sensitivity_df[sensitivity_df['combination_id'] == combo_id]
        g32a = combo_data[combo_data['mutation'] == 'G32A'].iloc[0]
        r403c = combo_data[combo_data['mutation'] == 'R403C'].iloc[0]

        # Calculate non-Complex totals for more meaningful "dominance"
        g32a_classifiable = g32a['n_pathways'] - g32a['Complex'] - g32a.get('Insufficient_data', 0)
        r403c_classifiable = r403c['n_pathways'] - r403c['Complex'] - r403c.get('Insufficient_data', 0)

        # Compensation fraction of classifiable pathways
        g32a_comp_of_classifiable = g32a['Compensation'] / g32a_classifiable * 100 if g32a_classifiable > 0 else 0
        r403c_comp_of_classifiable = r403c['Compensation'] / r403c_classifiable * 100 if r403c_classifiable > 0 else 0

        claims.append({
            'combination_id': combo_id,
            'NES_EFFECT': g32a['NES_EFFECT'],
            'NES_STRONG': g32a['NES_STRONG'],
            'IMPROVEMENT_RATIO': g32a['IMPROVEMENT_RATIO'],
            'WORSENING_RATIO': g32a['WORSENING_RATIO'],
            'is_default': g32a['is_default'],
            # Actual counts and percentages
            'comp_count_G32A': g32a['Compensation'],
            'comp_count_R403C': r403c['Compensation'],
            'comp_pct_G32A': g32a['Compensation_pct'],
            'comp_pct_R403C': r403c['Compensation_pct'],
            'classifiable_G32A': g32a_classifiable,
            'classifiable_R403C': r403c_classifiable,
            'comp_of_classifiable_G32A': g32a_comp_of_classifiable,
            'comp_of_classifiable_R403C': r403c_comp_of_classifiable,
            # Claim 1: Compensation is largest classifiable pattern (> 50% of non-Complex)
            'comp_dominates_classifiable_G32A': g32a_comp_of_classifiable > 50,
            'comp_dominates_classifiable_R403C': r403c_comp_of_classifiable > 50,
            # Claim 2: Which mutation shows more compensation
            'R403C_more_compensation': r403c['Compensation'] > g32a['Compensation'],
            'compensation_diff': r403c['Compensation'] - g32a['Compensation'],
            # Claim 3: Progressive is rare (<5% of total)
            'progressive_rare_G32A': g32a['Progressive_pct'] < 5,
            'progressive_rare_R403C': r403c['Progressive_pct'] < 5,
            'prog_pct_G32A': g32a['Progressive_pct'],
            'prog_pct_R403C': r403c['Progressive_pct'],
            # Claim 4: Compensation > Passive
            'comp_exceeds_passive_G32A': g32a['Compensation'] > g32a['Passive_total'],
            'comp_exceeds_passive_R403C': r403c['Compensation'] > r403c['Passive_total'],
        })

    return pd.DataFrame(claims)


def create_sensitivity_heatmap(sensitivity_df: pd.DataFrame, mutation: str, output_dir: Path):
    """
    Create heatmap showing pattern percentages across threshold combinations.
    """
    # Filter for this mutation
    df = sensitivity_df[sensitivity_df['mutation'] == mutation].copy()

    # Focus on NES_EFFECT vs IMPROVEMENT_RATIO (fixing others at default)
    df_subset = df[
        (df['NES_STRONG'] == 1.0) &
        (df['WORSENING_RATIO'] == 1.3)
    ].copy()

    if len(df_subset) == 0:
        print(f"No data for heatmap subset for {mutation}")
        return

    # Create pivot table for Compensation percentage
    pivot = df_subset.pivot_table(
        values='Compensation_pct',
        index='IMPROVEMENT_RATIO',
        columns='NES_EFFECT',
        aggfunc='mean'
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        pivot,
        annot=True,
        fmt='.1f',
        cmap='YlGnBu',
        ax=ax,
        cbar_kws={'label': 'Compensation %'}
    )
    ax.set_title(f'{mutation}: Compensation % by Threshold\n(NES_STRONG=1.0, WORSENING_RATIO=1.3)')
    ax.set_xlabel('NES_EFFECT threshold')
    ax.set_ylabel('IMPROVEMENT_RATIO threshold')

    # Mark default
    default_row = list(pivot.index).index(0.7)
    default_col = list(pivot.columns).index(0.5)
    ax.add_patch(plt.Rectangle(
        (default_col, default_row), 1, 1,
        fill=False, edgecolor='red', linewidth=3
    ))

    plt.tight_layout()
    plt.savefig(output_dir / f'sensitivity_heatmap_{mutation}.pdf', dpi=150)
    plt.savefig(output_dir / f'sensitivity_heatmap_{mutation}.png', dpi=150)
    plt.close()
    print(f"  Saved heatmap for {mutation}")


def create_summary_figure(sensitivity_df: pd.DataFrame, claims_df: pd.DataFrame, output_dir: Path):
    """
    Three-panel journal-style summary of threshold-robustness.

    A: pattern distribution at default thresholds, on the 12,221-pathway full
       background (% of all tested pathways).
    B: Compensation % across all 81 threshold combinations, full background.
    C: Compensation % of classifiable pathways across all 81 combinations on
       the 5,267 ever-significantly enriched universe (RESULTS denominator),
       with a dashed reference at 50 % (majority of classifiable).
    """
    import matplotlib as mpl

    # Okabe-Ito MUTATION_COLORS (matches pattern_summary_normalized.pdf)
    g32a_color = MUTATION_COLORS['G32A']    # #0072B2 blue
    r403c_color = MUTATION_COLORS['R403C']  # #D55E00 vermillion

    default_data = sensitivity_df[sensitivity_df['is_default']]
    g32a_default = default_data[default_data['mutation'] == 'G32A'].iloc[0]
    r403c_default = default_data[default_data['mutation'] == 'R403C'].iloc[0]
    n_combos = sensitivity_df['combination_id'].nunique()
    n_full = int(g32a_default['n_pathways'])  # 12,221

    # Load 5,267-universe sensitivity table (Comp_class_pct lives here)
    universe_5267_path = (
        PROJECT_ROOT / "03_Results" / "02_Analysis" / "Supplementary"
        / "sensitivity_5267universe.csv"
    )
    has_5267 = universe_5267_path.exists()
    df_5267 = pd.read_csv(universe_5267_path) if has_5267 else None
    n_5267 = 5267

    rc = {
        'font.family': 'sans-serif',
        'font.size': 9.5,
        'axes.titlesize': 10.5,
        'axes.labelsize': 9.5,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    }

    with mpl.rc_context(rc):
        fig, (ax_a, ax_b, ax_c) = plt.subplots(
            1, 3, figsize=(10.8, 3.5),
            gridspec_kw={'width_ratios': [0.85, 1.0, 1.05]}
        )

        # ---- Panel A: pattern distribution at default thresholds ----
        show_patterns = ['Compensation', 'Sign_reversal', 'Natural_improvement', 'Complex']
        pretty = ['Compensation', 'Sign reversal', 'Natural improvement', 'Complex']

        n_g32a = g32a_default['n_pathways']
        n_r403c = r403c_default['n_pathways']
        g32a_pct = [g32a_default[p] / n_g32a * 100 for p in show_patterns]
        r403c_pct = [r403c_default[p] / n_r403c * 100 for p in show_patterns]

        x = np.arange(len(show_patterns))
        w = 0.38
        ax_a.bar(x - w/2, g32a_pct, w, color=g32a_color, label='G32A', edgecolor='none')
        ax_a.bar(x + w/2, r403c_pct, w, color=r403c_color, label='R403C', edgecolor='none')
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(pretty, rotation=25, ha='right', rotation_mode='anchor')
        ax_a.set_ylabel('% of pathways')
        ax_a.set_ylim(0, max(max(g32a_pct), max(r403c_pct)) * 1.18)
        ax_a.legend(frameon=False, loc='upper right',
                    handlelength=1.0, handletextpad=0.5)
        ax_a.text(-0.22, 1.10, 'A', transform=ax_a.transAxes,
                  fontsize=12, fontweight='bold', va='top', ha='left')
        ax_a.text(0.5, 1.02,
                  f'default thresholds  ·  N = {n_full:,} (full background)',
                  transform=ax_a.transAxes, fontsize=7.8, color='#555555',
                  ha='center', va='bottom')

        # ---- Strip-plot helper (used for B and C) ----
        def strip(ax, values_by_mut, default_by_mut, xlabel, denom_label,
                  reference_x=None, x_pad_right=3.5):
            y_positions = {'G32A': 0, 'R403C': 1}
            rng = np.random.RandomState(0)
            all_vals = []
            for mut, color in [('G32A', g32a_color), ('R403C', r403c_color)]:
                vals = np.asarray(values_by_mut[mut], dtype=float)
                all_vals.extend(vals.tolist())
                default_val = float(default_by_mut[mut])
                y = y_positions[mut]
                jitter = rng.uniform(-0.32, 0.32, size=len(vals))
                ax.scatter(vals, np.full_like(vals, y, dtype=float) + jitter,
                           s=14, c=color, alpha=0.35, edgecolors='none', zorder=3)
                ax.hlines(y, vals.min(), vals.max(),
                          color=color, linewidth=1.2, zorder=2, alpha=0.9)
                ax.scatter([default_val], [y], marker='|', s=220,
                           c='black', linewidth=1.8, zorder=5)
                ax.text(vals.max() + 0.25, y,
                        f'{vals.min():.1f}–{vals.max():.1f}%',
                        va='center', fontsize=8.5, color='#333333')
            if reference_x is not None:
                ax.axvline(reference_x, color='#888888', linewidth=0.8,
                           linestyle='--', zorder=1)
            ax.set_yticks([0, 1])
            ax.set_yticklabels(['G32A', 'R403C'])
            ax.invert_yaxis()
            ax.set_ylim(1.9, -0.9)
            ax.set_xlabel(xlabel)
            x_lo = min(all_vals) - 0.4
            x_hi = max(all_vals) + x_pad_right
            ax.set_xlim(x_lo, x_hi)
            ax.text(0.5, 1.02, denom_label,
                    transform=ax.transAxes, fontsize=7.8, color='#555555',
                    ha='center', va='bottom')

        # ---- Panel B: Compensation % on full background (12,221) ----
        values_B = {
            mut: sensitivity_df.loc[sensitivity_df['mutation'] == mut,
                                    'Compensation_pct'].values
            for mut in ['G32A', 'R403C']
        }
        defaults_B = {
            'G32A': g32a_default['Compensation_pct'],
            'R403C': r403c_default['Compensation_pct'],
        }
        strip(ax_b, values_B, defaults_B,
              xlabel=f'Compensation %  ·  {n_combos} thresholds',
              denom_label=f'n / N  ·  N = {n_full:,} (full background)')
        ax_b.text(-0.20, 1.10, 'B', transform=ax_b.transAxes,
                  fontsize=12, fontweight='bold', va='top', ha='left')
        ax_b.text(0.98, 0.06, 'tick = default thresholds',
                  transform=ax_b.transAxes, fontsize=7.5, color='#777777',
                  ha='right', va='bottom')

        # ---- Panel C: Compensation % of classifiable on 5,267 universe ----
        if has_5267:
            values_C = {
                mut: df_5267.loc[df_5267['mut'] == mut, 'Comp_class_pct'].values
                for mut in ['G32A', 'R403C']
            }
            def_5267 = df_5267[(df_5267['ne'] == 0.5) & (df_5267['ns'] == 1.0) &
                               (df_5267['ir'] == 0.7) & (df_5267['wr'] == 1.3)]
            defaults_C = {
                'G32A': def_5267[def_5267['mut'] == 'G32A']['Comp_class_pct'].values[0],
                'R403C': def_5267[def_5267['mut'] == 'R403C']['Comp_class_pct'].values[0],
            }
            strip(ax_c, values_C, defaults_C,
                  xlabel=f'Compensation % of classifiable  ·  {n_combos} thresholds',
                  denom_label=f'n / classifiable  ·  ever-sig. universe = {n_5267:,}',
                  reference_x=50, x_pad_right=5.5)
            ax_c.text(-0.20, 1.10, 'C', transform=ax_c.transAxes,
                      fontsize=12, fontweight='bold', va='top', ha='left')
            # Place "50 %" label just above the dashed reference line
            ax_c.text(50, -0.78, '50 %',
                      fontsize=7.5, color='#888888', va='bottom', ha='center')
        else:
            ax_c.text(0.5, 0.5, '5,267-universe table not found',
                      transform=ax_c.transAxes, ha='center', va='center',
                      color='#999999')
            ax_c.axis('off')

        plt.tight_layout()
        plt.savefig(output_dir / 'sensitivity_analysis_summary.pdf', bbox_inches='tight')
        plt.savefig(output_dir / 'sensitivity_analysis_summary.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved summary figure (3 panels)")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("SENSITIVITY ANALYSIS FOR PATTERN CLASSIFICATION")
    print("=" * 70)

    # Load data
    print("\n1. Loading GSEA results (wide format)...")
    wide_file = DATA_DIR / "Python_exports" / "gsea_results_wide.csv"
    df = pd.read_csv(wide_file)
    print(f"   Loaded {len(df)} unique pathways")

    # Run sensitivity analysis
    print("\n2. Running sensitivity analysis...")
    sensitivity_df = run_sensitivity_analysis(df)
    print(f"   Generated {len(sensitivity_df)} classification results")

    # Analyze key claims
    print("\n3. Analyzing claim stability...")
    claims_df = analyze_key_claims(sensitivity_df)

    # Save results
    print("\n4. Saving results...")
    sensitivity_df.to_csv(OUTPUT_DIR / 'sensitivity_results.csv', index=False)
    claims_df.to_csv(OUTPUT_DIR / 'claim_stability.csv', index=False)
    print(f"   Saved to {OUTPUT_DIR}")

    # Create visualizations
    print("\n5. Creating visualizations...")
    for mutation in ['G32A', 'R403C']:
        create_sensitivity_heatmap(sensitivity_df, mutation, OUTPUT_DIR)
    create_summary_figure(sensitivity_df, claims_df, OUTPUT_DIR)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    default_data = sensitivity_df[sensitivity_df['is_default']]
    print("\nDefault threshold results:")
    for mutation in ['G32A', 'R403C']:
        row = default_data[default_data['mutation'] == mutation].iloc[0]
        print(f"\n  {mutation}:")
        print(f"    Compensation: {row['Compensation']:.0f} ({row['Compensation_pct']:.1f}%)")
        print(f"    Progressive:  {row['Progressive']:.0f} ({row['Progressive_pct']:.1f}%)")
        print(f"    Passive:      {row['Passive_total']:.0f} ({row['Passive_pct']:.1f}%)")
        print(f"    Late_onset:   {row['Late_onset']:.0f} ({row['Late_onset_pct']:.1f}%)")
        print(f"    Transient:    {row['Transient']:.0f} ({row['Transient_pct']:.1f}%)")
        print(f"    Complex:      {row['Complex']:.0f} ({row['Complex_pct']:.1f}%)")

    print("\nClaim stability across threshold variations:")
    n_combos = len(claims_df)

    # Key claims
    claims_summary = {
        'R403C shows more compensation than G32A': claims_df['R403C_more_compensation'].sum(),
        'Progressive patterns rare for G32A (<5%)': claims_df['progressive_rare_G32A'].sum(),
        'Progressive patterns rare for R403C (<5%)': claims_df['progressive_rare_R403C'].sum(),
        'Compensation > Passive for G32A': claims_df['comp_exceeds_passive_G32A'].sum(),
        'Compensation > Passive for R403C': claims_df['comp_exceeds_passive_R403C'].sum(),
    }

    for claim, count in claims_summary.items():
        pct = count / n_combos * 100
        status = "✓ STABLE" if pct == 100 else ("⚠ MOSTLY STABLE" if pct > 80 else "✗ UNSTABLE")
        print(f"  {status} ({pct:.0f}%): {claim}")

    # Compensation percentage range
    print("\n  Compensation percentage range across thresholds:")
    for mutation in ['G32A', 'R403C']:
        mut_data = sensitivity_df[sensitivity_df['mutation'] == mutation]
        min_pct = mut_data['Compensation_pct'].min()
        max_pct = mut_data['Compensation_pct'].max()
        median_pct = mut_data['Compensation_pct'].median()
        print(f"    {mutation}: {min_pct:.1f}% - {max_pct:.1f}% (median: {median_pct:.1f}%)")

    # Compensation as fraction of classifiable pathways
    print("\n  Compensation as % of classifiable (non-Complex) pathways:")
    for mutation in ['G32A', 'R403C']:
        mut_claims = claims_df[[f'comp_of_classifiable_{mutation}']].values.flatten()
        print(f"    {mutation}: {mut_claims.min():.1f}% - {mut_claims.max():.1f}% (median: {np.median(mut_claims):.1f}%)")

    print("\n" + "=" * 70)
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
