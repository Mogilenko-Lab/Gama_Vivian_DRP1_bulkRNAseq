#!/usr/bin/env python3
"""
GOChord-style Chord Diagrams: Gene-Pathway Leading Edge Membership

Generates THREE chord-diagram variants per mutation (G32A, R403C):

  1. focused        - SynGO synaptic-localized pathways (original scope; refreshed
                      with D35-outer / D65-inner ring layout and journal-grade
                      typography)
  2. expanded_clean - Adds GO:CC cytosolic ribosome to the SynGO set, so the
                      synaptic-ribosome arcs visibly nest inside the broader
                      cytoplasmic structural ribosome. Pathway arcs are grouped
                      under two outer brackets: 'SynGO synaptic-localized' and
                      'Broader cytoplasmic ribosome'.
  3. expanded_full  - Adds GO:CC cytosolic ribosome, KEGG translation initiation,
                      REACTOME eukaryotic translation elongation, and (counter-
                      direction) GO:BP ribosome biogenesis. Three outer brackets.

All diagrams now use:
  * D35 fold change on the OUTER gene-ring, D65 on the INNER (concentric
    left-to-right temporal reading)
  * NES annotations (D35 / D65) printed under each pathway label
  * 600 dpi, large canvas, journal-print-grade font sizing
  * Adaptive gene-arc span and font scaling to keep labels overlap-free

Usage:
    python3 02_Analysis/3.7.viz_chord_diagrams.py

Output:
    03_Results/02_Analysis/Plots/Chord_Diagrams/
        chord_diagram_<MUT>.{pdf,png}                  (focused)
        chord_diagram_<MUT>_expanded_clean.{pdf,png}
        chord_diagram_<MUT>_expanded_full.{pdf,png}
        chord_{genes,pathways,connections}_<MUT>_<variant>.csv
        gsea_pathway_data.csv
"""

import sys
import subprocess
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / '01_Scripts'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Wedge, Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings('ignore')

from Python.color_config import create_diverging_cmap, MUTATION_COLORS

# =============================================================================
# CONFIGURATION
# =============================================================================

# Per-variant pathway definitions.
# Each entry: (display_name, source_category, group_label)
#   source_category: which GSEA result list/category to pull from.
#     'syngo' -> syngo_gsea_results.rds
#     'all:<cat>' -> all_gsea_results.rds[[contrast]][[<cat>]]
PATHWAY_DEFS = {
    'SYNGO:presyn_ribosome':                          ('Presyn. Ribosome',           'syngo',          'SynGO synaptic-localized'),
    'SYNGO:postsyn_ribosome':                         ('Postsyn. Ribosome',          'syngo',          'SynGO synaptic-localized'),
    'GO:0099523':                                     ('Presyn. Cytosol',            'syngo',          'SynGO synaptic-localized'),
    'GO:0099524':                                     ('Postsyn. Cytosol',           'syngo',          'SynGO synaptic-localized'),
    'GO:0014069':                                     ('Postsyn. Density',           'syngo',          'SynGO synaptic-localized'),
    'GO:0045211':                                     ('Postsyn. Membrane',          'syngo',          'SynGO synaptic-localized'),
    'GOCC_CYTOSOLIC_RIBOSOME':                        ('GO:CC Cytosolic Ribosome',   'all:gocc',       'Broader cytoplasmic ribosome / translation'),
    'GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT':         ('GO:CC Cyto. Large Subunit',  'all:gocc',       'Broader cytoplasmic ribosome / translation'),
    'GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT':         ('GO:CC Cyto. Small Subunit',  'all:gocc',       'Broader cytoplasmic ribosome / translation'),
    'WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS':              ('WP Cyto. Ribosomal Proteins','all:wiki',       'Broader cytoplasmic ribosome / translation'),
    'KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION':  ('KEGG Translation Initiation','all:kegg',       'Broader cytoplasmic ribosome / translation'),
    'REACTOME_EUKARYOTIC_TRANSLATION_INITIATION':     ('REACTOME Transl. Initiation','all:reactome',   'Broader cytoplasmic ribosome / translation'),
    'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION':     ('REACTOME Transl. Elongation','all:reactome',   'Broader cytoplasmic ribosome / translation'),
    'GOBP_RIBOSOME_BIOGENESIS':                       ('GO:BP Ribosome Biogenesis ↑','all:gobp',       'Biogenesis (counter-direction)'),
}

VARIANTS = {
    'focused': [
        'SYNGO:presyn_ribosome', 'SYNGO:postsyn_ribosome',
        'GO:0099523', 'GO:0099524', 'GO:0014069', 'GO:0045211',
    ],
    'expanded_clean': [
        'SYNGO:presyn_ribosome', 'SYNGO:postsyn_ribosome',
        'GO:0099523', 'GO:0099524', 'GO:0014069', 'GO:0045211',
        'GOCC_CYTOSOLIC_RIBOSOME',
    ],
    'expanded_full': [
        'SYNGO:presyn_ribosome', 'SYNGO:postsyn_ribosome',
        'GO:0099523', 'GO:0099524', 'GO:0014069', 'GO:0045211',
        'GOCC_CYTOSOLIC_RIBOSOME',
        'GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT',
        'GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT',
        'WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS',
        'KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION',
        'REACTOME_EUKARYOTIC_TRANSLATION_INITIATION',
        'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION',
        'GOBP_RIBOSOME_BIOGENESIS',
    ],
}

# =============================================================================
# Three-tier colorblind-safe palette (Wong 2011 + Paul Tol vibrant/muted)
#
# Visual hierarchy encodes the manuscript message: synaptic ribosomes (focal,
# saturated cool jewel tones) sit nested INSIDE a warm cytoplasmic ribosome /
# translation envelope (low-alpha gold→vermilion); other synaptic compartments
# are present but recessive (muted cool); ribosome biogenesis is a counter-
# direction reference (desaturated steel blue).
#
# Tiers control three things at once: pathway-arc color, arc alpha, ribbon
# alpha, and z-order. Envelope and counter are drawn FIRST (under), recessive
# synaptic compartments NEXT, focal synaptic ribosomes LAST (on top).
# =============================================================================

TIER_ENVELOPE  = 'envelope'
TIER_RECESSIVE = 'recessive'
TIER_FOCAL     = 'focal'
TIER_COUNTER   = 'counter'

# Per-pathway: (hex color, tier)
PATHWAY_STYLE = {
    # Focal pair - mid-value cool jewel tones (distinguished by hue, not darkness,
    # so they do not overpower the warm envelope by sheer value contrast).
    'SYNGO:presyn_ribosome':                          ('#CC6677', TIER_FOCAL),     # Muted rose (Tol)
    'SYNGO:postsyn_ribosome':                         ('#AA4499', TIER_FOCAL),     # Muted purple (Tol)
    # Other synaptic compartments - recessive (muted cool)
    'GO:0099523':                                     ('#117733', TIER_RECESSIVE), # Forest green
    'GO:0099524':                                     ('#44AA99', TIER_RECESSIVE), # Teal
    'GO:0014069':                                     ('#999933', TIER_RECESSIVE), # Olive
    'GO:0045211':                                     ('#88CCEE', TIER_RECESSIVE), # Sky blue
    # Broader cytoplasmic ribosome / translation envelope - rich warm earth gradient.
    # Values bumped up so the envelope reads as a substantial container, not a wash.
    'GOCC_CYTOSOLIC_RIBOSOME':                        ('#B58029', TIER_ENVELOPE),  # Rich ochre
    'GOCC_CYTOSOLIC_LARGE_RIBOSOMAL_SUBUNIT':         ('#8B6914', TIER_ENVELOPE),  # Dark gold
    'GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT':         ('#C9A227', TIER_ENVELOPE),  # Mustard
    'WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS':              ('#A0522D', TIER_ENVELOPE),  # Sienna
    'KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION':  ('#E69F00', TIER_ENVELOPE),  # Orange (Wong)
    'REACTOME_EUKARYOTIC_TRANSLATION_INITIATION':     ('#CC5500', TIER_ENVELOPE),  # Burnt orange
    'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION':     ('#D55E00', TIER_ENVELOPE),  # Vermilion (Wong)
    # Counter-direction (biogenesis) - desaturated steel
    'GOBP_RIBOSOME_BIOGENESIS':                       ('#6699CC', TIER_COUNTER),
}

# Per-tier rendering parameters
# Envelope alpha raised so warm envelope is visually substantial; focal alpha
# softened slightly so cool jewel tones sit ON the envelope, not crush it.
TIER_PARAMS = {
    TIER_ENVELOPE:  {'arc_alpha': 0.78, 'ribbon_alpha': 0.42, 'zorder': 1},
    TIER_COUNTER:   {'arc_alpha': 0.65, 'ribbon_alpha': 0.28, 'zorder': 1},
    TIER_RECESSIVE: {'arc_alpha': 0.55, 'ribbon_alpha': 0.32, 'zorder': 2},
    TIER_FOCAL:     {'arc_alpha': 0.80, 'ribbon_alpha': 0.62, 'zorder': 3},
}

PATHWAY_COLOR_MAP = {pid: c for pid, (c, _t) in PATHWAY_STYLE.items()}

GROUP_COLORS = {
    'SynGO synaptic-localized':                   '#555555',
    'Broader cytoplasmic ribosome / translation': '#8B4500',
    'Biogenesis (counter-direction)':             '#2E5D8B',
}

# Short labels used on the on-figure bracket arc (full label lives in legend)
GROUP_SHORT_LABELS = {
    'SynGO synaptic-localized':                   'SynGO',
    'Broader cytoplasmic ribosome / translation': 'Broader cyto.',
    'Biogenesis (counter-direction)':             'Biogenesis',
}

OUTPUT_DIR = project_root / '03_Results/02_Analysis/Plots/Chord_Diagrams'

# =============================================================================
# DATA LOADING
# =============================================================================

def load_gsea_data():
    """Pull NES + core_enrichment for every pathway across all variants/contrasts.

    Returns a long-form DataFrame with columns:
        ID, NES, p.adjust, core_enrichment, contrast
    """
    all_ids = sorted(PATHWAY_DEFS.keys())

    # Split by source
    syngo_ids = [pid for pid in all_ids if PATHWAY_DEFS[pid][1] == 'syngo']
    all_cat_ids = {}  # category -> [ids]
    for pid in all_ids:
        src = PATHWAY_DEFS[pid][1]
        if src.startswith('all:'):
            cat = src.split(':', 1)[1]
            all_cat_ids.setdefault(cat, []).append(pid)

    def r_vec(ids):
        return 'c("' + '","'.join(ids) + '")' if ids else 'character(0)'

    cat_blocks = []
    for cat, ids in all_cat_ids.items():
        cat_blocks.append(f'''
        ids_{cat} <- {r_vec(ids)}
        for (ct in contrasts) {{
            obj <- tryCatch(all_results[[ct]][["{cat}"]], error=function(e) NULL)
            if (is.null(obj)) next
            df <- obj@result
            df <- df[df$ID %in% ids_{cat}, c("ID","NES","pvalue","p.adjust","core_enrichment"), drop=FALSE]
            if (nrow(df) > 0) {{ df$contrast <- ct; out_list[[length(out_list)+1]] <- df }}
        }}
        ''')

    r_script = f'''
    syngo_results <- readRDS("03_Results/02_Analysis/checkpoints/syngo_gsea_results.rds")
    all_results <- readRDS("03_Results/02_Analysis/checkpoints/all_gsea_results.rds")

    contrasts <- c("G32A_vs_Ctrl_D35","G32A_vs_Ctrl_D65",
                   "R403C_vs_Ctrl_D35","R403C_vs_Ctrl_D65")

    syngo_ids <- {r_vec(syngo_ids)}
    out_list <- list()

    for (ct in contrasts) {{
        obj <- tryCatch(syngo_results[[ct]], error=function(e) NULL)
        if (is.null(obj)) next
        df <- obj@result
        df <- df[df$ID %in% syngo_ids, c("ID","NES","pvalue","p.adjust","core_enrichment"), drop=FALSE]
        if (nrow(df) > 0) {{ df$contrast <- ct; out_list[[length(out_list)+1]] <- df }}
    }}

    {''.join(cat_blocks)}

    combined <- do.call(rbind, out_list)
    write.csv(combined, "{OUTPUT_DIR}/gsea_pathway_data.csv", row.names=FALSE)
    cat("Rows: ", nrow(combined), "\\n")
    '''

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(['Rscript', '-e', r_script],
                            cwd=str(project_root), capture_output=True, text=True)
    if result.returncode != 0:
        print(f"R script error:\n{result.stderr}")
        raise RuntimeError("Failed to extract GSEA data")

    return pd.read_csv(OUTPUT_DIR / 'gsea_pathway_data.csv')


def load_gene_set_memberships():
    """Pull FULL gene-set membership (annotation, not GSEA leading edge) for
    every pathway in PATHWAY_DEFS. Returns DataFrame with columns: pathway_id,
    gene_symbol.

    SynGO ids come from syngo_lists.rds$T2G; MSigDB ids (GO:CC, GO:BP, KEGG,
    REACTOME) come from msigdbr (the same source the GSEA was run against).
    """
    cache = OUTPUT_DIR / 'pathway_gene_set_memberships.csv'

    all_ids = sorted(PATHWAY_DEFS.keys())
    syngo_ids   = [p for p in all_ids if PATHWAY_DEFS[p][1] == 'syngo']
    msigdb_cats = {}  # category -> list of pathway IDs
    for p in all_ids:
        src = PATHWAY_DEFS[p][1]
        if src.startswith('all:'):
            msigdb_cats.setdefault(src.split(':', 1)[1], []).append(p)

    cat_to_msigdb = {
        'gocc':     ('C5', 'GO:CC'),
        'gobp':     ('C5', 'GO:BP'),
        'kegg':     ('C2', 'CP:KEGG_MEDICUS'),
        'reactome': ('C2', 'CP:REACTOME'),
        'wiki':     ('C2', 'CP:WIKIPATHWAYS'),
    }

    def r_vec(ids):
        return 'c("' + '","'.join(ids) + '")' if ids else 'character(0)'

    msigdb_blocks = []
    for cat, ids in msigdb_cats.items():
        gs_cat, gs_sub = cat_to_msigdb[cat]
        msigdb_blocks.append(f'''
        m <- tryCatch(msigdbr(species="Homo sapiens", category="{gs_cat}", subcategory="{gs_sub}"),
                      error=function(e) NULL)
        if (!is.null(m)) {{
            ids <- {r_vec(ids)}
            sub <- m[m$gs_name %in% ids, c("gs_name","gene_symbol")]
            colnames(sub) <- c("pathway_id","gene_symbol")
            out_list[[length(out_list)+1]] <- sub
        }}
        ''')

    r_script = f'''
    suppressPackageStartupMessages({{
        library(msigdbr); library(dplyr)
    }})

    syngo <- readRDS("03_Results/02_Analysis/checkpoints/syngo_lists.rds")$T2G
    syngo_ids <- {r_vec(syngo_ids)}
    syngo_sub <- syngo[syngo$gs_name %in% syngo_ids, c("gs_name","gene_symbol")]
    colnames(syngo_sub) <- c("pathway_id","gene_symbol")

    out_list <- list(syngo_sub)
    {''.join(msigdb_blocks)}

    combined <- do.call(rbind, lapply(out_list, as.data.frame))
    write.csv(combined, "{cache}", row.names=FALSE)
    cat("Rows:", nrow(combined), "\\n")
    '''

    result = subprocess.run(['Rscript', '-e', r_script],
                            cwd=str(project_root), capture_output=True, text=True)
    if result.returncode != 0:
        print(f"R script error:\n{result.stderr}")
        raise RuntimeError("Failed to extract gene-set memberships")
    print(f"    {result.stdout.strip()}")
    return pd.read_csv(cache)


def load_de_results():
    de_dir = project_root / '03_Results/02_Analysis/DE_results'
    contrasts = ['G32A_vs_Ctrl_D35', 'G32A_vs_Ctrl_D65',
                 'R403C_vs_Ctrl_D35', 'R403C_vs_Ctrl_D65']
    return {c: pd.read_csv(de_dir / f'{c}_DE_results.csv', index_col=0)
            for c in contrasts if (de_dir / f'{c}_DE_results.csv').exists()}


# =============================================================================
# DATA PREP
# =============================================================================

def prepare_chord_data(gsea_df, de_data, mutation, pathway_ids, max_genes=60,
                       gs_membership_df=None):
    d35_ct = f'{mutation}_vs_Ctrl_D35'
    d65_ct = f'{mutation}_vs_Ctrl_D65'
    mdf = gsea_df[gsea_df['contrast'].isin([d35_ct, d65_ct])].copy()

    pathways = []
    gene_pathway_map = {}
    all_genes = set()

    for pid in pathway_ids:
        rows = mdf[mdf['ID'] == pid]
        if rows.empty:
            print(f"   WARNING: no rows for {pid}")
            continue

        def grab(ct, col):
            v = rows[rows['contrast'] == ct][col].values
            return v[0] if len(v) else np.nan

        nes_d35 = grab(d35_ct, 'NES')
        nes_d65 = grab(d65_ct, 'NES')
        padj_d35 = grab(d35_ct, 'p.adjust')
        padj_d65 = grab(d65_ct, 'p.adjust')

        # Prefer D35 core_enrichment, fall back to D65
        ce = grab(d35_ct, 'core_enrichment')
        if not isinstance(ce, str) or pd.isna(ce):
            ce = grab(d65_ct, 'core_enrichment')
        if not isinstance(ce, str):
            continue

        core = ce.split('/')
        display_name, _src, group = PATHWAY_DEFS[pid]
        pathways.append({
            'id': pid,
            'display_name': display_name,
            'group': group,
            'nes_d35': nes_d35, 'nes_d65': nes_d65,
            'padj_d35': padj_d35, 'padj_d65': padj_d65,
            'core_genes': core,
            'n_genes': len(core),
        })
        all_genes.update(core)
        for g in core:
            gene_pathway_map.setdefault(g, []).append(pid)

    # Rank genes: ribosomal first, then by # pathways
    def gene_key(g):
        ribo = g.startswith(('RPL', 'RPS', 'MRPL', 'MRPS'))
        return (0 if ribo else 1, -len(gene_pathway_map[g]), g)
    sorted_genes = sorted(all_genes, key=gene_key)
    if len(sorted_genes) > max_genes:
        print(f"   Filtering {len(sorted_genes)} -> {max_genes} genes")
    selected = set(sorted_genes[:max_genes])

    connections = [(g, pid) for g in selected for pid in gene_pathway_map[g]]

    # Build gene-set-only connections: gene is annotated to the pathway but
    # was NOT in the leading edge for this contrast. These render as faint
    # background ribbons so the diagram shows true gene-set membership.
    connections_gs_only = []
    if gs_membership_df is not None:
        gs_pathway_to_genes = {}
        for pid in pathway_ids:
            members = set(gs_membership_df.loc[
                gs_membership_df['pathway_id'] == pid, 'gene_symbol'])
            gs_pathway_to_genes[pid] = members

        le_set = set(connections)
        for pid in pathway_ids:
            members = gs_pathway_to_genes.get(pid, set())
            for g in selected:
                if g in members and (g, pid) not in le_set:
                    connections_gs_only.append((g, pid))

    de_d35 = de_data.get(d35_ct)
    de_d65 = de_data.get(d65_ct)
    genes = []
    for g in selected:
        def safe(df, gene, col, default):
            return df.loc[gene, col] if df is not None and gene in df.index else default
        genes.append({
            'name': g,
            'logfc_d35': safe(de_d35, g, 'logFC', np.nan),
            'logfc_d65': safe(de_d65, g, 'logFC', np.nan),
            'padj_d35':  safe(de_d35, g, 'adj.P.Val', 1.0),
            'padj_d65':  safe(de_d65, g, 'adj.P.Val', 1.0),
            'n_pathways': len(gene_pathway_map[g]),
            'is_ribosomal': g.startswith(('RPL', 'RPS', 'MRPL', 'MRPS')),
        })
    genes.sort(key=lambda x: (0 if x['is_ribosomal'] else 1, -x['n_pathways'], x['name']))

    return {
        'mutation': mutation,
        'pathways': pathways,
        'genes': genes,
        'connections': connections,
        'connections_gs_only': connections_gs_only,
        'total_genes': len(all_genes),
        'filtered_genes': len(selected),
    }


# =============================================================================
# DRAWING
# =============================================================================

def draw_ribbon(ax, gene_angle, gene_width, p_start, p_end, color, alpha,
                gene_radius, pathway_radius, zorder=2):
    n_pts = 36
    g_a = np.radians(gene_angle - gene_width / 2)
    g_b = np.radians(gene_angle + gene_width / 2)
    p_a = np.radians(p_start)
    p_b = np.radians(p_end)
    ctrl_r = 0.25

    p_theta = np.linspace(p_a, p_b, n_pts)
    p_x = pathway_radius * np.cos(p_theta)
    p_y = pathway_radius * np.sin(p_theta)

    verts, codes = [], []
    # Start: gene left
    verts.append((gene_radius * np.cos(g_a), gene_radius * np.sin(g_a)))
    codes.append(MplPath.MOVETO)
    # Cubic Bezier to pathway start
    verts += [
        (ctrl_r * np.cos(g_a), ctrl_r * np.sin(g_a)),
        (ctrl_r * np.cos(p_a), ctrl_r * np.sin(p_a)),
        (p_x[0], p_y[0]),
    ]
    codes += [MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    # Arc along pathway
    for i in range(1, len(p_x)):
        verts.append((p_x[i], p_y[i]))
        codes.append(MplPath.LINETO)
    # Cubic Bezier back to gene right
    verts += [
        (ctrl_r * np.cos(p_b), ctrl_r * np.sin(p_b)),
        (ctrl_r * np.cos(g_b), ctrl_r * np.sin(g_b)),
        (gene_radius * np.cos(g_b), gene_radius * np.sin(g_b)),
    ]
    codes += [MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    # Close
    verts.append((gene_radius * np.cos(g_a), gene_radius * np.sin(g_a)))
    codes.append(MplPath.CLOSEPOLY)

    ax.add_patch(PathPatch(MplPath(verts, codes),
                           facecolor=color, edgecolor=color,
                           alpha=alpha, linewidth=0.25, zorder=zorder))


def draw_gene_rings(ax, gene_angle, gene_width, logfc_d35, logfc_d65, cmap, norm,
                    inner_r=0.86, ring_h=0.07):
    """D35 on OUTER ring (r2), D65 on INNER ring (r1).

    Concentric left-to-right reading: as the eye moves from circle interior
    outward to the gene label, the gene rings progress D65 -> D35 (temporal
    early). This matches the user's request to put D35 on the outside.
    """
    r_inner_inner = inner_r              # 0.86
    r_inner_outer = inner_r + ring_h     # 0.93
    r_outer_inner = inner_r + ring_h     # 0.93
    r_outer_outer = inner_r + 2 * ring_h # 1.00

    color_d65 = '#d3d3d3' if pd.isna(logfc_d65) else cmap(norm(logfc_d65))
    color_d35 = '#d3d3d3' if pd.isna(logfc_d35) else cmap(norm(logfc_d35))

    t1 = gene_angle - gene_width / 2
    t2 = gene_angle + gene_width / 2

    # Inner ring = D65
    ax.add_patch(Wedge((0, 0), r_inner_outer, t1, t2,
                       width=ring_h, facecolor=color_d65,
                       edgecolor='white', linewidth=0.25))
    # Outer ring = D35
    ax.add_patch(Wedge((0, 0), r_outer_outer, t1, t2,
                       width=ring_h, facecolor=color_d35,
                       edgecolor='white', linewidth=0.25))


def draw_pathway_arc(ax, t_start, t_end, color, label,
                     fontsize_label, arc_alpha=0.70,
                     inner_r=0.86, outer_r=1.0, label_r_offset=0.07,
                     zorder=4):
    ax.add_patch(Wedge((0, 0), outer_r, t_start, t_end,
                       width=outer_r - inner_r,
                       facecolor=color, edgecolor='white',
                       linewidth=1.5, alpha=arc_alpha, zorder=zorder))

    mid = (t_start + t_end) / 2
    rad = np.radians(mid)

    if -90 <= mid <= 90:
        rotation = mid
        ha_outer = 'left'
    else:
        rotation = mid + 180
        ha_outer = 'right'

    label_r = outer_r + label_r_offset
    lx, ly = label_r * np.cos(rad), label_r * np.sin(rad)
    ax.text(lx, ly, label, fontsize=fontsize_label, ha=ha_outer, va='center',
            rotation=rotation, rotation_mode='anchor', fontweight='bold',
            color='#222222', zorder=10)


def draw_group_bracket(ax, t_start, t_end, label, color,
                       bracket_r=1.40, label_r_offset=0.10, fontsize=11):
    """Thin outer arc that brackets a group of pathway arcs, plus a label."""
    if t_end - t_start < 0.2:
        return
    arc_theta = np.linspace(np.radians(t_start), np.radians(t_end), 60)
    ax.plot(bracket_r * np.cos(arc_theta), bracket_r * np.sin(arc_theta),
            color=color, linewidth=2.2, solid_capstyle='round', zorder=5)
    # Small inward ticks at the ends
    for t_deg in (t_start, t_end):
        t_rad = np.radians(t_deg)
        ax.plot([bracket_r * np.cos(t_rad), (bracket_r - 0.04) * np.cos(t_rad)],
                [bracket_r * np.sin(t_rad), (bracket_r - 0.04) * np.sin(t_rad)],
                color=color, linewidth=2.2, solid_capstyle='round', zorder=5)

    mid = (t_start + t_end) / 2
    rad = np.radians(mid)
    lr = bracket_r + label_r_offset
    lx, ly = lr * np.cos(rad), lr * np.sin(rad)
    if -90 <= mid <= 90:
        rotation = mid
        ha = 'left'
    else:
        rotation = mid + 180
        ha = 'right'
    ax.text(lx, ly, label, fontsize=fontsize, ha=ha, va='center',
            rotation=rotation, rotation_mode='anchor',
            fontweight='bold', color=color)


def draw_chord_diagram(data, output_path, variant, figsize=(15, 11)):
    pathways = data['pathways']
    genes = data['genes']
    connections = data['connections']
    connections_gs_only = data.get('connections_gs_only', [])
    mutation = data['mutation']

    n_genes = len(genes)
    n_pathways = len(pathways)

    # Adaptive font sizing for genes (boosted for projector/small-print legibility)
    if n_genes <= 40:
        gene_fontsize = 10.0
    elif n_genes <= 55:
        gene_fontsize = 8.5
    else:
        gene_fontsize = 7.5

    # Adaptive pathway font sizing
    if n_pathways <= 6:
        pw_fontsize = 11.5
    elif n_pathways <= 8:
        pw_fontsize = 10.5
    elif n_pathways <= 11:
        pw_fontsize = 9.5
    else:
        pw_fontsize = 8.5

    # Asymmetric xlim: extra room on the right for long pathway labels and
    # (on expanded_full) the outer group-bracket arc, which sits outside the
    # radial label extent.
    if variant == 'expanded_full':
        xlim_left, xlim_right = -1.85, 3.05
    else:
        xlim_left, xlim_right = -1.85, 2.55
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'aspect': 'equal'})
    ax.set_xlim(xlim_left, xlim_right)
    ax.set_ylim(-1.85, 1.85)
    ax.axis('off')

    cmap = create_diverging_cmap()
    norm = Normalize(vmin=-3, vmax=3)

    # Symmetric layout: pathways centered on mid=0 (right), genes centered on
    # mid=180 (left). Equal top/bottom gaps ensure no collision between the
    # pathway-arc ends and the gene-arc ends.
    if n_pathways <= 6:
        gene_arc_total    = 180
        pathway_arc_total = 110
    elif n_pathways <= 7:
        gene_arc_total    = 180
        pathway_arc_total = 120
    elif n_pathways <= 10:
        gene_arc_total    = 170
        pathway_arc_total = 140
    else:
        gene_arc_total    = 160
        pathway_arc_total = 160

    pathway_start_angle = -pathway_arc_total / 2          # symmetric about 0°
    gene_start_angle    = 180 - gene_arc_total / 2        # symmetric about 180°

    gene_gap    = 0.4
    pathway_gap = 4.0 if n_pathways <= 10 else 2.8

    # Gene layout
    total_gene_gap = gene_gap * max(n_genes - 1, 0)
    gene_width = (gene_arc_total - total_gene_gap) / max(n_genes, 1)
    gene_positions = {}
    current = gene_start_angle
    for g in genes:
        mid = current + gene_width / 2
        gene_positions[g['name']] = {
            'start': current, 'end': current + gene_width, 'mid': mid,
            'logfc_d35': g['logfc_d35'], 'logfc_d65': g['logfc_d65'],
        }
        current += gene_width + gene_gap

    # Pathway layout
    total_pw_gap = pathway_gap * max(n_pathways - 1, 0)
    pw_width = (pathway_arc_total - total_pw_gap) / max(n_pathways, 1)
    pathway_positions = {}
    current = pathway_start_angle
    for p in pathways:
        color, tier = PATHWAY_STYLE.get(p['id'], ('#888888', TIER_RECESSIVE))
        tparams = TIER_PARAMS[tier]
        pathway_positions[p['id']] = {
            'start': current, 'end': current + pw_width,
            'mid': current + pw_width / 2,
            'color': color, 'tier': tier,
            'arc_alpha': tparams['arc_alpha'],
            'ribbon_alpha': tparams['ribbon_alpha'],
            'zorder': tparams['zorder'],
            'group': p['group'],
            **p,
        }
        current += pw_width + pathway_gap

    gene_radius = 0.86
    pathway_radius = 0.86

    # ---- Layer 1 (background): gene-set-only ribbons --------------------
    # Faint ribbons connecting a displayed gene to a pathway when the gene is
    # ANNOTATED to that pathway (gene-set membership) but did NOT make the
    # GSEA leading edge for this contrast. Makes the true biological
    # containment visible without overstating GSEA-derived membership.
    gs_only_alpha = 0.13  # faint background layer
    for gene_name, pid in connections_gs_only:
        if gene_name not in gene_positions or pid not in pathway_positions:
            continue
        gp = gene_positions[gene_name]
        pp = pathway_positions[pid]
        draw_ribbon(ax, gp['mid'], gene_width,
                    pp['start'], pp['end'],
                    pp['color'], alpha=gs_only_alpha,
                    gene_radius=gene_radius, pathway_radius=pathway_radius,
                    zorder=0)

    # ---- Layer 2 (foreground): leading-edge ribbons ---------------------
    # Drawn in tier z-order: envelope/counter (under), recessive (middle),
    # focal synaptic ribosomes (on top). Cool focal pair above warm envelope.
    tier_order = [TIER_ENVELOPE, TIER_COUNTER, TIER_RECESSIVE, TIER_FOCAL]
    for tier in tier_order:
        for gene_name, pid in connections:
            if gene_name not in gene_positions or pid not in pathway_positions:
                continue
            pp = pathway_positions[pid]
            if pp['tier'] != tier:
                continue
            gp = gene_positions[gene_name]
            draw_ribbon(ax, gp['mid'], gene_width,
                        pp['start'], pp['end'],
                        pp['color'], alpha=pp['ribbon_alpha'],
                        gene_radius=gene_radius, pathway_radius=pathway_radius,
                        zorder=pp['zorder'])

    # Gene rings + labels
    for g in genes:
        pos = gene_positions[g['name']]
        draw_gene_rings(ax, pos['mid'], gene_width,
                        pos['logfc_d35'], pos['logfc_d65'],
                        cmap, norm)
        label_r = 1.04
        rad = np.radians(pos['mid'])
        lx, ly = label_r * np.cos(rad), label_r * np.sin(rad)
        if 90 <= pos['mid'] <= 270:
            rotation = pos['mid'] - 180
            ha = 'right'
        else:
            rotation = pos['mid']
            ha = 'left'
        ax.text(lx, ly, g['name'], fontsize=gene_fontsize, ha=ha, va='center',
                rotation=rotation, rotation_mode='anchor',
                fontweight='medium' if g['is_ribosomal'] else 'normal',
                color='#222222' if g['is_ribosomal'] else '#444444')

    # Pathway arcs + labels
    for pid, pp in pathway_positions.items():
        draw_pathway_arc(ax, pp['start'], pp['end'], pp['color'],
                         pp['display_name'],
                         fontsize_label=pw_fontsize,
                         arc_alpha=pp['arc_alpha'],
                         zorder=pp['zorder'] + 3)

    # Group brackets: kept on expanded_full only (warm-vs-cool palette already
    # conveys grouping at a glance on expanded_clean; focused has one group).
    if variant == 'expanded_full':
        groups = {}
        for pid, pp in pathway_positions.items():
            groups.setdefault(pp['group'], []).append((pp['start'], pp['end']))
        for grp_name, spans in groups.items():
            if len(spans) < 2:
                continue
            t0 = min(s[0] for s in spans)
            t1 = max(s[1] for s in spans)
            draw_group_bracket(ax, t0, t1,
                               GROUP_SHORT_LABELS.get(grp_name, grp_name),
                               GROUP_COLORS.get(grp_name, '#444444'),
                               bracket_r=2.55, label_r_offset=0.12,
                               fontsize=12.0)

    # Title
    mutation_color = MUTATION_COLORS.get(mutation, '#222222')
    variant_label = {
        'focused':         'SynGO synaptic-localized pathways',
        'expanded_clean':  'Synaptic ribosomes nested in broader cytoplasmic ribosome',
        'expanded_full':   'Synaptic ribosomes within broader cytoplasmic translation program',
    }[variant]
    fig.suptitle(f'{mutation} mutation  ·  {variant_label}',
                 fontsize=18, fontweight='bold', color=mutation_color,
                 y=0.97)

    # logFC colorbar
    cbar_ax = fig.add_axes([0.04, 0.30, 0.020, 0.32])
    sm = ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('log₂ fold change (gene rings)', fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Ring + ribbon explanation (top-left legend box)
    leg = fig.add_axes([0.005, 0.63, 0.14, 0.30]); leg.axis('off')
    leg.set_xlim(0, 1); leg.set_ylim(0, 1)
    leg.text(0.5, 0.97, 'Gene rings', fontsize=12, fontweight='bold',
             ha='center', va='top')
    leg.text(0.5, 0.85, 'Outer ring = D35', fontsize=11, ha='center', va='center')
    leg.text(0.5, 0.76, 'Inner ring = D65', fontsize=11, ha='center', va='center')
    leg.text(0.5, 0.64, 'Concentric reading:\nouter → inner = early → late',
             fontsize=9.5, ha='center', va='center', color='#555555',
             fontstyle='italic')

    leg.text(0.5, 0.44, 'Ribbons', fontsize=12, fontweight='bold',
             ha='center', va='top')
    leg.add_patch(Rectangle((0.05, 0.28), 0.20, 0.06,
                            facecolor='#AA4499', edgecolor='none', alpha=0.62))
    leg.text(0.30, 0.31, 'GSEA leading edge\n(this contrast)',
             fontsize=9.0, va='center', ha='left')
    leg.add_patch(Rectangle((0.05, 0.13), 0.20, 0.06,
                            facecolor='#AA4499', edgecolor='none', alpha=0.13))
    leg.text(0.30, 0.16, 'Gene-set annotation\n(not in leading edge)',
             fontsize=9.0, va='center', ha='left')

    # Pathway legend (bottom-right), grouped by category
    if variant == 'expanded_full':
        pl_x, pl_y, pl_w, pl_h = 0.855, 0.04, 0.145, 0.62
        grp_fs, name_fs, swatch_w = 9.5, 8.8, 0.07
        row_dy_grp, row_dy_item = 0.038, 0.044
    else:
        pl_x, pl_y, pl_w, pl_h = 0.855, 0.08, 0.145, 0.54
        grp_fs, name_fs, swatch_w = 10.0, 9.5, 0.08
        row_dy_grp, row_dy_item = 0.050, 0.060

    pl = fig.add_axes([pl_x, pl_y, pl_w, pl_h]); pl.axis('off')
    pl.set_xlim(0, 1); pl.set_ylim(0, 1)
    pl.text(0.5, 0.99, 'Pathways', fontsize=12, fontweight='bold',
            ha='center', va='top')

    grouped = {}
    for p in pathways:
        grouped.setdefault(p['group'], []).append(p)

    y = 0.93
    for grp_name, plist in grouped.items():
        grp_color = GROUP_COLORS.get(grp_name, '#444444')
        pl.text(0.0, y, grp_name, fontsize=grp_fs, fontweight='bold',
                va='center', color=grp_color)
        y -= row_dy_grp
        for p in plist:
            c = PATHWAY_COLOR_MAP.get(p['id'], '#888888')
            pl.add_patch(Rectangle((0.03, y - 0.020), swatch_w, 0.035,
                                   facecolor=c, edgecolor='white', linewidth=0.5))
            name = p['display_name'].replace('\n', ' ')
            if len(name) > 32:
                name = name[:30] + '…'
            pl.text(0.03 + swatch_w + 0.03, y, name, fontsize=name_fs, va='center')
            y -= row_dy_item
        y -= 0.012
        if y < 0.02:
            break

    if variant == 'expanded_full':
        plt.subplots_adjust(left=0.12, right=0.85, top=0.94, bottom=0.04)
    else:
        plt.subplots_adjust(left=0.13, right=0.85, top=0.94, bottom=0.04)

    pdf_path = output_path.with_suffix('.pdf')
    png_path = output_path.with_suffix('.png')
    fig.savefig(pdf_path, dpi=600, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    fig.savefig(png_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"   Saved: {pdf_path.name}, {png_path.name}")


def save_chord_data(data, variant, mutation):
    suffix = '' if variant == 'focused' else f'_{variant}'
    gene_df = pd.DataFrame(data['genes']); gene_df['mutation'] = mutation
    path_df = pd.DataFrame(data['pathways']); path_df['mutation'] = mutation
    conn_df = pd.DataFrame(data['connections'], columns=['gene', 'pathway'])
    conn_df['mutation'] = mutation

    gene_df.to_csv(OUTPUT_DIR / f'chord_genes_{mutation}{suffix}.csv', index=False)
    path_df.to_csv(OUTPUT_DIR / f'chord_pathways_{mutation}{suffix}.csv', index=False)
    conn_df.to_csv(OUTPUT_DIR / f'chord_connections_{mutation}{suffix}.csv', index=False)

    gs_only = data.get('connections_gs_only', [])
    if gs_only:
        gs_df = pd.DataFrame(gs_only, columns=['gene', 'pathway'])
        gs_df['mutation'] = mutation
        gs_df.to_csv(OUTPUT_DIR / f'chord_connections_geneset_only_{mutation}{suffix}.csv',
                     index=False)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 64)
    print("Chord-diagram visualization (focused + expanded variants)")
    print("=" * 64)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n[1] Loading GSEA data (syngo + all categories)...")
    gsea_df = load_gsea_data()
    print(f"    {len(gsea_df)} pathway-contrast rows")

    print("\n[2] Loading DE results...")
    de_data = load_de_results()
    print(f"    {len(de_data)} contrasts")

    print("\n[2b] Loading gene-set memberships (full annotation, not LE)...")
    gs_membership_df = load_gene_set_memberships()
    print(f"    {len(gs_membership_df)} pathway-gene annotations")

    for mutation in ['G32A', 'R403C']:
        for variant, pids in VARIANTS.items():
            print(f"\n[3] {mutation} :: {variant}")
            data = prepare_chord_data(gsea_df, de_data, mutation, pids,
                                       max_genes=60 if variant == 'focused' else 70,
                                       gs_membership_df=gs_membership_df)
            print(f"    pathways={len(data['pathways'])}  genes={len(data['genes'])}  "
                  f"LE ribbons={len(data['connections'])}  "
                  f"gene-set-only ribbons={len(data['connections_gs_only'])}")
            suffix = '' if variant == 'focused' else f'_{variant}'
            output_path = OUTPUT_DIR / f'chord_diagram_{mutation}{suffix}'
            draw_chord_diagram(data, output_path, variant)
            save_chord_data(data, variant, mutation)

    print("\n" + "=" * 64)
    print(f"Done. Output dir: {OUTPUT_DIR}")
    print("=" * 64)


if __name__ == '__main__':
    main()
