# Critical Period Trajectories — GSVA Temporal Trajectory Analysis

This directory holds GSVA-based trajectory visualisations for seven functional modules tracked across the D35 → D65 neuronal maturation window in control and DRP1 mutant neurons. The figures document the paradoxical decoupling between cytoplasmic ribosome production and synaptic ribosome programs during the critical period for cortical synaptogenesis. These are supplementary figures that contextualise the GSEA-level findings in `../Synaptic_ribosomes/` (Fig 5F) and `../Trajectory_Flow/` (Fig 5B, 5C).

## Directory structure

```
Critical_period_trajectories/
└── gsva/
    ├── combined/     # Multi-panel and overlay figures
    ├── data/         # Exported CSV tables
    ├── divergence/   # Per-module divergence-from-control panels
    └── trajectory/   # Per-module absolute GSVA trajectory panels
```

## File inventory

### gsva/combined/

| File | Description |
|---|---|
| `Trajectory_7panel_grid.pdf` | 3 × 3 grid of all 7 module trajectory panels |
| `Trajectory_7panel_grid-1.png` | PNG version of the 7-panel grid |
| `Divergence_overlay_all.pdf` | All 7 modules overlaid as divergence from control (multi-color) |
| `Divergence_overlay_mitocarta.pdf` | 4 mitochondrial modules only (Panels D–G) overlaid |

### gsva/trajectory/ (per-module GSVA panels — absolute scores)

Filenames carry an explicit `_trajectory` suffix to disambiguate from the parallel divergence views below.

| File | Module | Key pattern |
|---|---|---|
| `Panel_A_Ribosome_Biogenesis_trajectory.pdf` | Ribosome Biogenesis | Control decreases; mutants resist downregulation (compensatory maintenance) |
| `Panel_B_Cytoplasmic_Translation_trajectory.pdf` | Cytoplasmic Translation | Control increases; mutants fail to increase (functional failure) |
| `Panel_C_Synaptic_Ribosomes_trajectory.pdf` | Synaptic Ribosomes | Control increases; mutants decline (critical-period crisis) |
| `Panel_D_Mitochondrial_Ribosome_trajectory.pdf` | Mitochondrial Ribosome | Shared compensation in both mutations |
| `Panel_E_Mito_Ribosome_Assembly_trajectory.pdf` | Mito Ribosome Assembly | Shared compensation in both mutations |
| `Panel_F_ATP_Hydrolysis_trajectory.pdf` | ATP Hydrolysis | Shared compensation in both mutations |
| `Panel_F_mtDNA_Maintenance_trajectory.pdf` | mtDNA Maintenance | Shared compensation in both mutations |
| `Panel_G_OXPHOS_trajectory.pdf` | OXPHOS | Shared compensation in both mutations |

### gsva/divergence/ (per-module GSVA panels — mutation-specific divergence from control)

Same set of panels as `trajectory/` above; y-axis is the divergence from within-timepoint control mean (control trajectory removed). Filenames carry a `_divergence` suffix.

| File | Module |
|---|---|
| `Panel_A_Ribosome_Biogenesis_divergence.pdf` | Ribosome Biogenesis |
| `Panel_B_Cytoplasmic_Translation_divergence.pdf` | Cytoplasmic Translation |
| `Panel_C_Synaptic_Ribosomes_divergence.pdf` | Synaptic Ribosomes |
| `Panel_D_Mitochondrial_Ribosome_divergence.pdf` | Mitochondrial Ribosome |
| `Panel_E_Mito_Ribosome_Assembly_divergence.pdf` | Mito Ribosome Assembly |
| `Panel_F_ATP_Hydrolysis_divergence.pdf` | ATP Hydrolysis |
| `Panel_F_mtDNA_Maintenance_divergence.pdf` | mtDNA Maintenance |
| `Panel_G_OXPHOS_divergence.pdf` | OXPHOS |

### gsva/data/

| File | Description |
|---|---|
| `gsva_scores_matrix.csv` | GSVA enrichment scores for all 7 modules across 6 groups |
| `trajectory_data.csv` | Mean trajectory values (relative to Ctrl D35 baseline) |
| `divergence_data.csv` | Mutation-specific divergence from control at each timepoint |
| `module_gene_counts.csv` | Gene counts per module |

## Generating script

`02_Analysis/2.4.viz_critical_period_trajectories_gsva.R`

```bash
Rscript 02_Analysis/2.4.viz_critical_period_trajectories_gsva.R
```

Requires checkpoints from `02_Analysis/1.1.main_pipeline.R`. GSVA scores are cached in `checkpoints/gsva_7modules_checkpoint.rds`.

## Reading guide

**Trajectory panels (absolute scale)**: Y-axis = GSVA enrichment score relative to Control D35 baseline. The Control trajectory shows normal developmental change; mutant trajectories show genotype-specific deviation. Individual sample points (semi-transparent dots, n = 3–6 per group) show biological variability. Mean lines connect group means at D35 and D65. Color: gray = Control, orange = G32A, green = R403C.

**Divergence panels (relative scale)**: Y-axis = mutation effect at each timepoint after removing developmental baseline (Control is zero by definition). These isolate genotype effects and are more interpretable for comparing mutation severity at each stage.

**The paradox in these panels**: Panels A–B show opposing patterns — ribosome biogenesis goes up (Panel A, Compensation pattern) while cytoplasmic translation goes down (Panel B, Sign_reversal). Panel C (Synaptic Ribosomes) shows the strongest failure, consistent with the GSEA-level NES ≈ −3.0 for the TrajDev contrast. Panels D–G show coordinated mitochondrial compensatory upregulation shared by both mutations.

**Note on biogenesis driver**: Cytoplasmic ribosome biogenesis (Panel A) is classified as Compensation at the contrast level but is `ctrl_driven` at the replicate level — controls descend toward mutant baseline during maturation rather than mutants actively rebounding. Panel C (Synaptic Ribosomes) shows a genuine sample-level crossover (Sign_reversal). See `../Supplementary_10/` for the replicate-level view.

## Manuscript supplementary caption (`Trajectory_7panel_grid.pdf`)

**GSVA trajectories for seven curated DRP1-relevant modules across the D35 → D65 critical period.** Each panel shows one functional module: (A) Ribosome Biogenesis (GO:BP, 158 genes), (B) Cytoplasmic Translation (GO:BP, 76 genes), (C) Synaptic Ribosomes (SynGO postsynaptic ribosome, 65 genes detected), (D) Mitochondrial Ribosome (MitoCarta, 77 genes), (E) Mito Ribosome Assembly (MitoCarta, 24 genes), (F) mtDNA Maintenance (MitoCarta, 29 genes), (G) OXPHOS (MitoCarta, 139 genes). Y-axis = per-sample GSVA enrichment score (Gaussian kernel; pathway-size filter 10–500 genes); x-axis = developmental timepoint (D35, D65). Individual sample points (semi-transparent dots, n = 3–6 per group) overlay genotype mean lines (Ctrl, grey; G32A, orange; R403C, green). Significance markers on each panel are FDR-corrected two-sample t-tests against time-matched controls (BH-corrected; p_adj < 0.05). The figure captures three replicate-level patterns: Panel A is a `ctrl_driven` Compensation (controls descend toward mutant baseline rather than mutants rebounding), Panel C is a genuine sample-level crossover (`both_moving` Sign_reversal: controls ascend while mutants descend), and Panels D–G are coordinated `both_moving` compensatory upregulations shared by both mutations. The dissociation between ribosome biogenesis upregulation (Panel A) and assembled cytoplasmic / synaptic translation collapse (Panels B–C) is the GSVA-level signature of the synaptic translation paradox; the coordinated mitochondrial compensation (Panels D–G) is the parallel mitochondrial arm. The accompanying `Divergence_overlay_*.pdf` panels remove the control developmental trajectory to isolate genotype effects at each timepoint.

---

**Last Updated**: 2026-05-15
**Generating script**: `02_Analysis/2.4.viz_critical_period_trajectories_gsva.R`
