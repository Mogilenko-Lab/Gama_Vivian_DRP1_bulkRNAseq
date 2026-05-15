# Trajectory_Flow/alluvial — Sankey / Alluvial Flows (Exploratory)

Exploratory alluvial / Sankey diagrams that visualise the **flow** of pathway dynamics from Early defect status → mechanism (Active TrajDev vs Passive buffering) → Late outcome. These are exploratory companions to the bump charts in `../` and `../bump/`; they are not selected as main-text figures but are linked from the supplementary text discussion of pattern-class flow.

## File inventory

| File | Description |
|---|---|
| `alluvial_binary_G32A.html` / `.png` | Binary flow (Has Early Defect vs No Defect) for G32A; interactive Plotly + static PNG |
| `alluvial_binary_R403C.html` / `.png` | Same, R403C |
| `alluvial_graded_G32A.html` / `.png` | Graded flow (Strong / Moderate / No defect) for G32A |
| `alluvial_graded_R403C.html` / `.png` | Same, R403C |
| `alluvial_ggalluvial_G32A.pdf` / `.png` | Classical R ggalluvial rendering for G32A |
| `alluvial_ggalluvial_R403C.pdf` / `.png` | Same, R403C |
| `alluvial_combined.pdf` / `.png` | Side-by-side combined ggalluvial composite |

## Flow structure

```
   Early status            Mechanism                   Late outcome
   ────────────            ─────────                   ────────────
  ┌─────────────┐         ┌──────────┐               ┌──────────┐
  │ Early defect│  ─────► │ Active   │  ──────────►  │ Improved │
  └─────────────┘         │ (TrajDev)│               └──────────┘
        │                 └──────────┘                    ▲
        │                       │                         │
        │                 ┌──────────┐                    │
        └───────────────► │ Passive  │  ──────────────────┘
                          │ (Buffer) │
                          └──────────┘

  ┌─────────────┐         ┌──────────┐               ┌──────────┐
  │ No early    │  ─────► │Late-onset│  ──────────►  │New defect│
  │   defect    │         └──────────┘               └──────────┘
  └─────────────┘
```

- **LEFT** = Early defect status / severity
- **MIDDLE** = Mechanism (Active TrajDev vs Passive buffering)
- **RIGHT** = Late outcome (Improved, Resolved, Worsened, New defect)
- **Late_onset** branches form a separate stream with no left-side input

## Generating scripts

| Script | Outputs |
|---|---|
| `02_Analysis/3.5.viz_trajectory_flow.py` | `alluvial_binary_*.html/.png`, `alluvial_graded_*.html/.png` (interactive Plotly) |
| `02_Analysis/3.6.viz_alluvial_ggalluvial.R` | `alluvial_ggalluvial_*.pdf/.png`, `alluvial_combined.*` (classical R) |

```bash
python3 02_Analysis/3.5.viz_trajectory_flow.py
Rscript 02_Analysis/3.6.viz_alluvial_ggalluvial.R
```

Both scripts read `03_Results/02_Analysis/master_gsea_table.csv` and pattern definitions from `01_Scripts/Python/pattern_definitions.py`.

## How to read this folder

For an at-a-glance flow summary, open `alluvial_combined.pdf`. For interactive exploration (hover for stream-level pathway counts and pattern composition), open the `*.html` files in a browser. The graded variants split the Early defect column by severity (Strong / Moderate / None) and are useful for showing that Compensation predominantly originates from Strong Early defects, while Late_onset and Natural_improvement populate the No-Early-Defect stream.

## See also

- `../README.md` — parent folder, Fig 5B / 5C bump charts
- `../bump/README.md` — bump chart variant inventory

---

**Last Updated**: 2026-05-15
**Generating scripts**: `02_Analysis/3.5.viz_trajectory_flow.py`, `02_Analysis/3.6.viz_alluvial_ggalluvial.R`
