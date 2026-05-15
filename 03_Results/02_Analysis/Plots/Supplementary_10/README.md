# Supplementary_10 — Replicate-Level GSVA Analysis (Supp Fig S10)

Replicate-level GSVA scores for the focused ribosome and mitochondrial modules, backing **Supplementary Fig S10**. These panels are essential for distinguishing the biological meaning behind contrast-level GSEA pattern labels: a contrast-level "Compensation" can be `ctrl_driven` (controls descend toward mutant baseline) rather than a genuine mutant rebound, and a contrast-level "Sign_reversal" can be backed by a genuine sample-level crossover or just a noisy slope flip.

## File inventory

| File | Description |
|---|---|
| `replicate_level_gsva.pdf` | Supp Fig S10: strip plots of per-replicate GSVA scores for focused modules with median lines and IQR bands |
| `replicate_level_gsva.png` | Raster version (300 dpi) |
| `replicate_level_gsva_per_sample.csv` | Raw GSVA enrichment scores for each biological replicate × module |
| `replicate_level_gsva_group_summary.csv` | Summary statistics (median, IQR) per genotype × timepoint × module cell |

## Generating scripts

- `02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py` — main visualisation
- `02_Analysis/revision/supplements/Supp10a.export_gsva_modules.R` — exports per-sample GSVA scores (input CSV for the Python step)

```bash
Rscript 02_Analysis/revision/supplements/Supp10a.export_gsva_modules.R
python3 02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py
```

## How to read this folder

`replicate_level_gsva.pdf` is the panel that goes into Supp Fig S10 — open this first. The two CSVs back the figure: `_per_sample.csv` is the raw per-replicate matrix (one row per biological replicate × module) used to draw the dots, and `_group_summary.csv` carries the medians and IQRs used for the white rings and shaded bands. To regenerate, rerun the R export first, then the Python plotting step.

## Manuscript figure caption (Supp Fig S10)

**Replicate-level GSVA scores reveal that GSEA pattern labels can carry different biological meanings across pathway buckets.** Each dot is one biological-replicate sample, coloured by genotype (Ctrl, grey; G32A, blue; R403C, vermillion); white rings mark the within-cell median and shaded bands span the IQR. White trajectory lines connect the D35 and D65 medians for each genotype, so the underlying replicate-level dynamics behind each module's GSEA pattern label are directly visible. Panels show the focused structural-vs-biogenesis ribosome modules together with a combined mitochondrial-programme context panel (ATP Hydrolysis ○ + OXPHOS ◇). Three pattern interpretations are critical to read together. (1) **Synaptic Ribosome** and **Cytoplasmic Translation** (top row, left and middle) show a genuine **sample-level crossover** between control and mutants — controls rise during maturation while G32A and R403C lines drop — which is the meaning of the *Sign_reversal* GSEA label for these modules. (2) **Cytoplasmic Ribosome — biogenesis** (top-right) is classified as *Compensation* in the contrast-level GSEA view, but the per-replicate panel shows that this label is driven by **control developmental descent closing the gap** rather than by a mutant rebound: mutant medians remain near baseline while controls fall from +0.4 toward −0.4. We refer to this driver classification as `ctrl_driven`. (3) The **mitochondrial-ribosome modules** (bottom row) show a third pattern, where both controls and mutants descend together with mutants offset more negatively — a coherent down-trajectory rather than a divergence (driver classification: `both_moving`). Approximate 95% confidence intervals are shown for qualitative comparison only (n = 3–6 per cell, GSVA Gaussian kernel; pathway-size filter 10–500 genes). The same per-pathway replicate-level view is available for every pathway in the master GSEA table via `../Trajectory_Flow/interactive_bump_dashboard.html` (Supplementary Data File 1).

---

**Last Updated**: 2026-05-15
**Generating scripts**: `02_Analysis/revision/supplements/Supp10.replicate_level_gsva.py`, `02_Analysis/revision/supplements/Supp10a.export_gsva_modules.R`
