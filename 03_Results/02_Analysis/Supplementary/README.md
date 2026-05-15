# 03_Results/02_Analysis/Supplementary — reviewer-response digests and expression-filter audits

## Overview

Compact CSV / TSV artifacts produced during reviewer-response work, intended to be cited verbatim from supplementary text or figure legends. The 6a sensitivity digests confirm that the key finding — Compensation > passive recovery in 81/81 threshold combinations — holds regardless of classifier parameterisation. The 5a files audit the `AveExpr > 0` expression filter impact on calcium gene visibility. The 5b, 6c files support UpSet/Euler and ribosome-scope claims respectively.

## Files (concern 6a — fully documented)

### Trajectory-pattern sensitivity grid

The classifier from `01_Scripts/Python/pattern_definitions.py` was re-run across an
81-combination threshold grid (NES_EFFECT × NES_STRONG × IMPROVEMENT_RATIO × WORSENING_RATIO).
Two digests are kept side-by-side because they answer different questions:

| File | Universe | Generator | When to cite |
|---|---|---|---|
| `6a_sensitivity_stability_digest.tsv` | **12,221 universe** (all GSEA pathways, including those not significant in any contrast) | `02_Analysis/Supp4.sensitivity_analysis.py` (Apr 23) | When the question is "is the classifier itself robust to thresholds, irrespective of which pathways are ever significant?" — does NOT match RESULTS L11/13/17 percentages |
| `6a_sensitivity_stability_digest_5267universe.tsv` | **5,267 universe** (FDR<0.05 in any of 9 GSEA contrasts; matches manuscript Results section) | `02_Analysis/6a.sensitivity_5267universe.py` (Apr 24, external audit) | **Default for paper citations.** Numbers in this file align with RESULTS L11/13/17 and with for-the-paper.md Methods Edit 1 |
| `sensitivity_5267universe.csv` | 5,267 universe (full grid: 81 combos × 2 mutations × 8 patterns + denominator columns) | same generator | Reproducibility / audit; not directly cited |

**Critical distinction:** the two digests use *different denominators* and therefore quote
*different percentage ranges* for the same claim. Cross-reference table:

| Quantity | 12,221-universe digest | 5,267-universe digest (paper-aligned) |
|---|---|---|
| G32A Comp/classifiable range | 52.3–55.3% (81/81 majority) | **54.5–58.5%** (81/81 majority) |
| R403C Comp/classifiable range | 43.6–49.4% (0/81 majority) | **46.9–54.5%** (54/81 majority, plurality in 27) |
| G32A NI/classifiable range | 25–36% | **21.9–32.7%** |
| R403C NI/classifiable range | 30–43% | **25.0–37.3%** |
| G32A Compensation count range | 1190–1590 | 1190–1590 (denominator-invariant) |
| R403C Compensation count range | 1246–1832 | 1246–1832 (denominator-invariant) |
| `comp_exceeds_passive` / `progressive_rare` / `R403C_more_compensation` | 81/81 | 81/81 (all unchanged) |

See [`../Tables/README_pattern_summary_denominators.md`](../Tables/README_pattern_summary_denominators.md) for the canonical denominator vocabulary.

## Files (other concerns — listed for completeness)

| File | Concern | Generator |
|---|---|---|
| `5a_aveexpr_diagnostics.csv` | 5a — AveExpr filter impact per contrast | `02_Analysis/revision/supplements/5a.filtered_volcano_supplement.R` |
| `5a_aveexpr_calcium_genes.csv` | 5a — per-gene AveExpr values for the 13 calcium genes | same |
| `5a_aveexpr_filter_impact.csv` | 5a — gene counts before/after AveExpr > 0 filter | same |
| `5b_intersection_counts.tsv` | 5b — DEG intersection counts for UpSet/Euler plots | `02_Analysis/revision/supplements/5b.highlighted_upset.R` / `5b.maturation_euler.R` |
| `6c_cytoplasmic_ribo_nes.csv` | 6c — NES values for cytoplasmic ribosome pathways | `02_Analysis/revision/supplements/6c.extract_cyto_ribo_nes.py` |
| `6c_ribosome_jaccard.csv` | 6c — pairwise Jaccard between ribosome compartments | `02_Analysis/revision/supplements/6c.compute_jaccard.py` |

### Column dictionaries (5a files)

**`5a_aveexpr_diagnostics.csv`** — one row per contrast:
`contrast`, `n_total` (genes entering DE), `n_AveExpr_lt_0` (genes with AveExpr < 0), `pct_AveExpr_lt_0`, `n_AveExpr_lt_0p5`, `n_AveExpr_lt_1`, `AveExpr_min/max`, `n_sig_FDR05` (significant genes at FDR 0.05), `n_sig_FDR05_AveExpr_le_0` (significant genes that would be removed by the filter).

**`5a_aveexpr_calcium_genes.csv`** — one row per (contrast × calcium gene):
`contrast`, `gene`, `AveExpr`, `logFC`, `P.Value`, `adj.P.Val`, `passes_AveExpr_gt_0` (TRUE/FALSE — whether the gene survives the expression filter).

**`5a_aveexpr_filter_impact.csv`** — one row per contrast:
`contrast`, `n_before`, `n_after`, `n_removed`, `pct_removed` (fraction of genes removed by the AveExpr > 0 filter; ~28% of genes on average have AveExpr < 0 and would be excluded).

## Cross-references

- **Manuscript anchors that cite the 6a files:**
  - `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/for-the-paper.md`:
    - Header denominator-vocabulary table (lists both digest files)
    - Edit 1 Methods (cites `6a_sensitivity_stability_digest_5267universe.tsv` ranges)
    - Edit 5 Supp Fig S8 legend (cites both digests)
  - `Manuscript/current_submission/docs/6a_interaction_gsea_vs_trajectories/status.md` Q6 / Q7 audit-log entries
- **Upstream data:** `../master_gsea_table.csv`; `../Sensitivity_Analysis/sensitivity_results.csv` (the 12,221-universe sensitivity table from which the 12,221 digest is summarised).

## How to regenerate (concern 6a only)

```bash
# 12,221-universe digest (legacy)
python3 02_Analysis/Supp4.sensitivity_analysis.py
# (also produces sensitivity_results.csv and the heatmap PDFs in ../Sensitivity_Analysis/)

# 5,267-universe digest + full grid (paper-aligned; preferred)
python3 02_Analysis/6a.sensitivity_5267universe.py
```

Both are deterministic, read-only on `master_gsea_table.csv`, and run in 1–3 minutes.

## Read-only constraints

- The 12,221-universe digest (`6a_sensitivity_stability_digest.tsv`) is **untouched** and remains valid for the question it answers. The 5,267-universe digest is **strictly additive**.
- No upstream artifact (master GSEA table, classifier code, sensitivity script) is modified.
