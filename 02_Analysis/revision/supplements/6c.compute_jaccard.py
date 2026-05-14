#!/usr/bin/env python3
"""
6c.compute_jaccard.py
=====================

Helper script for the ribosome-scope concern.

Computes Jaccard indices and raw overlap counts between the three ribosomal
gene sets used in the paper:

  - Synaptic    (SynGO Cellular Component "ribosome" terms)
  - Mitochondrial (MitoCarta Mitochondrial_ribosome + assembly + translation factors)
  - Cytoplasmic (curated GO cytoplasmic ribosome + biogenesis factors)

Gene-set loaders are reused by importing from 3.3.ribosome_upset_plot.py (the
same loaders that feed Fig1b_Ribosome_Gene_Overlap_UpSet.pdf). We do NOT
re-render the UpSet plot — we only extract the Jaccard matrix that is logged
to stdout inside that script (lines 237-245 of 3.3) and persist it to CSV for auditable manifest.

Outputs:
  - stdout log (Jaccard values + raw overlaps + MRPL/MRPS sanity check)
  - 03_Results/02_Analysis/Supplementary/6c_ribosome_jaccard.csv

Run from project root:
  python 02_Analysis/6c.compute_jaccard.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


# ------------------------------------------------------------------- paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
UPSET_SCRIPT = PROJECT_ROOT / "02_Analysis" / "3.3.ribosome_upset_plot.py"
OUTPUT_DIR = PROJECT_ROOT / "03_Results" / "02_Analysis" / "Supplementary"
OUTPUT_CSV = OUTPUT_DIR / "6c_ribosome_jaccard.csv"

# Ensure config module is importable (3.3 imports Python.config)
sys.path.insert(0, str(PROJECT_ROOT / "01_Scripts"))


def _load_upset_loaders():
    """Dynamically import the 3.3.ribosome_upset_plot.py loaders by spec."""
    spec = importlib.util.spec_from_file_location("ribosome_upset_plot", UPSET_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("6c.compute_jaccard — ribosome gene-set overlap")
    print("=" * 70)

    upset = _load_upset_loaders()

    # Build gene sets via the same loaders that feed the UpSet plot.
    print("\n[1/3] Loading MitoCarta mitochondrial ribosome gene set ...")
    mito_pathways = [
        "Mitochondrial_ribosome",
        "Mitochondrial_ribosome_assembly",
        "Translation_factors",
    ]
    mito_dict = upset.load_mitocarta_genes(mito_pathways)
    mito_all = set()
    for gs in mito_dict.values():
        mito_all |= gs

    print("\n[2/3] Loading SynGO synaptic ribosome gene set ...")
    syn_all = upset.load_syngo_genes()

    print("\n[3/3] Loading curated cytoplasmic GO ribosome gene set ...")
    cyto_all = upset.load_go_ribosome_genes()

    gene_sets: dict[str, set] = {
        "Synaptic_SynGO": syn_all,
        "Mitochondrial_MitoCarta": mito_all,
        "Cytoplasmic_GO": cyto_all,
    }

    # -------------------------------------------------------------- pairwise
    records: list[dict] = []
    names = list(gene_sets.keys())
    print("\n" + "-" * 70)
    print("PAIRWISE JACCARD & OVERLAP")
    print("-" * 70)
    for i, n1 in enumerate(names):
        for n2 in names[i + 1 :]:
            inter = gene_sets[n1] & gene_sets[n2]
            union = gene_sets[n1] | gene_sets[n2]
            j = len(inter) / len(union) if union else 0.0
            rec = {
                "set_1": n1,
                "set_2": n2,
                "size_set_1": len(gene_sets[n1]),
                "size_set_2": len(gene_sets[n2]),
                "intersection_n": len(inter),
                "union_n": len(union),
                "jaccard": round(j, 4),
                "jaccard_pct": round(j * 100, 2),
                "shared_genes_sample": ",".join(sorted(inter)[:20]),
            }
            records.append(rec)
            print(
                f"  Jaccard({n1}, {n2}) = {j:.4f}  "
                f"({len(inter)} shared of {len(union)} union)"
            )
            if inter:
                sample = sorted(inter)
                print(
                    f"    Shared (up to 20): "
                    f"{', '.join(sample[:20])}"
                    f"{'...' if len(sample) > 20 else ''}"
                )

    # -------------------------------------------------------------- triple
    triple = gene_sets[names[0]] & gene_sets[names[1]] & gene_sets[names[2]]
    union_all = gene_sets[names[0]] | gene_sets[names[1]] | gene_sets[names[2]]
    j_triple = len(triple) / len(union_all) if union_all else 0.0
    records.append(
        {
            "set_1": "ALL_THREE",
            "set_2": "ALL_THREE",
            "size_set_1": len(gene_sets[names[0]]),
            "size_set_2": len(gene_sets[names[1]]),
            "intersection_n": len(triple),
            "union_n": len(union_all),
            "jaccard": round(j_triple, 4),
            "jaccard_pct": round(j_triple * 100, 2),
            "shared_genes_sample": ",".join(sorted(triple)),
        }
    )
    print(
        f"\n  Triple intersection (all 3 sets): {len(triple)} genes "
        f"of {len(union_all)} union (Jaccard = {j_triple:.4f})"
    )
    if triple:
        print(f"    Genes: {', '.join(sorted(triple))}")

    # -------------------------------------------------------------- MRPL/MRPS sanity
    print("\n" + "-" * 70)
    print("SANITY CHECK: MRPL/MRPS ∩ SynGO synaptic (plan.md §8.5)")
    print("-" * 70)
    mrpl_mrps = {g for g in mito_all if isinstance(g, str) and g.startswith(("MRPL", "MRPS"))}
    mrpl_syn_overlap = mrpl_mrps & syn_all
    print(f"  MRPL/MRPS in MitoCarta ribosome set: {len(mrpl_mrps)}")
    print(f"  SynGO synaptic ribosome size:       {len(syn_all)}")
    print(f"  MRPL/MRPS ∩ SynGO:                  {len(mrpl_syn_overlap)}")
    if mrpl_syn_overlap:
        print(f"    Genes: {', '.join(sorted(mrpl_syn_overlap))}")
    else:
        print("    (empty — safe to state 'largely non-overlapping')")

    records.append(
        {
            "set_1": "MRPL_MRPS_subunits",
            "set_2": "Synaptic_SynGO",
            "size_set_1": len(mrpl_mrps),
            "size_set_2": len(syn_all),
            "intersection_n": len(mrpl_syn_overlap),
            "union_n": len(mrpl_mrps | syn_all),
            "jaccard": round(
                len(mrpl_syn_overlap) / max(1, len(mrpl_mrps | syn_all)), 4
            ),
            "jaccard_pct": round(
                100 * len(mrpl_syn_overlap) / max(1, len(mrpl_mrps | syn_all)), 2
            ),
            "shared_genes_sample": ",".join(sorted(mrpl_syn_overlap)),
        }
    )

    # -------------------------------------------------------------- write CSV
    out = pd.DataFrame(records)
    out.to_csv(OUTPUT_CSV, index=False)
    print("\n" + "-" * 70)
    print(f"Wrote: {OUTPUT_CSV}")
    print("-" * 70)


if __name__ == "__main__":
    main()
