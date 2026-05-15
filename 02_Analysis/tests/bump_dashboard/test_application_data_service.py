"""
Tests for application.data_service
=====================================

Uses a minimal stub DataFrame to exercise the data service without touching
the real master_gsea_table.csv.

The stub mirrors the exact column contract expected by ``DashboardDataService``
so that future schema changes in the real data are caught here first.
"""

import math
import pytest
import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.schema import DashboardConfig, PathwayRecord, DashboardMetadata
from Python.bump_dashboard.application.data_service import DashboardDataService


# ---------------------------------------------------------------------------
# Minimal stub data factory
# ---------------------------------------------------------------------------

def _make_stub_df(n: int = 5) -> pd.DataFrame:
    """
    Create a minimal stub DataFrame with the column contract of master_gsea_table.csv.

    Each row is duplicated for several contrasts (the upstream table is long-format:
    one row per pathway × contrast).
    """
    contrasts = [
        "G32A_vs_Ctrl_D35",
        "G32A_vs_Ctrl_D65",
        "Maturation_G32A_specific",
        "R403C_vs_Ctrl_D35",
        "R403C_vs_Ctrl_D65",
        "Maturation_R403C_specific",
        "Time_Ctrl",
    ]
    rows = []
    for i in range(n):
        pid = f"DB::pathway_{i:03d}"
        for contrast in contrasts:
            rows.append({
                "pathway_id":  pid,
                "ID":          f"PW{i}",
                "Description": f"Pathway {i} description text",
                "database":    "MitoCarta" if i % 2 == 0 else "SynGO",
                "contrast":    contrast,
                "new_name":    contrast,
                "category":    "Early" if "D35" in contrast else "Late",
                "mutation":    "G32A" if "G32A" in contrast else "R403C",
                "NES":         0.8 + i * 0.1,
                "pvalue":      0.01,
                "p.adjust":    0.02 + i * 0.01,
                "qvalue":      0.03,
                "enrichmentScore": 0.3,
                "setSize":     50,
                "Pattern_G32A":    "Compensation" if i < 3 else "Complex",
                "Confidence_G32A": "High",
                "Super_Category_G32A": "Active_Compensation",
                "Pattern_R403C":    "Natural_improvement" if i < 2 else "Progressive",
                "Confidence_R403C": "High",
                "Super_Category_R403C": "Passive",
                "change_consistency": "Consistent",
                "NES_Early_G32A":    -1.5 - i * 0.1,
                "NES_Early_R403C":   -1.2 - i * 0.1,
                "NES_TrajDev_G32A":   1.2 + i * 0.05,
                "NES_TrajDev_R403C":  0.8 - i * 0.05,
                "NES_Late_G32A":     -0.4 + i * 0.05,
                "NES_Late_R403C":    -1.8 + i * 0.1,
                "NES_Time_Ctrl":      0.3,
                "NES_Time_G32A":      0.2,
                "NES_Time_R403C":     0.1,
                "p.adjust_Early_G32A":    0.01,
                "p.adjust_Early_R403C":   0.02,
                "p.adjust_TrajDev_G32A":  0.03,
                "p.adjust_TrajDev_R403C": 0.9,
                "p.adjust_Late_G32A":     0.9,
                "p.adjust_Late_R403C":    0.01,
                "p.adjust_Time_Ctrl":     0.5,
                "p.adjust_Time_G32A":     0.4,
                "p.adjust_Time_R403C":    0.3,
                "ever_significant":             True,
                "ever_significant_trajectory":  True,
                "Sig_TrajDev_G32A":             True,
                "Sig_TrajDev_R403C":            False,
            })

    return pd.DataFrame(rows)


def _focused_filter(df, scope):
    """Stub scope filter that always returns the full deduplicated set."""
    return df.drop_duplicates(subset=["pathway_id"]).copy()


def _make_service(n: int = 5) -> DashboardDataService:
    stub_df = _make_stub_df(n)
    cfg = DashboardConfig(scopes=["focused"])
    return DashboardDataService(
        config=cfg,
        raw_loader=lambda: stub_df,
        scope_filter=_focused_filter,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestDashboardDataService:

    def test_build_returns_tuple(self):
        records, meta = _make_service().build()
        assert isinstance(records, list)
        assert isinstance(meta, DashboardMetadata)

    def test_record_count_matches_unique_pathways(self):
        n = 5
        records, _ = _make_service(n).build()
        assert len(records) == n

    def test_records_are_pathway_record_instances(self):
        records, _ = _make_service().build()
        for r in records:
            assert isinstance(r, PathwayRecord)

    def test_pathway_ids_unique(self):
        records, _ = _make_service(5).build()
        ids = [r.pathway_id for r in records]
        assert len(ids) == len(set(ids))

    def test_mutation_stats_populated(self):
        records, _ = _make_service().build()
        for r in records:
            # G32A should have NES data
            assert r.g32a.nes_early is not None
            assert r.g32a.nes_late is not None

    def test_padj_columns_populated(self):
        records, _ = _make_service().build()
        for r in records:
            assert r.g32a.padj_early is not None

    def test_rank_early_populated_and_positive(self):
        records, _ = _make_service(5).build()
        for r in records:
            assert r.g32a.rank_early is not None
            assert r.g32a.rank_early >= 1.0

    def test_ranks_are_unique_when_nes_values_unique(self):
        """All NES values in stub are distinct → all ranks should differ."""
        records, _ = _make_service(5).build()
        early_ranks = [r.g32a.rank_early for r in records]
        assert len(set(early_ranks)) == len(early_ranks)

    def test_metadata_has_databases(self):
        _, meta = _make_service().build()
        assert len(meta.databases) > 0
        assert all(isinstance(d, str) for d in meta.databases)

    def test_metadata_databases_sorted(self):
        _, meta = _make_service().build()
        assert meta.databases == sorted(meta.databases)

    def test_metadata_weight_categories_non_empty(self):
        _, meta = _make_service().build()
        assert len(meta.weight_categories) > 0

    def test_metadata_weight_categories_valid_values(self):
        from Python.bump_dashboard.domain.schema import WeightCategory
        _, meta = _make_service().build()
        valid = {"dominant", "common", "uncommon", "rare"}
        for cat in meta.weight_categories.values():
            assert cat in valid

    def test_metadata_pattern_colors_present(self):
        _, meta = _make_service().build()
        # Colors must come from the canonical source (pattern_definitions.py)
        assert len(meta.pattern_colors) > 0

    def test_sig_trajdev_set_correctly(self):
        records, _ = _make_service(5).build()
        # Stub sets Sig_TrajDev_G32A = True for all
        for r in records:
            assert r.g32a.sig_trajdev is True
            assert r.r403c.sig_trajdev is False

    def test_single_pathway_service(self):
        records, meta = _make_service(1).build()
        assert len(records) == 1

    def test_missing_pattern_pathway_skipped(self):
        """Pathways where both patterns are NaN must be skipped."""
        df = _make_stub_df(3)
        # Corrupt one pathway by removing patterns for both mutations
        pid_to_corrupt = "DB::pathway_000"
        df.loc[df["pathway_id"] == pid_to_corrupt, ["Pattern_G32A", "Pattern_R403C"]] = np.nan

        svc = DashboardDataService(
            config=DashboardConfig(scopes=["focused"]),
            raw_loader=lambda: df,
            scope_filter=_focused_filter,
        )
        records, _ = svc.build()
        ids = [r.pathway_id for r in records]
        assert pid_to_corrupt not in ids


# ---------------------------------------------------------------------------
# Driver-label loader — unit tests
# ---------------------------------------------------------------------------
#
# Driver classification now lives in 02_Analysis/1.6.gsva_analysis.R and is
# materialised as Driver_G32A / Driver_R403C columns on master_gsva_all_table.csv.
# The dashboard service merely reads them via _default_driver_loader().

class TestDriverLabelLoader:
    """Verify that _load_driver_labels dedupes, coerces, and tolerates absence."""

    def _service(self, driver_df):
        return DashboardDataService(
            config=DashboardConfig(scopes=["focused"]),
            raw_loader=lambda: _make_stub_df(2),
            scope_filter=_focused_filter,
            driver_loader=lambda: driver_df,
        )

    def test_dedupes_per_pathway(self):
        """Master table is cell-level; loader collapses to one entry per pathway."""
        df = pd.DataFrame({
            "pathway_id":   ["P::a", "P::a", "P::b", "P::b"],
            "Driver_G32A":  ["ctrl_driven", "ctrl_driven", "mutant_driven", "mutant_driven"],
            "Driver_R403C": ["both_moving", "both_moving", "neither_moving", "neither_moving"],
        })
        svc = self._service(df)
        out = svc._load_driver_labels()
        assert out == {
            "P::a": {"G32A": "ctrl_driven",   "R403C": "both_moving"},
            "P::b": {"G32A": "mutant_driven", "R403C": "neither_moving"},
        }

    def test_nan_coerced_to_none(self):
        """NaN driver values (missing arm in upstream) must yield None, not 'nan'."""
        df = pd.DataFrame({
            "pathway_id":   ["P::a"],
            "Driver_G32A":  [np.nan],
            "Driver_R403C": ["ctrl_driven"],
        })
        out = self._service(df)._load_driver_labels()
        assert out["P::a"] == {"G32A": None, "R403C": "ctrl_driven"}

    def test_loader_returning_none_yields_empty_dict(self):
        """Missing master CSV → empty dict, dashboard degrades gracefully."""
        assert self._service(None)._load_driver_labels() == {}

    def test_missing_columns_logs_and_returns_empty(self, caplog):
        """If master table lacks Driver_* columns, return {} with a warning."""
        df = pd.DataFrame({"pathway_id": ["P::a"]})  # no Driver_* columns
        with caplog.at_level("WARNING"):
            out = self._service(df)._load_driver_labels()
        assert out == {}
        assert any("Driver_" in r.message for r in caplog.records)
