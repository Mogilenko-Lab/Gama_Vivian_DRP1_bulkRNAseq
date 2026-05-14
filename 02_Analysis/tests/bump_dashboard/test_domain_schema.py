"""
Tests for domain.schema
=======================

Covers Pydantic validation on ``PathwayRecord``, ``MutationStats``, and
``DashboardConfig``.  No I/O; no dependency on the real data files.
"""

import pytest
from pydantic import ValidationError

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.schema import (
    GsvaDriver,
    MutationStats,
    PathwayRecord,
    DashboardConfig,
    DashboardMetadata,
    SCHEMA_VERSION,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_stats(**overrides) -> MutationStats:
    defaults = dict(
        pattern="Compensation",
        nes_early=-1.5,
        nes_trajdev=1.2,
        nes_late=-0.4,
        padj_early=0.01,
        padj_trajdev=0.03,
        padj_late=0.9,
        sig_trajdev=True,
        rank_early=5.0,
        rank_late=12.0,
    )
    defaults.update(overrides)
    return MutationStats(**defaults)


def _make_record(**overrides) -> PathwayRecord:
    defaults = dict(
        pathway_id="TEST::pathway_001",
        description="Test pathway description",
        database="MitoCarta",
        g32a=_make_stats(),
        r403c=_make_stats(pattern="Progressive"),
    )
    defaults.update(overrides)
    return PathwayRecord(**defaults)


# ---------------------------------------------------------------------------
# MutationStats
# ---------------------------------------------------------------------------

class TestMutationStats:

    def test_valid_stats_constructs(self):
        stats = _make_stats()
        assert stats.pattern == "Compensation"
        assert stats.nes_early == pytest.approx(-1.5)
        assert stats.sig_trajdev is True

    def test_all_none_values_allowed(self):
        """A stats with no NES data is valid — it represents missing data."""
        stats = MutationStats(
            pattern=None, nes_early=None, nes_trajdev=None, nes_late=None,
            padj_early=None, padj_trajdev=None, padj_late=None,
            sig_trajdev=False, rank_early=None, rank_late=None,
        )
        assert stats.pattern is None

    def test_invalid_pattern_rejected(self):
        with pytest.raises(ValidationError, match="pattern"):
            _make_stats(pattern="NotAPattern")

    def test_padj_out_of_range_rejected(self):
        with pytest.raises(ValidationError):
            _make_stats(padj_early=1.5)  # > 1.0

    def test_padj_negative_rejected(self):
        with pytest.raises(ValidationError):
            _make_stats(padj_early=-0.01)  # < 0.0

    def test_frozen_mutation_stats(self):
        stats = _make_stats()
        with pytest.raises(Exception):  # frozen model → AttributeError / ValidationError
            stats.nes_early = 0.0  # type: ignore[misc]


# ---------------------------------------------------------------------------
# GsvaDriver field on MutationStats
# ---------------------------------------------------------------------------

class TestGsvaDriverField:
    """gsva_driver is an optional controlled-vocabulary field on MutationStats."""

    @pytest.mark.parametrize("label", ["mutant_driven", "ctrl_driven", "both_moving", "neither_moving"])
    def test_valid_driver_label_round_trips(self, label):
        """Each of the four valid labels must be accepted and round-trip correctly."""
        stats = _make_stats(gsva_driver=label)
        assert stats.gsva_driver == label

    def test_driver_none_is_valid(self):
        """None (absent GSVA data) must be accepted without error."""
        stats = _make_stats(gsva_driver=None)
        assert stats.gsva_driver is None

    def test_invalid_driver_label_rejected(self):
        """An unrecognised string must raise a ValidationError."""
        with pytest.raises(ValidationError):
            _make_stats(gsva_driver="both_flat")  # not in the Literal


# ---------------------------------------------------------------------------
# PathwayRecord
# ---------------------------------------------------------------------------

class TestPathwayRecord:

    def test_valid_record_constructs(self):
        record = _make_record()
        assert record.pathway_id == "TEST::pathway_001"
        assert record.g32a.pattern == "Compensation"

    def test_empty_pathway_id_rejected(self):
        with pytest.raises(ValidationError):
            _make_record(pathway_id="")

    def test_whitespace_only_description_rejected(self):
        with pytest.raises(ValidationError):
            _make_record(description="   ")

    def test_empty_database_rejected(self):
        with pytest.raises(ValidationError):
            _make_record(database="")

    def test_both_patterns_none_rejected(self):
        """A pathway with no pattern in either mutation is not usable."""
        with pytest.raises(ValidationError, match="no pattern"):
            _make_record(
                g32a=_make_stats(pattern=None),
                r403c=_make_stats(pattern=None),
            )

    def test_one_pattern_none_allowed(self):
        """One mutation can have no pattern as long as the other does."""
        record = _make_record(g32a=_make_stats(pattern=None))
        assert record.g32a.pattern is None
        assert record.r403c.pattern == "Progressive"

    def test_pathway_id_stripped(self):
        record = _make_record(pathway_id="  TEST::pathway  ")
        assert record.pathway_id == "TEST::pathway"

    def test_frozen_record(self):
        record = _make_record()
        with pytest.raises(Exception):
            record.database = "SynGO"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# DashboardConfig
# ---------------------------------------------------------------------------

class TestDashboardConfig:

    def test_default_config(self):
        cfg = DashboardConfig()
        assert "focused" in cfg.scopes
        assert "significant" in cfg.scopes
        assert "plotly" in cfg.plotly_js_url

    def test_custom_scopes_accepted(self):
        cfg = DashboardConfig(scopes=["focused"])
        assert cfg.scopes == ["focused"]

    def test_invalid_scope_rejected(self):
        with pytest.raises(ValidationError):
            DashboardConfig(scopes=["nonexistent"])  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# DashboardMetadata
# ---------------------------------------------------------------------------

class TestDashboardMetadata:

    def test_schema_version_present(self):
        meta = DashboardMetadata()
        assert meta.schema_version == SCHEMA_VERSION

    def test_empty_metadata_valid(self):
        meta = DashboardMetadata()
        assert meta.databases == []
        assert meta.patterns == []
