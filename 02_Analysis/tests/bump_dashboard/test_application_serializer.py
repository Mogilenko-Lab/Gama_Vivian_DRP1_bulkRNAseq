"""
Tests for application.serializer
=================================

Validates the shape and content of the JSON-ready payloads.
No I/O; uses in-memory domain objects.
"""

import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.schema import (
    DashboardMetadata,
    MutationStats,
    PathwayRecord,
    SCHEMA_VERSION,
)
from Python.bump_dashboard.application.serializer import DashboardSerializer


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _stats(pattern="Compensation", **kw) -> MutationStats:
    defaults = dict(
        pattern=pattern, nes_early=-1.5, nes_trajdev=1.2, nes_late=-0.4,
        padj_early=0.01, padj_trajdev=0.03, padj_late=0.9,
        sig_trajdev=True, rank_early=5.0, rank_late=12.0,
    )
    defaults.update(kw)
    return MutationStats(**defaults)


def _record(pid="PW001", **kw) -> PathwayRecord:
    defaults = dict(
        pathway_id=pid,
        description="Oxidative phosphorylation",
        database="MitoCarta",
        g32a=_stats("Compensation"),
        r403c=_stats("Progressive"),
    )
    defaults.update(kw)
    return PathwayRecord(**defaults)


def _metadata() -> DashboardMetadata:
    return DashboardMetadata(
        weight_categories={"Compensation": "common", "Progressive": "rare"},
        pattern_colors={"Compensation": "#009E73", "Progressive": "#D55E00"},
        pattern_definitions={"Compensation": "Active adaptive response"},
        databases=["MitoCarta", "SynGO"],
        patterns=["Compensation", "Progressive"],
    )


# ---------------------------------------------------------------------------
# Serializer tests
# ---------------------------------------------------------------------------

class TestDashboardSerializer:

    def setup_method(self):
        self.s = DashboardSerializer()

    # --- serialize_pathways ------------------------------------------------

    def test_returns_list(self):
        result = self.s.serialize_pathways([_record()])
        assert isinstance(result, list)

    def test_one_dict_per_record(self):
        records = [_record("A"), _record("B"), _record("C")]
        result = self.s.serialize_pathways(records)
        assert len(result) == 3

    def test_base_fields_present(self):
        row = self.s.serialize_pathways([_record()])[0]
        assert "pathway_id" in row
        assert "description" in row
        assert "database" in row

    def test_mutation_fields_present_for_both_mutations(self):
        row = self.s.serialize_pathways([_record()])[0]
        for mut in ("G32A", "R403C"):
            assert f"Pattern_{mut}" in row
            assert f"NES_Early_{mut}" in row
            assert f"NES_TrajDev_{mut}" in row
            assert f"NES_Late_{mut}" in row
            assert f"padj_Early_{mut}" in row
            assert f"padj_TrajDev_{mut}" in row
            assert f"padj_Late_{mut}" in row
            assert f"Rank_Early_{mut}" in row
            assert f"Rank_Late_{mut}" in row
            assert f"Sig_TrajDev_{mut}" in row

    def test_total_fields_count(self):
        """3 base + 11 per mutation × 2 + 1 gsva_scores = 26 fields (schema 1.2: +Driver_*)."""
        row = self.s.serialize_pathways([_record()])[0]
        assert len(row) == 26
        assert "gsva_scores" in row

    def test_driver_fields_present_for_both_mutations(self):
        """Driver_{mut} fields must be present and default to None when no GSVA scores."""
        row = self.s.serialize_pathways([_record()])[0]
        for mut in ("G32A", "R403C"):
            assert f"Driver_{mut}" in row
            # No GSVA scores on the stub record → driver must be None.
            assert row[f"Driver_{mut}"] is None

    def test_none_pattern_serialised_as_none(self):
        row = self.s.serialize_pathways([_record(g32a=_stats(pattern=None))])[0]
        assert row["Pattern_G32A"] is None

    def test_sig_trajdev_is_bool(self):
        row = self.s.serialize_pathways([_record()])[0]
        assert isinstance(row["Sig_TrajDev_G32A"], bool)

    def test_nes_values_are_floats(self):
        row = self.s.serialize_pathways([_record()])[0]
        assert isinstance(row["NES_Early_G32A"], float)

    def test_pathway_id_matches_input(self):
        row = self.s.serialize_pathways([_record("MY::PATH")])[0]
        assert row["pathway_id"] == "MY::PATH"

    def test_empty_list_returns_empty_list(self):
        assert self.s.serialize_pathways([]) == []

    # --- serialize_metadata ------------------------------------------------

    def test_metadata_schema_version_present(self):
        result = self.s.serialize_metadata(_metadata())
        assert result["schema_version"] == SCHEMA_VERSION

    def test_metadata_databases_list(self):
        result = self.s.serialize_metadata(_metadata())
        assert result["databases"] == ["MitoCarta", "SynGO"]

    def test_metadata_pattern_colors_dict(self):
        result = self.s.serialize_metadata(_metadata())
        assert result["pattern_colors"]["Compensation"] == "#009E73"

    def test_metadata_is_json_serialisable(self):
        import json
        result = self.s.serialize_metadata(_metadata())
        # Should not raise
        json.dumps(result)

    def test_pathway_payload_is_json_serialisable(self):
        import json
        result = self.s.serialize_pathways([_record()])
        json.dumps(result)
