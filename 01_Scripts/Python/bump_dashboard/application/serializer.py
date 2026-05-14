"""
application.serializer
======================

``DashboardSerializer`` converts validated domain objects into plain Python
dicts / lists that can be JSON-serialised by the renderer.

This is a *pure transformation* step — no I/O, no validation, no rendering.
Keeping serialisation here (rather than inside the renderer) makes it trivial
to unit-test the output shape and to swap the renderer without touching data
logic.

Output shape (single pathway dict)
-----------------------------------
::

    {
        "pathway_id": "...",
        "description": "...",
        "database": "...",
        # --- per-mutation fields (G32A example) ---
        "Pattern_G32A":         "Compensation" | null,
        "NES_Early_G32A":       1.23 | null,
        "NES_TrajDev_G32A":     0.45 | null,
        "NES_Late_G32A":        0.11 | null,
        "padj_Early_G32A":      0.001 | null,
        "padj_TrajDev_G32A":    0.032 | null,
        "padj_Late_G32A":       0.98  | null,
        "Rank_Early_G32A":      3.0   | null,
        "Rank_Late_G32A":       7.0   | null,
        "Sig_TrajDev_G32A":     true | false,
        # ... same for R403C ...
    }

Metadata dict shape
-------------------
See ``DashboardMetadata`` model — the serialiser emits it verbatim via
``.model_dump()``.
"""

from __future__ import annotations

from typing import Any

from ..domain.schema import DashboardMetadata, MutationStats, PathwayRecord

# Mutation keys in the order the JS client expects them
_MUTATIONS = ("G32A", "R403C")


class DashboardSerializer:
    """
    Converts ``PathwayRecord`` instances and ``DashboardMetadata`` into plain
    JSON-serialisable Python structures.

    Usage
    -----
    ::

        serializer = DashboardSerializer()
        pathways_json_payload = serializer.serialize_pathways(records)
        metadata_json_payload = serializer.serialize_metadata(metadata)
    """

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def serialize_pathways(
        self,
        records: list[PathwayRecord],
    ) -> list[dict[str, Any]]:
        """
        Serialize a list of ``PathwayRecord`` objects.

        Parameters
        ----------
        records:
            Validated domain pathway records.

        Returns
        -------
        list[dict[str, Any]]
            One flat dict per record; all values are JSON-safe.
        """
        return [self._serialize_record(r) for r in records]

    def serialize_metadata(
        self,
        metadata: DashboardMetadata,
    ) -> dict[str, Any]:
        """
        Serialize ``DashboardMetadata`` to a plain dict.

        Returns
        -------
        dict[str, Any]
            JSON-safe metadata dict.
        """
        # model_dump() is Pydantic v2's preferred method
        return metadata.model_dump()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _mut_stats_to_dict(
        stats: MutationStats,
        mutation: str,
    ) -> dict[str, Any]:
        """
        Flatten one ``MutationStats`` value object into dashboard field names.

        The key naming convention (``Pattern_G32A``, ``NES_Early_G32A``, …)
        deliberately mirrors the JS client's expectations — changing it here
        requires a corresponding update in the HTML template.
        """
        return {
            f"Pattern_{mutation}": stats.pattern,
            f"NES_Early_{mutation}": stats.nes_early,
            f"NES_TrajDev_{mutation}": stats.nes_trajdev,
            f"NES_Late_{mutation}": stats.nes_late,
            f"padj_Early_{mutation}": stats.padj_early,
            f"padj_TrajDev_{mutation}": stats.padj_trajdev,
            f"padj_Late_{mutation}": stats.padj_late,
            f"Rank_Early_{mutation}": stats.rank_early,
            f"Rank_Late_{mutation}": stats.rank_late,
            f"Sig_TrajDev_{mutation}": stats.sig_trajdev,
            f"Driver_{mutation}": stats.gsva_driver,
        }

    def _serialize_record(self, record: PathwayRecord) -> dict[str, Any]:
        """Flatten a single ``PathwayRecord`` to a JSON-ready dict."""
        base: dict[str, Any] = {
            "pathway_id": record.pathway_id,
            "description": record.description,
            "database": record.database,
        }
        base.update(self._mut_stats_to_dict(record.g32a, "G32A"))
        base.update(self._mut_stats_to_dict(record.r403c, "R403C"))
        # Per-sample GSVA scores (None when pathway has no GSVA data).
        # Aligned positionally with metadata.gsva_sample_index.
        base["gsva_scores"] = record.gsva_scores
        return base
