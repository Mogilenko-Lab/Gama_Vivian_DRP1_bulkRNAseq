"""
application.data_service
========================

``DashboardDataService`` is the application-layer boundary between the
upstream data pipeline (viz_bump_charts / master_gsea_table) and the
dashboard domain objects.

Responsibilities
----------------
1. Load raw data via the existing infrastructure adapter.
2. Enrich the raw frame with per-contrast p-values (pivot step).
3. Filter to the configured scope superset.
4. Compute per-pattern weight categories.
5. Compute per-mutation Early/Late ranks.
6. Validate and assemble ``PathwayRecord`` objects.
7. Assemble ``DashboardMetadata``.

This class is intentionally *not* a god-object: each responsibility is
isolated in a private method so it can be tested or replaced independently.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

from ..domain.schema import (
    DashboardConfig,
    DashboardMetadata,
    GsvaSampleMeta,
    MutationStats,
    PathwayRecord,
    SCHEMA_VERSION,
)
from ..domain.rules import (
    PADJ_THRESHOLD,
    classify_weight_category,
    compute_ranks_descending,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Contrast → column name mapping (single source of truth inside this module)
# ---------------------------------------------------------------------------

_CONTRAST_TO_PADJ_COL: dict[str, str] = {
    "G32A_vs_Ctrl_D35": "padj_Early_G32A",
    "G32A_vs_Ctrl_D65": "padj_Late_G32A",
    "Maturation_G32A_specific": "padj_TrajDev_G32A",
    "R403C_vs_Ctrl_D35": "padj_Early_R403C",
    "R403C_vs_Ctrl_D65": "padj_Late_R403C",
    "Maturation_R403C_specific": "padj_TrajDev_R403C",
}

_MUTATIONS: tuple[str, ...] = ("G32A", "R403C")


class DashboardDataService:
    """
    Loads, validates, and transforms data into dashboard domain objects.

    Parameters
    ----------
    config:
        Dashboard configuration.  Defaults to ``DashboardConfig()`` if omitted.
    raw_loader:
        Callable that returns the raw wide-format DataFrame.  Defaults to the
        project's ``viz_bump_charts.load_data()``.  Inject a custom callable
        in tests to avoid hitting the filesystem.
    scope_filter:
        Callable ``(df, scope) → df`` that filters rows by scope.  Defaults to
        ``viz_bump_charts.filter_by_scope()``.
    """

    def __init__(
        self,
        config: Optional[DashboardConfig] = None,
        raw_loader=None,
        scope_filter=None,
        gsva_long_loader=None,
    ) -> None:
        self._config = config or DashboardConfig()
        self._raw_loader = raw_loader or self._default_loader()
        self._scope_filter = scope_filter or self._default_scope_filter()
        # gsva_long_loader returns a long-format DataFrame with columns
        # (pathway_id, sample_id, genotype, day, gsva_score) — None disables.
        self._gsva_long_loader = (
            gsva_long_loader if gsva_long_loader is not None
            else self._default_gsva_long_loader
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> tuple[list[PathwayRecord], DashboardMetadata]:
        """
        Execute the full data preparation pipeline.

        Returns
        -------
        tuple[list[PathwayRecord], DashboardMetadata]
            Validated pathway records and accompanying metadata.
        """
        logger.info("DashboardDataService: loading raw data…")
        raw_df = self._raw_loader()

        logger.info("DashboardDataService: enriching with per-contrast p-values…")
        enriched_df = self._enrich_pvalues(raw_df)

        logger.info("DashboardDataService: filtering by scope(s) %s…", self._config.scopes)
        filtered_df = self._filter_scope_union(enriched_df)

        logger.info("DashboardDataService: %d pathways after scope filter.", len(filtered_df))

        logger.info("DashboardDataService: computing weight categories…")
        weight_cats = self._compute_weight_categories(filtered_df)

        logger.info("DashboardDataService: computing ranks…")
        ranked_df = self._add_ranks(filtered_df)

        logger.info("DashboardDataService: attaching per-sample GSVA scores…")
        gsva_sample_index, gsva_lookup = self._load_gsva()

        logger.info("DashboardDataService: assembling domain objects…")
        records = self._assemble_records(ranked_df, gsva_lookup, gsva_sample_index)

        metadata = self._assemble_metadata(ranked_df, weight_cats, gsva_sample_index)

        logger.info(
            "DashboardDataService: done. %d valid records, schema_version=%s.",
            len(records),
            SCHEMA_VERSION,
        )
        return records, metadata

    # ------------------------------------------------------------------
    # Private helpers — each maps to one discrete responsibility
    # ------------------------------------------------------------------

    @staticmethod
    def _default_loader():
        """Return the project's standard raw-data loader function."""
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Python.viz_bump_charts import load_data  # type: ignore[import]

        return load_data

    @staticmethod
    def _default_scope_filter():
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Python.viz_bump_charts import filter_by_scope  # type: ignore[import]

        return filter_by_scope

    @staticmethod
    def _default_gsva_long_loader() -> Optional[pd.DataFrame]:
        """Default loader for the per-sample GSVA long-format CSV.

        Returns the parsed DataFrame, or None if the file is missing — in
        which case the dashboard will simply lack the per-replicate panel
        (degrades gracefully rather than failing the whole pipeline).
        """
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Python.config import resolve_path  # type: ignore[import]

        csv_path = Path(resolve_path("03_Results/02_Analysis/replicate_level_gsva_long.csv"))
        if not csv_path.exists():
            logger.warning(
                "Per-sample GSVA long CSV not found at %s — dashboard will be "
                "rendered without the replicate-level click-through panel.",
                csv_path,
            )
            return None
        return pd.read_csv(csv_path)

    def _load_gsva(self) -> tuple[list[GsvaSampleMeta], dict[str, list[Optional[float]]]]:
        """Load per-sample GSVA scores and build the dashboard payload.

        Returns
        -------
        gsva_sample_index : ordered list of (sample_id, genotype, day) tuples
            shared across all pathways. The order is:
                day asc (D35, D65), genotype canonical (Ctrl, G32A, R403C),
                sample_id asc.
        gsva_lookup : dict mapping pathway_id → list of GSVA scores aligned
            positionally with gsva_sample_index. Missing pathways are absent
            from the dict, and missing per-sample scores are encoded as None.
        """
        try:
            df = self._gsva_long_loader() if callable(self._gsva_long_loader) else self._gsva_long_loader
        except Exception as exc:  # noqa: BLE001 — soft-fail keeps pipeline usable
            logger.warning("GSVA long-format loader failed (%s); skipping per-sample panel.", exc)
            return [], {}

        if df is None or len(df) == 0:
            return [], {}

        required = {"pathway_id", "sample_id", "genotype", "day", "gsva_score"}
        missing = required - set(df.columns)
        if missing:
            logger.warning(
                "GSVA long CSV is missing columns %s; skipping per-sample panel.",
                missing,
            )
            return [], {}

        genotype_rank = {"Ctrl": 0, "G32A": 1, "R403C": 2}
        day_rank      = {"D35": 0, "D65": 1}

        samples_df = (
            df[["sample_id", "genotype", "day"]]
            .drop_duplicates()
            .assign(
                _g=lambda d: d["genotype"].map(genotype_rank).fillna(99),
                _d=lambda d: d["day"].map(day_rank).fillna(99),
            )
            .sort_values(["_d", "_g", "sample_id"])
            .drop(columns=["_g", "_d"])
            .reset_index(drop=True)
        )

        sample_index: list[GsvaSampleMeta] = [
            GsvaSampleMeta(sample_id=row.sample_id, genotype=row.genotype, day=row.day)
            for row in samples_df.itertuples(index=False)
        ]
        column_order = [s.sample_id for s in sample_index]

        wide = df.pivot_table(
            index="pathway_id",
            columns="sample_id",
            values="gsva_score",
            aggfunc="first",
        ).reindex(columns=column_order)

        import math
        lookup: dict[str, list[Optional[float]]] = {}
        for pid, row in wide.iterrows():
            scores: list[Optional[float]] = []
            for v in row.values.tolist():
                if v is None:
                    scores.append(None)
                    continue
                try:
                    f = float(v)
                except (TypeError, ValueError):
                    scores.append(None)
                    continue
                scores.append(None if math.isnan(f) else f)
            lookup[str(pid)] = scores

        logger.info(
            "GSVA: %d pathways × %d samples loaded for click-through panel.",
            len(lookup),
            len(sample_index),
        )
        return sample_index, lookup

    def _enrich_pvalues(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Pivot the per-contrast p.adjust values into wide format and join back.

        The master table stores one row per (pathway_id, contrast).  We pivot
        only the six contrasts we care about and join to a deduplicated base.
        """
        relevant = df[df["contrast"].isin(_CONTRAST_TO_PADJ_COL.keys())].copy()

        pivoted = relevant.pivot(
            index="pathway_id",
            columns="contrast",
            values="p.adjust",
        ).rename(columns=_CONTRAST_TO_PADJ_COL)

        # Deduplicated base (all per-pathway wide columns already present)
        base = df.drop_duplicates(subset=["pathway_id"]).copy()
        result = base.set_index("pathway_id").join(pivoted).reset_index()
        return result

    def _filter_scope_union(self, df: pd.DataFrame) -> pd.DataFrame:
        """Union-merge rows from all configured scopes, deduplicate."""
        frames = [
            self._scope_filter(df, scope) for scope in self._config.scopes
        ]
        merged = pd.concat(frames).drop_duplicates(subset=["pathway_id"])
        return merged.reset_index(drop=True)

    @staticmethod
    def _compute_weight_categories(df: pd.DataFrame) -> dict[str, str]:
        """
        Compute frequency-based weight category for each pattern globally.

        Counts are accumulated across both mutations so that the categories
        reflect the full dataset, not a per-mutation view.
        """
        pattern_counts: dict[str, int] = {}
        for mut in _MUTATIONS:
            col = f"Pattern_{mut}"
            if col not in df.columns:
                continue
            nes_cols = [f"NES_Early_{mut}", f"NES_Late_{mut}"]
            valid = df[df[nes_cols].notna().all(axis=1)]
            for pat, cnt in valid[col].value_counts().items():
                pattern_counts[pat] = pattern_counts.get(pat, 0) + int(cnt)

        total = sum(pattern_counts.values())
        if total == 0:
            return {}

        return {
            pat: classify_weight_category(cnt / total)
            for pat, cnt in pattern_counts.items()
        }

    @staticmethod
    def _add_ranks(df: pd.DataFrame) -> pd.DataFrame:
        """Append Early/Late rank columns for each mutation."""
        df = df.copy()
        for mut in _MUTATIONS:
            for stage in ("Early", "Late"):
                nes_col = f"NES_{stage}_{mut}"
                rank_col = f"Rank_{stage}_{mut}"
                if nes_col in df.columns:
                    ranks = compute_ranks_descending(df[nes_col].tolist())
                    df[rank_col] = ranks
        return df

    @staticmethod
    def _safe_float(row: pd.Series, col: str) -> float | None:
        """Extract a float from *row[col]*, returning None on NaN / missing."""
        import math

        if col not in row.index:
            return None
        v = row[col]
        if v is None:
            return None
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return None if math.isnan(f) else f

    GSVA_DRIVER_EPSILON: float = 0.10

    @staticmethod
    def _classify_gsva_driver(
        scores: Optional[list[Optional[float]]],
        sample_index: list[GsvaSampleMeta],
        mutation: str,
        eps: float = 0.10,
    ) -> Optional[str]:
        """Classify a pathway's GSVA trajectory driver for one mutation.

        Compares per-arm Δ_GSVA = median(D65) − median(D35) for Ctrl and the
        named mutation. Returns one of {"mutant_driven", "ctrl_driven",
        "both_moving", "neither_moving"} when all four (Ctrl/M × D35/D65)
        cells have at least one valid score; otherwise None.
        """
        import math
        from statistics import median

        if scores is None or not sample_index:
            return None
        if len(scores) != len(sample_index):
            return None

        cells: dict[tuple[str, str], list[float]] = {}
        for s, v in zip(sample_index, scores):
            if v is None:
                continue
            try:
                f = float(v)
            except (TypeError, ValueError):
                continue
            if math.isnan(f):
                continue
            cells.setdefault((s.genotype, s.day), []).append(f)

        needed = [
            ("Ctrl", "D35"), ("Ctrl", "D65"),
            (mutation, "D35"), (mutation, "D65"),
        ]
        if any(not cells.get(k) for k in needed):
            return None

        delta_ctrl = median(cells[("Ctrl", "D65")]) - median(cells[("Ctrl", "D35")])
        delta_mut  = median(cells[(mutation, "D65")]) - median(cells[(mutation, "D35")])

        moving_ctrl = abs(delta_ctrl) >= eps
        moving_mut  = abs(delta_mut)  >= eps
        if moving_ctrl and moving_mut:
            return "both_moving"
        if moving_mut:
            return "mutant_driven"
        if moving_ctrl:
            return "ctrl_driven"
        return "neither_moving"

    def _build_mutation_stats(
        self,
        row: pd.Series,
        mut: str,
        gsva_scores: Optional[list[Optional[float]]] = None,
        sample_index: Optional[list[GsvaSampleMeta]] = None,
    ) -> MutationStats:
        """Construct a ``MutationStats`` from one DataFrame row and mutation."""
        sf = self._safe_float

        raw_pattern = row.get(f"Pattern_{mut}")
        pattern = raw_pattern if isinstance(raw_pattern, str) and raw_pattern else None

        driver = self._classify_gsva_driver(
            gsva_scores, sample_index or [], mut, self.GSVA_DRIVER_EPSILON
        )

        return MutationStats(
            pattern=pattern,  # type: ignore[arg-type]
            nes_early=sf(row, f"NES_Early_{mut}"),
            nes_trajdev=sf(row, f"NES_TrajDev_{mut}"),
            nes_late=sf(row, f"NES_Late_{mut}"),
            padj_early=sf(row, f"padj_Early_{mut}"),
            padj_trajdev=sf(row, f"padj_TrajDev_{mut}"),
            padj_late=sf(row, f"padj_Late_{mut}"),
            sig_trajdev=bool(row.get(f"Sig_TrajDev_{mut}", False)),
            rank_early=sf(row, f"Rank_Early_{mut}"),
            rank_late=sf(row, f"Rank_Late_{mut}"),
            gsva_driver=driver,  # type: ignore[arg-type]
        )

    def _assemble_records(
        self,
        df: pd.DataFrame,
        gsva_lookup: Optional[dict[str, list[Optional[float]]]] = None,
        sample_index: Optional[list[GsvaSampleMeta]] = None,
    ) -> list[PathwayRecord]:
        """Build and validate ``PathwayRecord`` objects from each DataFrame row."""
        records: list[PathwayRecord] = []
        skipped = 0
        gsva_lookup = gsva_lookup or {}
        sample_index = sample_index or []

        for _, row in df.iterrows():
            try:
                pid = str(row["pathway_id"])
                scores = gsva_lookup.get(pid)
                record = PathwayRecord(
                    pathway_id=pid,
                    description=str(row.get("Description", row["pathway_id"])),
                    database=str(row.get("database", "unknown")),
                    g32a=self._build_mutation_stats(row, "G32A", scores, sample_index),
                    r403c=self._build_mutation_stats(row, "R403C", scores, sample_index),
                    gsva_scores=scores,
                )
                records.append(record)
            except Exception as exc:  # noqa: BLE001
                skipped += 1
                logger.debug("Skipping pathway %r: %s", row.get("pathway_id"), exc)

        if skipped:
            logger.warning(
                "Skipped %d pathways that failed validation (see DEBUG logs).",
                skipped,
            )

        attached = sum(1 for r in records if r.gsva_scores is not None)
        logger.info(
            "GSVA attached to %d / %d records (%.1f%%).",
            attached,
            len(records),
            (100.0 * attached / max(1, len(records))),
        )

        # Driver-label distribution (per-mutation tally across the dashboard scope).
        from collections import Counter
        for mut, attr in (("G32A", "g32a"), ("R403C", "r403c")):
            tally = Counter(
                (getattr(r, attr).gsva_driver or "unavailable") for r in records
            )
            ordered = ", ".join(
                f"{k}={tally.get(k, 0)}"
                for k in ("mutant_driven", "ctrl_driven", "both_moving",
                          "neither_moving", "unavailable")
            )
            logger.info("GSVA driver labels [%s]: %s", mut, ordered)
        return records

    @staticmethod
    def _assemble_metadata(
        df: pd.DataFrame,
        weight_cats: dict[str, str],
        gsva_sample_index: Optional[list[GsvaSampleMeta]] = None,
    ) -> DashboardMetadata:
        """Build ``DashboardMetadata`` from the filtered DataFrame."""
        # Import pattern colors from the canonical source
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Python.pattern_definitions import (  # type: ignore[import]
            get_pattern_colors,
            PATTERN_DEFINITIONS,
        )

        pattern_colors = get_pattern_colors()
        pattern_defs = {
            k: v.get("interpretation", "")
            for k, v in PATTERN_DEFINITIONS.items()
        }
        databases = sorted(df["database"].dropna().unique().tolist())

        # Use G32A patterns as the reference list for UI checkboxes
        patterns = sorted(
            p
            for p in df["Pattern_G32A"].dropna().unique()
            if isinstance(p, str) and p
        )

        return DashboardMetadata(
            weight_categories=weight_cats,
            pattern_colors=pattern_colors,
            pattern_definitions=pattern_defs,
            databases=databases,
            patterns=patterns,
            gsva_sample_index=gsva_sample_index or [],
        )
