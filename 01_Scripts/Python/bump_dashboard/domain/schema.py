"""
domain.schema
=============

Canonical data contracts for the bump-dashboard module.

These Pydantic models serve as the validated boundary between upstream
data loaders and the dashboard rendering pipeline.  Every row that enters
the pipeline is validated against ``PathwayRecord``; every config dict is
validated against ``DashboardConfig``.

Schema versioning
-----------------
``SCHEMA_VERSION`` must be bumped when a *breaking* field is added or
removed.  Additive-only changes (new optional fields with defaults) do
**not** require a version bump.

Backward-compatibility note
---------------------------
The field names deliberately mirror the column names that already exist in
the upstream ``master_gsea_table.csv`` / the wide-format pathway table
produced by ``viz_bump_charts.load_data()``.  This avoids a translation
layer for the current data shape while still enforcing types and allowing
future schema evolution.
"""

from __future__ import annotations

from typing import Literal, Optional
from pydantic import BaseModel, Field, field_validator, model_validator

# ---------------------------------------------------------------------------
# Schema version – increment on breaking changes
# ---------------------------------------------------------------------------

SCHEMA_VERSION: str = "1.1.0"

# ---------------------------------------------------------------------------
# Controlled vocabularies
# ---------------------------------------------------------------------------

MutationId = Literal["G32A", "R403C"]

PatternName = Literal[
    "Compensation",
    "Sign_reversal",
    "Progressive",
    "Natural_worsening",
    "Natural_improvement",
    "Late_onset",
    "Transient",
    "Complex",
    "Insufficient_data",
]

GsvaDriver = Literal[
    "mutant_driven",
    "ctrl_driven",
    "both_moving",
    "neither_moving",
]

WeightCategory = Literal["dominant", "common", "uncommon", "rare"]

ColorByMode = Literal["pattern", "nes_early", "nes_late", "nes_trajdev"]

YAxisMetric = Literal["nes", "rank"]

LineDisplayMode = Literal["uniform", "weighted"]


# ---------------------------------------------------------------------------
# Per-mutation trajectory statistics (one per mutation per pathway)
# ---------------------------------------------------------------------------

class MutationStats(BaseModel):
    """
    All trajectory statistics for a single mutation × pathway combination.

    Fields map 1-to-1 with the wide-format columns already present in the
    master table; no renaming is required by callers.
    """

    pattern: Optional[PatternName] = Field(
        None,
        description="Trajectory pattern classification for this mutation.",
    )
    nes_early: Optional[float] = Field(
        None, 
        description="NES at Early stage (Mutation_vs_Ctrl_D35)."
        )
    nes_trajdev: Optional[float] = Field(
        None, 
        description="NES for TrajDev (Maturation_Mutation_specific)."
        )
    nes_late: Optional[float] = Field(
        None, 
        description="NES at Late stage (Mutation_vs_Ctrl_D65)."
        )
    padj_early: Optional[float] = Field(
        None, ge=0.0, le=1.0, 
        description="p.adjust for Early stage.")
    padj_trajdev: Optional[float] = Field(
        None, ge=0.0, le=1.0, 
        description="p.adjust for TrajDev."
        )
    padj_late: Optional[float] = Field(
        None, ge=0.0, le=1.0, 
        description="p.adjust for Late stage."
        )
    sig_trajdev: bool = Field(
        False, 
        description="True when TrajDev is significant (padj < 0.05).")
    rank_early: Optional[float] = Field(
        None, 
        description="NES-based rank at Early stage (1 = highest NES)."
        )
    rank_late: Optional[float] = Field(
        None,
        description="NES-based rank at Late stage."
        )
    gsva_driver: Optional[GsvaDriver] = Field(
        None,
        description=(
            "GSVA-derived driver classification. Compares per-arm Δ_GSVA "
            "(median D65 − D35) for Ctrl vs this mutation; tells whether "
            "the contrast-level pattern is carried by the mutant arm, the "
            "Ctrl arm, both, or neither. None when GSVA data is unavailable."
        ),
    )

    model_config = {"frozen": True}


# ---------------------------------------------------------------------------
# Per-sample GSVA score metadata (Phase 2b: replicate-level click-through)
# ---------------------------------------------------------------------------

class GsvaSampleMeta(BaseModel):
    """Identity of one GSVA sample column.

    Sample columns are stored once in ``DashboardMetadata.gsva_sample_index``;
    each ``PathwayRecord`` then carries a positionally-aligned float array of
    GSVA scores under ``gsva_scores``. This avoids repeating sample metadata
    across thousands of pathway records in the embedded JSON.
    """

    sample_id: str = Field(..., min_length=1, description="Original library sample id.")
    genotype: Literal["Ctrl", "G32A", "R403C"] = Field(..., description="Genotype group.")
    day: Literal["D35", "D65"] = Field(..., description="Timepoint.")

    model_config = {"frozen": True}


# ---------------------------------------------------------------------------
# Core pathway record (one row in the dashboard data)
# ---------------------------------------------------------------------------

class PathwayRecord(BaseModel):
    """
    Single pathway entry as consumed by the dashboard.

    Validation rules
    ----------------
    * ``pathway_id`` must be non-empty.
    * ``description`` must be non-empty.
    * ``database`` must be non-empty.
    * At least one of ``g32a`` or ``r403c`` must have a non-None ``pattern``.
    """

    pathway_id: str = Field(..., min_length=1, description="Unique pathway identifier.")
    description: str = Field(..., min_length=1, description="Human-readable pathway name.")
    database: str = Field(..., min_length=1, description="Source database (e.g. MitoCarta, SynGO, gobp).")
    g32a: MutationStats = Field(..., description="Trajectory statistics for the G32A mutation.")
    r403c: MutationStats = Field(..., description="Trajectory statistics for the R403C mutation.")
    gsva_scores: Optional[list[Optional[float]]] = Field(
        None,
        description=(
            "Per-sample GSVA enrichment scores, positionally aligned with "
            "DashboardMetadata.gsva_sample_index. None when the pathway is "
            "absent from the GSVA universe (e.g. failed size filter)."
        ),
    )

    model_config = {"frozen": True}

    @field_validator("pathway_id", "description", "database", mode="before")
    @classmethod
    def strip_and_nonempty(cls, v: object) -> str:
        if not isinstance(v, str):
            v = str(v)
        v = v.strip()
        if not v:
            raise ValueError("Field must be a non-empty string after stripping whitespace.")
        return v

    @model_validator(mode="after")
    def at_least_one_mutation_has_pattern(self) -> "PathwayRecord":
        if self.g32a.pattern is None and self.r403c.pattern is None:
            raise ValueError(
                f"Pathway '{self.pathway_id}' has no pattern for either mutation. "
                "It should not be included in the dashboard data."
            )
        return self


# ---------------------------------------------------------------------------
# Dashboard metadata (supporting lookup tables, not per-pathway)
# ---------------------------------------------------------------------------

class DashboardMetadata(BaseModel):
    """
    Metadata passed alongside pathway records to the renderer.

    All dicts use string keys for safe JSON serialisation.
    """

    schema_version: str = Field(SCHEMA_VERSION, description="Data-contract version.")
    weight_categories: dict[str, WeightCategory] = Field(
        default_factory=dict,
        description="Mapping pattern_name → weight category.",
    )
    pattern_colors: dict[str, str] = Field(
        default_factory=dict,
        description="Mapping pattern_name → hex color string.",
    )
    pattern_definitions: dict[str, str] = Field(
        default_factory=dict,
        description="Mapping pattern_name → short interpretation text.",
    )
    databases: list[str] = Field(
        default_factory=list,
        description="Sorted list of database names present in the data.",
    )
    patterns: list[str] = Field(
        default_factory=list,
        description="Sorted list of pattern names present in the data (G32A as reference).",
    )
    gsva_sample_index: list[GsvaSampleMeta] = Field(
        default_factory=list,
        description=(
            "Ordered list of GSVA sample columns (sample_id, genotype, day). "
            "Each PathwayRecord.gsva_scores array is positionally aligned with "
            "this index. Empty list means GSVA data is unavailable."
        ),
    )

    model_config = {"frozen": True}


# ---------------------------------------------------------------------------
# Dashboard render configuration
# ---------------------------------------------------------------------------

class DashboardConfig(BaseModel):
    """
    Runtime configuration for the dashboard generator.

    All values have sensible defaults so callers only need to override what
    they care about.
    """

    plotly_js_url: str = Field(
        "https://cdn.plot.ly/plotly-2.27.0.min.js",
        description="CDN URL for Plotly JS.",
    )
    output_dir: str = Field(
        "03_Results/02_Analysis/Plots/Trajectory_Flow",
        description="Output directory for the HTML file (relative to project root).",
    )
    output_filename: str = Field(
        "interactive_bump_dashboard.html",
        description="Output filename.",
    )
    # Data scope: which scopes to include in the superset
    scopes: list[Literal["focused", "significant"]] = Field(
        default_factory=lambda: ["focused", "significant"],
        description="Data scopes to union-merge into the dashboard data.",
    )

    model_config = {"frozen": True}
