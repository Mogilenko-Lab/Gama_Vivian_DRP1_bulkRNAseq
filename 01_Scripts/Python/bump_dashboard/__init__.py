"""
bump_dashboard
==============

Self-contained module for generating the interactive pathway-trajectory bump-chart
dashboard (HTML / Plotly).

Public surface
--------------
    DashboardPipeline          – top-level orchestrator (entry point)
    DashboardDataService       – data loading + enrichment
    DashboardSerializer        – produces JSON-ready plain dicts
    DashboardHtmlRenderer      – assembles the final HTML artefact
    DashboardOutputWriter      – writes the HTML file

Sub-packages
------------
    domain/     – pure-Python value objects and domain rules (no I/O)
    application/ – use-case orchestration (thin, coordinates domain + infra)
    infrastructure/ – I/O adapters (data loading, file writing)
    presentation/  – HTML / JS generation

Design notes
------------
* Domain objects are *anemic* value objects (dataclasses / NamedTuples) that
  carry data; all business rules live in explicit service/validator functions.
* No circular imports: domain ← application ← infrastructure/presentation.
* Framework: pydantic v2 for validated value objects; standard-library
  dataclasses for lightweight containers where pydantic overhead is unnecessary.
"""

from .application.pipeline import DashboardPipeline
from .application.data_service import DashboardDataService
from .application.serializer import DashboardSerializer
from .presentation.html_renderer import DashboardHtmlRenderer
from .infrastructure.output_writer import DashboardOutputWriter

__all__ = [
    "DashboardPipeline",
    "DashboardDataService",
    "DashboardSerializer",
    "DashboardHtmlRenderer",
    "DashboardOutputWriter",
]
