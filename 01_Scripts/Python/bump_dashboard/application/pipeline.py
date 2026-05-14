"""
application.pipeline
====================

``DashboardPipeline`` is the single entry point for the bump-dashboard
generation workflow.

It coordinates three application services:
  * ``DashboardDataService``   – data loading + domain assembly
  * ``DashboardSerializer``    – domain → JSON-safe dicts
  * ``DashboardHtmlRenderer``  – dicts → HTML string
  * ``DashboardOutputWriter``  – HTML string → file on disk

Dependency injection is used for every collaborator so the pipeline can be
run in tests with stubs and driven by the real implementations in production.

Execution trace
---------------
::

    pipeline.run()
        → data_service.build()
            → loader() → enrich → filter → weight_cats → ranks → validate
        → serializer.serialize_pathways() + serialize_metadata()
        → renderer.render()
        → output_writer.write()
"""

from __future__ import annotations

import logging
from typing import Optional

from ..domain.schema import DashboardConfig
from .data_service import DashboardDataService
from .serializer import DashboardSerializer
from ..presentation.html_renderer import DashboardHtmlRenderer
from ..infrastructure.output_writer import DashboardOutputWriter

logger = logging.getLogger(__name__)


class DashboardPipeline:
    """
    Orchestrates the end-to-end bump-dashboard generation.

    Parameters
    ----------
    config:
        Dashboard configuration.  Defaults to ``DashboardConfig()``.
    data_service:
        Optional override for the data-preparation step.
    serializer:
        Optional override for the serialisation step.
    renderer:
        Optional override for the HTML-rendering step.
    output_writer:
        Optional override for the file-write step.

    Examples
    --------
    Typical production use::

        from bump_dashboard import DashboardPipeline
        pipeline = DashboardPipeline()
        output_path = pipeline.run()

    With a custom config::

        from bump_dashboard import DashboardPipeline
        from bump_dashboard.domain.schema import DashboardConfig
        cfg = DashboardConfig(scopes=["focused"])
        pipeline = DashboardPipeline(config=cfg)
        output_path = pipeline.run()
    """

    def __init__(
        self,
        config: Optional[DashboardConfig] = None,
        data_service: Optional[DashboardDataService] = None,
        serializer: Optional[DashboardSerializer] = None,
        renderer: Optional[DashboardHtmlRenderer] = None,
        output_writer: Optional[DashboardOutputWriter] = None,
    ) -> None:
        self._config = config or DashboardConfig()
        self._data_service = data_service or DashboardDataService(config=self._config)
        self._serializer = serializer or DashboardSerializer()
        self._renderer = renderer or DashboardHtmlRenderer(config=self._config)
        self._output_writer = output_writer or DashboardOutputWriter(config=self._config)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> str:
        """
        Execute the full pipeline.

        Returns
        -------
        str
            Absolute path to the written HTML file.
        """
        logger.info("=" * 70)
        logger.info("BUMP DASHBOARD PIPELINE — START")
        logger.info("=" * 70)

        # 1. Build domain objects
        records, metadata = self._data_service.build()
        logger.info("Pipeline: %d pathway records loaded.", len(records))

        # 2. Serialise
        pathways_payload = self._serializer.serialize_pathways(records)
        metadata_payload = self._serializer.serialize_metadata(metadata)
        logger.info("Pipeline: serialisation complete.")

        # 3. Render HTML
        html = self._renderer.render(pathways_payload, metadata_payload)
        logger.info("Pipeline: HTML rendered (%d bytes).", len(html))

        # 4. Write to disk
        output_path = self._output_writer.write(html)
        logger.info("Pipeline: output written → %s", output_path)
        logger.info("=" * 70)
        logger.info("BUMP DASHBOARD PIPELINE — DONE")
        logger.info("=" * 70)

        return output_path
