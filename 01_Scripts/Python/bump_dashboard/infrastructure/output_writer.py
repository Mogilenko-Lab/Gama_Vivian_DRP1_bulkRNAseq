"""
infrastructure.output_writer
============================

``DashboardOutputWriter`` writes the rendered HTML string to the filesystem.

This is the sole I/O write-concern in the module: all other layers are
pure-Python transformations.  Isolating writes here makes it trivial to:
  * Test the renderer and data pipeline without touching the disk.
  * Swap the write target (e.g. S3, stdout) by replacing this class.

Path resolution delegates to the project's ``config.resolve_path()`` so the
module remains independent of the current working directory.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from ..domain.schema import DashboardConfig

logger = logging.getLogger(__name__)


class DashboardOutputWriter:
    """
    Writes the dashboard HTML to the configured output path.

    Parameters
    ----------
    config:
        Dashboard configuration.  Uses ``output_dir`` and ``output_filename``.
    project_root:
        Optional override for the project root used in path resolution.
        When ``None``, deferred to ``config.resolve_path()`` at write time.
    """

    def __init__(
        self,
        config: Optional[DashboardConfig] = None,
        project_root: Optional[Path] = None,
    ) -> None:
        self._config = config or DashboardConfig()
        self._project_root = project_root  # None → resolve at write time

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write(self, html: str) -> str:
        """
        Write *html* to the configured output file.

        Parameters
        ----------
        html:
            Rendered HTML string.

        Returns
        -------
        str
            Absolute path string of the written file.

        Raises
        ------
        OSError
            On any filesystem error during directory creation or file write.
        """
        output_path = self._resolve_output_path()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")
        logger.info("OutputWriter: wrote %d bytes → %s", len(html), output_path)
        return str(output_path)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _resolve_output_path(self) -> Path:
        """Resolve the full output path from config."""
        if self._project_root is not None:
            root = self._project_root
        else:
            root = self._get_project_root()

        return root / self._config.output_dir / self._config.output_filename

    @staticmethod
    def _get_project_root() -> Path:
        """Locate the project root by delegating to the project config."""
        import sys

        sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
        from Python.config import get_project_root  # type: ignore[import]

        return get_project_root()
