"""
presentation.html_renderer
==========================

``DashboardHtmlRenderer`` assembles the final self-contained HTML file from:
  * Static fragments defined in ``html_fragments``
  * Serialised data payloads (JSON strings) produced by ``DashboardSerializer``

Design principles
-----------------
* **No data transformation** here — all data is received already serialised.
* **No file I/O** here — the renderer returns a string; writing is the
  ``DashboardOutputWriter``'s responsibility.
* **Safe token substitution** replaces ``%%…%%`` sentinels; there is no
  risk of f-string / JS syntax collision.
* **Testable in isolation**: inject any payloads, inspect the HTML string.

Token contract
--------------
The following tokens must appear (exactly once each) in the assembled
template:

  ``%%PLOTLY_JS_URL%%``   – CDN URL for Plotly
  ``%%PATHWAYS_JSON%%``   – JSON array of pathway objects
  ``%%METADATA_JSON%%``   – JSON object with metadata

These are guaranteed to be safe inside a ``<script>`` block because ``%%``
is not a valid JS token prefix.
"""

from __future__ import annotations

import json
import logging
from typing import Any, Optional

from ..domain.schema import DashboardConfig
from . import html_fragments as _frags

# Module-level aliases — tests can monkey-patch these names to inject
# corrupted fragments without touching the original module.
HEAD   = _frags.HEAD
BODY   = _frags.BODY
SCRIPT = _frags.SCRIPT

logger = logging.getLogger(__name__)

# Sentinel tokens — must match their counterparts in html_fragments.SCRIPT
_TOKEN_PLOTLY_URL    = "%%PLOTLY_JS_URL%%"
_TOKEN_PATHWAYS_JSON = "%%PATHWAYS_JSON%%"
_TOKEN_METADATA_JSON = "%%METADATA_JSON%%"


class DashboardHtmlRenderer:
    """
    Renders the dashboard HTML by combining static fragments with live data.

    Parameters
    ----------
    config:
        Dashboard configuration.  Uses ``plotly_js_url``.
    """

    def __init__(self, config: Optional[DashboardConfig] = None) -> None:
        self._config = config or DashboardConfig()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def render(
        self,
        pathways_payload: list[dict[str, Any]],
        metadata_payload: dict[str, Any],
    ) -> str:
        """
        Assemble and return the full self-contained HTML string.

        Parameters
        ----------
        pathways_payload:
            JSON-safe list of pathway dicts (from ``DashboardSerializer``).
        metadata_payload:
            JSON-safe metadata dict (from ``DashboardSerializer``).

        Returns
        -------
        str
            Complete HTML document.

        Raises
        ------
        ValueError
            When a sentinel token is missing from the assembled template.
        """
        pathways_json = self._to_json(pathways_payload)
        metadata_json = self._to_json(metadata_payload)

        # Read module-level names so test monkey-patching takes effect
        import Python.bump_dashboard.presentation.html_renderer as _self_mod
        template = _self_mod.HEAD + _self_mod.BODY + _self_mod.SCRIPT
        self._validate_tokens(template)

        html = template
        html = html.replace(_TOKEN_PLOTLY_URL,    self._config.plotly_js_url)
        html = html.replace(_TOKEN_PATHWAYS_JSON, pathways_json)
        html = html.replace(_TOKEN_METADATA_JSON, metadata_json)

        logger.debug("Renderer: produced %d-byte HTML document.", len(html))
        return html

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_json(obj: Any) -> str:
        """Serialise *obj* to a compact JSON string safe for inline ``<script>`` embedding."""
        # Use ``default=str`` so that any stragglers (e.g. numpy floats, Paths)
        # are stringified rather than raising TypeError.
        return json.dumps(obj, default=str, ensure_ascii=False)

    @staticmethod
    def _validate_tokens(template: str) -> None:
        """Assert all required sentinel tokens are present in *template*."""
        missing = [
            tok for tok in (_TOKEN_PLOTLY_URL, _TOKEN_PATHWAYS_JSON, _TOKEN_METADATA_JSON)
            if tok not in template
        ]
        if missing:
            raise ValueError(
                f"HTML template is missing sentinel tokens: {missing}. "
                "Check html_fragments.py for consistency."
            )
