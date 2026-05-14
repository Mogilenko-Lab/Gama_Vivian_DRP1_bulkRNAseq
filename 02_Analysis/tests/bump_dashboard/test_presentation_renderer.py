"""
Tests for presentation.html_renderer
======================================

Validates HTML assembly: token substitution, completeness, and error handling.
No file I/O; no data pipeline dependency.
"""

import json
import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.schema import DashboardConfig
from Python.bump_dashboard.presentation.html_renderer import DashboardHtmlRenderer


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_SAMPLE_PATHWAYS = [
    {
        "pathway_id": "TEST::001",
        "description": "Test pathway",
        "database": "MitoCarta",
        "Pattern_G32A": "Compensation",
        "NES_Early_G32A": -1.5,
        "NES_TrajDev_G32A": 1.2,
        "NES_Late_G32A": -0.3,
        "padj_Early_G32A": 0.01,
        "padj_TrajDev_G32A": 0.03,
        "padj_Late_G32A": 0.9,
        "Rank_Early_G32A": 3.0,
        "Rank_Late_G32A": 7.0,
        "Sig_TrajDev_G32A": True,
        "Pattern_R403C": "Progressive",
        "NES_Early_R403C": -1.5,
        "NES_TrajDev_R403C": -0.8,
        "NES_Late_R403C": -2.1,
        "padj_Early_R403C": 0.02,
        "padj_TrajDev_R403C": 0.04,
        "padj_Late_R403C": 0.001,
        "Rank_Early_R403C": 5.0,
        "Rank_Late_R403C": 2.0,
        "Sig_TrajDev_R403C": False,
        "Driver_G32A": "ctrl_driven",
        "Driver_R403C": None,
    }
]

_SAMPLE_METADATA = {
    "schema_version": "1.0.0",
    "weight_categories": {"Compensation": "common"},
    "pattern_colors": {"Compensation": "#009E73"},
    "pattern_definitions": {"Compensation": "Active adaptive response"},
    "databases": ["MitoCarta"],
    "patterns": ["Compensation"],
}


# ---------------------------------------------------------------------------
# Renderer tests
# ---------------------------------------------------------------------------

class TestDashboardHtmlRenderer:

    def _renderer(self, url="https://cdn.plot.ly/plotly-test.min.js") -> DashboardHtmlRenderer:
        cfg = DashboardConfig(plotly_js_url=url)
        return DashboardHtmlRenderer(config=cfg)

    def _render(self, renderer=None):
        r = renderer or self._renderer()
        return r.render(_SAMPLE_PATHWAYS, _SAMPLE_METADATA)

    # ── Basic structure ──────────────────────────────────────────────────────

    def test_returns_string(self):
        assert isinstance(self._render(), str)

    def test_html_non_empty(self):
        assert len(self._render()) > 1000  # minimal sanity on size

    def test_starts_with_doctype(self):
        assert self._render().strip().startswith("<!DOCTYPE html>")

    def test_ends_with_html_close(self):
        html = self._render().strip()
        assert html.endswith("</html>")

    def test_has_head_and_body(self):
        html = self._render()
        assert "<head>" in html
        assert "<body>" in html

    # ── Token substitution ───────────────────────────────────────────────────

    def test_plotly_url_substituted(self):
        url = "https://cdn.plot.ly/plotly-test.min.js"
        html = self._render(self._renderer(url=url))
        assert url in html
        assert "%%PLOTLY_JS_URL%%" not in html

    def test_pathways_json_embedded(self):
        html = self._render()
        assert "%%PATHWAYS_JSON%%" not in html
        # The pathway_id should appear in the rendered JSON
        assert "TEST::001" in html

    def test_metadata_json_embedded(self):
        html = self._render()
        assert "%%METADATA_JSON%%" not in html
        assert "Compensation" in html

    def test_no_sentinel_tokens_remain(self):
        html = self._render()
        for token in ("%%PLOTLY_JS_URL%%", "%%PATHWAYS_JSON%%", "%%METADATA_JSON%%"):
            assert token not in html, f"Sentinel token {token!r} was not replaced"

    # ── JSON validity ────────────────────────────────────────────────────────

    def test_pathways_json_parseable(self):
        """The embedded pathways JSON must be syntactically valid."""
        html = self._render()
        # Extract the JSON array from the known assignment pattern
        start = html.index("const RAW_DATA = ") + len("const RAW_DATA = ")
        end = html.index(";", start)
        extracted = html[start:end]
        parsed = json.loads(extracted)
        assert isinstance(parsed, list)
        assert len(parsed) == 1

    def test_metadata_json_parseable(self):
        html = self._render()
        start = html.index("const METADATA  = ") + len("const METADATA  = ")
        end = html.index(";", start)
        extracted = html[start:end]
        parsed = json.loads(extracted)
        assert parsed["schema_version"] == "1.0.0"

    # ── Content checks ───────────────────────────────────────────────────────

    def test_g32a_chart_div_present(self):
        assert 'id="chart-g32a"' in self._render()

    def test_r403c_chart_div_present(self):
        assert 'id="chart-r403c"' in self._render()

    def test_reset_button_present(self):
        assert "resetView" in self._render()

    def test_plotly_script_referenced(self):
        html = self._render()
        assert "<script src=" in html

    # ── Error handling ───────────────────────────────────────────────────────

    def test_missing_token_raises_value_error(self):
        """If a fragment is broken (missing sentinel), renderer must raise."""
        import Python.bump_dashboard.presentation.html_renderer as html_renderer_mod
        from Python.bump_dashboard.presentation.html_renderer import _TOKEN_PATHWAYS_JSON

        renderer = DashboardHtmlRenderer()
        # Patch the SCRIPT name *inside the renderer module* (where it is used)
        original = html_renderer_mod.SCRIPT
        html_renderer_mod.SCRIPT = html_renderer_mod.SCRIPT.replace(_TOKEN_PATHWAYS_JSON, "")
        try:
            with pytest.raises(ValueError, match="sentinel"):
                renderer.render(_SAMPLE_PATHWAYS, _SAMPLE_METADATA)
        finally:
            html_renderer_mod.SCRIPT = original  # restore

    def test_empty_pathways_list_renders_valid_html(self):
        html = self._renderer().render([], _SAMPLE_METADATA)
        assert "<!DOCTYPE html>" in html
        assert "[]" in html  # empty array embedded
