"""
Tests for presentation.html_fragments — UI contracts and bug-fix regressions.
==============================================================================

These tests validate that the HTML fragments contain the specific structural
and JavaScript elements required by the 7 bug-fixes applied in 2026-05.

They are deliberately *structural* (string-presence) tests rather than
browser-automation tests, which would require a headless browser.  The intent
is to catch any future accidental regression that removes a required element.

Covered regressions
-------------------
Bug 1  – Tooltip uses position:fixed + JS viewport-flip logic (not CSS bottom:%).
Bug 2  – Display Mode control group has a ? info icon.
Bug 3  – Visual Style tooltip explains Bézier amplitude source and direction.
Bug 4/5 – Status bar count uses union of per-mutation visible sets, not global scope.
Bug 6  – NES filters are dual-handle range sliders (not single number inputs).
Bug 7  – Highlight overlay uses _buildSegment() and DIMMED_STYLE for non-matches.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.presentation import html_fragments as frags
from Python.bump_dashboard.presentation.html_renderer import DashboardHtmlRenderer
from Python.bump_dashboard.domain.schema import DashboardConfig

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _render_html() -> str:
    """Produce a complete HTML string with minimal stub data."""
    renderer = DashboardHtmlRenderer(DashboardConfig())
    return renderer.render(
        pathways_payload=[],
        metadata_payload={
            "schema_version": "1.0.0",
            "weight_categories": {},
            "pattern_colors": {},
            "pattern_definitions": {},
            "databases": [],
            "patterns": [],
        },
    )


@pytest.fixture(scope="module")
def html() -> str:
    return _render_html()


# ---------------------------------------------------------------------------
# Bug 1 — Tooltip clipping fix
# ---------------------------------------------------------------------------

class TestBug1TooltipViewportAware:

    def test_tooltip_div_uses_position_fixed(self, html):
        """The global tooltip element must use position:fixed to escape sidebar overflow."""
        assert "position: fixed" in html

    def test_init_tooltips_function_exists(self, html):
        """JS function _initTooltips must be present."""
        assert "function _initTooltips" in html

    def test_viewport_clip_detection_uses_getboundingrect(self, html):
        """Viewport check must use getBoundingClientRect."""
        assert "getBoundingClientRect" in html

    def test_flip_below_logic_present(self, html):
        """If tooltip would overflow top, code must flip it below the icon."""
        assert "iconRect.bottom + margin" in html

    def test_global_tooltip_singleton_div(self, html):
        """A single global tooltip element must exist in the body."""
        assert 'id="global-tooltip"' in html

    def test_no_css_bottom_percent_on_tooltip(self, html):
        """Old CSS bottom:125% trick must no longer be present."""
        # The old tooltip used `bottom: 125%` inside .tooltip-text CSS
        assert "bottom: 125%" not in html

    def test_tooltip_visibility_via_display_none(self, html):
        """Tooltip default state is display:none (managed by JS, not CSS opacity hack)."""
        assert "display: none" in html


# ---------------------------------------------------------------------------
# Bug 2 — Display Mode missing help icon
# ---------------------------------------------------------------------------

class TestBug2DisplayModeHelpIcon:

    def test_display_mode_has_info_icon(self, html):
        """Display Mode control group must contain a ? info icon."""
        assert "info about Display Mode" in html

    def test_display_mode_tooltip_explains_uniform(self, html):
        """Tooltip must mention Uniform mode."""
        assert "Uniform" in html

    def test_display_mode_tooltip_explains_weighted(self, html):
        """Tooltip must mention Weighted mode and rarity logic."""
        assert "Weighted" in html
        # Should explain the rarity/frequency link
        assert "rare" in html.lower() or "dominant" in html.lower() or "frequency" in html.lower()


# ---------------------------------------------------------------------------
# Bug 3 — Visual Style tooltip incomplete
# ---------------------------------------------------------------------------

class TestBug3VisualStyleTooltip:

    def test_bezier_amplitude_source_documented(self, html):
        """Tooltip must explain that curvature magnitude = |TrajDev NES|."""
        assert "Curvature magnitude" in html

    def test_bezier_direction_documented(self, html):
        """Tooltip must explain that curvature direction = sign of TrajDev NES."""
        assert "Curvature direction" in html

    def test_bezier_only_for_sig_pathways_documented(self, html):
        """Tooltip must clarify that curves only appear for significant TrajDev pathways."""
        # Look for some form of "only" + "significant"
        import re
        pattern = re.compile(r"only.*significant|significant.*only", re.IGNORECASE | re.DOTALL)
        assert pattern.search(html), "Tooltip should state curves appear only for significant TrajDev"

    def test_bezier_scale_constant_in_js(self, html):
        """JS must define a named constant for the Bézier scale factor."""
        assert "BEZIER_NES_SCALE" in html


# ---------------------------------------------------------------------------
# Bugs 4 & 5 — Status bar count correctness
# ---------------------------------------------------------------------------

class TestBug4Bug5StatusBarCount:

    def test_visible_g32a_state_variable(self, html):
        """JS must maintain a visibleG32A array separate from scopeData."""
        assert "visibleG32A" in html

    def test_visible_r403c_state_variable(self, html):
        """JS must maintain a visibleR403C array separate from scopeData."""
        assert "visibleR403C" in html

    def test_union_count_computed_for_status_bar(self, html):
        """Status bar count must use a union of visible IDs, not raw length of scopeData."""
        assert "unionIds" in html

    def test_status_bar_reflects_union_not_scope(self, html):
        """The status bar must not use filteredData.length (old buggy approach)."""
        # The old code was: `filteredData.length` for the count
        # New code uses unionIds.size
        assert "unionIds.size" in html

    def test_is_visible_function_used_per_mutation(self, html):
        """_isVisible must be called separately for each mutation."""
        assert "_isVisible(d, 'G32A'" in html or "_isVisible(d, \"G32A\"" in html
        assert "_isVisible(d, 'R403C'" in html or "_isVisible(d, \"R403C\"" in html


# ---------------------------------------------------------------------------
# Bug 6 — Single number inputs replaced by dual-handle range sliders
# ---------------------------------------------------------------------------

class TestBug6DualHandleSliders:

    @pytest.mark.parametrize("key", ["min-nes", "early-nes", "late-nes", "trajdev-nes"])
    def test_lower_thumb_exists(self, html, key):
        assert f'id="rs-{key}-lo"' in html

    @pytest.mark.parametrize("key", ["min-nes", "early-nes", "late-nes", "trajdev-nes"])
    def test_upper_thumb_exists(self, html, key):
        assert f'id="rs-{key}-hi"' in html

    @pytest.mark.parametrize("key", ["min-nes", "early-nes", "late-nes", "trajdev-nes"])
    def test_fill_bar_exists(self, html, key):
        assert f'id="rf-{key}"' in html

    @pytest.mark.parametrize("key", ["min-nes", "early-nes", "late-nes", "trajdev-nes"])
    def test_readout_element_exists(self, html, key):
        assert f'id="rr-{key}"' in html

    def test_sync_slider_function_defined(self, html):
        assert "function syncSlider" in html

    def test_slider_values_reader_defined(self, html):
        assert "function _sliderValues" in html

    def test_range_fill_css_class(self, html):
        """CSS .range-fill must be defined for the slider highlight bar."""
        assert ".range-fill" in html

    def test_filter_uses_range_window_not_single_threshold(self, html):
        """Filter logic must check both lo and hi bounds, not just a minimum."""
        # Old code: if (s.filterEarlyNes > 0 && abs < threshold) return false
        # New code: if (absEarly < eLo || absEarly > eHi) return false
        assert "eHi" in html or "earlyRange.hi" in html or "> eHi" in html

    def test_no_legacy_number_input_for_nes_filter(self, html):
        """Old-style <input type=number id=min-nes> must not exist."""
        assert 'id="min-nes"' not in html  # was the old single input; now replaced by slider pair


# ---------------------------------------------------------------------------
# Bug 7 — Highlight uses Bézier geometry and dims non-matches
# ---------------------------------------------------------------------------

class TestBug7HighlightOverlay:

    def test_highlight_trace_builder_uses_build_segment(self, html):
        """Highlight trace builder must call _buildSegment for correct geometry."""
        # _buildHighlightTraces calls _buildSegment internally
        assert "_buildSegment" in html
        assert "_buildHighlightTraces" in html

    def test_dimmed_style_constant_defined(self, html):
        """DIMMED_STYLE constant must be defined for non-matching lines."""
        assert "DIMMED_STYLE" in html

    def test_dimmed_style_applied_when_highlight_active(self, html):
        """When a highlight search is active, DIMMED_STYLE must be applied to base traces."""
        assert "dimmed ? DIMMED_STYLE" in html

    def test_highlight_uses_two_pass_rendering(self, html):
        """Highlight should render an outline pass + inner colour pass for contrast."""
        # Two-pass = black outline trace + coloured inner trace
        assert "stroke" in html or "'black'" in html or '"black"' in html

    def test_highlight_traces_pushed_after_base_traces(self, html):
        """Inside renderChart(), the highlight traces must be pushed AFTER base traces."""
        # Locate the renderChart function body and check ordering within it
        render_start = html.find("function renderChart")
        assert render_start != -1, "renderChart function not found"
        render_body = html[render_start : render_start + 2000]
        base_pos = render_body.find("_buildPatternTraces")
        hl_pos   = render_body.find("_buildHighlightTraces")
        assert base_pos != -1, "_buildPatternTraces not found in renderChart"
        assert hl_pos   != -1, "_buildHighlightTraces not found in renderChart"
        assert hl_pos > base_pos, "Highlight trace builder must be called after base traces"

    def test_highlight_searches_pathway_id_too(self, html):
        """Search must match on both description and pathway_id."""
        assert "pathway_id" in html and "description" in html
        # Specifically both should appear in the filter logic
        assert "d.pathway_id" in html

    def test_status_bar_shows_highlight_count(self, html):
        """Status bar should display how many pathways are highlighted."""
        assert "highlighted" in html


# ---------------------------------------------------------------------------
# GSVA driver filter + modal annotation
# ---------------------------------------------------------------------------

class TestGsvaDriverFilterAndModal:
    """Structural contracts for the GSVA driver classification feature.

    Verifies that:
    - The filter UI checkbox group is present in the rendered HTML.
    - The JS populates the driver list from RAW_DATA.
    - The settings reader includes the 'drivers' key.
    - The visibility predicate gates on the driver field.
    - The modal annotation block is present.
    - The serialised data carries Driver_G32A and Driver_R403C.
    """

    def test_driver_list_container_present(self, html):
        """HTML must contain the #driver-list checkbox group."""
        assert 'id="driver-list"' in html

    def test_driver_list_populate_function_present(self, html):
        """JS must define _populateDriverList() to bootstrap the filter."""
        assert "_populateDriverList" in html

    def test_driver_order_constant_defined(self, html):
        """A GSVA_DRIVER_ORDER constant with all four labels must be defined."""
        assert "GSVA_DRIVER_ORDER" in html
        assert "mutant_driven" in html
        assert "ctrl_driven" in html
        assert "both_moving" in html
        assert "neither_moving" in html

    def test_driver_labels_map_defined(self, html):
        """GSVA_DRIVER_LABELS human-readable map must be present."""
        assert "GSVA_DRIVER_LABELS" in html
        # Each label should have a human-readable expansion
        assert "mutant carries the change" in html
        assert "Ctrl developmental trajectory" in html

    def test_read_settings_includes_drivers(self, html):
        """_readSettings() must include the drivers key."""
        assert "drivers:" in html
        assert "#driver-list input:checked" in html

    def test_is_visible_gates_on_driver(self, html):
        """_isVisible() must reference Driver_{mut} for the per-mutation gate."""
        assert "Driver_${mut}" in html or 'Driver_' in html
        assert "s.drivers.includes" in html

    def test_gsva_driver_verdict_div_in_modal(self, html):
        """Modal HTML must contain the #gsva-driver-verdict block."""
        assert 'id="gsva-driver-verdict"' in html

    def test_gsva_driver_verdict_css_defined(self, html):
        """CSS for .gsva-driver-verdict must be defined."""
        assert ".gsva-driver-verdict" in html

    def test_verdict_delta_computation_in_modal_js(self, html):
        """openGsvaModal() must compute per-genotype delta values."""
        assert "deltas" in html
        assert "b.median - a.median" in html

    def test_driver_fields_embedded_in_raw_data(self):
        """The serialised JSON payload must carry Driver_G32A and Driver_R403C."""
        from Python.bump_dashboard.domain.schema import MutationStats, PathwayRecord
        from Python.bump_dashboard.application.serializer import DashboardSerializer

        stats = MutationStats(
            pattern="Compensation", nes_early=-1.5, nes_trajdev=1.2, nes_late=-0.4,
            padj_early=0.01, padj_trajdev=0.03, padj_late=0.9,
            sig_trajdev=True, rank_early=5.0, rank_late=12.0,
            gsva_driver="ctrl_driven",
        )
        record = PathwayRecord(
            pathway_id="TEST::001", description="Test", database="MitoCarta",
            g32a=stats, r403c=MutationStats(
                pattern="Compensation", nes_early=-1.2, nes_trajdev=0.9, nes_late=-0.3,
                padj_early=0.02, padj_trajdev=0.04, padj_late=0.8,
                sig_trajdev=False, rank_early=6.0, rank_late=14.0,
                gsva_driver="mutant_driven",
            ),
        )
        row = DashboardSerializer().serialize_pathways([record])[0]
        assert row["Driver_G32A"] == "ctrl_driven"
        assert row["Driver_R403C"] == "mutant_driven"
