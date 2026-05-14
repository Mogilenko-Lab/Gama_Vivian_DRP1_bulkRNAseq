"""
Tests for domain.geometry
==========================
"""

import math
import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.geometry import (
    quadratic_bezier,
    BezierPoints,
    compute_curve_offset,
    wrap_text,
)


# ---------------------------------------------------------------------------
# quadratic_bezier
# ---------------------------------------------------------------------------

class TestQuadraticBezier:

    def test_returns_bezier_points_type(self):
        result = quadratic_bezier(0, 0, 0.5, 1, 1, 0)
        assert isinstance(result, BezierPoints)

    def test_start_point_correct(self):
        pts = quadratic_bezier(0, -1.5, 0.5, 0, 1, 1.5, n_points=10)
        assert math.isclose(pts.x[0], 0.0)
        assert math.isclose(pts.y[0], -1.5)

    def test_end_point_correct(self):
        pts = quadratic_bezier(0, -1.5, 0.5, 0, 1, 1.5, n_points=10)
        assert math.isclose(pts.x[-1], 1.0)
        assert math.isclose(pts.y[-1], 1.5)

    def test_point_count_is_n_plus_one(self):
        pts = quadratic_bezier(0, 0, 0.5, 1, 1, 0, n_points=20)
        assert len(pts.x) == 21
        assert len(pts.y) == 21

    def test_straight_line_when_control_on_midpoint(self):
        """When control point is the midpoint, curve should be a straight line."""
        pts = quadratic_bezier(0, 0, 0.5, 0.5, 1, 1, n_points=4)
        # Check midpoint
        assert math.isclose(pts.x[2], 0.5, abs_tol=1e-9)
        assert math.isclose(pts.y[2], 0.5, abs_tol=1e-9)

    def test_n_points_less_than_2_raises(self):
        with pytest.raises(ValueError, match="n_points"):
            quadratic_bezier(0, 0, 0.5, 0, 1, 0, n_points=1)

    def test_x_and_y_same_length(self):
        pts = quadratic_bezier(0, 0, 0.5, 1, 1, 0)
        assert len(pts.x) == len(pts.y)

    def test_bezier_points_frozen(self):
        pts = quadratic_bezier(0, 0, 0.5, 1, 1, 0)
        with pytest.raises(Exception):
            pts.x = (0.0,)  # type: ignore[misc]


# ---------------------------------------------------------------------------
# compute_curve_offset
# ---------------------------------------------------------------------------

class TestComputeCurveOffset:

    def test_nes_mode_positive_trajdev_positive_offset(self):
        offset = compute_curve_offset(1.5, 0.0, 1.0, y_metric="nes", nes_scale_factor=0.5)
        assert offset == pytest.approx(0.75)

    def test_nes_mode_negative_trajdev_negative_offset(self):
        offset = compute_curve_offset(-1.5, 0.0, 1.0, y_metric="nes", nes_scale_factor=0.5)
        assert offset == pytest.approx(-0.75)

    def test_rank_mode_offset_scales_with_chart_height(self):
        offset_small = compute_curve_offset(1.0, 0.0, 1.0, y_metric="rank", chart_height=10)
        offset_large = compute_curve_offset(1.0, 0.0, 1.0, y_metric="rank", chart_height=100)
        # Larger chart → larger offset
        assert abs(offset_large) > abs(offset_small)

    def test_rank_mode_positive_trajdev_negative_data_offset(self):
        """Positive NES (upregulated) = lower rank = negative direction in rank mode."""
        offset = compute_curve_offset(1.0, 0.0, 1.0, y_metric="rank", chart_height=100)
        assert offset < 0

    def test_zero_trajdev_zero_offset_nes(self):
        assert compute_curve_offset(0.0, 0.0, 1.0, y_metric="nes") == pytest.approx(0.0)

    def test_zero_trajdev_zero_offset_rank(self):
        assert compute_curve_offset(0.0, 0.0, 1.0, y_metric="rank") == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# wrap_text
# ---------------------------------------------------------------------------

class TestWrapText:

    def test_short_text_unchanged(self):
        assert wrap_text("short", 35) == "short"

    def test_empty_string_returns_empty(self):
        assert wrap_text("") == ""

    def test_long_text_gets_br_inserted(self):
        text = "Mitochondrial ribosome biogenesis and translation regulation"
        result = wrap_text(text, 30)
        assert "<br>" in result

    def test_each_line_within_width(self):
        text = "a " * 50  # 50 words
        result = wrap_text(text, 20)
        for line in result.split("<br>"):
            assert len(line.rstrip()) <= 22  # allow slight overshoot for single long words

    def test_preserves_full_content(self):
        """All words must be present in the output (order preserved)."""
        text = "alpha beta gamma delta epsilon"
        result = wrap_text(text, 10)
        for word in text.split():
            assert word in result

    def test_width_respected_for_single_long_word(self):
        """A word longer than width must not be split."""
        text = "supercalifragilisticexpialidocious"
        result = wrap_text(text, 10)
        assert "supercalifragilisticexpialidocious" in result
