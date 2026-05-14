"""
Tests for domain.rules
=======================

Pure unit tests — no I/O, no data files required.
Each test exercises one business rule in isolation.
"""

import math
import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "01_Scripts"))

from Python.bump_dashboard.domain.rules import (
    classify_weight_category,
    is_significant,
    significance_marker,
    nes_to_diverging_color,
    compute_ranks_descending,
    interpolate_hex_color,
    PADJ_THRESHOLD,
    NES_DIVERGING_MAX,
)


# ---------------------------------------------------------------------------
# classify_weight_category
# ---------------------------------------------------------------------------

class TestClassifyWeightCategory:

    @pytest.mark.parametrize("freq,expected", [
        (0.31,  "dominant"),
        (0.30,  "common"),     # boundary: strictly > 0.30 required; 0.30 is NOT dominant
        (0.29,  "common"),
        (0.11,  "common"),
        (0.10,  "uncommon"),  # boundary: strictly > 0.10 required; 0.10 is NOT common
        (0.09,  "uncommon"),
        (0.011, "uncommon"),
        (0.01,  "rare"),      # boundary: strictly > 0.01 required; 0.01 is NOT uncommon
        (0.009, "rare"),
        (0.0,   "rare"),
    ])
    def test_category_boundaries(self, freq, expected):
        assert classify_weight_category(freq) == expected

    def test_frequency_above_one_raises(self):
        with pytest.raises(ValueError, match="frequency"):
            classify_weight_category(1.1)

    def test_frequency_below_zero_raises(self):
        with pytest.raises(ValueError):
            classify_weight_category(-0.01)

    def test_exact_one_is_dominant(self):
        """100% frequency = dominant (edge case)."""
        assert classify_weight_category(1.0) == "dominant"


# ---------------------------------------------------------------------------
# is_significant / significance_marker
# ---------------------------------------------------------------------------

class TestSignificance:

    def test_below_threshold_is_significant(self):
        assert is_significant(0.04) is True

    def test_at_threshold_is_not_significant(self):
        assert is_significant(0.05) is False

    def test_none_is_not_significant(self):
        assert is_significant(None) is False

    def test_custom_threshold(self):
        assert is_significant(0.09, threshold=0.10) is True
        assert is_significant(0.10, threshold=0.10) is False

    def test_marker_star_for_significant(self):
        assert significance_marker(0.001) == "*"

    def test_marker_empty_for_non_significant(self):
        assert significance_marker(0.99) == ""

    def test_marker_empty_for_none(self):
        assert significance_marker(None) == ""


# ---------------------------------------------------------------------------
# nes_to_diverging_color
# ---------------------------------------------------------------------------

class TestNesToDivergingColor:

    def test_zero_nes_returns_neutral_white(self):
        color = nes_to_diverging_color(0.0)
        assert color.lower() == "#f7f7f7"

    def test_positive_nes_at_max_returns_orange(self):
        color = nes_to_diverging_color(NES_DIVERGING_MAX)
        assert color.lower() == "#b35806"

    def test_negative_nes_at_max_returns_blue(self):
        color = nes_to_diverging_color(-NES_DIVERGING_MAX)
        assert color.lower() == "#2166ac"

    def test_none_returns_grey(self):
        assert nes_to_diverging_color(None) == "#999999"

    def test_nan_returns_grey(self):
        assert nes_to_diverging_color(float("nan")) == "#999999"

    def test_nes_clamped_above_max(self):
        """Values beyond NES_MAX should be clamped to the extremity color."""
        assert nes_to_diverging_color(100.0) == nes_to_diverging_color(NES_DIVERGING_MAX)

    def test_returns_hex_string(self):
        color = nes_to_diverging_color(1.5)
        assert color.startswith("#")
        assert len(color) == 7

    def test_positive_is_warmer_than_zero(self):
        """Positive NES should produce a color with more red/orange."""
        neutral = nes_to_diverging_color(0.0)
        positive = nes_to_diverging_color(2.0)
        assert neutral != positive

    def test_negative_is_cooler_than_zero(self):
        neutral = nes_to_diverging_color(0.0)
        negative = nes_to_diverging_color(-2.0)
        assert neutral != negative


# ---------------------------------------------------------------------------
# interpolate_hex_color
# ---------------------------------------------------------------------------

class TestInterpolateHexColor:

    def test_t_zero_returns_first_color(self):
        result = interpolate_hex_color("#000000", "#ffffff", 0.0)
        assert result == "#000000"

    def test_t_one_returns_second_color(self):
        result = interpolate_hex_color("#000000", "#ffffff", 1.0)
        assert result == "#ffffff"

    def test_t_half_midpoint(self):
        result = interpolate_hex_color("#000000", "#ffffff", 0.5)
        # 0 + (255 - 0) * 0.5 = 127.5; Python round() uses banker's rounding → 128 = 0x80
        assert result == "#808080"


# ---------------------------------------------------------------------------
# compute_ranks_descending
# ---------------------------------------------------------------------------

class TestComputeRanksDescending:

    def test_simple_descending_order(self):
        ranks = compute_ranks_descending([3.0, 1.0, 2.0])
        assert ranks == [1.0, 3.0, 2.0]

    def test_none_values_get_none_rank(self):
        ranks = compute_ranks_descending([2.0, None, 1.0])
        assert ranks[0] == 1.0  # highest value
        assert ranks[1] is None
        assert ranks[2] == 2.0

    def test_all_none_returns_all_none(self):
        ranks = compute_ranks_descending([None, None])
        assert ranks == [None, None]

    def test_single_value_gets_rank_1(self):
        assert compute_ranks_descending([5.0]) == [1.0]

    def test_ties_get_minimum_rank(self):
        """Two equal values both get rank 1 (min method)."""
        ranks = compute_ranks_descending([2.0, 2.0, 1.0])
        assert ranks[0] == 1.0
        assert ranks[1] == 1.0
        assert ranks[2] == 3.0  # next rank after two-way tie at 1

    def test_empty_list(self):
        assert compute_ranks_descending([]) == []

    def test_rank_count_matches_input_length(self):
        values = [1.5, None, 3.0, 2.0, None]
        ranks = compute_ranks_descending(values)
        assert len(ranks) == len(values)

    def test_nan_treated_as_missing(self):
        ranks = compute_ranks_descending([2.0, float("nan"), 1.0])
        assert ranks[1] is None
