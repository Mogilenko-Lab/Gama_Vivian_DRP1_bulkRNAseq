"""
domain.rules
============

Deterministic domain rules — no I/O, no framework dependencies.

Each function is a standalone predicate or computation on domain primitives.
They are imported by the application layer but tested independently.

Significance threshold
----------------------
All significance decisions use ``PADJ_THRESHOLD``.  It is defined once here
and nowhere else inside this module; callers that need it import it explicitly.
"""

from __future__ import annotations

from typing import Optional

from .schema import WeightCategory

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PADJ_THRESHOLD: float = 0.05
"""FDR-adjusted p-value threshold for statistical significance."""

NES_DIVERGING_MAX: float = 3.5
"""Absolute NES ceiling used when mapping NES → diverging colour scale."""

# Weight category frequency boundaries
_DOMINANT_FREQ: float = 0.30   # > 30 % → dominant
_COMMON_FREQ: float = 0.10     # > 10 % → common
_UNCOMMON_FREQ: float = 0.01   # >  1 % → uncommon
# < 1 % → rare


# ---------------------------------------------------------------------------
# Weight category rules
# ---------------------------------------------------------------------------

def classify_weight_category(frequency: float) -> WeightCategory:
    """
    Map a relative frequency to a named weight category.

    Parameters
    ----------
    frequency:
        The pattern's share of total pathway-mutation observations (0.0–1.0).

    Returns
    -------
    WeightCategory
        One of ``"dominant"``, ``"common"``, ``"uncommon"``, ``"rare"``.

    Raises
    ------
    ValueError
        When ``frequency`` is outside [0, 1].
    """
    if not (0.0 <= frequency <= 1.0):
        raise ValueError(f"frequency must be in [0, 1], got {frequency!r}")

    if frequency > _DOMINANT_FREQ:
        return "dominant"
    if frequency > _COMMON_FREQ:
        return "common"
    if frequency > _UNCOMMON_FREQ:
        return "uncommon"
    return "rare"


# ---------------------------------------------------------------------------
# Significance rules
# ---------------------------------------------------------------------------

def is_significant(padj: Optional[float], threshold: float = PADJ_THRESHOLD) -> bool:
    """Return True iff *padj* is not None and strictly below *threshold*."""
    return padj is not None and padj < threshold


def significance_marker(padj: Optional[float], threshold: float = PADJ_THRESHOLD) -> str:
    """Return ``"*"`` for significant results, ``""`` otherwise."""
    return "*" if is_significant(padj, threshold) else ""


# ---------------------------------------------------------------------------
# NES → colour interpolation rules
# ---------------------------------------------------------------------------

# Colorblind-safe diverging palette (matches R ``color_config.R`` and
# ``color_config.py``).
_COLOR_NEGATIVE: str = "#2166AC"   # Blue  (negative NES)
_COLOR_NEUTRAL: str = "#F7F7F7"    # White (zero NES)
_COLOR_POSITIVE: str = "#B35806"   # Orange (positive NES)


def _hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    h = hex_color.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_to_hex(r: int, g: int, b: int) -> str:
    return f"#{r:02x}{g:02x}{b:02x}"


def interpolate_hex_color(c1: str, c2: str, t: float) -> str:
    """
    Linearly interpolate between two hex colors.

    Parameters
    ----------
    c1, c2:
        Hex colour strings (``"#rrggbb"``).
    t:
        Interpolation factor in [0, 1] where 0 → c1 and 1 → c2.
    """
    r1, g1, b1 = _hex_to_rgb(c1)
    r2, g2, b2 = _hex_to_rgb(c2)
    r = round(r1 + (r2 - r1) * t)
    g = round(g1 + (g2 - g1) * t)
    b = round(b1 + (b2 - b1) * t)
    return _rgb_to_hex(r, g, b)


def nes_to_diverging_color(
    nes: Optional[float],
    max_nes: float = NES_DIVERGING_MAX,
) -> str:
    """
    Map an NES value to a hex color on the Blue-White-Orange diverging scale.

    Parameters
    ----------
    nes:
        Normalized enrichment score.  ``None`` / ``NaN`` returns ``"#999999"``.
    max_nes:
        Absolute value used to saturate the scale.

    Returns
    -------
    str
        Hex color string.
    """
    import math

    if nes is None or (isinstance(nes, float) and math.isnan(nes)):
        return "#999999"

    t = max(-1.0, min(1.0, nes / max_nes))
    if t <= 0.0:
        # negative: blue → white
        return interpolate_hex_color(_COLOR_NEGATIVE, _COLOR_NEUTRAL, 1.0 + t)
    else:
        # positive: white → orange
        return interpolate_hex_color(_COLOR_NEUTRAL, _COLOR_POSITIVE, t)


# ---------------------------------------------------------------------------
# Rank computation rule
# ---------------------------------------------------------------------------

def compute_ranks_descending(values: list[Optional[float]]) -> list[Optional[float]]:
    """
    Assign integer ranks to a list of values (descending: highest value = rank 1).

    ``None`` values receive rank ``None``.  Ties use *min* method (equal values
    share the lowest rank in their group).

    Parameters
    ----------
    values:
        Flat list of NES values (or any numeric), possibly containing ``None``.

    Returns
    -------
    list[Optional[float]]
        Ranks in the same order as *values*.
    """
    import math

    indexed = [(i, v) for i, v in enumerate(values) if v is not None and not math.isnan(v)]
    indexed.sort(key=lambda x: -x[1])  # descending

    rank_map: dict[int, int] = {}
    # Min-rank ties: give all tied values the rank of the first occurrence
    i = 0
    while i < len(indexed):
        j = i
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        for k in range(i, j):
            rank_map[indexed[k][0]] = i + 1  # 1-based
        i = j

    return [float(rank_map[i]) if i in rank_map else None for i in range(len(values))]
