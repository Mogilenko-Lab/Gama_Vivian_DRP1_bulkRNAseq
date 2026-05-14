"""
domain.geometry
===============

Pure geometric helpers used by the chart-rendering pipeline.

All functions are stateless and dependency-free (standard library only).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class BezierPoints:
    """Sampled x / y coordinates of a quadratic Bézier curve."""

    x: tuple[float, ...]
    y: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have equal length.")


def quadratic_bezier(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    n_points: int = 20,
) -> BezierPoints:
    """
    Sample *n_points* along a quadratic Bézier curve defined by three control
    points ``(x0,y0) → (x1,y1) → (x2,y2)``.

    Parameters
    ----------
    x0, y0 : float
        Start point.
    x1, y1 : float
        Control point (the "pull" point — not on the curve).
    x2, y2 : float
        End point.
    n_points : int
        Number of sample points (≥ 2).

    Returns
    -------
    BezierPoints

    Raises
    ------
    ValueError
        When ``n_points < 2``.
    """
    if n_points < 2:
        raise ValueError(f"n_points must be ≥ 2, got {n_points}.")

    xs: list[float] = []
    ys: list[float] = []
    for i in range(n_points + 1):
        t = i / n_points
        a = (1.0 - t) ** 2
        b = 2.0 * (1.0 - t) * t
        c = t ** 2
        xs.append(a * x0 + b * x1 + c * x2)
        ys.append(a * y0 + b * y1 + c * y2)

    return BezierPoints(x=tuple(xs), y=tuple(ys))


def compute_curve_offset(
    trajdev_nes: float,
    y1: float,
    y2: float,
    y_metric: str,
    chart_height: float = 100.0,
    nes_scale_factor: float = 0.5,
) -> float:
    """
    Compute the vertical offset of the Bézier control point from the
    midpoint of the two endpoints.

    The offset direction and magnitude encode the trajectory deviation:
    * Positive NES (up-regulated TrajDev) → curve bulges *upward*.
    * Negative NES (down-regulated TrajDev) → curve bulges *downward*.

    Parameters
    ----------
    trajdev_nes : float
        NES of the trajectory deviation contrast.
    y1, y2 : float
        y coordinates of the start and end points (Early / Late).
    y_metric : str
        ``"nes"`` or ``"rank"``.  Rank mode uses chart-height-relative scaling.
    chart_height : float
        Total height of the y-axis (used only in rank mode).
    nes_scale_factor : float
        Fraction of ``trajdev_nes`` applied as offset in NES mode.

    Returns
    -------
    float
        Vertical offset from the midpoint.
    """
    if y_metric == "rank":
        # Positive NES (up-regulated) → lower rank value → curve should go up
        # (up in data space = smaller rank number = negative direction in a
        # reversed y-axis chart, but here we return data-space offset)
        return -1.0 * trajdev_nes * (chart_height * 0.1)
    else:
        return trajdev_nes * nes_scale_factor


def wrap_text(text: str, width: int = 35) -> str:
    """
    Wrap a string at *width* characters, inserting ``<br>`` line breaks.

    Preserves word boundaries.  Trailing space on each line is included for
    simplicity (browsers collapse it).

    Parameters
    ----------
    text : str
        Input string.
    width : int
        Maximum characters per line (≥ 1).

    Returns
    -------
    str
        HTML-wrapped string.
    """
    if not text:
        return ""
    words = text.split()
    current: list[str] = []
    lines: list[str] = []
    current_len = 0

    for word in words:
        # +1 for the space that would precede the word
        test_len = current_len + len(word) + (1 if current else 0)
        if current and test_len > width:
            lines.append(" ".join(current))
            current = [word]
            current_len = len(word)
        else:
            current.append(word)
            current_len = test_len

    if current:
        lines.append(" ".join(current))

    return "<br>".join(lines)
