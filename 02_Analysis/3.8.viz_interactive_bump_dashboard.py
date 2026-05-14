#!/usr/bin/env python3
"""
3.8.viz_interactive_bump_dashboard.py
======================================

Entry point for the interactive bump-chart dashboard generator.

This script is intentionally thin: all logic lives inside the
``bump_dashboard`` module (``01_Scripts/Python/bump_dashboard/``).

Usage
-----
    python3 02_Analysis/3.8.viz_interactive_bump_dashboard.py

Output
------
    03_Results/02_Analysis/Plots/Trajectory_Flow/interactive_bump_dashboard.html
"""

import logging
import sys
from pathlib import Path

_SCRIPTS_ROOT = Path(__file__).resolve().parent.parent / "01_Scripts"
sys.path.insert(0, str(_SCRIPTS_ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)

from Python.bump_dashboard import DashboardPipeline  # noqa: E402


def main() -> None:
    pipeline = DashboardPipeline()
    output_path = pipeline.run()
    print(f"\nDashboard written to:\n  {output_path}")
    print("\nOpen in a browser to explore pathway trajectories interactively.")


if __name__ == "__main__":
    main()
