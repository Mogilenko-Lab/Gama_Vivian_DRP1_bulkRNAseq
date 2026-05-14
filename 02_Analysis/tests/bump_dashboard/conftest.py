"""
conftest.py — path bootstrap for bump_dashboard tests.

Ensures ``01_Scripts/`` is on ``sys.path`` before any test module is
imported, regardless of the working directory pytest is invoked from.
"""
import sys
from pathlib import Path

_SCRIPTS_ROOT = Path(__file__).resolve().parents[4] / "01_Scripts"
if str(_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_ROOT))
