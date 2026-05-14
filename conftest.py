"""
Root conftest.py — ensures 01_Scripts/ is on sys.path for all tests.
"""
import sys
from pathlib import Path

_SCRIPTS_ROOT = Path(__file__).resolve().parent / "01_Scripts"
if str(_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_ROOT))
