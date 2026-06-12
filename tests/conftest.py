import os
import sys
from pathlib import Path

# Headless plotting for all tests (figures are created but never shown).
os.environ.setdefault("MPLBACKEND", "Agg")

# Make the repo-root modules importable when pytest is run from anywhere.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
