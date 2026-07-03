"""Golden regression test for KPointsModifier.interactive_modify.

The interactive Step 0-5 driver has no other coverage (the rest of the suite
tests the engine methods).  This drives the full altermagnetic flow on case 12
(MnF2, tetragonal P4/mmm, d-wave) with canned keyboard input and asserts the
produced KPOINTS matches a stored reference byte-for-byte.  It guards the phase-5
extraction of interactive_modify into helper methods (and any future edits).
"""
import io
import sys
from pathlib import Path

import pytest

# The full interactive flow runs Step 0 (FindSpinGroup) and reads the structure
# with ASE; skip cleanly (e.g. CI) when those optional deps are not installed.
pytest.importorskip("findspingroup")
pytest.importorskip("ase")

REPO = Path(__file__).resolve().parents[1]
POSCAR = REPO / "data" / "case_12_t-P4mmm-2121" / "POSCAR"
REFERENCE = Path(__file__).parent / "references" / "case12_golden_kpoints.txt"


@pytest.mark.skipif(not POSCAR.exists(), reason="case_12 POSCAR not present")
def test_interactive_modify_case12_golden(tmp_path, monkeypatch):
    try:
        from alterseek.kpoints import KPointsModifier
    except Exception as exc:  # pragma: no cover - deps missing
        pytest.skip(f"kpoints/deps unavailable: {exc}")

    # Canned answers: structure file, spin axis, moments, [Enter]=auto path,
    # [Enter]=Option 1 op, [Enter]=vasp. Output filename is fixed
    # (KPOINTS_alter), no longer prompted.
    answers = "\n".join([str(POSCAR), "0 0 1", "5 -5", "", "", ""]) + "\n"
    monkeypatch.chdir(tmp_path)               # side-effect files land in tmp
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    KPointsModifier().interactive_modify()

    produced = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    expected = REFERENCE.read_text(encoding="utf-8")
    assert produced.splitlines() == expected.splitlines()
