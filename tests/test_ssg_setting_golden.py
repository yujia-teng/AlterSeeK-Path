"""Golden regression test for the experimental --ssg-setting workflow.

--ssg-setting (ssg_setting.py + KPointsModifier(magnetic_setting=True) +
FindSpinGroup's acc-primitive setting) had no automated coverage. This drives it
end-to-end on a GdAuGe 2x1x1 supercell (SUPERCELL_211.vasp): the structural cell
is low-symmetry (oC1/mP1) but --ssg-setting recovers the hP2 magnetic primitive
and writes the hP2 butterfly path. Verified byte-identical between `main` and the
`restructure` branch when this test was written.

Guards the shared helpers in ssg_setting.py / find_sf_operations.py (the deferred
de-duplication) and the SSG branch of interactive_modify against future breakage.
"""
import io
import sys
from pathlib import Path

import pytest

# --ssg-setting needs FindSpinGroup (acc-primitive) and ASE; skip cleanly (e.g. CI)
# when they are not installed.
pytest.importorskip("findspingroup")
pytest.importorskip("ase")

REF_DIR = Path(__file__).parent / "references"
POSCAR = REF_DIR / "SUPERCELL_211.vasp"
GOLDEN = REF_DIR / "ssg_supercell211_golden_kpoints.txt"


@pytest.mark.skipif(not POSCAR.exists(), reason="SSG test input not present")
def test_ssg_setting_supercell211_golden(tmp_path, monkeypatch):
    try:
        from kpoints import KPointsModifier
        import ssg_setting
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"kpoints/ssg_setting unavailable: {exc}")
    if ssg_setting.find_spin_group_acc_primitive is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    # struct file, spin axis, moments (4 Gd: c-axis AFM + + - -), [Enter]=auto
    # path, [Enter]=Option 1 op, [Enter]=vasp, output filename.
    answers = "\n".join([str(POSCAR), "0 0 1", "5 5 -5 -5", "", "", "", "OUT"]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    KPointsModifier(magnetic_setting=True).interactive_modify()

    produced = (tmp_path / "OUT").read_text(encoding="utf-8")
    expected = GOLDEN.read_text(encoding="utf-8")
    assert produced.splitlines() == expected.splitlines()
