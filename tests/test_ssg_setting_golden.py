"""Golden regression test for the experimental --ssg-setting workflow.

--ssg-setting (ssg_setting.py + KPointsModifier(magnetic_setting=True) +
FindSpinGroup's acc-primitive setting) had no automated coverage. This drives it
end-to-end on a GdAuGe 2x1x1 supercell (SUPERCELL_211.vasp) with the correct
AFM magnetic order (moments 1 -1 1 -1 along c, alternating by the a1-doubling
translation): --ssg-setting recovers the oC1 orthorhombic magnetic primitive
(MSG P_Cmc2_1, BNS 26.76) and writes its default (non-altermagnetic) path.

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
        from alterseek.kpoints import KPointsModifier
        from alterseek import ssg_setting
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"kpoints/ssg_setting unavailable: {exc}")
    if ssg_setting.find_spin_group_acc_primitive is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    # struct file, spin axis, moments (4 Gd: AFM 1 -1 1 -1), [Enter]=auto path
    # (not altermagnetic, so no operation-choice prompt), [Enter]=vasp,
    # output filename.
    answers = "\n".join([str(POSCAR), "0 0 1", "1 -1 1 -1", "", "", "OUT"]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    KPointsModifier(magnetic_setting=True).interactive_modify()

    produced = (tmp_path / "OUT").read_text(encoding="utf-8")
    expected = GOLDEN.read_text(encoding="utf-8")
    assert produced.splitlines() == expected.splitlines()
