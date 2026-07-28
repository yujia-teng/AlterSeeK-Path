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

import numpy as np
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
    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    # struct file, spin axis, moments (4 Gd: AFM 1 -1 1 -1), [Enter]=auto path
    # (not altermagnetic, so no operation-choice prompt), [Enter]=vasp.
    # Output filename is fixed (KPOINTS_alter), no longer prompted.
    answers = "\n".join([str(POSCAR), "0 0 1", "1 -1 1 -1", "", ""]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    KPointsModifier(magnetic_setting=True).interactive_modify()

    produced = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    expected = GOLDEN.read_text(encoding="utf-8")
    assert produced.splitlines() == expected.splitlines()


def test_prepare_mcif_uses_refined_from_data_route(tmp_path, monkeypatch):
    from alterseek import ssg_setting

    structure_path = tmp_path / "rounded-special-position.mcif"
    structure_path.write_text("# loader is mocked in this focused routing test\n", encoding="utf-8")

    lattice = np.eye(3)
    positions = np.array([[1.0 / 3.0, 2.0 / 3.0, 0.0]])
    elements = ["Mn"]
    moments = np.array([[0.0, 0.0, 1.0]])
    calls = {}

    def fake_load(structure_file, moments_str, spin_axis_cart):
        calls["load"] = (structure_file, moments_str, spin_axis_cart)
        return lattice, positions, elements, moments, "cartesian"

    identity_operation = {
        "index": 1,
        "spin_rotation": np.eye(3).tolist(),
        "real_rotation": np.eye(3).tolist(),
        "translation": [0.0, 0.0, 0.0],
    }

    def fake_from_data(
        source_name,
        lattice_factors,
        input_positions,
        input_elements,
        occupancies,
        input_moments,
        input_spin_setting,
    ):
        calls["from_data"] = {
            "source_name": source_name,
            "lattice": np.asarray(lattice_factors),
            "positions": np.asarray(input_positions),
            "elements": list(input_elements),
            "occupancies": list(occupancies),
            "moments": np.asarray(input_moments),
            "spin_setting": input_spin_setting,
        }
        return {
            "index": "test",
            "acc_symbol": "test",
            "acc_primitive_cell_setting": "acc_primitive",
            "acc_primitive_cell_detail": {
                "lattice": lattice.tolist(),
                "positions": positions.tolist(),
                "elements": elements,
                "moments": moments.tolist(),
            },
            "operation_views": {
                "magnetic_primitive_cartesian": {
                    "views": {"nssg": {"ops": [identity_operation]}}
                }
            },
        }

    monkeypatch.setattr(ssg_setting, "_load_magnetic_input_data", fake_load)
    monkeypatch.setattr(
        ssg_setting,
        "find_spin_group_acc_primitive_from_data",
        fake_from_data,
    )

    result = ssg_setting.prepare_magnetic_setting_files(
        str(structure_path),
        output_dir=str(tmp_path),
    )

    assert calls["load"] == (str(structure_path), "", None)
    assert calls["from_data"]["source_name"] == str(structure_path)
    assert np.allclose(calls["from_data"]["lattice"], lattice)
    assert np.allclose(calls["from_data"]["positions"], positions)
    assert calls["from_data"]["elements"] == elements
    assert calls["from_data"]["occupancies"] == [1.0]
    assert np.allclose(calls["from_data"]["moments"], moments)
    assert calls["from_data"]["spin_setting"] == "cartesian"
    assert Path(result["real_poscar_path"]).exists()
    assert Path(result["helper_path"]).exists()
