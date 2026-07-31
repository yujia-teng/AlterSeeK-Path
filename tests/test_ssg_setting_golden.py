"""Golden regression test for magnetic-primitive-cell path construction.

This is the default path-construction route (ssg_setting.py +
KPointsModifier(magnetic_setting=True) + FindSpinGroup's acc-primitive
setting); --parent-setting opts out of it. This drives it end-to-end on a
GdAuGe 2x1x1 supercell (SUPERCELL_211.vasp) with the correct AFM magnetic
order (moments 1 -1 1 -1 along c, alternating by the a1-doubling
translation): the magnetic primitive cell is recovered as oC1 orthorhombic
(MSG P_Cmc2_1, BNS 26.76) and its default (non-altermagnetic) path is
written.

magnetic_setting=True is passed explicitly even though it is now the
default, so this test keeps covering the magnetic-primitive route if the
default is ever changed again.

Guards the shared helpers in ssg_setting.py / find_sf_operations.py (the deferred
de-duplication) and the magnetic-cell branch of interactive_modify against
future breakage.
"""
import io
import sys
from pathlib import Path

import numpy as np
import pytest

# Magnetic-primitive construction needs FindSpinGroup (acc-primitive) and ASE;
# skip cleanly (e.g. CI) when they are not installed.
pytest.importorskip("findspingroup")
pytest.importorskip("ase")

REF_DIR = Path(__file__).parent / "references"
POSCAR = REF_DIR / "SUPERCELL_211.vasp"
POSCAR_221 = REF_DIR / "SUPERCELL_221.vasp"
GOLDEN = REF_DIR / "ssg_supercell211_golden_kpoints.txt"


@pytest.mark.skipif(not POSCAR.exists(), reason="SSG test input not present")
def test_ssg_setting_supercell211_golden(tmp_path, monkeypatch):
    try:
        from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
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

    # The basis mapping is the only record of how the path's reciprocal basis
    # relates to the submitted cell, and this route is where it matters most:
    # the standardized cell is one the user never supplied. It used to be
    # deleted here while the ordinary route kept it.
    out = tmp_path / OUTPUT_DIR
    mapping = out / f"{POSCAR.stem}_seekpath_basis_mapping.txt"
    assert mapping.exists(), "magnetic route discarded the SeeK-path basis mapping"
    mapping_text = mapping.read_text(encoding="utf-8")
    header = mapping_text.splitlines()[1]
    # It must name the marker-bearing magnetic cell actually submitted to
    # SeeK-path, not either the user's submitted structure or the final
    # marker-stripped standardized diagnostic.
    assert "_magnetic_marker_input.vasp" in header, header
    assert "analysis_input_lattice:" in mapping_text
    assert "seekpath_standard_primitive_lattice:" in mapping_text
    assert "seekpath_standard_conventional_lattice:" in mapping_text
    assert "kpoints_output_lattice:" in mapping_text
    assert "SUPERCELL_211_magnetic_primitive.vasp" in mapping_text
    assert (out / f"{POSCAR.stem}_seekpath_standard.vasp").exists()
    assert (out / f"{POSCAR.stem}_magnetic_primitive.mcif").exists()
    assert (out / "spin_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0] == (
        "# Basis: submitted structure 'SUPERCELL_211.vasp' real-space "
        "fractional basis (a1, a2, a3)."
    )
    assert (out / "spin_flip_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0] == (
        "# Basis: magnetic primitive cell "
        "'SUPERCELL_211_magnetic_primitive.mcif' real-space fractional "
        "basis (a1, a2, a3)."
    )

    # Only files the user acts on stay at the top level: KPOINTS_alter feeds
    # the band calculation, alterband.toml is read from the working directory
    # by the band plotter, and the magnetic primitive POSCAR plus its matching
    # reordered moments are the inputs the calculation must use. The last two
    # appear only on this route -- when the submitted cell already carries the
    # path, there is no new calculation cell to hand over.
    calculation_poscar = tmp_path / f"{POSCAR.stem}_magnetic_primitive.vasp"
    calculation_magmom = \
        tmp_path / f"{POSCAR.stem}_magnetic_primitive_MAGMOM.txt"
    top_level = {p.name for p in tmp_path.iterdir()}
    assert top_level == {
        "KPOINTS_alter",
        "alterband.toml",
        calculation_poscar.name,
        calculation_magmom.name,
        OUTPUT_DIR,
    }, top_level
    assert "moments: SUPERCELL_211_magnetic_primitive_MAGMOM.txt" in \
        calculation_poscar.read_text(encoding="utf-8").splitlines()[0]
    magmom_text = calculation_magmom.read_text(encoding="utf-8")
    assert "# SUPERCELL_211_magnetic_primitive.vasp" in magmom_text
    assert "# Species order: Au Gd Ge" in magmom_text
    assert "# Counts: 4 4 4" in magmom_text
    assert "# Collinear spin axis (Cartesian): 0 0 1" in magmom_text
    assert "MAGMOM = 0 0 0 0 1 1 -1 -1 0 0 0 0" in magmom_text
    assert "# Vector moments" not in magmom_text
    assert ".alterseek_magnetic_tmp" not in magmom_text

    generated_names = {
        path.name for path in tmp_path.rglob("*") if path.is_file()
    }
    assert not any("ssgprim" in name or "ssgstd" in name
                   for name in generated_names)


@pytest.mark.skipif(not POSCAR_221.exists(), reason="221 test input not present")
def test_ssg_setting_keeps_submitted_221_calculation_supercell(
    tmp_path,
    monkeypatch,
):
    """Analyze G0 in its primitive cell but write k-points in the submitted
    2x2x1 calculation-supercell basis.

    The A-type order is four +c Gd moments at z=0 followed by four -c
    moments at z=1/2; Au and Ge are nonmagnetic.  FindSpinGroup reduces this
    24-site input to the six-site magnetic primitive cell.  That reduction
    must not force the user to replace an intentional calculation supercell.
    """
    from ase.io import read
    from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
    from alterseek import ssg_setting

    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    answers = "\n".join([
        str(POSCAR_221),
        "0 0 1",
        "4*8 4*-8 16*0",
        "",
        "1",
        "",
    ]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    assert KPointsModifier(magnetic_setting=True).interactive_modify() is True

    out = tmp_path / OUTPUT_DIR
    mapping = out / f"{POSCAR_221.stem}_seekpath_basis_mapping.txt"
    mapping_text = mapping.read_text(encoding="utf-8")
    mapping_lines = mapping_text.splitlines()
    output_start = mapping_lines.index("kpoints_output_lattice:") + 1
    output_lattice = np.array([
        [float(value) for value in line.split()]
        for line in mapping_lines[output_start:output_start + 3]
    ])
    submitted_lattice = np.asarray(read(POSCAR_221).cell, dtype=float)
    assert np.allclose(output_lattice, submitted_lattice)

    top_level = {path.name for path in tmp_path.iterdir()}
    assert top_level == {
        "KPOINTS_alter",
        "alterband.toml",
        OUTPUT_DIR,
    }, top_level
    assert not (tmp_path / f"{POSCAR_221.stem}_magnetic_primitive.vasp").exists()
    assert not (
        tmp_path / f"{POSCAR_221.stem}_magnetic_primitive_MAGMOM.txt"
    ).exists()

    kpoints_text = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    k_rows = [
        line.split()
        for line in kpoints_text.splitlines()
        if line.split() and line.split()[-1] == "K"
    ]
    assert any(
        np.allclose(
            [float(value) for value in row[:3]],
            [2.0 / 3.0, -4.0 / 3.0, 0.0],
        )
        for row in k_rows
    ), k_rows
    assert (out / f"{POSCAR_221.stem}_magnetic_primitive.mcif").exists()
    assert (out / "spin_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0] == (
        "# Basis: submitted structure 'SUPERCELL_221.vasp' real-space "
        "fractional basis (a1, a2, a3)."
    )
    assert (out / "spin_flip_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0] == (
        "# Left basis: magnetic primitive cell "
        "'SUPERCELL_221_magnetic_primitive.mcif' real-space fractional "
        "basis (a1, a2, a3)."
    )


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
