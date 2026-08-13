"""Golden regression test for magnetic-primitive-cell path construction.

This drives the current magnetic path-construction route (ssg_setting.py +
KPointsModifier + FindSpinGroup's acc-primitive setting) end-to-end on a
GdAuGe 2x1x1 supercell (SUPERCELL_211.vasp) with the correct AFM magnetic
order (moments 1 -1 1 -1 along c, alternating by the a1-doubling
translation): the magnetic primitive cell is recovered as oC1 orthorhombic
(MSG P_Cmc2_1, BNS 26.76) and its default (non-altermagnetic) path is
written.

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
CHANGED_CELL_AM = REF_DIR / "case15_changed_cell_altermagnet.vasp"
GOLDEN = REF_DIR / "ssg_supercell211_golden_kpoints.txt"


def test_marker_species_is_deterministic_and_absent_from_real_atoms():
    from alterseek import ssg_setting

    operation = {
        "real_rotation": np.eye(3),
        "translation": [0.0, 0.0, 0.0],
    }
    common = {
        "lattice": np.diag([4.0, 4.0, 4.0]),
        "counts": [1],
        "positions": [np.zeros(3)],
        "operations": [operation],
    }

    ordinary = ssg_setting._build_marker_helper(
        symbols=["Fe"], **common
    )
    helium_structure = ssg_setting._build_marker_helper(
        symbols=["He"], **common
    )

    assert ordinary["marker_species"] == "He"
    assert ordinary["symbols"] == ["Fe", "He"]
    assert helium_structure["marker_species"] == "Ne"
    assert helium_structure["symbols"] == ["He", "Ne"]
    assert ssg_setting._select_marker_species(["He", "Ne"]) == "Ar"


def test_marker_cleanup_preserves_physical_helium(tmp_path):
    from alterseek import ssg_setting
    from alterseek.io import _read_grouped_poscar

    lattice = np.diag([4.0, 4.0, 8.0])
    helper = tmp_path / "standardized_marker_helper.vasp"
    ssg_setting._write_poscar(
        helper,
        "physical He plus artificial Ne marker",
        lattice,
        ["He", "Ne"],
        [1, 1],
        [np.zeros(3), np.array([0.21, 0.17, 0.13])],
    )

    result = ssg_setting.finalize_magnetic_setting_outputs(
        {
            "basename": "helium_case",
            "marker_species": "Ne",
            "magnetic_primitive_lattice": lattice,
            "submitted_lattice": lattice.copy(),
            "magnetic_symbols": ["He"],
        },
        {"standardized_structure_path": str(helper), "symprec": 1e-3},
        output_dir=str(tmp_path),
        calculation_cell_dir=str(tmp_path),
    )

    _, elements, _ = _read_grouped_poscar(result["standard_real_path"])
    assert elements == ["He"]


def _assert_structural_standard_exists(out, stem):
    standard_vasp = out / f"{stem}_seekpath_standard.vasp"
    assert standard_vasp.exists()


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

    KPointsModifier().interactive_modify()

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
    _assert_structural_standard_exists(out, POSCAR.stem)
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

    assert KPointsModifier().interactive_modify() is True

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
    _assert_structural_standard_exists(out, POSCAR_221.stem)

    kpoints_text = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    # The recorded matrix is the raw Step-3 selection, which on this route is
    # expressed in the magnetic primitive basis -- not the submitted 2x2x1
    # calculation cell. The header must name that basis explicitly.
    assert kpoints_text.splitlines()[0].startswith(
        "Selected spin-flip operation (Option 1) in magnetic primitive cell "
        "'SUPERCELL_221_magnetic_primitive.mcif' real-space fractional basis:"
    ), kpoints_text.splitlines()[0]
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


def test_custom_path_in_magnetic_route_reads_submitted_basis(
    tmp_path,
    monkeypatch,
):
    """A custom Step-1 path is defined in the submitted structure's basis.

    The 221 supercell run analyzes symmetry in the six-site magnetic
    primitive cell while keeping the submitted 2x2x1 cell as the output
    basis. A custom endpoint given in the submitted basis must therefore
    come back out numerically unchanged in KPOINTS_alter. Interpreting the
    file in the magnetic analysis basis instead (the old behavior) rescales
    the in-plane coordinates and breaks this round trip.
    """
    from alterseek.kpoints import KPointsModifier
    from alterseek import ssg_setting

    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    custom_path = tmp_path / "KPATH_custom.in"
    custom_path.write_text(
        "Custom path in the submitted SUPERCELL_221 basis\n"
        "   30\n"
        "Line-Mode\n"
        "Reciprocal\n"
        "   0.0000000000   0.0000000000   0.0000000000     GAMMA\n"
        "   0.6666666667  -1.3333333333   0.0000000000     K\n",
        encoding="utf-8",
    )

    answers = "\n".join([
        str(POSCAR_221),
        "0 0 1",
        "4*8 4*-8 16*0",
        str(custom_path),
        "1",
        "",
    ]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    assert KPointsModifier().interactive_modify() is True

    kpoints_text = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    rows = {}
    for line in kpoints_text.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[-1] in {"GAMMA", "K"}:
            rows[parts[-1]] = [float(value) for value in parts[:3]]
    assert rows["GAMMA"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-8)
    assert rows["K"] == pytest.approx(
        [0.6666666667, -1.3333333333, 0.0], abs=1e-8
    )


def test_changed_calculation_cell_builds_butterfly_path(
    tmp_path,
    monkeypatch,
    capsys,
):
    """Exercise cell replacement and butterfly construction in one run.

    The structure is the tracked Case-15 Ca(Al2Fe)4 crystal, but this test uses
    the balanced Fe order ++-- rather than the paper case's original +--+
    order. FindSpinGroup identifies this validation configuration as an
    altermagnet with G0 Immm (71). Its same-volume magnetic primitive setting
    is a genuine nontrivial basis change, so the calculation POSCAR must be
    replaced before the oI3 butterfly path is written.
    """
    from ase.io import read
    from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
    from alterseek import ssg_setting

    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    answers = "\n".join([
        str(CHANGED_CELL_AM),
        "0 0 1",
        "0 2*1 2*-1 8*0",
        "",
        "",
        "",
    ]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    assert KPointsModifier().interactive_modify() is True

    stdout = capsys.readouterr().out
    assert "Magnetic primitive cell (G0): SG Immm (71)" in stdout
    assert "Phase: AFM(Altermagnet)" in stdout
    assert "11 original segments -> 28 generated segments" in stdout

    stem = CHANGED_CELL_AM.stem
    calculation_poscar = tmp_path / f"{stem}_magnetic_primitive.vasp"
    calculation_magmom = tmp_path / f"{stem}_magnetic_primitive_MAGMOM.txt"
    assert {path.name for path in tmp_path.iterdir()} == {
        "KPOINTS_alter",
        "alterband.toml",
        calculation_poscar.name,
        calculation_magmom.name,
        OUTPUT_DIR,
    }

    submitted_lattice = np.asarray(read(CHANGED_CELL_AM).cell, dtype=float)
    calculation_lattice = np.asarray(read(calculation_poscar).cell, dtype=float)
    assert not np.allclose(calculation_lattice, submitted_lattice)
    assert np.isclose(
        abs(np.linalg.det(calculation_lattice)),
        abs(np.linalg.det(submitted_lattice)),
    )

    magmom_text = calculation_magmom.read_text(encoding="utf-8")
    assert f"# {calculation_poscar.name}" in magmom_text
    assert "# Species order: Al Ca Fe" in magmom_text
    assert "# Counts: 8 1 4" in magmom_text
    assert "MAGMOM = 0 0 0 0 0 0 0 0 0 1 -1 -1 1" in magmom_text
    _assert_structural_standard_exists(tmp_path / OUTPUT_DIR, stem)

    coordinate_rows = []
    for line in (tmp_path / "KPOINTS_alter").read_text(
        encoding="utf-8"
    ).splitlines()[4:]:
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            coords = [float(value) for value in fields[:3]]
        except ValueError:
            continue
        coordinate_rows.append((coords, fields[3]))

    assert len(coordinate_rows) == 56
    labels = [label for _, label in coordinate_rows]
    assert "k" in labels and "k'" in labels
    assert any(label.endswith("'") and label != "k'" for label in labels)
    assert any(
        label == "k" and np.allclose(coords, [0.1920785823, 0.1273278402, 0.1273278402])
        for coords, label in coordinate_rows
    )
    assert any(
        label == "k'" and np.allclose(coords, [0.4467337620, -0.1273278402, -0.1273278402])
        for coords, label in coordinate_rows
    )

    out = tmp_path / OUTPUT_DIR
    mapping = out / f"{stem}_seekpath_basis_mapping.txt"
    mapping_text = mapping.read_text(encoding="utf-8")
    mapping_lines = mapping_text.splitlines()
    output_start = mapping_lines.index("kpoints_output_lattice:") + 1
    mapped_output_lattice = np.array([
        [float(value) for value in line.split()]
        for line in mapping_lines[output_start:output_start + 3]
    ])
    assert np.allclose(mapped_output_lattice, calculation_lattice)
    assert (out / "spin_flip_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0] == (
        "# Left basis: magnetic primitive cell "
        f"'{stem}_magnetic_primitive.mcif' real-space fractional "
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
    assert result["marker_species"] == "He"
    assert result["summary"]["marker_species"] == "He"
