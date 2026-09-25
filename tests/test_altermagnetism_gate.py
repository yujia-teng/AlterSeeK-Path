"""Regression tests for the altermagnetism Laue-group gate.

The gate decides whether a structure can host altermagnetic splitting at all;
when it fires, the ordinary structural path is written and no butterfly path is
built. Getting it right needs two things that MnSe2 (MAGNDATA 1.0.47) shows are
independent:

1. The submitted translation lattice defines the BZ and output basis.

2. The point symmetry inside that lattice must still be the magnetic G0.
   MnSe2's moment-free coordinates recover the cubic parent and would therefore
   reproduce the forbidding m-3. At spglib's default symprec the MCIF's
   5-decimal rounding can instead make it look orthorhombic, which is an
   accident, not a magnetic-symmetry result.

G0, the spatial part of the spin space group (Pbca 61 -> Laue mmm here), is the
symmetry that actually holds, and FindSpinGroup reports it directly.
"""
import io
import sys
from pathlib import Path

import numpy as np
import pytest

from alterseek.kpoints import _altermagnetism_gate, _g0_symmetry
from alterseek.symmetry import laue_group_from_spacegroup_number

REF_DIR = Path(__file__).parent / "references"
POSCAR = REF_DIR / "SUPERCELL_211.vasp"

CUBIC_PARENT = {'point_group': 'm-3', 'laue_group': 'm-3'}
# MnSe2's G0: orthorhombic Pbca, the symmetry left once the moments are counted.
ORTHORHOMBIC_G0 = {'label': 'Pbca (61)', 'spacegroup_number': 61, 'laue_group': 'mmm'}
CUBIC_G0 = {'label': 'Pa-3 (205)', 'spacegroup_number': 205, 'laue_group': 'm-3'}


def _synthetic_analysis_preparation():
    return {
        "analysis_cell": (
            np.eye(3).tolist(),
            [[0.11, 0.12, 0.15], [0.31, 0.27, 0.19]],
            [119, 119],
        ),
        "analysis_has_markers": True,
        "submitted_lattice": np.eye(3),
        "submitted_sites": 2,
        "magnetic_primitive_sites": 2,
        "uses_conventional_supercell_bz": False,
        "summary": {
            "submitted_to_primitive_volume_index": 1,
        },
        "input_cell_symmetry": {
            "number": 47,
            "symbol": "Pmmm",
            "point_group": "mmm",
            "seekpath_bravais": "oP1",
        },
        "nonmagnetic_primitive_symmetry": {
            "number": 47,
            "symbol": "Pmmm",
            "point_group": "mmm",
            "seekpath_bravais": "oP1",
            "sites": 2,
        },
        "analysis_symmetry": {
            "number": 47,
            "symbol": "Pmmm",
            "point_group": "mmm",
            "seekpath_bravais": "oP1",
        },
    }


def test_submitted_cell_decides_when_no_magnetic_cell_was_built():
    """The submitted-cell symmetry applies when no helper cell is needed."""
    assert _altermagnetism_gate(CUBIC_PARENT) is not None
    assert _altermagnetism_gate({'point_group': 'mmm', 'laue_group': 'mmm'}) is None


def test_g0_overrides_a_forbidding_parent():
    """The MnSe2 case: cubic m-3 parent, orthorhombic mmm G0."""
    assert _altermagnetism_gate(CUBIC_PARENT, working_cell_symmetry=ORTHORHOMBIC_G0) is None


def test_g0_can_also_forbid_where_the_parent_permitted():
    """Not a one-way override -- the working cell's own symmetry decides."""
    reason = _altermagnetism_gate(
        {'point_group': 'mmm', 'laue_group': 'mmm'}, working_cell_symmetry=CUBIC_G0
    )
    assert reason is not None
    assert reason['laue_group'] == 'm-3'


@pytest.mark.parametrize("spacegroup_number,laue", [(1, '-1'), (146, '-3'), (205, 'm-3')])
def test_forbidden_laue_groups_fire_via_g0(spacegroup_number, laue):
    reason = _altermagnetism_gate(
        CUBIC_PARENT,
        working_cell_symmetry={'label': 'x', 'spacegroup_number': spacegroup_number,
                               'laue_group': laue},
    )
    assert reason is not None
    assert reason['laue_group'] == laue


def test_g0_symmetry_is_read_from_findspingroup_not_from_coordinates():
    """The whole point: no structure file, no tolerance, no spglib call."""
    sf_result = {'g0_symbol': 'Pbca', 'g0_number': 61}
    assert _g0_symmetry(sf_result) == {
        'label': 'Pbca (61)', 'spacegroup_number': 61,
        'point_group': 'mmm', 'laue_group': 'mmm', 'sites': None,
    }
    # MnSe2's real FindSpinGroup output must permit altermagnetism.
    assert _altermagnetism_gate(CUBIC_PARENT, _g0_symmetry(sf_result)) is None


def test_g0_symmetry_is_none_when_findspingroup_reported_no_g0():
    assert _g0_symmetry({}) is None
    assert _g0_symmetry({'g0_symbol': 'Pbca', 'g0_number': 999}) is None


@pytest.mark.parametrize("number,expected", [
    (1, '-1'), (2, '-1'), (3, '2/m'), (15, '2/m'), (16, 'mmm'), (61, 'mmm'),
    (74, 'mmm'), (75, '4/m'), (88, '4/m'), (89, '4/mmm'), (142, '4/mmm'),
    (143, '-3'), (148, '-3'), (149, '-3m'), (167, '-3m'), (168, '6/m'),
    (176, '6/m'), (177, '6/mmm'), (194, '6/mmm'), (195, 'm-3'), (205, 'm-3'),
    (206, 'm-3'), (207, 'm-3m'), (230, 'm-3m'),
])
def test_laue_group_covers_every_crystal_class_boundary(number, expected):
    assert laue_group_from_spacegroup_number(number) == expected


@pytest.mark.parametrize("bad", [0, 231, -1, None, "x", 1.5e300])
def test_laue_group_rejects_non_spacegroup_numbers(bad):
    assert laue_group_from_spacegroup_number(bad) is None


def test_point_group_table_matches_spglib_for_all_230_groups():
    """The table is a hardcoded range map, so pin it to spglib's own answer."""
    spglib = pytest.importorskip("spglib")
    from alterseek.symmetry import point_group_from_spacegroup_number

    truth = {}
    for hall in range(1, 531):
        entry = spglib.get_spacegroup_type(hall)
        truth.setdefault(int(entry.number), entry.pointgroup_international)
    assert len(truth) == 230
    mismatches = {n: (point_group_from_spacegroup_number(n), pg)
                  for n, pg in truth.items()
                  if point_group_from_spacegroup_number(n) != pg}
    assert not mismatches


@pytest.mark.parametrize("bad", [0, 231, -1, None, "x"])
def test_point_group_rejects_non_spacegroup_numbers(bad):
    from alterseek.symmetry import point_group_from_spacegroup_number
    assert point_group_from_spacegroup_number(bad) is None


def test_g0_symmetry_reports_point_group_and_site_count():
    """Both cell lines carry the same fields so they can be compared directly."""
    g0 = _g0_symmetry({'g0_symbol': 'Cmc2_1', 'g0_number': 36}, sites=12)
    assert g0['label'] == 'Cmc2_1 (36)'
    assert g0['point_group'] == 'mm2'
    assert g0['laue_group'] == 'mmm'
    assert g0['sites'] == 12


def test_nonmagnetic_report_describes_the_primitive_cell_not_the_input():
    """The reported space group belongs to the primitive cell.

    A magnetic supercell input still reports its parent's space group, so the
    site count printed beside it must be the primitive cell's, not the
    submitted cell's -- otherwise the line misstates what the group describes.
    """
    pytest.importorskip("spglib")
    from alterseek.find_sf_operations import _non_magnetic_symmetry

    # Simple cubic, then the same crystal given as a 3x1x1 supercell.
    lattice = np.eye(3) * 4.0
    single = _non_magnetic_symmetry(
        "POSCAR", lattice, np.array([[0.0, 0.0, 0.0]]), np.array([84]), False)
    supercell = _non_magnetic_symmetry(
        "POSCAR", np.diag([12.0, 4.0, 4.0]),
        np.array([[0.0, 0.0, 0.0], [1 / 3, 0.0, 0.0], [2 / 3, 0.0, 0.0]]),
        np.array([84, 84, 84]), False)

    # Same crystal -> same space group and same primitive site count, even
    # though one input had three times as many atoms.
    assert single['spacegroup_number'] == supercell['spacegroup_number']
    assert single['sites'] == supercell['sites'] == 1


def test_cell_rows_align_every_field_into_columns(capsys):
    """Symbols vary in width, so unpadded fields make the reader hunt for the
    difference between the two cells."""
    from alterseek.kpoints import _print_cell_rows

    _print_cell_rows([
        ('Nonmagnetic primitive cell:', 'P6_3mc (186)', '6mm', '6/mmm', '[6 atoms, hP2]'),
        ('Magnetic primitive cell (G0):', 'Cmc2_1 (36)', 'mm2', 'mmm', '[12 atoms, oC1]'),
    ])
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    for token in ('SG ', 'PG ', 'Laue ', '['):
        assert lines[0].index(token) == lines[1].index(token), token


def test_cell_rows_place_a_note_under_the_selected_row(capsys):
    from alterseek.kpoints import _print_cell_rows

    _print_cell_rows([
        ('Input cell:', 'Pbca (61)', 'mmm', 'mmm', '[24 atoms, oP1]'),
        ('Nonmagnetic primitive cell:', 'Pa-3 (205)', 'm-3', 'm-3', '[12 atoms]'),
    ], note="recovered from the input cell", note_after_index=1)
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 3
    # Indented into the field column, not hanging off the label.
    assert lines[2].startswith(" " * 30)
    assert lines[2].strip() == "recovered from the input cell"


def test_cell_rows_leave_no_trailing_whitespace(capsys):
    """A row with no suffix must not emit padding to nowhere."""
    from alterseek.kpoints import _print_cell_rows

    _print_cell_rows([('Nonmagnetic primitive cell:', 'Pa-3 (205)', 'm-3', 'm-3', '')])
    line = capsys.readouterr().out.splitlines()[0]
    assert line == line.rstrip()


def test_cell_suffix_only_describes_what_is_known():
    from alterseek.kpoints import _cell_suffix

    assert _cell_suffix(12, "oC1") == "[12 atoms, oC1]"
    assert _cell_suffix(None, "cP1") == "[cP1]"
    assert _cell_suffix(12, None) == "[12 atoms]"
    assert _cell_suffix(None, None) == ""
    # 'unknown' is a placeholder, not a lattice type worth printing.
    assert _cell_suffix(None, 'unknown') == ""


def test_submitted_helper_failure_aborts_before_centroid(
    tmp_path, monkeypatch, capsys
):
    """A failed submitted-cell helper must abort without another-cell fallback."""
    from alterseek import kpoints as kpoints_module

    structure = tmp_path / "POSCAR"
    structure.write_text("test structure placeholder\n", encoding="utf-8")
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = "1 -1"\n',
        encoding="utf-8",
    )
    sf_result = {
        "g0_number": 61,
        "nonmagnetic_spacegroup_number": 205,
        "nonmagnetic_sites": 2,
        "num_atoms": 2,
    }

    def fail_construction(*args, **kwargs):
        raise RuntimeError("synthetic marker validation failure")

    def forbid_centroid(*args, **kwargs):
        raise AssertionError("centroid generation must not run")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        kpoints_module, "find_sf_run", lambda *args, **kwargs: sf_result
    )
    monkeypatch.setattr(
        kpoints_module, "prepare_submitted_cell_analysis", fail_construction
    )
    monkeypatch.setattr(kpoints_module, "compute_centroid", forbid_centroid)

    assert kpoints_module.KPathBuilder().interactive_build() is False
    output = capsys.readouterr().out
    assert "Submitted-cell analysis helper construction failed" in output
    assert "synthetic marker validation failure" in output
    assert "Aborting" in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_centroid_failure_is_reported_once_under_its_own_headline(
    tmp_path, monkeypatch, capsys
):
    """A required centroid/basis failure is reported once and aborts."""
    from alterseek import kpoints as kpoints_module

    (tmp_path / "POSCAR").write_text("test structure placeholder\n", encoding="utf-8")
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = "1 -1"\n'
        'output_code = "vasp"\n',
        encoding="utf-8",
    )
    # mmm passes the Laue gate, and G0 == the nonmagnetic group with equal site
    # counts keeps the run on the ordinary route (no magnetic cell to build).
    sf_result = {
        "structure_file": "POSCAR",
        "g0_number": 61,
        "g0_symbol": "Pbca",
        "nonmagnetic_spacegroup_number": 61,
        "nonmagnetic_sites": 2,
        "nonmagnetic_lattice": "oP1",
        "num_atoms": 2,
        "space_group": "Pbca (61)",
        "point_group": "mmm",
        "laue_group": "mmm",
        "magnetic_phase": "AFM(Altermagnet)",
        "ssg_index": "61.1.1.1.L",
        "ssg_symbol": "test",
        "magnetic_space_group_without_soc": "Pb'c'a (BNS 61.436)",
        "actual_spin_flip_point_operations": 4,
        "actual_spin_preserve_point_operations": 4,
        "spin_flip_operations": 1,
    }

    calls = []

    def failing_centroid(*args, **kwargs):
        calls.append(args[0] if args else None)
        raise RuntimeError("synthetic seekpath failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))
    monkeypatch.setattr(kpoints_module, "find_sf_run", lambda *a, **k: sf_result)
    monkeypatch.setattr(
        kpoints_module,
        "prepare_submitted_cell_analysis",
        lambda *a, **k: _synthetic_analysis_preparation(),
    )
    monkeypatch.setattr(kpoints_module, "compute_centroid", failing_centroid)

    result = kpoints_module.KPathBuilder().interactive_build()
    output = capsys.readouterr().out

    assert result is False
    assert len(calls) == 1, f"compute_centroid called {len(calls)} times"
    assert "IBZ centroid construction failed" in output
    assert "synthetic seekpath failure" in output
    assert "Aborting" in output
    assert "Falling back to manual file input" not in output
    assert "Auto path generation failed" not in output
    assert "Centroid computation failed" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_altermagnet_with_no_available_spin_flip_operation_aborts(
    tmp_path, monkeypatch, capsys
):
    """An internally inconsistent altermagnetic result must not look successful."""
    from alterseek import kpoints as kpoints_module

    (tmp_path / "POSCAR").write_text(
        "test structure placeholder\n", encoding="utf-8"
    )
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = "1 -1"\n'
        'output_code = "vasp"\n',
        encoding="utf-8",
    )
    sf_result = {
        "structure_file": "POSCAR",
        "g0_number": 61,
        "g0_symbol": "Pbca",
        "nonmagnetic_spacegroup_number": 61,
        "nonmagnetic_sites": 2,
        "nonmagnetic_lattice": "oP1",
        "num_atoms": 2,
        "space_group": "Pbca (61)",
        "point_group": "mmm",
        "laue_group": "mmm",
        "magnetic_phase": "AFM(Altermagnet)",
        "ssg_index": "61.1.1.1.L",
        "ssg_symbol": "test",
        "magnetic_space_group_without_soc": "Pb'c'a (BNS 61.436)",
        "actual_spin_flip_point_operations": 0,
        "actual_spin_preserve_point_operations": 8,
        "spin_flip_operations": 0,
    }
    centroid_result = {
        "display_figures": [],
        "sp_path": [("GAMMA", "X")],
        "sp_point_coords": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
        },
        "b_matrix": np.eye(3),
        "b_matrix_input": np.eye(3),
        "centroid_frac": [0.1, 0.2, 0.3],
        "sc_type": "oP1",
        "point_group": "mmm",
        "spacegroup": 61,
    }

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(kpoints_module, "find_sf_run", lambda *a, **k: sf_result)
    monkeypatch.setattr(
        kpoints_module,
        "prepare_submitted_cell_analysis",
        lambda *a, **k: _synthetic_analysis_preparation(),
    )
    monkeypatch.setattr(
        kpoints_module, "compute_centroid", lambda *a, **k: centroid_result
    )

    assert (
        kpoints_module.KPathBuilder().interactive_build()
        is False
    )
    output = capsys.readouterr().out
    assert "Phase: AFM(Altermagnet)" in output
    assert "no detected spin-flip point operation is available" in output
    assert "inconsistent" in output
    assert "Aborting" in output
    assert "Writing ordinary k-path" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_centroid_failure_without_spin_analysis_also_aborts(
    tmp_path, monkeypatch, capsys
):
    """The ordinary no-moments route has the same required-analysis policy."""
    from alterseek import kpoints as kpoints_module

    (tmp_path / "POSCAR").write_text(
        "test structure placeholder\n", encoding="utf-8"
    )
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = ""\n',
        encoding="utf-8",
    )

    calls = []

    def failing_centroid(*args, **kwargs):
        calls.append(args[0] if args else None)
        raise RuntimeError("synthetic ordinary seekpath failure")

    def forbid_spin_analysis(*args, **kwargs):
        raise AssertionError("no-moments route must not run spin analysis")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))
    monkeypatch.setattr(kpoints_module, "find_sf_run", forbid_spin_analysis)
    monkeypatch.setattr(
        kpoints_module,
        "prepare_submitted_cell_analysis",
        lambda *a, **k: _synthetic_analysis_preparation(),
    )
    monkeypatch.setattr(kpoints_module, "compute_centroid", failing_centroid)

    result = kpoints_module.KPathBuilder().interactive_build()
    output = capsys.readouterr().out

    assert result is False
    assert calls == ["POSCAR"]
    assert "IBZ centroid construction failed" in output
    assert "synthetic ordinary seekpath failure" in output
    assert "Aborting" in output
    assert "Falling back to manual file input" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_step1_internal_error_propagates_without_manual_fallback(
    tmp_path, monkeypatch, capsys
):
    """Malformed required path data reaches the final workflow boundary."""
    from alterseek import kpoints as kpoints_module

    (tmp_path / "POSCAR").write_text(
        "test structure placeholder\n", encoding="utf-8"
    )
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = ""\n',
        encoding="utf-8",
    )

    def forbid_spin_analysis(*args, **kwargs):
        raise AssertionError("no-moments route must not run spin analysis")

    input_stream = io.StringIO("")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", input_stream)
    monkeypatch.setattr(kpoints_module, "find_sf_run", forbid_spin_analysis)
    monkeypatch.setattr(
        kpoints_module,
        "prepare_submitted_cell_analysis",
        lambda *a, **k: _synthetic_analysis_preparation(),
    )
    monkeypatch.setattr(
        kpoints_module,
        "compute_centroid",
        lambda *args, **kwargs: {
            "display_figures": [],
            "b_matrix_input": np.eye(3),
        },
    )

    with pytest.raises(KeyError, match="sp_path"):
        kpoints_module.KPathBuilder().interactive_build()
    output = capsys.readouterr().out

    assert input_stream.tell() == 0
    assert "Auto path generation failed" not in output
    assert "Falling back to manual file input" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_missing_centroid_aborts_without_manual_k_prompt(
    tmp_path, monkeypatch, capsys
):
    """A malformed centroid result must not trigger basis-ambiguous input."""
    from alterseek import kpoints as kpoints_module

    (tmp_path / "POSCAR").write_text(
        "test structure placeholder\n", encoding="utf-8"
    )
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        'moments = ""\n',
        encoding="utf-8",
    )
    centroid_result = {
        "display_figures": [],
        "sp_path": [("GAMMA", "X")],
        "sp_point_coords": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
        },
        "b_matrix": np.eye(3),
        "b_matrix_input": np.eye(3),
        "sc_type": "oP1",
        "point_group": "mmm",
        "spacegroup": 61,
    }

    def forbid_spin_analysis(*args, **kwargs):
        raise AssertionError("no-moments route must not run spin analysis")

    input_stream = io.StringIO("")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", input_stream)
    monkeypatch.setattr(kpoints_module, "find_sf_run", forbid_spin_analysis)
    monkeypatch.setattr(
        kpoints_module,
        "prepare_submitted_cell_analysis",
        lambda *a, **k: _synthetic_analysis_preparation(),
    )
    monkeypatch.setattr(
        kpoints_module, "compute_centroid", lambda *a, **k: centroid_result
    )

    assert kpoints_module.KPathBuilder().interactive_build() is False
    output = capsys.readouterr().out

    assert input_stream.tell() == 0
    assert "IBZ centroid is unavailable. Aborting." in output
    assert "Enter k-point" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_workflow_gates_on_findspingroup_g0(tmp_path, monkeypatch):
    """The gate uses G0 while the path remains in the submitted lattice.

    No assertion on the gate helper alone catches a workflow that computes the
    right thing and then feeds the gate something else.
    """
    pytest.importorskip("findspingroup")
    pytest.importorskip("ase")
    from alterseek import kpoints as kpoints_module
    from alterseek import submitted_cell_analysis

    if submitted_cell_analysis.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    seen = []
    real_gate = kpoints_module._altermagnetism_gate

    def spy_gate(sf_result, working_cell_symmetry=None):
        seen.append(working_cell_symmetry)
        return real_gate(sf_result, working_cell_symmetry)

    monkeypatch.setattr(kpoints_module, "_altermagnetism_gate", spy_gate)

    answers = "\n".join([str(POSCAR), "0 0 1", "1 -1 1 -1", "", ""]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    kpoints_module.KPathBuilder().interactive_build()

    assert seen, "the gate was never consulted"
    working = seen[0]
    assert working is not None, "the gate ignored FindSpinGroup G0"
    # A G0 description, not a re-detection from coordinates.
    assert set(working) == {'label', 'spacegroup_number', 'point_group',
                            'laue_group', 'sites'}
    assert laue_group_from_spacegroup_number(working['spacegroup_number']) == \
        working['laue_group']
    # This row describes the retained FindSpinGroup magnetic reference.
    assert working['sites'] == 12


def test_figure_basename_comes_from_the_submitted_structure():
    """All figures use the submitted structure's basename."""
    from alterseek.kpoints import _figure_basename

    assert _figure_basename("POSCAR") == "POSCAR"
    assert _figure_basename("1.0.47_MnSe2.mcif") == "1.0.47_MnSe2"
    assert _figure_basename("/tmp/case/SUPERCELL_211.vasp") == "SUPERCELL_211"
    assert _figure_basename(None) is None


def _run_workflow_on(tmp_path, monkeypatch, structure, moments, mode_2d=False):
    """Run the whole workflow from a settings file in a fresh folder."""
    pytest.importorskip("findspingroup")
    pytest.importorskip("ase")
    from alterseek import run_workflow

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "POSCAR").write_bytes(Path(structure).read_bytes())
    settings = (
        'structure = "POSCAR"\n'
        'spin_axis = "0 0 1"\n'
        f'moments = "{moments}"\n'
        'flip_option = 1\n'
        'output_code = "vasp"\n'
    )
    if mode_2d:
        settings += 'vacuum_axis = "c"\n'
    (run_dir / "alterseek_input.toml").write_text(settings, encoding="utf-8")
    monkeypatch.chdir(run_dir)
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))
    return run_workflow(mode_2d=mode_2d), run_dir


@pytest.mark.parametrize(
    "mode_2d, repeat", [(True, "2 x 2 in the plane"), (False, "2 x 2 x 2")]
)
def test_workflow_stops_when_no_spin_flip_operation_fits_the_cell(
    tmp_path, monkeypatch, capsys, mode_2d, repeat
):
    """Square d-wave 2x1: the symbol and MSG describe the magnet, the counts the cell."""
    success, run_dir = _run_workflow_on(
        tmp_path, monkeypatch, REF_DIR / "square_dwave_2x1_POSCAR",
        "1 1 -1 -1 1 1 -1 -1", mode_2d=mode_2d,
    )

    assert success is False
    output = capsys.readouterr().out
    lines = output.splitlines()
    start = lines.index("Phase: AFM(Altermagnet)")
    mcif = Path("alterseek_output") / "POSCAR_magnetic_primitive.mcif"
    assert lines[start:start + 10] == [
        "Phase: AFM(Altermagnet)",
        "Oriented SSG: 123.47.1.1.L",
        "SSG Symbol (Chen-Liu): P -1|4/ 1|m 1|m -1|m infinity_{001}m|1",
        "MSG without SOC: P4'/mmm' (BNS 123.343), Type III",
        "Spin operations in the input cell: 0 flip, 8 preserve",
        "[Error] The structure is an altermagnet, but none of its 8 spin-flip "
        "point operations maps the input cell onto itself:",
        "  C2 [1 -1 0], C2 [1 1 0], C4+ [0 0 1], C4- [0 0 1], S4+ [0 0 1], "
        "S4- [0 0 1], mirror m (1 -1 0), mirror m (1 1 0).",
        "  Rotation axis/mirror plane indices are in the reciprocal basis of "
        "the magnetic primitive cell.",
        "Only an operation that maps the cell onto itself relates all bands "
        "folded into k and k' in this cell.",
        f"Use the magnetic primitive cell (4 atoms, written to {mcif}), or a "
        f"supercell of it that keeps these operations, such as {repeat}. "
        "Aborting.",
    ]
    assert "Step 1" not in output
    assert "inconsistent" not in output
    assert (run_dir / mcif).exists()
    assert not (run_dir / "KPOINTS_alter").exists()


def test_workflow_labels_the_input_cell_counts_and_writes_the_path(
    tmp_path, monkeypatch, capsys
):
    """MnF2 2x1x1 keeps four of its eight spin-flip operations."""
    success, run_dir = _run_workflow_on(
        tmp_path, monkeypatch, REF_DIR / "case12_211_supercell.vasp",
        "5 -5 5 -5",
    )

    assert success is True
    output = capsys.readouterr().out
    assert (
        "SSG Symbol (Chen-Liu): P -1|4_{2}/ 1|m -1|n 1|m infinity_{001}m|1"
        in output
    )
    assert "MSG without SOC: P4_2'/mn'm (BNS 136.498), Type III\n" in output
    assert (
        "Spin operations in the input cell: 4 flip, 4 preserve\n"
        "[Note] 4 of the 8 spin-flip operations of the magnetic primitive "
        "cell are not symmetries of this supercell, so they are not used for "
        "k': C4+ [0 0 1], C4- [0 0 1], S4+ [0 0 1], S4- [0 0 1].\n"
    ) in output
    assert "Pn'n'm" not in output
    assert (run_dir / "KPOINTS_alter").exists()


@pytest.mark.parametrize(
    "outside, in_cell, expected",
    [
        (
            ["C2 [0 1 0]"], 1,
            "[Note] 1 of the 2 spin-flip operations of the magnetic primitive "
            "cell is not a symmetry of this supercell, so it is not used for "
            "k': C2 [0 1 0].",
        ),
        # None fits: the stop message explains it instead.
        (["C2 [0 1 0]"], 0, None),
        ([], 2, None),
    ],
)
def test_unused_flip_operations_note(outside, in_cell, expected):
    from alterseek.kpoints import _unused_flip_operations_note

    analysis = {"magnet_operations": {
        "spin_flip_outside_cell": outside,
        "spin_flip_point_operations": len(outside) + in_cell,
        "spin_flip_point_operations_in_cell": in_cell,
    }}
    assert _unused_flip_operations_note(analysis) == expected


def test_workflow_prints_a_primitive_cell_as_before(
    tmp_path, monkeypatch, capsys
):
    success, _ = _run_workflow_on(
        tmp_path, monkeypatch, REF_DIR / "case12_POSCAR", "5 -5"
    )

    assert success is True
    output = capsys.readouterr().out
    assert "MSG without SOC: P4_2'/mn'm (BNS 136.498), Type III\n" in output
    assert "\nSpin operations: 8 flip, 8 preserve\n" in output
    assert "in the input cell" not in output
    assert "axes of the magnetic primitive cell" not in output
    assert "not symmetries of this supercell" not in output
