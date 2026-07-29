"""Regression tests for the altermagnetism Laue-group gate.

The gate decides whether a structure can host altermagnetic splitting at all;
when it fires, the ordinary structural path is written and no butterfly path is
built. Getting it right needs two things that MnSe2 (MAGNDATA 1.0.47) shows are
independent:

1. It must judge the cell the path is built in. MnSe2 is deposited in a cubic
   Pa-3 (205) parent whose Laue group m-3 forbids altermagnetism outright, so
   judging the submitted cell discards a genuine altermagnet.

2. It must judge that cell by its *magnetic* symmetry. MnSe2's magnetic
   primitive cell is a 3x1x1 supercell of the same cubic crystal -- strip the
   moments and spglib still finds the 12-atom cubic primitive. Re-detecting
   symmetry from its coordinates therefore reproduces the forbidding m-3, and
   does so tolerance-dependently: at spglib's default symprec the mcif's
   5-decimal rounding hides the 3-fold axis and it looks orthorhombic, which
   is an accident, not a result.

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


def test_submitted_cell_decides_when_no_magnetic_cell_was_built():
    """Under --parent-setting the submitted cell is the only cell."""
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


def test_cell_rows_place_a_note_under_the_first_row(capsys):
    from alterseek.kpoints import _print_cell_rows

    _print_cell_rows([('Nonmagnetic primitive cell:', 'Pa-3 (205)', 'm-3', 'm-3', '[12 atoms]')],
                     note_after_first="recovered from the input cell")
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2
    # Indented into the field column, not hanging off the label.
    assert lines[1].startswith(" " * 30)
    assert lines[1].strip() == "recovered from the input cell"


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


@pytest.mark.skipif(not POSCAR.exists(), reason="SSG test input not present")
def test_workflow_gates_on_the_constructed_cells_g0(tmp_path, monkeypatch):
    """The wiring: the gate's symmetry must come from G0 of the built cell.

    No assertion on the gate helper alone catches a workflow that computes the
    right thing and then feeds the gate something else.
    """
    pytest.importorskip("findspingroup")
    pytest.importorskip("ase")
    from alterseek import kpoints as kpoints_module
    from alterseek import ssg_setting

    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
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

    kpoints_module.KPointsModifier(magnetic_setting=True).interactive_modify()

    assert seen, "the gate was never consulted"
    working = seen[0]
    assert working is not None, "the gate judged the submitted cell, not the built one"
    # A G0 description, not a re-detection from coordinates.
    assert set(working) == {'label', 'spacegroup_number', 'point_group',
                            'laue_group', 'sites'}
    assert laue_group_from_spacegroup_number(working['spacegroup_number']) == \
        working['laue_group']
    # The site count must come from the constructed cell, not the input cell.
    assert working['sites'] == 12


def test_figure_basename_comes_from_the_submitted_structure():
    """Figure 1 must not be named after the internal helper cell.

    The magnetic route computes the centroid from a helper structure written
    under a derived name, which used to leak into Figure 1's filename and
    title while Figures 2-4 used the submitted name.
    """
    from alterseek.kpoints import _figure_basename

    assert _figure_basename("POSCAR") == "POSCAR"
    assert _figure_basename("1.0.47_MnSe2.mcif") == "1.0.47_MnSe2"
    assert _figure_basename("/tmp/case/SUPERCELL_211.vasp") == "SUPERCELL_211"
    assert _figure_basename(None) is None
