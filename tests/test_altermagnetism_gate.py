"""Regression tests for the altermagnetism Laue-group gate.

The gate decides whether a structure can host altermagnetic splitting at all;
when it fires, the ordinary structural path is written and no butterfly path is
built. It must judge the cell the path is actually constructed in.

MnSe2 (MAGNDATA 1.0.47) is the case that exposed the difference. It is
deposited in a cubic Pa-3 (No. 205) parent whose Laue group m-3 forbids
altermagnetism outright, but its magnetic primitive cell -- the cell
AlterSeeK-Path builds the path in -- is orthorhombic with Laue group mmm, which
permits it. Reading the submitted cell there discards a genuine altermagnet.
"""
import io
import sys
from pathlib import Path

import numpy as np
import pytest

from alterseek.kpoints import _altermagnetism_gate

REF_DIR = Path(__file__).parent / "references"
POSCAR = REF_DIR / "SUPERCELL_211.vasp"


CUBIC_PARENT = {'point_group': 'm-3', 'laue_group': 'm-3'}
ORTHORHOMBIC_MAGNETIC_CELL = {'point_group': 'mmm', 'laue_group': 'mmm'}


def test_submitted_cell_decides_when_no_conversion_happened():
    """Without a magnetic primitive cell, the submitted cell is the only cell."""
    assert _altermagnetism_gate(CUBIC_PARENT) is not None
    assert _altermagnetism_gate(ORTHORHOMBIC_MAGNETIC_CELL) is None


def test_magnetic_primitive_cell_overrides_a_forbidding_parent():
    """The MnSe2 case: cubic m-3 parent, orthorhombic mmm magnetic cell.

    Judging the parent would report 'no altermagnetism' and fall back to the
    ordinary path, losing the altermagnetic result entirely.
    """
    assert _altermagnetism_gate(
        CUBIC_PARENT, working_cell_symmetry=ORTHORHOMBIC_MAGNETIC_CELL
    ) is None


def test_magnetic_primitive_cell_can_also_forbid_a_permitting_parent():
    """The gate is not a one-way override -- the working cell simply decides."""
    reason = _altermagnetism_gate(
        ORTHORHOMBIC_MAGNETIC_CELL, working_cell_symmetry=CUBIC_PARENT
    )
    assert reason is not None
    assert reason['laue_group'] == 'm-3'


@pytest.mark.parametrize("laue_group", ['-1', '-3', 'm-3'])
def test_all_forbidden_laue_groups_fire(laue_group):
    reason = _altermagnetism_gate({'point_group': laue_group, 'laue_group': laue_group})
    assert reason is not None
    assert reason['laue_group'] == laue_group


def _write_poscar(path, lattice, symbol, positions):
    lines = [symbol, "1.0"]
    lines += [" ".join(f"{value:.10f}" for value in row) for row in lattice]
    lines += [symbol, str(len(positions)), "Direct"]
    lines += [" ".join(f"{value:.10f}" for value in row) for row in positions]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_non_magnetic_symmetry_of_file_reads_real_laue_groups(tmp_path):
    """The helper feeding the gate must report the cell's own Laue group."""
    pytest.importorskip("ase")
    from alterseek.find_sf_operations import non_magnetic_symmetry_of_file

    # Simple cubic: Laue m-3m.
    cubic = tmp_path / "cubic.vasp"
    _write_poscar(cubic, np.eye(3) * 4.0, "Po", [[0.0, 0.0, 0.0]])
    assert non_magnetic_symmetry_of_file(str(cubic))['laue_group'] == 'm-3m'

    # Same cell stretched along c: Laue drops to 4/mmm, not cubic any more.
    tetragonal = tmp_path / "tetragonal.vasp"
    _write_poscar(tetragonal, np.diag([4.0, 4.0, 6.0]), "Po", [[0.0, 0.0, 0.0]])
    assert non_magnetic_symmetry_of_file(str(tetragonal))['laue_group'] == '4/mmm'

    # Three distinct axes: orthorhombic, Laue mmm -- the MnSe2 magnetic-cell class.
    orthorhombic = tmp_path / "orthorhombic.vasp"
    _write_poscar(orthorhombic, np.diag([4.0, 5.0, 6.0]), "Po", [[0.0, 0.0, 0.0]])
    assert non_magnetic_symmetry_of_file(str(orthorhombic))['laue_group'] == 'mmm'


def test_non_magnetic_symmetry_of_file_returns_none_on_unreadable_input(tmp_path):
    from alterseek.find_sf_operations import non_magnetic_symmetry_of_file

    missing = tmp_path / "does-not-exist.vasp"
    assert non_magnetic_symmetry_of_file(str(missing)) is None


@pytest.mark.skipif(not POSCAR.exists(), reason="SSG test input not present")
def test_workflow_gates_on_the_magnetic_primitive_cell(tmp_path, monkeypatch):
    """The wiring: the gate's symmetry must be read from the constructed cell.

    Reading it from the submitted structure instead is the MnSe2 bug, and no
    assertion on the gate helper alone can catch that -- only checking which
    file the workflow hands it does.
    """
    pytest.importorskip("findspingroup")
    pytest.importorskip("ase")
    from alterseek import kpoints as kpoints_module
    from alterseek import ssg_setting

    if ssg_setting.find_spin_group_acc_primitive_from_data is None:
        pytest.skip("findspingroup acc-primitive setting unavailable")

    prepared = {}
    real_prepare = kpoints_module.prepare_magnetic_setting_files

    def spy_prepare(*args, **kwargs):
        result = real_prepare(*args, **kwargs)
        prepared["real_poscar_path"] = result["real_poscar_path"]
        return result

    inspected = []
    real_symmetry = kpoints_module.non_magnetic_symmetry_of_file

    def spy_symmetry(structure_file):
        inspected.append(structure_file)
        return real_symmetry(structure_file)

    monkeypatch.setattr(kpoints_module, "prepare_magnetic_setting_files", spy_prepare)
    monkeypatch.setattr(kpoints_module, "non_magnetic_symmetry_of_file", spy_symmetry)

    answers = "\n".join([str(POSCAR), "0 0 1", "1 -1 1 -1", "", ""]) + "\n"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    kpoints_module.KPointsModifier(magnetic_setting=True).interactive_modify()

    assert inspected, "the gate never inspected any cell's symmetry"
    assert inspected == [prepared["real_poscar_path"]], (
        "the gate must read the magnetic primitive cell, not the submitted structure"
    )
    assert str(POSCAR) not in inspected
