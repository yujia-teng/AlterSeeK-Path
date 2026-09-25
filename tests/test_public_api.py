"""The public entry points: names, the vacuum-axis letters, output locations."""
import shutil
from pathlib import Path

import pytest

pytest.importorskip("findspingroup")
pytest.importorskip("ase")

from alterseek import workflow
from alterseek.constants import _normalize_vacuum_axis

POSCAR = Path(__file__).parent / "references" / "case12_POSCAR"
MOMENTS = "5 -5"


# --- the vacuum-axis helper ------------------------------------------------

@pytest.mark.parametrize(
    "value, index",
    [("a", 0), ("b", 1), ("c", 2), ("C", 2), (" c ", 2), ("A", 0)],
)
def test_axis_letters_become_cell_vector_indices(value, index):
    assert _normalize_vacuum_axis(value) == index


@pytest.mark.parametrize("value", ["d", "", "ab", 0, 1, 2, 3, -1, True, 1.0, None])
def test_axis_rejects_anything_that_is_not_a_letter(value):
    with pytest.raises(ValueError, match='vacuum_axis must be "a", "b", or "c"'):
        _normalize_vacuum_axis(value)


def test_axis_none_is_accepted_only_where_it_is_allowed():
    assert _normalize_vacuum_axis(None, allow_none=True) is None
    with pytest.raises(ValueError):
        _normalize_vacuum_axis(None)
    with pytest.raises(ValueError):
        _normalize_vacuum_axis(None, allow_index=True)


@pytest.mark.parametrize("value", [0, 1, 2])
def test_axis_index_is_accepted_only_where_it_is_allowed(value):
    assert _normalize_vacuum_axis(value, allow_index=True) == value
    with pytest.raises(ValueError):
        _normalize_vacuum_axis(value)


@pytest.mark.parametrize("value", [3, -1, True, 1.0, "d", None])
def test_axis_index_form_still_rejects_bad_values(value):
    with pytest.raises(ValueError):
        _normalize_vacuum_axis(value, allow_index=True)


def test_axis_letters_are_still_accepted_in_the_index_form():
    assert _normalize_vacuum_axis("b", allow_index=True) == 1


# --- run_workflow ----------------------------------------------------------

@pytest.fixture
def builders(monkeypatch, tmp_path):
    """Replace both builders so no calculation runs."""
    monkeypatch.chdir(tmp_path)
    made = []

    class _Recorder:
        tag = None
        result = True
        error = None

        def __init__(self, **kwargs):
            made.append((type(self).tag, kwargs))

        def interactive_build(self):
            print("workflow transcript line")
            if _Recorder.error is not None:
                raise _Recorder.error
            return _Recorder.result

    class _Recorder3D(_Recorder):
        tag = "3d"

    class _Recorder2D(_Recorder):
        tag = "2d"

    _Recorder.result = True
    _Recorder.error = None
    monkeypatch.setattr(workflow, "KPathBuilder", _Recorder3D)
    monkeypatch.setattr(workflow, "KPathBuilder2D", _Recorder2D)
    return made, _Recorder


def test_run_workflow_uses_the_three_dimensional_builder(builders):
    made, _ = builders

    assert workflow.run_workflow() is True
    assert made == [("3d", {"input_vacuum_axis": None})]


def test_run_workflow_uses_the_two_dimensional_builder(builders):
    made, _ = builders

    assert workflow.run_workflow(mode_2d=True) is True
    assert made == [("2d", {"input_vacuum_axis": None})]


@pytest.mark.parametrize("axis, index", [("a", 0), ("b", 1), ("c", 2)])
def test_run_workflow_hands_the_axis_index_to_the_builder(builders, axis, index):
    made, _ = builders

    workflow.run_workflow(mode_2d=True, vacuum_axis=axis)
    assert made == [("2d", {"input_vacuum_axis": index})]


def test_run_workflow_rejects_an_invalid_axis(builders):
    made, _ = builders

    with pytest.raises(ValueError, match='vacuum_axis must be'):
        workflow.run_workflow(vacuum_axis="d")
    assert made == []


def test_no_axis_leaves_the_toml_setting_in_charge(builders):
    from alterseek.kpoints import KPathBuilder

    made, _ = builders
    workflow.run_workflow()

    assert made[0][1]["input_vacuum_axis"] is None
    assert KPathBuilder(input_vacuum_axis=None)._vacuum_axis_from_cli is False


def test_a_successful_run_saves_the_run_log(builders, tmp_path):
    workflow.run_workflow()

    log = tmp_path / "alterseek_output" / "alterseek_run.log"
    assert log.exists()
    assert "workflow transcript line" in log.read_text(encoding="utf-8")


def test_a_failed_run_leaves_the_previous_run_log_alone(builders, tmp_path):
    _, recorder = builders
    log = tmp_path / "alterseek_output" / "alterseek_run.log"
    log.parent.mkdir()
    log.write_text("previous successful run\n", encoding="utf-8")
    recorder.result = False

    assert workflow.run_workflow() is False
    assert log.read_text(encoding="utf-8") == "previous successful run\n"


def test_an_unexpected_failure_propagates(builders, tmp_path):
    _, recorder = builders
    recorder.error = RuntimeError("synthetic workflow failure")

    with pytest.raises(RuntimeError, match="synthetic workflow failure"):
        workflow.run_workflow()
    assert not (tmp_path / "alterseek_output" / "alterseek_run.log").exists()


# --- general_k_path axis ---------------------------------------------------

@pytest.mark.parametrize(
    "axis, index", [("a", 0), ("b", 1), ("c", 2), ("C", 2), (" c ", 2)]
)
def test_general_k_path_accepts_a_letter_and_forwards_its_index(
    tmp_path, monkeypatch, axis, index
):
    from alterseek import kpoints as kpoints_module

    seen = {}

    class _Stop(RuntimeError):
        pass

    def _record(*args, **kwargs):
        seen.update(kwargs)
        raise _Stop("stop before any calculation")

    monkeypatch.setattr(
        kpoints_module, "prepare_submitted_cell_analysis", _record
    )

    with pytest.raises(_Stop):
        kpoints_module.general_k_path(
            str(POSCAR),
            mode_2d=True,
            vacuum_axis=axis,
            output_dir=str(tmp_path / "out"),
        )

    assert seen["input_vacuum_axis"] == index


@pytest.mark.parametrize("axis", ["d", "", 0, 2, 3, True, 1.0, None])
def test_general_k_path_rejects_an_invalid_axis(tmp_path, axis):
    from alterseek import general_k_path

    with pytest.raises(ValueError, match='vacuum_axis must be "a", "b", or "c"'):
        general_k_path(
            str(POSCAR), vacuum_axis=axis, output_dir=str(tmp_path / "out")
        )


# --- output locations ------------------------------------------------------

def _working_folder_entries(tmp_path):
    return sorted(entry.name for entry in tmp_path.iterdir())


def test_spin_symmetry_writes_into_the_shared_output_folder(
    tmp_path, monkeypatch
):
    from alterseek import spin_symmetry

    monkeypatch.chdir(tmp_path)
    shutil.copy(POSCAR, tmp_path / "POSCAR")

    spin_symmetry("POSCAR", MOMENTS)

    assert (tmp_path / "alterseek_output" / "spin_operations.txt").exists()
    assert _working_folder_entries(tmp_path) == ["POSCAR", "alterseek_output"]


def test_brillouin_zone_writes_into_the_shared_output_folder(
    tmp_path, monkeypatch
):
    from alterseek import brillouin_zone

    monkeypatch.chdir(tmp_path)
    shutil.copy(POSCAR, tmp_path / "POSCAR")

    brillouin_zone("POSCAR", show_plot=False)

    written = sorted(
        entry.name for entry in (tmp_path / "alterseek_output").iterdir()
    )
    assert "POSCAR_seekpath_basis_mapping.txt" in written
    assert _working_folder_entries(tmp_path) == ["POSCAR", "alterseek_output"]


def test_brillouin_zone_uses_the_magnetic_zone_of_an_mcif(tmp_path):
    """La2NiO4: the crystal alone is tP1, the magnetic structure is oP1."""
    from alterseek import brillouin_zone

    structure = Path(__file__).parent / "references" / "La2NiO4.mcif"
    zone = brillouin_zone(
        str(structure), show_plot=False, output_dir=str(tmp_path / "out")
    )
    assert zone["sc_type"] == "oP1"
    assert not list((tmp_path / "out").glob("*_tP1.png"))


def test_general_k_path_writes_into_the_shared_output_folder(
    tmp_path, monkeypatch
):
    from alterseek import general_k_path

    monkeypatch.chdir(tmp_path)
    shutil.copy(POSCAR, tmp_path / "POSCAR")

    general_k_path("POSCAR", moments=MOMENTS)

    assert (tmp_path / "alterseek_output" / "spin_operations.txt").exists()
    assert _working_folder_entries(tmp_path) == ["POSCAR", "alterseek_output"]


# --- package surface -------------------------------------------------------

_PUBLIC_NAMES = {
    "__version__",
    "spin_symmetry",
    "brillouin_zone",
    "general_k_path",
    "run_workflow",
    "SpinSymmetryError",
}


def test_the_package_exports_exactly_the_public_names():
    import alterseek

    assert set(alterseek.__all__) == _PUBLIC_NAMES
    assert dir(alterseek) == sorted(_PUBLIC_NAMES)


def test_every_public_name_can_be_imported():
    from alterseek import (  # noqa: F401
        SpinSymmetryError,
        brillouin_zone,
        general_k_path,
        run_workflow,
        spin_symmetry,
    )

    assert run_workflow is workflow.run_workflow


def test_star_import_brings_in_the_public_names():
    namespace = {}
    exec("from alterseek import *", namespace)

    assert _PUBLIC_NAMES <= set(namespace)
    assert "KPathBuilder" not in namespace
    assert "KPathBuilder2D" not in namespace


@pytest.mark.parametrize("name", ["KPathBuilder", "KPathBuilder2D"])
def test_the_builders_are_not_top_level_names(name):
    import alterseek

    assert name not in alterseek.__all__
    with pytest.raises(AttributeError):
        getattr(alterseek, name)
