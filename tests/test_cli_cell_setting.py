"""The CLI contract for which cell the k-path is built in.

The magnetic primitive cell is the default: it is the cell that reflects the
symmetry of the magnetic state, and defaulting to the nonmagnetic parent
silently produces a path for a cell the magnetic structure does not have. The
parent cell remains reachable through --parent-setting, which is useful for
comparing several magnetic orders against one fixed reference path.
"""
import sys

import pytest

import alterseek_path


class _Recorder:
    """Stands in for KPointsModifier so no workflow actually runs."""

    instances = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _Recorder.instances.append(self)

    def interactive_modify(self):
        return True


@pytest.fixture
def recorded(monkeypatch):
    _Recorder.instances = []
    monkeypatch.setattr(alterseek_path, "KPointsModifier", _Recorder)
    return _Recorder.instances


def _run(monkeypatch, recorded, argv):
    monkeypatch.setattr(sys, "argv", ["alterseek-path", *argv])
    assert alterseek_path.main() == 0
    assert len(recorded) == 1
    return recorded[0].kwargs


def test_magnetic_primitive_cell_is_the_default(monkeypatch, recorded):
    assert _run(monkeypatch, recorded, [])["magnetic_setting"] is True


def test_parent_setting_opts_out(monkeypatch, recorded):
    assert _run(monkeypatch, recorded, ["--parent-setting"])["magnetic_setting"] is False


def test_ssg_setting_flag_no_longer_exists(monkeypatch, recorded):
    """It named the behavior that is now the default, so it was removed
    outright rather than kept as a silent alias."""
    monkeypatch.setattr(sys, "argv", ["alterseek-path", "--ssg-setting"])
    with pytest.raises(SystemExit) as excinfo:
        alterseek_path.main()
    assert excinfo.value.code == 2
    assert not recorded


def test_default_constructor_matches_the_default_cli(monkeypatch, recorded):
    """Library callers get the same default as CLI users."""
    from alterseek.kpoints import KPointsModifier
    import inspect

    default = inspect.signature(KPointsModifier).parameters["magnetic_setting"].default
    assert default is True
    assert default is _run(monkeypatch, recorded, [])["magnetic_setting"]


def test_other_flags_still_route_through(monkeypatch, recorded):
    kwargs = _run(monkeypatch, recorded, ["--2d", "--vacuum-axis", "a", "--output", "verbose"])
    assert kwargs["mode_2d"] is True
    assert kwargs["input_vacuum_axis"] == 0
    assert kwargs["output_verbose"] is True
    assert kwargs["magnetic_setting"] is True
