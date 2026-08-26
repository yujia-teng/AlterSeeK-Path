"""The submitted-cell workflow has no alternate cell-setting CLI mode."""
import sys

import pytest

import alterseek_path


class _Recorder:
    """Stands in for KPointsModifier so no workflow actually runs."""

    instances = []
    result = True
    error = None

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _Recorder.instances.append(self)

    def interactive_modify(self):
        if self.error is not None:
            raise self.error
        return self.result


@pytest.fixture
def recorded(monkeypatch):
    _Recorder.instances = []
    _Recorder.result = True
    _Recorder.error = None
    monkeypatch.setattr(alterseek_path, "KPointsModifier", _Recorder)
    monkeypatch.setattr(alterseek_path, "KPointsModifier2D", _Recorder)
    return _Recorder.instances


def _run(monkeypatch, recorded, argv):
    monkeypatch.setattr(sys, "argv", ["alterseek-path", *argv])
    assert alterseek_path.main() == 0
    assert len(recorded) == 1
    return recorded[0].kwargs


def test_ssg_setting_flag_no_longer_exists(monkeypatch, recorded):
    """It named the behavior that is now the default, so it was removed
    outright rather than kept as a silent alias."""
    monkeypatch.setattr(sys, "argv", ["alterseek-path", "--ssg-setting"])
    with pytest.raises(SystemExit) as excinfo:
        alterseek_path.main()
    assert excinfo.value.code == 2
    assert not recorded


def test_other_flags_still_route_through(monkeypatch, recorded):
    kwargs = _run(monkeypatch, recorded, ["--2d", "--vacuum-axis", "a"])
    assert kwargs["input_vacuum_axis"] == 0
    assert "mode_2d" not in kwargs
    assert "magnetic_setting" not in kwargs


def test_cli_returns_failure_when_workflow_reports_failure(monkeypatch, recorded):
    _Recorder.result = False
    monkeypatch.setattr(sys, "argv", ["alterseek-path"])

    assert alterseek_path.main() == 1
    assert len(recorded) == 1


def test_cli_reports_unexpected_workflow_failure_once(
    monkeypatch, recorded, capsys
):
    _Recorder.error = RuntimeError("synthetic workflow failure")
    monkeypatch.setattr(sys, "argv", ["alterseek-path"])

    assert alterseek_path.main() == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.count("synthetic workflow failure") == 1
    assert "[Error] AlterSeeK-Path failed:" in captured.err
    assert "Traceback" not in captured.err


def test_cli_does_not_catch_keyboard_interrupt(monkeypatch, recorded):
    _Recorder.error = KeyboardInterrupt()
    monkeypatch.setattr(sys, "argv", ["alterseek-path"])

    with pytest.raises(KeyboardInterrupt):
        alterseek_path.main()
