"""Tests for the unified alterseek-plot command."""
import pytest

import alterseek_plot


def test_explicit_code_forwards_remaining_arguments(monkeypatch):
    seen = []

    def record(code, argv):
        seen.append((code, argv))
        return 0

    monkeypatch.setattr(alterseek_plot, "_run_plotter", record)

    assert alterseek_plot.main(["qe", "--config", "custom.toml"]) == 0
    assert seen == [("qe", ["--config", "custom.toml"])]


@pytest.mark.parametrize(
    ("filename", "expected_code"),
    [
        ("alterseek_plot_vasp.toml", "vasp"),
        ("alterseek_plot_qe.toml", "qe"),
        ("alterseek_plot_abinit.toml", "abinit"),
    ],
)
def test_bare_command_detects_one_plot_config(
    tmp_path, monkeypatch, filename, expected_code
):
    (tmp_path / filename).write_text("", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    seen = []
    monkeypatch.setattr(
        alterseek_plot,
        "_run_plotter",
        lambda code, argv: seen.append((code, argv)) or 0,
    )

    assert alterseek_plot.main([]) == 0
    assert seen == [(expected_code, [])]


def test_detected_code_receives_bare_command_options(tmp_path, monkeypatch):
    (tmp_path / "alterseek_plot_vasp.toml").write_text("", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    seen = []
    monkeypatch.setattr(
        alterseek_plot,
        "_run_plotter",
        lambda code, argv: seen.append((code, argv)) or 0,
    )

    assert alterseek_plot.main(["-o", "bands.pdf"]) == 0
    assert seen == [("vasp", ["-o", "bands.pdf"])]


def test_bare_command_shows_help_when_no_config_exists(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)

    assert alterseek_plot.main([]) == 0
    captured = capsys.readouterr()
    assert "usage: alterseek-plot {vasp,qe,abinit}" in captured.out
    assert captured.err == ""


def test_bare_command_rejects_multiple_code_configs(
    tmp_path, monkeypatch, capsys
):
    (tmp_path / "alterseek_plot_vasp.toml").write_text("", encoding="utf-8")
    (tmp_path / "alterseek_plot_qe.toml").write_text("", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    assert alterseek_plot.main([]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Multiple plot configurations found (vasp, qe)" in captured.err
    assert "specify a code" in captured.err


def test_unknown_code_reports_error(capsys):
    assert alterseek_plot.main(["unknown"]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "[Error] Unknown code: unknown" in captured.err
    assert "usage: alterseek-plot" in captured.err


def test_options_without_detectable_config_require_code(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.chdir(tmp_path)

    assert alterseek_plot.main(["-o", "bands.pdf"]) == 2
    assert "No plot configuration found; specify vasp, qe, or abinit" in (
        capsys.readouterr().err
    )


@pytest.mark.parametrize("code", ["vasp", "qe", "abinit"])
def test_backend_help_keeps_subcommand_in_usage(code, capsys):
    with pytest.raises(SystemExit) as exc_info:
        alterseek_plot.main([code, "--help"])

    assert exc_info.value.code == 0
    assert f"usage: alterseek-plot {code}" in capsys.readouterr().out
