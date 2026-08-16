"""Validation and exit-status tests for alterseek_input.toml and the CLI."""

import io
import os
import sys

import pytest

from alterseek.kpoints import KPointsModifier, _read_input_config
from alterseek.find_sf_operations import SpinSymmetryError


@pytest.mark.parametrize(
    "content, message",
    [
        ('output_code = "wien2k"\n', "output_code"),
        ('flip_option = 0\n', "flip_option"),
        ('save_pdf = "yes"\n', "save_pdf"),
        ('view_elev = 25\n', "supplied together"),
        ('view_elev = true\nview_azim = 20\n', "view_elev must be a number"),
        ('view_elev = "10.5"\nview_azim = 20\n', "view_elev must be a number"),
        ('structure = "POSCAR"\nstrcture = "typo"\n', "unknown setting"),
        ('structure = ["POSCAR"]\n', "structure"),
    ],
)
def test_invalid_input_config_is_rejected(tmp_path, content, message):
    path = tmp_path / "alterseek_input.toml"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        _read_input_config(str(path))


def test_malformed_toml_is_rejected(tmp_path):
    path = tmp_path / "alterseek_input.toml"
    path.write_text('structure = "unterminated\n', encoding="utf-8")
    with pytest.raises(ValueError, match="Failed to read"):
        _read_input_config(str(path))


def test_spin_symmetry_failure_aborts_workflow(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    structure = tmp_path / "POSCAR"
    structure.write_text("placeholder\n", encoding="utf-8")
    answers = f"{structure}\n0 0 1\n1 -1\n"
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))
    def fail_spin_symmetry(*args, **kwargs):
        raise SpinSymmetryError("synthetic spin-symmetry failure")

    monkeypatch.setattr("alterseek.kpoints.find_sf_run", fail_spin_symmetry)
    assert KPointsModifier().interactive_modify() is False
    output = capsys.readouterr().out
    assert "Spin-symmetry analysis failed" in output
    assert "synthetic spin-symmetry failure" in output
    assert "Step 1" not in output


def test_cli_returns_nonzero_for_missing_structure(tmp_path, monkeypatch):
    from alterseek_path import main

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["alterseek-path"])
    monkeypatch.setattr(sys, "stdin", io.StringIO("missing.vasp\n"))
    assert main() == 1


def test_save_pdf_does_not_mutate_process_environment(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("ALTERSEEK_BZ_EXTRA_FORMATS", raising=False)
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "missing.vasp"\nsave_pdf = true\n',
        encoding="utf-8",
    )

    assert KPointsModifier().interactive_modify() is False
    assert "ALTERSEEK_BZ_EXTRA_FORMATS" not in os.environ
