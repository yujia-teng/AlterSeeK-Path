"""Failure-contract tests for spin-symmetry analysis."""

import pytest

from alterseek import find_sf_operations as spin_ops


def test_missing_structure_raises_spin_symmetry_error(tmp_path):
    missing = tmp_path / "missing.vasp"

    with pytest.raises(
        spin_ops.SpinSymmetryError,
        match="Structure file .* was not found",
    ):
        spin_ops.run(
            str(missing),
            "1 -1",
            verbose=False,
            spin_axis_cart="0 0 1",
        )


def test_structure_parser_failure_raises_spin_symmetry_error(
    tmp_path, monkeypatch
):
    structure = tmp_path / "POSCAR"
    structure.write_text("placeholder\n", encoding="utf-8")

    def fail_read(*args, **kwargs):
        raise ValueError("synthetic parser failure")

    monkeypatch.setattr(spin_ops, "read", fail_read)

    with pytest.raises(
        spin_ops.SpinSymmetryError,
        match="Could not read structure file.*synthetic parser failure",
    ):
        spin_ops.run(
            str(structure),
            "1 -1",
            verbose=False,
            spin_axis_cart="0 0 1",
        )


def test_unexpected_late_failure_is_wrapped_with_its_cause(monkeypatch):
    original = OSError("synthetic output failure")

    def fail_internal(*args, **kwargs):
        raise original

    monkeypatch.setattr(spin_ops, "_run", fail_internal)

    with pytest.raises(
        spin_ops.SpinSymmetryError,
        match="Unexpected spin-symmetry analysis failure.*synthetic output failure",
    ) as exc_info:
        spin_ops.run("POSCAR", "1 -1", verbose=False)

    assert exc_info.value.__cause__ is original
