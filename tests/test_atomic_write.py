"""Regressions for the atomic text writers.

Every generated file now goes through these. Before, ten writers built their
output incrementally with a plain ``open(..., "w")``, so a failure partway
through left a truncated POSCAR / MCIF / operation list on disk that still
looked loadable. The buffered form commits only on a clean exit.
"""
import os

import pytest

import alterseek.atomic_write as atomic_write
from alterseek.atomic_write import (
    _atomic_open_text,
    _atomic_write_text,
    _atomic_write_text_set,
)


def _temp_siblings(directory):
    """Leftover temp files use a leading dot and a .tmp suffix."""
    return [n for n in os.listdir(directory)
            if n.startswith(".") and n.endswith(".tmp")]


def _transaction_siblings(directory):
    return [
        name for name in os.listdir(directory)
        if name.startswith(".")
        and (name.endswith(".stage") or name.endswith(".backup"))
    ]


def test_buffered_write_commits_on_clean_exit(tmp_path):
    target = tmp_path / "POSCAR"
    with _atomic_open_text(str(target)) as f:
        f.write("line one\n")
        f.write("line two\n")
    # newline="\n" is preserved, so generated files stay byte-stable across
    # platforms (the golden tests compare bytes).
    assert target.read_bytes() == b"line one\nline two\n"
    assert _temp_siblings(tmp_path) == []


def test_failed_write_leaves_the_previous_file_untouched(tmp_path):
    target = tmp_path / "POSCAR"
    target.write_bytes(b"original contents\n")

    with pytest.raises(RuntimeError, match="halfway"):
        with _atomic_open_text(str(target)) as f:
            f.write("first half\n")
            raise RuntimeError("failed halfway through")

    # The old file survives intact -- not truncated, not half-rewritten.
    assert target.read_bytes() == b"original contents\n"
    assert _temp_siblings(tmp_path) == []


def test_failed_write_creates_no_file_when_there_was_none(tmp_path):
    target = tmp_path / "new_file.txt"

    with pytest.raises(ValueError):
        with _atomic_open_text(str(target)) as f:
            f.write("partial\n")
            raise ValueError("boom")

    assert not target.exists()
    assert _temp_siblings(tmp_path) == []


def test_a_real_writer_no_longer_leaves_a_truncated_file(tmp_path):
    """End-to-end on _write_poscar, not just on the helper.

    A malformed site (two coordinates instead of three) raises after the
    header and lattice rows have already been handed to the writer. With the
    old plain ``open(..., "w")`` that produced a POSCAR on disk containing a
    valid header and no atoms -- silently loadable, and wrong.
    """
    import numpy as np

    from alterseek.io import _write_poscar

    target = tmp_path / "POSCAR"
    with pytest.raises(IndexError):
        _write_poscar(
            str(target),
            "truncation probe",
            np.eye(3) * 4.0,
            ["Fe"],
            [2],
            [np.array([0.0, 0.0, 0.0]), np.array([0.5, 0.5])],  # malformed site
        )
    assert not target.exists()
    assert _temp_siblings(tmp_path) == []


def test_replacing_an_existing_file_is_all_or_nothing(tmp_path):
    target = tmp_path / "KPOINTS_alter"
    target.write_text("old\n", encoding="utf-8")
    _atomic_write_text(str(target), "new contents\n")
    assert target.read_bytes() == b"new contents\n"
    assert _temp_siblings(tmp_path) == []


def test_multi_file_commit_replaces_the_complete_set(tmp_path):
    targets = [
        tmp_path / "case_magnetic_primitive.vasp",
        tmp_path / "case_magnetic_primitive_MAGMOM.txt",
        tmp_path / "KPOINTS_alter",
    ]
    for index, target in enumerate(targets):
        target.write_bytes(f"old {index}\n".encode())

    _atomic_write_text_set({
        target: f"new {index}\n"
        for index, target in enumerate(targets)
    })

    assert [target.read_bytes() for target in targets] == [
        b"new 0\n", b"new 1\n", b"new 2\n"
    ]
    assert _transaction_siblings(tmp_path) == []


def test_multi_file_stage_collision_does_not_remove_sentinel(
    tmp_path, monkeypatch
):
    target = tmp_path / "KPOINTS_alter"
    token = "deadbeefcafebabe"
    colliding = tmp_path / f".alterseek-{token}-0.stage"
    colliding.write_bytes(b"previous transaction stage\n")
    monkeypatch.setattr(atomic_write, "_scratch_token", lambda: token)

    with pytest.raises(FileExistsError, match="scratch path already exists"):
        _atomic_write_text_set({target: "new k-points\n"})

    assert not target.exists()
    assert colliding.read_bytes() == b"previous transaction stage\n"
    assert _transaction_siblings(tmp_path) == [colliding.name]


def test_short_transaction_names_fit_legacy_windows_path_limit(
    tmp_path, monkeypatch
):
    parent = tmp_path / "long"
    parent.mkdir()
    target_length = 240
    basename_length = target_length - len(os.fspath(parent)) - 1
    target = parent / ("x" * basename_length)
    assert len(os.fspath(target)) == target_length
    target.write_bytes(b"old\n")

    real_open = os.open
    scratch_lengths = []

    def legacy_windows_open(path, flags, mode=0o777):
        path_length = len(os.fspath(path))
        scratch_lengths.append(path_length)
        if path_length >= 260:
            raise FileNotFoundError("synthetic legacy Windows MAX_PATH")
        return real_open(path, flags, mode)

    monkeypatch.setattr(os, "open", legacy_windows_open)
    _atomic_write_text_set({target: "new\n"})

    assert target.read_bytes() == b"new\n"
    assert max(scratch_lengths) < 260
    assert max(scratch_lengths) < len(os.fspath(target))
    assert _transaction_siblings(parent) == []


@pytest.mark.parametrize("failed_replace", [4, 5, 6])
def test_multi_file_promotion_failure_restores_existing_set(
    tmp_path, monkeypatch, failed_replace
):
    targets = [
        tmp_path / "case_magnetic_primitive.vasp",
        tmp_path / "case_magnetic_primitive_MAGMOM.txt",
        tmp_path / "KPOINTS_alter",
    ]
    originals = [b"old cell\n", b"old moments\n", b"old kpoints\n"]
    for target, contents in zip(targets, originals):
        target.write_bytes(contents)

    real_replace = os.replace
    calls = 0

    def fail_at_boundary(source, target):
        nonlocal calls
        calls += 1
        if calls == failed_replace:
            raise PermissionError("synthetic promotion failure")
        real_replace(source, target)

    monkeypatch.setattr(os, "replace", fail_at_boundary)
    with pytest.raises(PermissionError, match="synthetic promotion failure"):
        _atomic_write_text_set({
            target: f"new {index}\n"
            for index, target in enumerate(targets)
        })

    assert [target.read_bytes() for target in targets] == originals
    assert _transaction_siblings(tmp_path) == []


@pytest.mark.parametrize("failed_replace", [1, 2, 3])
def test_multi_file_promotion_failure_removes_new_set(
    tmp_path, monkeypatch, failed_replace
):
    targets = [
        tmp_path / "case_magnetic_primitive.vasp",
        tmp_path / "case_magnetic_primitive_MAGMOM.txt",
        tmp_path / "KPOINTS_alter",
    ]
    real_replace = os.replace
    calls = 0

    def fail_at_boundary(source, target):
        nonlocal calls
        calls += 1
        if calls == failed_replace:
            raise PermissionError("synthetic promotion failure")
        real_replace(source, target)

    monkeypatch.setattr(os, "replace", fail_at_boundary)
    with pytest.raises(PermissionError, match="synthetic promotion failure"):
        _atomic_write_text_set({
            target: f"new {index}\n"
            for index, target in enumerate(targets)
        })

    assert not any(target.exists() for target in targets)
    assert _transaction_siblings(tmp_path) == []
