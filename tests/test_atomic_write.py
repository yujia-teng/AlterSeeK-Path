"""Regressions for the atomic text writers.

Every generated file now goes through these. Before, ten writers built their
output incrementally with a plain ``open(..., "w")``, so a failure partway
through left a truncated POSCAR / MCIF / operation list on disk that still
looked loadable. The buffered form commits only on a clean exit.
"""
import os

import pytest

from alterseek.atomic_write import (
    _atomic_open_text,
    _atomic_write_text,
)


def _temp_siblings(directory):
    """Leftover temp files use a leading dot and a .tmp suffix."""
    return [n for n in os.listdir(directory)
            if n.startswith(".") and n.endswith(".tmp")]


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
