"""The library entry point returns the same path the interactive driver writes.

general_k_path is the scripted route to the k-path; this pins its segments to
the stored KPOINTS reference used by the interactive golden test.
"""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("findspingroup")
pytest.importorskip("ase")

POSCAR = Path(__file__).parent / "references" / "case12_POSCAR"
REFERENCE = Path(__file__).parent / "references" / "case12_golden_kpoints.txt"


def _reference_points():
    points = []
    for line in REFERENCE.read_text(encoding="utf-8").splitlines()[4:]:
        if line.strip():
            fields = line.split()
            points.append((fields[3], [float(x) for x in fields[:3]]))
    return points


def _kpoints_label(label):
    return "GAMMA" if label == "Γ" else label


def test_general_k_path_matches_the_written_kpoints(tmp_path):
    from alterseek import general_k_path

    result = general_k_path(
        str(POSCAR),
        moments="5 -5",
        spin_axis="0 0 1",
        flip_option=1,
        output_dir=str(tmp_path / "out"),
    )

    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    expected = _reference_points()
    assert len(produced) == len(expected)
    for (got_label, got_coords), (want_label, want_coords) in zip(produced, expected):
        assert got_label == want_label
        assert got_coords == pytest.approx(want_coords, abs=1e-9)


def test_general_k_path_reports_the_operation_and_the_general_point(tmp_path):
    from alterseek import general_k_path

    result = general_k_path(
        str(POSCAR),
        moments="5 -5",
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == "tP1"
    assert result["magnetic_phase"] == "AFM(Altermagnet)"
    assert result["no_splitting_reason"] is None
    assert result["flip_option_count"] == 8
    assert result["spin_flip_operation"].shape == (3, 3)
    assert result["k"] == pytest.approx([1 / 6, 1 / 3, 0.25], abs=1e-9)
    assert result["k_prime"] is not None
    assert result["path"].startswith("Γ-X-k | ")


def test_general_k_path_rejects_an_out_of_range_operation(tmp_path):
    from alterseek import general_k_path

    with pytest.raises(ValueError, match="out of range"):
        general_k_path(
            str(POSCAR),
            moments="5 -5",
            flip_option=99,
            output_dir=str(tmp_path / "out"),
        )


@pytest.mark.parametrize(
    "flip_option",
    [True, False, 0, -1, 1.0, "1", np.bool_(True), np.float64(1.0), np.int64(0)],
)
def test_general_k_path_requires_a_positive_integer_operation(
    tmp_path, flip_option
):
    from alterseek import general_k_path

    with pytest.raises(ValueError, match="positive integer"):
        general_k_path(
            str(POSCAR),
            moments="5 -5",
            flip_option=flip_option,
            output_dir=str(tmp_path / "out"),
        )


def test_general_k_path_accepts_a_numpy_integer_operation(tmp_path):
    """A numpy integer selects the same operation as the plain integer."""
    from alterseek import general_k_path

    expected = general_k_path(
        str(POSCAR),
        moments="5 -5",
        flip_option=2,
        output_dir=str(tmp_path / "plain"),
    )
    result = general_k_path(
        str(POSCAR),
        moments="5 -5",
        flip_option=np.int64(2),
        output_dir=str(tmp_path / "numpy"),
    )

    assert result["path"] == expected["path"]
    assert result["k_prime"] == pytest.approx(expected["k_prime"], abs=1e-12)
    assert result["spin_flip_operation"] == pytest.approx(
        expected["spin_flip_operation"]
    )


@pytest.mark.parametrize("moments", ["5 -5", None])
def test_general_k_path_writes_nothing_for_a_missing_structure(tmp_path, moments):
    """A missing structure file is refused before output_dir is created."""
    from alterseek import SpinSymmetryError, general_k_path

    output_dir = tmp_path / "out"
    with pytest.raises(SpinSymmetryError, match="was not found"):
        general_k_path(
            str(tmp_path / "MISSING_POSCAR"),
            moments=moments,
            output_dir=str(output_dir),
        )
    assert not output_dir.exists()


def _cli_points(tmp_path, monkeypatch, structure, answers):
    """Drive the interactive session and read back the KPOINTS it writes."""
    import io
    import sys

    from alterseek.kpoints import KPathBuilder

    tmp_path.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys, "stdin",
        io.StringIO("\n".join([str(structure)] + answers + [""] * 6) + "\n"),
    )
    assert KPathBuilder().interactive_build() is True

    points = []
    for line in (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8").splitlines()[4:]:
        fields = line.split()
        if len(fields) >= 4:
            points.append((fields[3], [float(x) for x in fields[:3]]))
    return points


def test_general_k_path_uses_the_marker_cell_for_a_supercell(tmp_path, monkeypatch):
    """A supercell's zone comes from the marker cell, not the raw structure."""
    from alterseek import general_k_path

    structure = Path(__file__).parent / "references" / "SUPERCELL_211.vasp"
    result = general_k_path(
        str(structure), moments="1 -1 1 -1", output_dir=str(tmp_path / "out")
    )
    assert result["lattice"] == "oC1"

    expected = _cli_points(tmp_path / "cli", monkeypatch, structure, ["0 0 1", "1 -1 1 -1", ""])
    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    assert len(produced) == len(expected)
    for (got_label, got_coords), (want_label, want_coords) in zip(produced, expected):
        assert got_label == want_label
        assert got_coords == pytest.approx(want_coords, abs=1e-9)


def test_general_k_path_reads_moments_from_an_mcif(tmp_path):
    """La2NiO4 in the SSG setting, against a path file written by the workflow."""
    from alterseek import general_k_path

    structure = Path(__file__).parent / "references" / "La2NiO4.mcif"
    reference = Path(__file__).parent / "references" / "la2nio4_golden_kpoints.txt"

    result = general_k_path(str(structure), output_dir=str(tmp_path / "out"))
    assert result["lattice"] == "oP1"
    assert result["no_splitting_reason"] is None

    expected = []
    for line in reference.read_text(encoding="utf-8").splitlines()[4:]:
        fields = line.split()
        if len(fields) >= 4:
            expected.append((fields[3], [float(x) for x in fields[:3]]))

    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    assert len(produced) == len(expected)
    for (got_label, got_coords), (want_label, want_coords) in zip(produced, expected):
        assert got_label == want_label
        assert got_coords == pytest.approx(want_coords, abs=1e-9)


def test_general_k_path_without_moments_returns_the_ordinary_path(tmp_path):
    from alterseek import general_k_path

    result = general_k_path(str(POSCAR), output_dir=str(tmp_path / "out"))

    assert result["no_splitting_reason"] == "No magnetic moments entered."
    assert result["spin_flip_operation"] is None
    assert result["k_prime"] is None
    assert result["magnetic_phase"] is None
    assert len(result["segments"]) == 15


@pytest.mark.parametrize(
    "stem, moments, lattice, extra",
    [
        ("case21_tP1_4m", "5 -5 -5 5", "tP1", ["X_A", "R_A"]),
        ("case07_hP2_6m", "5 -5", "hP2", ["L_A", "M_A"]),
        ("case23_tI2_4m", "8 -8", "tI2", ["R_A", "S_0A", "S_A"]),
    ],
)
def test_general_k_path_carries_the_doubled_ibz_vertices(
    tmp_path, stem, moments, lattice, extra
):
    """4/m and 6/m append copied vertices in the order of the published paths."""
    from alterseek import general_k_path

    references = Path(__file__).parent / "references"
    result = general_k_path(
        str(references / f"{stem}_POSCAR"),
        moments=moments,
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == lattice
    assert result["brillouin_zone"]["extra_general_vertices"] == extra

    expected = []
    for line in (references / f"{stem}_golden_kpoints.txt").read_text(
        encoding="utf-8"
    ).splitlines()[4:]:
        fields = line.split()
        if len(fields) >= 4:
            expected.append((fields[3], [float(x) for x in fields[:3]]))
    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    assert [label for label, _ in produced] == [label for label, _ in expected]
    for (_, got), (_, want) in zip(produced, expected):
        assert got == pytest.approx(want, abs=1e-9)


DET5_SUPERCELL = Path(__file__).parent / "references" / "case12_det5_supercell.vasp"
DET5_MOMENTS = " ".join(["5 -5 0 0 0 0"] * 5)


def test_general_k_path_honours_the_submitted_cell_laue_gate(tmp_path):
    """A generic supercell can forbid altermagnetism even when its magnetic
    primitive cell does not."""
    from alterseek import general_k_path

    result = general_k_path(
        str(DET5_SUPERCELL), moments=DET5_MOMENTS, output_dir=str(tmp_path / "out")
    )

    assert result["lattice"] == "aP3"
    assert result["no_splitting_reason"] == "Laue group -1: no altermagnetism."
    assert result["k_prime"] is None
    assert result["spin_flip_operation"] is None


DET3_SUPERCELL = Path(__file__).parent / "references" / "case12_det3_supercell.vasp"
DET3_MOMENTS = " ".join(["5 -5 0 0 0 0"] * 3)


def test_general_k_path_refuses_when_no_operation_was_found(tmp_path):
    """No current spin-flip operation and no reason to skip: refuse, as the
    session does, rather than return a path."""
    from alterseek import general_k_path

    with pytest.raises(ValueError, match="no detected spin-flip point operation"):
        general_k_path(
            str(DET3_SUPERCELL),
            moments=DET3_MOMENTS,
            output_dir=str(tmp_path / "clean"),
        )


def test_general_k_path_ignores_a_stale_operation_file(tmp_path):
    """An operation file this run did not write is never read.

    The structure has no spin-flip operation of its own and no reason to skip
    the search, so a stale file is the only thing that could supply one.
    """
    from alterseek import general_k_path

    source = tmp_path / "source"
    general_k_path(str(POSCAR), moments="5 -5", output_dir=str(source))
    stale = (source / "spin_flip_operations.txt").read_text(encoding="utf-8")

    seeded = tmp_path / "seeded"
    seeded.mkdir()
    (seeded / "spin_flip_operations.txt").write_text(stale, encoding="utf-8")

    with pytest.raises(ValueError, match="no detected spin-flip point operation"):
        general_k_path(
            str(DET3_SUPERCELL), moments=DET3_MOMENTS, output_dir=str(seeded)
        )


def test_general_k_path_treats_blank_moments_as_none(tmp_path):
    from alterseek import general_k_path

    result = general_k_path(str(POSCAR), moments="   ", output_dir=str(tmp_path / "out"))
    assert result["no_splitting_reason"] == "No magnetic moments entered."


@pytest.mark.parametrize(
    "stem, moments",
    [
        ("case2d06_square_4m", "1 -1 -1 1 6*0"),
        ("case2d01_hex_6mmm", "1 -1 6*0"),
    ],
)
def test_general_k_path_slab_mode(tmp_path, stem, moments):
    """2D mode against the stored path files for two slab cases."""
    from alterseek import general_k_path

    references = Path(__file__).parent / "references"
    result = general_k_path(
        str(references / f"{stem}_POSCAR"),
        moments=moments,
        mode_2d=True,
        vacuum_axis="c",
        output_dir=str(tmp_path / "out"),
    )

    expected = []
    for line in (references / f"{stem}_golden_kpoints.txt").read_text(
        encoding="utf-8"
    ).splitlines()[4:]:
        fields = line.split()
        if len(fields) >= 4:
            expected.append((fields[3], [float(x) for x in fields[:3]]))

    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    assert len(produced) == len(expected)
    for (got_label, got_coords), (want_label, want_coords) in zip(produced, expected):
        assert got_label == want_label
        assert got_coords == pytest.approx(want_coords, abs=1e-9)


def test_general_k_path_slab_mode_on_a_nearly_square_cell(tmp_path):
    """CrS2 (2D Case 05) in a sqrt2 x sqrt2 cell, square only within the
    symmetry tolerance."""
    from alterseek import general_k_path

    result = general_k_path(
        str(Path(__file__).parent / "references" / "case2d05_sqrt2_POSCAR"),
        moments="2.55 2.55 -2.55 -2.55 8*0",
        mode_2d=True,
        vacuum_axis="c",
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == "square"
    assert result["no_splitting_reason"] == (
        "PT symmetry detected, not altermagnet."
    )
    assert result["k"] == pytest.approx([0.25, 0.25, 0.0])
    produced = []
    for segment in result["segments"][:4]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))
    expected = [
        ("GAMMA", [0.0, 0.0, 0.0]), ("X", [0.5, 0.0, 0.0]),
        ("X", [0.5, 0.0, 0.0]), ("S", [0.5, 0.5, 0.0]),
        ("S", [0.5, 0.5, 0.0]), ("Y", [0.0, 0.5, 0.0]),
        ("Y", [0.0, 0.5, 0.0]), ("GAMMA", [0.0, 0.0, 0.0]),
    ]
    assert [label for label, _ in produced] == [label for label, _ in expected]
    for (_, got), (_, want) in zip(produced, expected):
        assert got == pytest.approx(want, abs=1e-9)


@pytest.mark.parametrize(
    "degeneracy_forcing, valid_in_plane, expected_reason",
    [
        (
            True,
            False,
            "C_2z T / U m_z symmetry detected, not a 2D altermagnet.",
        ),
        (
            False,
            False,
            "No in-plane spin-flip point operation available: not a 2D "
            "altermagnet.",
        ),
    ],
)
def test_general_k_path_distinguishes_2d_fallback_reasons(
    tmp_path, monkeypatch, degeneracy_forcing, valid_in_plane, expected_reason
):
    """The two 2D fallbacks return distinct reasons.

    Both predicates are stubbed to constants, because no shipped 2D case
    reaches either branch. This pins the returned reasons, not the symmetry
    tests behind them.
    """
    from alterseek import general_k_path
    from alterseek.kpoints import KPathBuilder

    monkeypatch.setattr(
        KPathBuilder,
        "_forces_2d_degeneracy",
        lambda self, operation, zone: degeneracy_forcing,
    )
    monkeypatch.setattr(
        KPathBuilder,
        "_is_valid_2d_operation",
        lambda self, operation, zone: valid_in_plane,
    )

    references = Path(__file__).parent / "references"
    result = general_k_path(
        str(references / "case2d01_hex_6mmm_POSCAR"),
        moments="1 -1 6*0",
        mode_2d=True,
        vacuum_axis="c",
        output_dir=str(tmp_path / "out"),
    )

    assert result["no_splitting_reason"] == expected_reason
    assert result["spin_flip_operation"] is None
    assert result["k_prime"] is None


@pytest.mark.parametrize(
    "name", ["spin_flip_operations.txt", "spin_preserve_operations.txt"]
)
def test_general_k_path_writes_the_operation_files_the_workflow_writes(
    tmp_path, monkeypatch, name
):
    """cF2: the submitted and SeeK-path cells differ, so both bases are listed."""
    import shutil

    from alterseek import general_k_path

    references = Path(__file__).parent / "references"
    monkeypatch.chdir(tmp_path)
    shutil.copy(references / "case02_cF2_POSCAR", tmp_path / "POSCAR")

    general_k_path("POSCAR", moments="5 -5 12*0", output_dir="out")

    stem = name.removesuffix(".txt")
    assert (tmp_path / "out" / name).read_text(encoding="utf-8") == (
        references / f"case02_cF2_golden_{stem}.txt"
    ).read_text(encoding="utf-8")
