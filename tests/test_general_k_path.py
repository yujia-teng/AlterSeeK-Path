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
    assert result["lattice"] == "oP1"

    expected = _cli_points(tmp_path / "cli", monkeypatch, structure, ["0 0 1", "1 -1 1 -1", ""])
    produced = []
    for segment in result["segments"]:
        produced.append((_kpoints_label(segment["start_label"]), segment["start"]))
        produced.append((_kpoints_label(segment["end_label"]), segment["end"]))

    assert len(produced) == len(expected)
    for (got_label, got_coords), (want_label, want_coords) in zip(produced, expected):
        assert got_label == want_label
        assert got_coords == pytest.approx(want_coords, abs=1e-9)


def test_general_k_path_gives_a_spin_flipping_translation_cell_its_own_zone(tmp_path):
    """GdAuGe 2x1x1 against the stored path file of its oP1 zone."""
    from alterseek import general_k_path

    references = Path(__file__).parent / "references"
    result = general_k_path(
        str(references / "SUPERCELL_211.vasp"),
        moments="1 -1 1 -1",
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == "oP1"
    assert result["no_splitting_reason"] == "Ut symmetry detected, not altermagnet."
    assert result["k"] == pytest.approx([0.0, 0.25, -0.25], abs=1e-9)

    expected = []
    for line in (references / "ssg_supercell211_golden_kpoints.txt").read_text(
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


@pytest.mark.parametrize(
    "mode_2d, lattice, corners",
    [
        (True, "rectangular", {"X": [0.5, 0.0, 0.0], "S": [0.5, 0.5, 0.0], "Y": [0.0, 0.5, 0.0]}),
        (False, "oP1", {"X": [0.0, -0.5, 0.0], "Y": [0.0, 0.0, -0.5], "Z": [0.5, 0.0, 0.0]}),
    ],
)
def test_general_k_path_handles_a_spin_flipping_translation_in_primitive_g0(
    tmp_path, mode_2d, lattice, corners
):
    """V2Se2O 2x1 with V moments 1 1 -1 -1: G0 Pmm2 repeats on the 1x1 cell."""
    from alterseek import general_k_path

    result = general_k_path(
        str(Path(__file__).parent / "references" / "case2d04_square_2x1_POSCAR"),
        moments="1 1 -1 -1",
        mode_2d=mode_2d,
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == lattice
    assert result["no_splitting_reason"] == "Ut symmetry detected, not altermagnet."
    assert result["magnetic_space_group_without_soc"] == "P_cmm2 (BNS 25.61), Type IV"
    points = {}
    for segment in result["segments"]:
        points[segment["start_label"]] = segment["start"]
        points[segment["end_label"]] = segment["end"]
    for label, coords in corners.items():
        assert points[label] == pytest.approx(coords, abs=1e-9)


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
DET3_SUPERCELL = Path(__file__).parent / "references" / "case12_det3_supercell.vasp"
DET3_MOMENTS = " ".join(["5 -5 0 0 0 0"] * 3)
SQUARE_DWAVE_2X1 = Path(__file__).parent / "references" / "square_dwave_2x1_POSCAR"


@pytest.mark.parametrize(
    "structure, moments",
    [(DET3_SUPERCELL, DET3_MOMENTS), (DET5_SUPERCELL, DET5_MOMENTS)],
)
def test_general_k_path_stops_when_a_skewed_supercell_keeps_no_spin_flip_operation(
    tmp_path, structure, moments
):
    """Case 12 in two skewed supercells keeps none of its spin-flip operations.

    The five-fold cell keeps only the identity and inversion, so its zone alone
    reports Laue group -1, although the structure is an altermagnet.
    """
    from alterseek import general_k_path

    with pytest.raises(
        ValueError,
        match=(
            "The structure is an altermagnet, but none of its 8 spin-flip "
            "point operations maps the input cell onto itself"
        ),
    ):
        general_k_path(
            str(structure), moments=moments, output_dir=str(tmp_path / "out")
        )


@pytest.mark.parametrize(
    "mode_2d, repeat", [(True, "2 x 2 in the plane"), (False, "2 x 2 x 2")]
)
def test_general_k_path_stops_when_no_spin_flip_operation_fits_a_2x1_cell(
    tmp_path, mode_2d, repeat
):
    from alterseek import general_k_path

    with pytest.raises(ValueError) as error:
        general_k_path(
            str(SQUARE_DWAVE_2X1),
            moments="1 1 -1 -1 1 1 -1 -1",
            mode_2d=mode_2d,
            output_dir=str(tmp_path / "out"),
        )

    message = str(error.value)
    assert message.startswith(
        "The structure is an altermagnet, but none of its 8 spin-flip point "
        "operations maps the input cell onto itself: C2 [1 -1 0], C2 [1 1 0], "
        "C4+ [0 0 1], C4- [0 0 1], S4+ [0 0 1], S4- [0 0 1], "
        "mirror m (1 -1 0), mirror m (1 1 0). "
    )
    mcif = tmp_path / "out" / "square_dwave_2x1_POSCAR_magnetic_primitive.mcif"
    assert (
        f"Use the magnetic primitive cell (4 atoms, written to {mcif})"
        in message
    )
    assert message.endswith(f"such as {repeat}.")
    assert mcif.exists()


def test_general_k_path_ignores_a_stale_operation_file(tmp_path, monkeypatch):
    """An operation file this run did not write is never read.

    The spin-symmetry result is made inconsistent on purpose: an altermagnet
    whose operations all fit the cell, but for which no spin-flip operation
    file was written. A file left by an earlier run is then the only thing
    that could supply an operation.
    """
    from alterseek import general_k_path
    from alterseek import kpoints as kpoints_module

    source = tmp_path / "source"
    general_k_path(str(POSCAR), moments="5 -5", output_dir=str(source))
    seeded = tmp_path / "seeded"
    seeded.mkdir()
    (seeded / "spin_flip_operations.txt").write_text(
        (source / "spin_flip_operations.txt").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    real_run = kpoints_module.find_sf_run

    def run_without_flip_file(structure_file, moments, **kwargs):
        kwargs["output_dir"] = str(tmp_path / "elsewhere")
        result = real_run(structure_file, moments, **kwargs)
        return {**result, "spin_flip_operations": 0}

    loaded = []
    real_load = kpoints_module.KPathBuilder.load_flip_operations

    def spy_load(self, *args, **kwargs):
        loaded.append(args)
        return real_load(self, *args, **kwargs)

    monkeypatch.setattr(kpoints_module, "find_sf_run", run_without_flip_file)
    monkeypatch.setattr(
        kpoints_module.KPathBuilder, "load_flip_operations", spy_load
    )

    with pytest.raises(ValueError, match="no detected spin-flip point operation"):
        general_k_path(str(POSCAR), moments="5 -5", output_dir=str(seeded))
    assert loaded == []


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


@pytest.mark.parametrize("thickness", [0.0, 0.001, 0.002, 0.003])
def test_general_k_path_slab_mode_on_a_slab_about_as_thin_as_the_tolerance(
    tmp_path, thickness
):
    """A P1 slab in a 2x1 cell, flattened to ``thickness`` angstrom (symprec 1e-3)."""
    import spglib
    from ase.io import read, write

    from alterseek import general_k_path

    atoms = read(Path(__file__).parent / "references" / "thin_p1_2x1_POSCAR")
    positions = atoms.get_scaled_positions()
    heights = positions[:, 2]
    spread = heights.max() - heights.min()
    positions[:, 2] = 0.5 + (heights - heights.mean()) / spread * thickness / 20.0
    atoms.set_scaled_positions(positions)
    structure = tmp_path / "POSCAR"
    write(structure, atoms, format="vasp", direct=True, sort=False)

    result = general_k_path(
        str(structure), mode_2d=True, vacuum_axis="c",
        output_dir=str(tmp_path / "out"),
    )

    assert result["lattice"] == "rectangular"
    real = spglib.get_symmetry_dataset(
        (atoms.cell[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
        symprec=1e-3,
    )
    assert result["brillouin_zone"]["layer_point_group"] == real.pointgroup


def test_general_k_path_slab_mode_on_a_sixty_degree_hexagonal_cell(tmp_path):
    """FeBr3 (2D Case 01) in its 60-degree cell against the 120-degree cell."""
    from alterseek import general_k_path

    references = Path(__file__).parent / "references"
    runs = {}
    for name in ("case2d01_hex_6mmm", "case2d01_hex_6mmm_60deg"):
        structure = references / f"{name}_POSCAR"
        result = general_k_path(
            str(structure),
            moments="1 -1 6*0",
            mode_2d=True,
            vacuum_axis="c",
            output_dir=str(tmp_path / name),
        )
        lattice = np.loadtxt(structure, skiprows=2, max_rows=3)
        reciprocal = 2.0 * np.pi * np.linalg.inv(lattice).T
        points = [(result["k"], "k"), (result["k_prime"], "k'")]
        for segment in result["segments"]:
            points.append((segment["start"], segment["start_label"]))
            points.append((segment["end"], segment["end_label"]))
        runs[name] = [
            (label, np.linalg.norm(np.asarray(point) @ reciprocal))
            for point, label in points
        ]

    sixty = runs["case2d01_hex_6mmm_60deg"]
    assert [label for label, _ in sixty] == [
        label for label, _ in runs["case2d01_hex_6mmm"]
    ]
    assert dict(sixty)["K"] == pytest.approx(
        4.0 * np.pi / (3.0 * np.linalg.norm(lattice[0])), abs=1e-9
    )
    for (_, got), (_, want) in zip(sixty, runs["case2d01_hex_6mmm"]):
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


def _left_handed_copy(tmp_path):
    """Case 12 with a and b swapped: the same crystal, left-handed vectors."""
    from ase.io import read, write

    atoms = read(POSCAR)
    positions = atoms.get_scaled_positions()[:, [1, 0, 2]]
    atoms.set_cell(atoms.cell[:][[1, 0, 2]])
    atoms.set_scaled_positions(positions)
    path = tmp_path / "POSCAR"
    write(path, atoms, format="vasp", direct=True)
    assert np.linalg.det(read(path).cell[:]) < 0.0
    return path


LEFT_HANDED = (
    "The lattice vectors in POSCAR are left-handed (a . (b x c) < 0). "
    "Swap two of them, or reverse one, and run again."
)


@pytest.mark.parametrize(
    "function, moments",
    [
        ("spin_symmetry", "5 -5"),
        ("brillouin_zone", "5 -5"),
        ("brillouin_zone", None),
        ("general_k_path", "5 -5"),
        ("general_k_path", None),
    ],
)
def test_left_handed_cell_stops_the_python_functions(tmp_path, function, moments):
    import alterseek
    from alterseek import SpinSymmetryError

    kwargs = {"output_dir": str(tmp_path / "out")}
    if function == "brillouin_zone":
        kwargs["show_plot"] = False
    with pytest.raises(SpinSymmetryError) as error:
        getattr(alterseek, function)(str(_left_handed_copy(tmp_path)), moments, **kwargs)
    assert str(error.value) == LEFT_HANDED


def test_left_handed_cell_stops_the_workflow(tmp_path, monkeypatch, capsys):
    import io
    import sys

    from alterseek import run_workflow

    _left_handed_copy(tmp_path)
    (tmp_path / "alterseek_input.toml").write_text(
        'structure = "POSCAR"\nspin_axis = "0 0 1"\nmoments = "5 -5"\n'
        'flip_option = 1\noutput_code = "vasp"\n',
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(""))

    assert run_workflow() is False
    assert f"[Error] {LEFT_HANDED} Aborting." in capsys.readouterr().out
    assert not (tmp_path / "KPOINTS_alter").exists()
