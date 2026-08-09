"""Aligning the magnetic standard cell onto the SeeK-path standard cell.

The optional `*_seekpath_standard.mcif` puts the input moments into the
SeeK-path standardized cell.  spglib's magnetic standardization and SeeK-path's
structural one describe the same crystal but can pick different axes and
different site representatives, so the writer has to find the basis change
between them before it may attach any moment to any atom.
"""

import numpy as np
import pytest

from alterseek.io import _basis_change_candidates, _write_seekpath_standard_mcif

# Real triclinic MnO2.  Both standardizations return this same lattice, but
# they list different oxygen representatives, and the transform relating them
# has a row that is a combination of all three axes rather than a relabelling.
MNO2_LATTICE = np.array([
    [2.9197337429464350, 0.0000000000000000, 0.0000000000000000],
    [-1.4594396286354010, 2.5288093565268150, 0.0000000000000000],
    [-1.3984483268394059, -0.8784852794219501, 4.7809860000000004],
])
MNO2_MAGNETIC_POSITIONS = np.array([
    [0.000000, 0.000000, 0.000000],
    [0.734748, 0.462431, 0.798869],
    [0.265252, 0.537569, 0.201131],
])
MNO2_SEEKPATH_POSITIONS = np.array([
    [0.000000, 0.000000, 0.000000],
    [0.663562, 0.935879, 0.201131],
    [0.336438, 0.064121, 0.798869],
])
MNO2_ALIGNING_BASIS_CHANGE = np.array([[0, 1, 0], [1, 0, 0], [-1, -1, -1]])


def _is_signed_permutation(matrix):
    matrix = np.asarray(matrix)
    return (
        all(np.count_nonzero(row) == 1 for row in matrix)
        and all(np.count_nonzero(column) == 1 for column in matrix.T)
        and set(np.abs(matrix[matrix != 0]).tolist()) <= {1}
    )


def _write_poscar(path, lattice, positions, elements, counts):
    lines = ["alignment fixture", "1.0"]
    lines += ["   %.16f   %.16f   %.16f" % tuple(row) for row in lattice]
    lines += [" ".join(elements), " ".join(str(count) for count in counts)]
    lines.append("Direct")
    lines += ["   %.16f   %.16f   %.16f" % tuple(row) for row in positions]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_signed_permutations_are_offered_first_in_the_historical_order():
    """Callers break ties by order, so already-aligned cases must not move.

    The 48 relabellings are emitted unconditionally and in exactly the order
    the previous nested loops produced, including the ones the caller then
    rejects on determinant or orthogonality. Anything else could silently
    re-point a case that already aligned.
    """
    import itertools

    expected = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            entry = np.zeros((3, 3))
            for row, column in enumerate(permutation):
                entry[row, column] = signs[row]
            expected.append(entry)

    candidates = _basis_change_candidates(MNO2_LATTICE, MNO2_LATTICE)

    assert len(candidates) > 48
    for offered, historical in zip(candidates[:48], expected):
        assert np.array_equal(offered, historical)
    assert all(_is_signed_permutation(entry) for entry in candidates[:48])


def test_candidates_include_the_combination_row_mno2_needs():
    """The transform MnO2 aligns under is not any relabelling of the axes."""
    candidates = _basis_change_candidates(MNO2_LATTICE, MNO2_LATTICE)
    matches = [
        entry for entry in candidates
        if np.array_equal(entry.astype(int), MNO2_ALIGNING_BASIS_CHANGE)
    ]

    assert matches, "the aligning basis change must be offered"
    assert not _is_signed_permutation(MNO2_ALIGNING_BASIS_CHANGE)


def test_generated_candidates_reproduce_the_target_lengths_and_angles():
    """Beyond the historical 48, every candidate must satisfy the metric.

    The relabellings stay unconditional so the caller keeps filtering them as
    before; the new ones are generated *from* the metric and would otherwise
    be free to propose a transform that reshapes the cell.
    """
    candidates = _basis_change_candidates(MNO2_LATTICE, MNO2_LATTICE)
    target_gram = MNO2_LATTICE @ MNO2_LATTICE.T

    for entry in candidates[48:]:
        reindexed = entry @ MNO2_LATTICE
        assert np.allclose(reindexed @ reindexed.T, target_gram, atol=1e-4)
        assert abs(abs(np.linalg.det(entry)) - 1.0) < 1e-9


def test_mno2_moment_reaches_the_right_atom_in_the_standard_cell(tmp_path):
    """End to end: the case that used to refuse with a 1.16 A mismatch."""
    target = tmp_path / "MnO2_seekpath_standard.vasp"
    _write_poscar(
        target, MNO2_LATTICE, MNO2_SEEKPATH_POSITIONS, ["Mn", "O"], [1, 2]
    )
    output = tmp_path / "MnO2_seekpath_standard.mcif"

    _write_seekpath_standard_mcif(
        str(target),
        str(output),
        "MnO2_seekpath_standard",
        MNO2_LATTICE,
        MNO2_MAGNETIC_POSITIONS,
        ["Mn", "O", "O"],
        np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        symprec=1e-3,
    )

    text = output.read_text(encoding="utf-8")
    moment_rows = [
        line.split()
        for line in text.splitlines()
        if line.startswith(("Mn", "O")) and len(line.split()) == 4
    ]
    moments = {row[0]: np.array([float(v) for v in row[1:]]) for row in moment_rows}

    # Only the manganese carries a moment, and it kept unit magnitude: the
    # crystal-axis components are written against unit-length axes.
    assert set(moments) == {"Mn1", "O1", "O2"}
    assert np.allclose(moments["O1"], 0.0)
    assert np.allclose(moments["O2"], 0.0)
    unit_axes = MNO2_LATTICE / np.linalg.norm(MNO2_LATTICE, axis=1)[:, None]
    assert np.linalg.norm(moments["Mn1"] @ unit_axes) == pytest.approx(1.0, abs=1e-6)


def test_unalignable_pair_still_refuses_to_write(tmp_path):
    """The fail-closed gate must survive the wider search.

    A uniform shift would be a legitimate origin choice that the matcher is
    meant to absorb, so this moves a single oxygen instead: that is a
    genuinely different arrangement and must never receive moments.
    """
    target = tmp_path / "bad_seekpath_standard.vasp"
    displaced = MNO2_SEEKPATH_POSITIONS.copy()
    displaced[2] += np.array([0.0, 0.0, 0.17])
    _write_poscar(target, MNO2_LATTICE, displaced, ["Mn", "O"], [1, 2])
    output = tmp_path / "bad_seekpath_standard.mcif"

    with pytest.raises(RuntimeError):
        _write_seekpath_standard_mcif(
            str(target),
            str(output),
            "bad",
            MNO2_LATTICE,
            MNO2_MAGNETIC_POSITIONS,
            ["Mn", "O", "O"],
            np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            symprec=1e-3,
        )
    assert not output.exists()
