"""Regression tests for KPOINTS label/path formatting in alterseek_path.py."""

import numpy as np
import pytest

from alterseek.kpoints import KPointsModifier, _fmt_coord


def test_signed_zero_is_normalized():
    """-0.0 and 0.0 are the same k-point and must render identically.

    Whether a zero component comes out signed depends on the basis conversion,
    so without this the same path emits different bytes depending on which cell
    it was built in.
    """
    assert _fmt_coord(-0.0) == _fmt_coord(0.0) == "0.0000000000"
    # A value that merely rounds to zero is also a zero component here.
    assert _fmt_coord(-1e-15) == "0.0000000000"


def test_nonzero_coordinates_keep_their_sign_and_precision():
    assert _fmt_coord(-0.5) == "-0.5000000000"
    assert _fmt_coord(0.5) == "0.5000000000"
    assert _fmt_coord(1.0 / 3.0) == "0.3333333333"
    # Small but genuinely nonzero at this precision: the sign must survive.
    assert _fmt_coord(-0.00000000006) == "-0.0000000001"


def test_gamma_label_is_vasp_safe():
    assert KPointsModifier._kpoints_label("Γ") == "GAMMA"
    assert KPointsModifier._kpoints_label("gamma") == "GAMMA"
    assert KPointsModifier._kpoints_label("GAMMA") == "GAMMA"
    # Ordinary labels pass through untouched, including primes and subscripts.
    assert KPointsModifier._kpoints_label("M_A'") == "M_A'"
    assert KPointsModifier._kpoints_label("X_1") == "X_1"


def test_display_label_matches_kpoints_label():
    # The console form and the VASP-file form must never drift apart.
    for label in ["Γ", "GAMMA", "K", "H_2", "M_A'", "k", "k'"]:
        assert (KPointsModifier._display_label(label)
                == KPointsModifier._kpoints_label(label))


def test_format_path_joins_continuous_and_breaks_discontinuous():
    segments = [("Γ", "X"), ("X", "M"), ("R", "Z")]
    assert KPointsModifier._format_path(segments) == "GAMMA-X-M | R-Z"


def test_dedupe_frac_positions_wraps_periodic_images():
    from alterseek.io import _dedupe_frac_positions
    unique = _dedupe_frac_positions([
        [0.25, 0.25, 0.25],
        [1.25, 0.25, -0.75],   # same site shifted by lattice vectors
        [0.75, 0.25, 0.25],
    ])
    assert len(unique) == 2


def test_reciprocal_basis_roundtrip():
    # k_input_frac = k_prim_frac @ B_prim @ inv(B_input) must be identity when
    # the bases coincide (the KPOINTS output-basis conversion contract).
    rng = np.random.default_rng(0)
    b = rng.normal(size=(3, 3)) + 3 * np.eye(3)
    k = np.array([0.386383, 0.363842, 0.351437])
    converted = k @ b @ np.linalg.inv(b)
    assert np.allclose(converted, k, atol=1e-12)


def _qe_waypoints(path):
    lines = path.read_text(encoding="utf-8").splitlines()[2:]
    return [(line.split("!", 1)[1].strip(), int(line.split("!", 1)[0].split()[3]))
            for line in lines]


def test_qe_writer_preserves_disconnected_chains(tmp_path):
    modifier = KPointsModifier()
    points = [
        [0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"], None,
        [0.0, 0.5, 0.0, "C"], [0.5, 0.5, 0.0, "D"],
    ]
    output = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_qe(points, str(output))
    assert _qe_waypoints(output) == [("A", 30), ("B", 1), ("C", 30), ("D", 1)]


def test_qe_writer_deduplicates_continuous_chain_boundary(tmp_path):
    modifier = KPointsModifier()
    points = [
        [0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"], None,
        [0.5, 0.0, 0.0, "B"], [0.5, 0.5, 0.0, "C"],
    ]
    output = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_qe(points, str(output))
    assert _qe_waypoints(output) == [("A", 30), ("B", 30), ("C", 1)]


def test_qe_writer_keeps_k_to_k_prime_gap_disconnected(tmp_path):
    modifier = KPointsModifier()
    points = [
        [0.0, 0.0, 0.0, "A"], [0.2, 0.2, 0.0, "k"],
        [0.2, 0.2, 0.0, "k'"], [0.5, 0.5, 0.0, "B"],
    ]
    output = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_qe(points, str(output))
    assert _qe_waypoints(output) == [("A", 30), ("k", 1), ("k'", 30), ("B", 1)]


def test_writers_record_the_operation_source_basis(tmp_path):
    """The recorded R is the raw Step-3 selection, so the header must name
    the operation-source basis (magnetic primitive cell in the magnetic
    route) instead of claiming a generic input cell."""
    modifier = KPointsModifier()
    modifier.header_lines = ["title", "20", "Line-Mode", "Reciprocal"]
    points = [[0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"]]
    matrix = np.eye(3)
    label = "magnetic primitive cell 'CASE_magnetic_primitive.mcif'"

    vasp_out = tmp_path / "KPOINTS_alter"
    assert modifier.write_kpoints_file(
        points, str(vasp_out), matrix, "Option 2", operation_basis_label=label
    )
    assert vasp_out.read_text(encoding="utf-8").splitlines()[0].startswith(
        f"Selected spin-flip operation (Option 2) in {label} "
        "real-space fractional basis:"
    )

    qe_out = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_qe(
        points, str(qe_out), matrix, "Option 2", operation_basis_label=label
    )
    assert qe_out.read_text(encoding="utf-8").splitlines()[-1].startswith(
        f"! Spin-flip operation (Option 2) in {label} "
        "real-space fractional basis:"
    )


def test_vasp_writer_does_not_damage_existing_file_on_conversion_failure(
        tmp_path, monkeypatch):
    modifier = KPointsModifier()
    output = tmp_path / "KPOINTS_alter"
    output.write_text("keep me\n", encoding="utf-8")
    calls = 0

    def fail_on_second_point(point):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("synthetic conversion failure")
        return point

    monkeypatch.setattr(modifier, "_kpoint_for_output_basis", fail_on_second_point)
    assert not modifier.write_kpoints_file(
        [[0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"]], str(output)
    )
    assert output.read_text(encoding="utf-8") == "keep me\n"


def test_qe_writer_does_not_damage_existing_file_on_conversion_failure(
        tmp_path, monkeypatch):
    modifier = KPointsModifier()
    output = tmp_path / "KPOINTS_alter_qe"
    output.write_text("keep me\n", encoding="utf-8")

    def fail_conversion(point):
        raise RuntimeError("synthetic conversion failure")

    monkeypatch.setattr(modifier, "_kpoint_for_output_basis", fail_conversion)
    assert not modifier.write_kpoints_file_qe(
        [[0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"]], str(output)
    )
    assert output.read_text(encoding="utf-8") == "keep me\n"


def test_custom_path_input_basis_is_converted_and_round_trips():
    modifier = KPointsModifier()
    modifier.kpoints_data = [[0.5, 0.25, 0.0, "X"], [0.0, 0.0, 0.0, "GAMMA"]]
    centroid_result = {
        "b_matrix": np.diag([2.0, 2.0, 3.0]),
        "b_matrix_input": np.diag([1.0, 2.0, 3.0]),
        "seekpath_rotation_matrix": np.eye(3),
    }
    modifier.convert_custom_path_from_input_basis(centroid_result)
    assert modifier.kpoints_data[0][:3] == pytest.approx([0.25, 0.25, 0.0])
    assert modifier._kpoint_for_output_basis(modifier.kpoints_data[0])[:3] == pytest.approx(
        [0.5, 0.25, 0.0]
    )


def test_custom_path_uses_submitted_basis_not_magnetic_analysis_basis():
    """Custom paths are defined in the submitted structure's basis.

    In the magnetic-cell route b_matrix_input belongs to the marker helper
    (the magnetic primitive cell), not to the submitted structure. Converting
    from it would silently reinterpret every custom coordinate; the conversion
    must read b_matrix_submitted whenever the two differ.
    """
    modifier = KPointsModifier()
    modifier.kpoints_data = [[0.5, 0.25, 0.0, "X"], [0.0, 0.0, 0.0, "GAMMA"]]
    centroid_result = {
        "b_matrix": np.diag([2.0, 2.0, 3.0]),
        # Magnetic analysis basis: equals b_internal, so using it by mistake
        # would leave the coordinates numerically unchanged.
        "b_matrix_input": np.diag([2.0, 2.0, 3.0]),
        "b_matrix_submitted": np.diag([1.0, 2.0, 3.0]),
        # Retained submitted calculation cell (supercell-style route).
        "b_matrix_output": np.diag([1.0, 2.0, 3.0]),
        "seekpath_rotation_matrix": np.eye(3),
    }
    modifier.convert_custom_path_from_input_basis(centroid_result)
    assert modifier.kpoints_data[0][:3] == pytest.approx([0.25, 0.25, 0.0])
    # Round trip back to the submitted (output) basis restores the file value.
    assert modifier._kpoint_for_output_basis(modifier.kpoints_data[0])[:3] == pytest.approx(
        [0.5, 0.25, 0.0]
    )
