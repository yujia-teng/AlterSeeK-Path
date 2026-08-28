"""Regression tests for KPOINTS label/path formatting in alterseek_path.py."""

import numpy as np
import pytest

from alterseek import kpoints as kpoints_module
from alterseek.kpoints import KPointsModifier, _fmt_coord


def test_manual_spin_flip_matrix_is_rejected(monkeypatch, capsys):
    modifier = KPointsModifier()
    detected_ops = [np.eye(3), np.diag([-1.0, -1.0, 1.0])]
    answers = iter(["manual", "2"])
    monkeypatch.setattr("builtins.input", lambda: next(answers))

    selected, label = modifier._select_spin_flip_operation(
        detected_ops, centroid_result=None
    )

    assert np.array_equal(selected, detected_ops[1])
    assert label == "Option 2"
    output = capsys.readouterr().out
    assert "Please choose 1-2 or 'list'." in output
    assert "custom transformation matrix" not in output
    assert "or 'manual'" not in output


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
    assert KPointsModifier._kpoints_label("\u0393") == "GAMMA"
    assert KPointsModifier._kpoints_label("gamma") == "GAMMA"
    assert KPointsModifier._kpoints_label("GAMMA") == "GAMMA"
    assert KPointsModifier._kpoints_label("M_A'") == "M_A'"
    assert KPointsModifier._kpoints_label("X_1") == "X_1"


def test_display_label_matches_kpoints_label():
    # The console form and the VASP-file form must never drift apart.
    for label in ["\u0393", "GAMMA", "K", "H_2", "M_A'", "k", "k'"]:
        assert (KPointsModifier._display_label(label)
                == KPointsModifier._kpoints_label(label))


def test_combined_gamma_label_is_vasp_safe():
    assert KPointsModifier._kpoints_label("\u0393/H_2") == "GAMMA/H_2"


def test_format_path_joins_continuous_and_breaks_discontinuous():
    segments = [("\u0393", "X"), ("X", "M"), ("R", "Z")]
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


def test_writers_combine_consecutive_coincident_path_labels(tmp_path):
    modifier = KPointsModifier()
    modifier.header_lines = ["Path", "30", "Line-Mode", "Reciprocal"]
    points = [
        [0.0, 0.0, 0.0, "A"], [0.0, 0.0, 0.0, "H"],
        [0.0, 0.0, 0.0, "H"], [0.5, 0.0, 0.0, "X"],
    ]

    vasp_output = tmp_path / "KPOINTS_alter"
    qe_output = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_vasp(points, str(vasp_output))
    assert modifier.write_kpoints_file_qe(points, str(qe_output))

    vasp_labels = [
        line.split()[-1]
        for line in vasp_output.read_text(encoding="utf-8").splitlines()[4:]
        if line.split()
    ]
    assert vasp_labels == ["A/H", "X"]
    assert _qe_waypoints(qe_output) == [("A/H", 30), ("X", 1)]


def test_ti1_boundary_propagates_alias_from_removed_zero_length_part(tmp_path):
    modifier = KPointsModifier()
    modifier.header_lines = ["Path", "30", "Line-Mode", "Reciprocal"]
    modifier.kpoints_data = [
        [0.0, 0.0, 0.0, "GAMMA"], [0.0, 0.0, 0.5, "X"],
        [0.0, 0.0, 0.5, "X"], [-0.5, 0.5, 0.5, "M"],
        [-0.5, 0.5, 0.5, "M"], [0.0, 0.0, 0.0, "GAMMA"],
        [0.0, 0.0, 0.0, "GAMMA"], [0.5, 0.5, -0.5, "Z"],
        [-0.5, 0.5, 0.5, "Z_0"], [-0.5, 0.5, 0.5, "M"],
        [0.0, 0.0, 0.5, "X"], [0.25, 0.25, 0.25, "P"],
        [0.25, 0.25, 0.25, "P"], [0.0, 0.5, 0.0, "N"],
        [0.0, 0.5, 0.0, "N"], [0.0, 0.0, 0.0, "GAMMA"],
    ]

    path = modifier.insert_general_kpoints(
        [0.2, 0.3, 0.4],
        np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    labels = [point[3] for point in path if point is not None]
    assert "M/Z_0" in labels
    assert "M'/Z_0'" in labels
    assert "M" not in labels
    assert "Z_0" not in labels

    vasp_output = tmp_path / "KPOINTS_alter"
    qe_output = tmp_path / "KPOINTS_alter_qe"
    assert modifier.write_kpoints_file_vasp(path, str(vasp_output))
    assert modifier.write_kpoints_file_qe(path, str(qe_output))

    vasp_text = vasp_output.read_text(encoding="utf-8")
    qe_labels = [label for label, _count in _qe_waypoints(qe_output)]
    assert "M/Z_0" in vasp_text
    assert "M'/Z_0'" in vasp_text
    assert "M/Z_0" in qe_labels
    assert "M'/Z_0'" in qe_labels


def test_butterfly_path_retains_primed_coincident_aliases():
    modifier = KPointsModifier()
    modifier.kpoints_data = [
        [0.5, 0.0, 0.0, "X"], [0.0, 0.0, 0.0, "A"],
        [0.0, 0.0, 0.0, "A"], [0.0, 0.0, 0.0, "H"],
    ]

    path = modifier.insert_general_kpoints(
        [0.2, 0.2, 0.2], np.diag([-1.0, 1.0, 1.0])
    )
    labels = [point[3] for point in path if point is not None]

    assert "A/H" in labels
    assert "A'/H'" in labels


def test_2d_4m_butterfly_uses_open_path_and_gamma_xa_connections():
    modifier = KPointsModifier()
    modifier.kpoints_basis_matrix = np.eye(3)
    modifier.output_basis_matrix = np.eye(3)
    butterfly = [
        [0.0, 0.0, 0.0, "GAMMA"], [0.0, 0.5, 0.0, "X"],
        [0.0, 0.5, 0.0, "X"], [0.5, 0.5, 0.0, "M"],
    ]
    connections = [
        [0.0, 0.0, 0.0, "GAMMA"],
        [0.5, 0.0, 0.0, "X_A"],
    ]
    c4 = np.array([
        [0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])

    path = modifier.insert_general_kpoints(
        [0.25, 0.25, 0.0],
        c4,
        connections,
        path_points=butterfly,
        report=False,
    )
    segments = [
        (start[3], end[3], break_before)
        for start, end, break_before, _index, _raw_label
        in modifier._prepare_output_segments(path)
    ]

    assert segments == [
        ("GAMMA", "X", False),
        ("X", "k", False),
        ("k'", "X'", True),
        ("X'", "M'", False),
        ("M'", "k'", False),
        ("k", "M", True),
        ("GAMMA", "k", True),
        ("k'", "GAMMA", True),
        ("X_A", "k", True),
        ("k'", "X_A'", True),
    ]
    assert ("M", "GAMMA", False) not in segments


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
    assert modifier.write_kpoints_file_vasp(
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
    modifier.header_lines = ["title", "20", "Line-Mode", "Reciprocal"]
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
    with pytest.raises(RuntimeError, match="synthetic conversion failure"):
        modifier.write_kpoints_file_vasp(
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
    with pytest.raises(RuntimeError, match="synthetic conversion failure"):
        modifier.write_kpoints_file_qe(
            [[0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"]], str(output)
        )
    assert output.read_text(encoding="utf-8") == "keep me\n"


@pytest.mark.parametrize(
    ("writer_name", "filename"),
    [
        ("write_kpoints_file_vasp", "KPOINTS_alter"),
        ("write_kpoints_file_qe", "KPOINTS_alter_qe"),
    ],
)
def test_writers_report_persistence_failure_without_damaging_existing_file(
    tmp_path, monkeypatch, capsys, writer_name, filename
):
    modifier = KPointsModifier()
    modifier.header_lines = ["title", "20", "Line-Mode", "Reciprocal"]
    output = tmp_path / filename
    output.write_text("keep me\n", encoding="utf-8")

    def fail_persistence(path, text):
        raise PermissionError("synthetic permission failure")

    monkeypatch.setattr(kpoints_module, "_atomic_write_text", fail_persistence)
    writer = getattr(modifier, writer_name)
    assert not writer(
        [[0.0, 0.0, 0.0, "A"], [0.5, 0.0, 0.0, "B"]], str(output)
    )
    assert "synthetic permission failure" in capsys.readouterr().out
    assert output.read_text(encoding="utf-8") == "keep me\n"
