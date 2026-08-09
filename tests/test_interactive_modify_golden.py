"""Golden regression test for KPointsModifier.interactive_modify.

The interactive Step 0-5 driver has no other coverage (the rest of the suite
tests the engine methods).  This drives the full altermagnetic flow on case 12
(MnF2, tetragonal P4/mmm, d-wave) with canned keyboard input and asserts the
produced KPOINTS matches a stored reference byte-for-byte.  It guards the phase-5
extraction of interactive_modify into helper methods (and any future edits).
"""
import io
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

# The full interactive flow runs Step 0 (FindSpinGroup) and reads the structure
# with ASE; skip cleanly (e.g. CI) when those optional deps are not installed.
pytest.importorskip("findspingroup")
pytest.importorskip("ase")

POSCAR = Path(__file__).parent / "references" / "case12_POSCAR"
REFERENCE = Path(__file__).parent / "references" / "case12_golden_kpoints.txt"


def _case12_answers():
    return "\n".join([str(POSCAR), "0 0 1", "5 -5", "", "", ""]) + "\n"


def test_interactive_modify_case12_golden(tmp_path, monkeypatch, capsys):
    try:
        from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
    except Exception as exc:  # pragma: no cover - deps missing
        pytest.skip(f"kpoints/deps unavailable: {exc}")

    # Canned answers: structure file, spin axis, moments, [Enter]=auto path,
    # [Enter]=Option 1 op, [Enter]=vasp. Output filename is fixed
    # (KPOINTS_alter), no longer prompted.
    answers = _case12_answers()
    monkeypatch.chdir(tmp_path)               # side-effect files land in tmp
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))

    KPointsModifier().interactive_modify()

    stdout = capsys.readouterr().out
    assert "Nonmagnetic primitive cell:   SG P4_2/mnm (136)" in stdout
    assert "Magnetic primitive cell (G0): SG P4_2/mnm (136)" in stdout

    produced = (tmp_path / "KPOINTS_alter").read_text(encoding="utf-8")
    expected = REFERENCE.read_text(encoding="utf-8")
    assert produced.splitlines() == expected.splitlines()

    output_dir = tmp_path / OUTPUT_DIR
    full_header = (output_dir / "spin_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[0]
    flip_headers = (output_dir / "spin_flip_operations.txt").read_text(
        encoding="utf-8"
    ).splitlines()[:3]
    assert full_header == (
        "# Basis: submitted structure 'case12_POSCAR' real-space "
        "fractional basis (a1, a2, a3)."
    )
    assert flip_headers == [
        "# Left basis: submitted structure 'case12_POSCAR' real-space "
        "fractional basis (a1, a2, a3).",
        "# Right basis: SeeK-path standardized primitive real-space "
        "fractional basis (a1, a2, a3).",
        "# k mapping: k' = R^(-T) k (mod G) in each corresponding reciprocal "
        "basis (b1, b2, b3).",
    ]

    standard_vasp = output_dir / "case12_POSCAR_seekpath_standard.vasp"
    standard_mcif = output_dir / "case12_POSCAR_seekpath_standard.mcif"
    assert standard_vasp.exists()
    assert standard_mcif.exists()
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifParser
    structural_standard = Structure.from_file(standard_vasp)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        magnetic_standard = CifParser(str(standard_mcif)).parse_structures(
            primitive=False
        )[0]
    assert np.allclose(
        magnetic_standard.lattice.matrix,
        structural_standard.lattice.matrix,
    )
    moment_norms = sorted(
        np.linalg.norm(site.properties["magmom"].moment)
        for site in magnetic_standard
        if "magmom" in site.properties
        and np.linalg.norm(site.properties["magmom"].moment) > 1e-8
    )
    assert np.allclose(moment_norms, [5.0, 5.0])


def test_step4_path_construction_error_reaches_workflow_boundary(
    tmp_path, monkeypatch, capsys
):
    from alterseek import kpoints as kpoints_module

    def fail_path_construction(*args, **kwargs):
        raise RuntimeError("synthetic butterfly construction failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(
        kpoints_module.KPointsModifier,
        "insert_general_kpoints",
        fail_path_construction,
    )

    with pytest.raises(RuntimeError, match="synthetic butterfly construction failure"):
        kpoints_module.KPointsModifier().interactive_modify()

    output = capsys.readouterr().out
    assert "Error processing k-points" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_optional_spin_figure_failure_does_not_block_kpoints(
    tmp_path, monkeypatch, capsys
):
    from alterseek import kpoints as kpoints_module

    def fail_spin_figures(*args, **kwargs):
        raise RuntimeError("synthetic spin-figure failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(
        kpoints_module.KPointsModifier,
        "_generate_spin_figures",
        fail_spin_figures,
    )
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: None)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "[Warning] Could not generate spin figures" in output
    assert "synthetic spin-figure failure" in output
    assert (tmp_path / "KPOINTS_alter").exists()


def test_optional_figure1_failure_does_not_block_kpoints(
    tmp_path, monkeypatch, capsys
):
    from alterseek import compute_centroid_hybrid as centroid_module
    from alterseek import kpoints as kpoints_module

    def fail_figure1(*args, **kwargs):
        raise RuntimeError("synthetic Figure 1 failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(centroid_module, "setup_3d_ax", fail_figure1)
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: None)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "[Warning] Could not generate Figure 1" in output
    assert "synthetic Figure 1 failure" in output
    assert "IBZ centroid construction failed" not in output
    assert (tmp_path / "KPOINTS_alter").exists()


def test_optional_bz_geometry_failure_does_not_block_kpoints(
    tmp_path, monkeypatch, capsys
):
    from alterseek import compute_centroid_hybrid as centroid_module
    from alterseek import kpoints as kpoints_module

    def fail_bz_geometry(*args, **kwargs):
        raise RuntimeError("synthetic BZ geometry failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(centroid_module, "get_bz_loops", fail_bz_geometry)
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: None)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "[Warning] Could not generate Figure 1" in output
    assert "synthetic BZ geometry failure" in output
    assert "IBZ centroid construction failed" not in output
    assert "Spin figures were skipped" in output
    assert (tmp_path / "KPOINTS_alter").exists()


def test_optional_standardization_diagnostic_failure_does_not_block_kpoints(
    tmp_path, monkeypatch, capsys
):
    from alterseek import compute_centroid_hybrid as centroid_module
    from alterseek import kpoints as kpoints_module

    def fail_standardization(*args, **kwargs):
        raise OSError("synthetic standardization diagnostic failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(
        centroid_module,
        "_write_seekpath_standard_poscar",
        fail_standardization,
    )
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: None)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "Could not write SeeK-path standardization files" in output
    assert "synthetic standardization diagnostic failure" in output
    assert "IBZ centroid construction failed" not in output
    assert (tmp_path / "KPOINTS_alter").exists()


def test_optional_symbolic_diagnostic_failure_does_not_block_kpoints(
    tmp_path, monkeypatch, capsys
):
    from alterseek import compute_centroid_hybrid as centroid_module
    from alterseek import kpoints as kpoints_module

    def fail_symbolic_centroid(*args, **kwargs):
        raise RuntimeError("synthetic symbolic diagnostic failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(
        centroid_module,
        "compute_symbolic_centroid",
        fail_symbolic_centroid,
    )
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: None)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "Symbolic IBZ centroid unavailable" in output
    assert "synthetic symbolic diagnostic failure" in output
    assert "IBZ centroid construction failed" not in output
    assert (tmp_path / "KPOINTS_alter").exists()


def test_deferred_figure_failure_does_not_skip_later_saves_or_cleanup(
    monkeypatch, capsys
):
    from alterseek import kpoints as kpoints_module

    events = []
    closed = []

    class DeferredFigure:
        def __init__(self, name, fail=False):
            self.name = name

            def save_after_show():
                events.append(name)
                if fail:
                    raise RuntimeError(f"synthetic {name} save failure")

            self._alterseek_save_after_show = save_after_show

    first = DeferredFigure("first", fail=True)
    second = DeferredFigure("second")
    monkeypatch.setattr(kpoints_module.plt, "show", lambda: events.append("show"))
    monkeypatch.setattr(kpoints_module.plt, "close", closed.append)

    kpoints_module._display_and_save_figures([first, second])

    assert events == ["show", "first", "second"]
    assert closed == [first, second]
    output = capsys.readouterr().out
    assert "synthetic first save failure" in output


def test_custom_path_conversion_error_reaches_workflow_boundary(
    tmp_path, monkeypatch, capsys
):
    from alterseek import kpoints as kpoints_module

    custom_path = tmp_path / "KPATH.in"
    custom_path.write_text(
        "Custom path\n30\nLine-Mode\nReciprocal\n"
        "0.0 0.0 0.0 GAMMA\n0.5 0.0 0.0 X\n",
        encoding="utf-8",
    )
    answers = "\n".join([
        str(POSCAR), "0 0 1", "5 -5", str(custom_path)
    ]) + "\n"

    def fail_conversion(*args, **kwargs):
        raise RuntimeError("synthetic custom-path conversion failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(answers))
    monkeypatch.setattr(
        kpoints_module.KPointsModifier,
        "convert_custom_path_from_input_basis",
        fail_conversion,
    )

    with pytest.raises(
        RuntimeError, match="synthetic custom-path conversion failure"
    ):
        kpoints_module.KPointsModifier().interactive_modify()

    output = capsys.readouterr().out
    assert "Successfully read 2 k-points" in output
    assert "[Error] synthetic custom-path conversion failure" not in output
    assert not (tmp_path / "KPOINTS_alter").exists()


def test_display_failure_after_kpoints_write_remains_successful(
    tmp_path, monkeypatch, capsys
):
    from alterseek import kpoints as kpoints_module

    def fail_display():
        raise RuntimeError("synthetic display failure")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(_case12_answers()))
    monkeypatch.setattr(kpoints_module.plt, "show", fail_display)

    assert kpoints_module.KPointsModifier().interactive_modify() is True
    output = capsys.readouterr().out
    assert "[Warning] Could not display/save generated figures" in output
    assert "synthetic display failure" in output
    assert (tmp_path / "KPOINTS_alter").exists()
