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


def test_interactive_modify_case12_golden(tmp_path, monkeypatch, capsys):
    try:
        from alterseek.kpoints import KPointsModifier, OUTPUT_DIR
    except Exception as exc:  # pragma: no cover - deps missing
        pytest.skip(f"kpoints/deps unavailable: {exc}")

    # Canned answers: structure file, spin axis, moments, [Enter]=auto path,
    # [Enter]=Option 1 op, [Enter]=vasp. Output filename is fixed
    # (KPOINTS_alter), no longer prompted.
    answers = "\n".join([str(POSCAR), "0 0 1", "5 -5", "", "", ""]) + "\n"
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
