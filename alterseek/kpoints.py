"""KPointsModifier: read a KPOINTS file, insert the IBZ-centroid general
k-point, drive the interactive Step 0-5 workflow, and write KPOINTS output.

Extracted from alterseek_path.py (restructuring phase 4).
"""
import os
from typing import List, Optional
import numpy as np

try:
    from .find_sf_operations import run as find_sf_run
    FIND_SF_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] find_sf_operations unavailable ({_exc}); "
          "Step 0 spin-operation detection is disabled.")
    find_sf_run = None
    FIND_SF_AVAILABLE = False

try:
    from .compute_centroid_hybrid import run as compute_centroid
    from .symmetry import (no_altermagnetism_reason,
                          laue_group_from_spacegroup_number,
                          point_group_from_spacegroup_number,
                          is_valid_2d_spin_flip,
                          is_trivial_2d_spin_flip,
                          describe_spinflip_op)
    from .plotting_3d import (plot_spin_flip_figure,
                             plot_spin_bz_figure,
                             plot_spin_bz_top_view_figure)
    from .plotting_2d import plot_2d_figures
    import matplotlib.pyplot as plt
    CENTROID_AVAILABLE = True
except ImportError as _exc:
    print(f"[Warning] compute_centroid_hybrid/matplotlib unavailable ({_exc}); "
          "centroid and figure generation are disabled.")
    CENTROID_AVAILABLE = False
    compute_centroid = None
    no_altermagnetism_reason = None
    laue_group_from_spacegroup_number = None
    point_group_from_spacegroup_number = None
    is_valid_2d_spin_flip = None
    is_trivial_2d_spin_flip = None
    plot_spin_flip_figure = None
    plot_spin_bz_figure = None
    plot_spin_bz_top_view_figure = None
    describe_spinflip_op = None
    plot_2d_figures = None
    plt = None

from .ssg_setting import (
    prepare_magnetic_setting_files,
    finalize_magnetic_setting_outputs,
)
from .io import (
    _atomic_write_text,
    _load_magnetic_input_data,
    _write_seekpath_standard_mcif,
    write_bandplot_lattice_config,
    write_qe_bandplot_config,
)


STEP0_VERBOSE_SUMMARY = False

INPUT_CONFIG_FILE = "alterseek_input.toml"

_INPUT_CONFIG_KEYS = {
    "structure", "spin_axis", "moments", "path", "flip_option",
    "output_code", "view_elev", "view_azim", "save_pdf", "symprec",
}


def _fmt_coord(value):
    """Format a k-point coordinate, collapsing signed zero to plain zero.

    -0.0 and 0.0 are the same k-point, but they render differently and so make
    otherwise identical KPOINTS files compare unequal. Which one comes out
    depends on sign carried through the basis conversion -- that is, on which
    cell the path was built in -- not on the physics.
    """
    text = f"{value:.10f}"
    return text[1:] if text.startswith("-") and float(text) == 0.0 else text


_NO_ALTERMAGNETISM_LAUE_GROUPS = {'-1', '-3', 'm-3'}

# Width of the "Nonmagnetic parent:" / "Magnetic primitive cell (G0):" labels,
# so the SG/PG/Laue fields line up and a difference between the two cells is
# visible at a glance rather than read word by word.
_CELL_LABEL_WIDTH = 30

# Everything AlterSeeK-Path generates goes here except the two files with a
# downstream consumer: KPOINTS_alter (feeds the band calculation) and
# alterband.toml (read from the working directory by the band plotter).
OUTPUT_DIR = "alterseek_output"



def _figure_basename(struct_file):
    """Name figures after the submitted structure, not an internal cell.

    The magnetic route computes the centroid from a helper structure written
    under a derived name, so without this Figure 1 alone picks up an extra
    filename token that Figures 2-4 do not carry.
    """
    if not struct_file:
        return None
    return os.path.splitext(os.path.basename(struct_file))[0]


def _cell_suffix(sites, lattice_tag):
    """Trailing size/lattice tag describing the cell on that line."""
    parts = []
    if sites:
        parts.append(f"{sites} atoms")
    if lattice_tag and lattice_tag != 'unknown':
        parts.append(lattice_tag)
    return f"[{', '.join(parts)}]" if parts else ""


def _print_cell_rows(rows, note_after_first=None):
    """Print the cell comparison with every field in its own column.

    Space-group symbols and numbers vary in width (``P6_3mc (186)`` against
    ``Cmc2_1 (36)``), so without padding the fields ragged and the reader has
    to hunt for what differs between the two cells -- which is the one thing
    the block exists to show. Widths are taken from the rows actually printed.
    """
    sg_w = max(len(row[1]) for row in rows)
    pg_w = max(len(row[2]) for row in rows)
    laue_w = max(len(row[3]) for row in rows)
    for index, (label, sg, pg, laue, suffix) in enumerate(rows):
        line = (f"{label:<{_CELL_LABEL_WIDTH}}"
                f"SG {sg:<{sg_w}}  PG {pg:<{pg_w}}  Laue {laue:<{laue_w}}  {suffix}")
        print(line.rstrip())
        if index == 0 and note_after_first:
            print(f"{'':<{_CELL_LABEL_WIDTH}}{note_after_first}")


def _g0_symmetry(sf_result, sites=None):
    """Describe the magnetic primitive cell by G0, the spatial part of its SSG.

    FindSpinGroup reports G0 directly, so this is the symmetry that actually
    holds once the magnetic order is accounted for, and reading it off the
    reported group needs no tolerance. Re-detecting symmetry from the magnetic
    cell's own coordinates would instead describe the moment-stripped crystal,
    which for a supercell altermagnet is still the higher-symmetry parent.
    """
    if laue_group_from_spacegroup_number is None:
        return None
    number = sf_result.get('g0_number')
    laue_group = laue_group_from_spacegroup_number(number)
    if laue_group is None:
        return None
    return {
        'label': f"{sf_result.get('g0_symbol')} ({int(number)})",
        'spacegroup_number': int(number),
        'point_group': point_group_from_spacegroup_number(number),
        'laue_group': laue_group,
        'sites': sites,
    }


def _magnetic_cell_needed(sf_result):
    """True when the band path has to be built in the magnetic primitive cell
    rather than in the cell the user submitted.

    Two things can make the submitted cell unusable, and both are read off
    quantities Step 0 already determined -- no extra symmetry detection, and no
    extra tolerance choice:

    1. The magnetic order changes the *space group*, i.e. G0 differs from the
       moment-free space group of the submitted cell (GdAuGe AFM5,
       P6_3mc -> Cmc2_1; MnSe2, Pa-3 -> Pbca). A changed Laue group is one way
       this happens, but not the only one: the order can also remove inversion
       or mirrors while leaving the Laue group intact, and the magnetic cell is
       then still the cell the path belongs in.
    2. The submitted cell is not the nonmagnetic primitive cell, i.e. the
       magnetic order enlarges the cell (BaMnO3, 10 -> 30 sites).
    """
    g0_number = sf_result.get('g0_number')
    nonmagnetic_number = sf_result.get('nonmagnetic_spacegroup_number')
    if g0_number is None or nonmagnetic_number is None:
        # Unknown on either side: keep the magnetic route, the conservative
        # choice -- it is the one that can represent a lowered symmetry.
        return True
    if int(g0_number) != int(nonmagnetic_number):
        return True
    nonmagnetic_sites = sf_result.get('nonmagnetic_sites')
    submitted_sites = sf_result.get('num_atoms')
    if nonmagnetic_sites is None or submitted_sites is None:
        return True
    return int(nonmagnetic_sites) != int(submitted_sites)


def _altermagnetism_gate(sf_result, working_cell_symmetry=None):
    """Return a reason dict when the working cell's Laue group forbids
    altermagnetism, or None when it permits it.

    The gate must judge the cell the k-path is actually built in, by the
    symmetry that cell actually has. Both halves matter for MnSe2 (MAGNDATA
    1.0.47). It is deposited in a cubic Pa-3 parent whose Laue group m-3
    forbids altermagnetism, so judging the submitted cell discards a real
    altermagnet -- but the magnetic cell replacing it is still the cubic Pa-3
    *crystal* once its moments are stripped, so judging that cell's bare
    coordinates discards it too (and does so tolerance-dependently). G0, the
    spatial part of the spin space group, is orthorhombic here (Pbca 61, Laue
    mmm), and that is the symmetry which permits altermagnetism.
    """
    if working_cell_symmetry is not None:
        return no_altermagnetism_reason(
            spacegroup=working_cell_symmetry['spacegroup_number'])
    point_group = sf_result.get('point_group')
    laue_group = sf_result.get('laue_group')
    if no_altermagnetism_reason is not None:
        return no_altermagnetism_reason(point_group)
    if laue_group in _NO_ALTERMAGNETISM_LAUE_GROUPS:
        return {'laue_group': laue_group, 'reason': 'No altermagnetism'}
    return None


def _validate_input_config(config):
    unknown = sorted(set(config) - _INPUT_CONFIG_KEYS)
    if unknown:
        raise ValueError(
            f"unknown setting{'s' if len(unknown) != 1 else ''}: {', '.join(unknown)}"
        )

    for key in ("structure", "spin_axis", "moments", "path", "output_code"):
        if key in config and not isinstance(config[key], str):
            raise ValueError(f"{key} must be a TOML string")

    if "output_code" in config:
        code = config["output_code"].strip().lower()
        if code not in {"vasp", "qe"}:
            raise ValueError("output_code must be \"vasp\" or \"qe\"")

    if "flip_option" in config:
        choice = config["flip_option"]
        if isinstance(choice, bool) or not isinstance(choice, int) or choice < 1:
            raise ValueError("flip_option must be a positive integer")

    has_elev = "view_elev" in config
    has_azim = "view_azim" in config
    if has_elev != has_azim:
        raise ValueError("view_elev and view_azim must be supplied together")
    for key in ("view_elev", "view_azim"):
        if key in config and (isinstance(config[key], bool) or
                              not isinstance(config[key], (int, float))):
            raise ValueError(f"{key} must be a number")

    if "save_pdf" in config and not isinstance(config["save_pdf"], bool):
        raise ValueError("save_pdf must be true or false")

    if "symprec" in config:
        value = config["symprec"]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("symprec must be a number")
        if not value > 0:
            raise ValueError("symprec must be positive")
    return config


def _read_input_config(path=INPUT_CONFIG_FILE):
    """Read optional per-run answers (structure, spin_axis, moments,
    flip_option, path, output_code) from a TOML file. Any key present skips
    its interactive prompt; missing keys (or a missing file) fall back to the
    normal prompt unchanged."""
    if not os.path.exists(path):
        return {}
    try:
        import tomllib
        with open(path, "rb") as f:
            config = tomllib.load(f)
        return _validate_input_config(config)
    except Exception as exc:
        raise ValueError(f"Failed to read {path}: {exc}") from exc


class KPointsModifier:
    def __init__(self, magnetic_setting: bool = True, output_verbose: bool = False,
                 mode_2d: bool = False, input_vacuum_axis: int = 2):
        """magnetic_setting: build the path in the magnetic primitive cell
        (the default). Set False to work in the nonmagnetic parent cell
        instead, which the CLI exposes as --parent-setting.
        """
        self.kpoints_data = []
        self.header_lines = []
        self.extra_general_points = []
        self.kpoints_basis_matrix = None
        self.output_basis_matrix = None
        self.kpoints_basis_rotation = None
        self.magnetic_setting = magnetic_setting
        self.output_verbose = output_verbose
        self.mode_2d = mode_2d
        self.input_vacuum_axis = input_vacuum_axis

    @staticmethod
    def _display_label(label: str) -> str:
        # Console display currently matches the VASP-safe form; delegate so the
        # two normalizations can never drift apart.
        return KPointsModifier._kpoints_label(label)

    @staticmethod
    def _kpoints_label(label: str) -> str:
        """Return labels in a VASP-safe form for KPOINTS files."""
        label = str(label)
        return 'GAMMA' if label.strip().upper() == 'GAMMA' or label == '\u0393' else label

    @classmethod
    def _format_path(cls, path_segments) -> str:
        parts = []
        prev_end = None
        for seg_start, seg_end in path_segments:
            start = cls._kpoints_label(seg_start)
            end = cls._kpoints_label(seg_end)
            if seg_start != prev_end:
                parts.append(f"| {start}-{end}" if prev_end else f"{start}-{end}")
            else:
                parts.append(end)
            prev_end = seg_end
        return "-".join(parts).replace("-| ", " | ")

    @staticmethod
    def _format_matrix(op: np.ndarray) -> str:
        rows = []
        for row in op:
            row_str = " ".join(
                f"{int(x): 3d}" if float(x) % 1 == 0 else f"{x: .2f}"
                for x in row
            )
            rows.append(f"    [ {row_str} ]")
        return "\n".join(rows)

    @staticmethod
    def _count_written_segments(kpoints) -> int:
        count = 0
        for i in range(len(kpoints) - 1):
            start_point = kpoints[i]
            end_point = kpoints[i + 1]
            if start_point is None or end_point is None:
                continue
            if start_point[3] == end_point[3]:
                continue
            if {start_point[3], end_point[3]} == {"k", "k'"}:
                continue
            count += 1
        return count
    
    def read_kpoints_file(self, filename: str = "KPOINTS") -> bool:
        """Read a line-mode KPOINTS file (e.g. generated by VASPKIT)."""
        self.kpoints_data = []
        self.header_lines = []
        self.kpoints_basis_matrix = None
        self.output_basis_matrix = None
        self.kpoints_basis_rotation = None
        try:
            # pymatgen validates the format; import locally to keep module
            # import light.
            from pymatgen.io.vasp.inputs import Kpoints

            kpoints = Kpoints.from_file(filename)
            if kpoints.style.name != "Line_mode":
                print(f"Error: {filename} is not a line-mode KPOINTS file "
                      f"(style: {kpoints.style.name}).")
                return False
            coord_type = kpoints.coord_type or "Reciprocal"
            if not coord_type.lower().startswith("r"):
                print(f"Error: {filename} uses {coord_type} coordinates; only "
                      "Reciprocal line-mode KPOINTS are supported.")
                return False

            with open(filename, 'r', encoding='utf-8-sig') as f:
                lines = f.readlines()
            if len(lines) < 4:
                print(f"Error: {filename} is missing the four-line KPOINTS header.")
                return False
            header_lines = [line.strip() for line in lines[:4]]

            coords_list = list(kpoints.kpts)
            labels = list(kpoints.labels or [])
            if len(coords_list) < 2:
                print(f"Error: {filename} must contain at least one complete path segment.")
                return False
            if len(coords_list) % 2:
                print(
                    f"Error: {filename} contains {len(coords_list)} endpoints; "
                    "line-mode paths require complete endpoint pairs."
                )
                return False
            if len(labels) != len(coords_list) or any(not str(label or "").strip()
                                                      for label in labels):
                print(
                    f"Error: every endpoint in {filename} must have a nonempty label."
                )
                return False

            loaded = []
            for coords, label in zip(coords_list, labels):
                label = (label or "").strip()
                xyz = np.asarray(coords[:3], dtype=float)
                if xyz.shape != (3,) or not np.all(np.isfinite(xyz)):
                    print(f"Error: endpoint '{label}' in {filename} has invalid coordinates.")
                    return False
                loaded.append([float(xyz[0]), float(xyz[1]), float(xyz[2]), label])

            self.header_lines = header_lines
            self.kpoints_data = loaded

            print(f"Successfully read {len(self.kpoints_data)} k-points from {filename}")
            return True

        except FileNotFoundError:
            print(f"Error: File {filename} not found!")
            return False
        except Exception as e:
            print(f"Error reading file: {e}")
            return False

    def convert_custom_path_from_input_basis(self, centroid_result) -> None:
        """Convert a submitted KPOINTS path from the input-cell reciprocal
        basis into the standardized primitive basis used internally.

        Custom path files are defined to use the reciprocal basis of the
        structure submitted at Step 0. Output uses the magnetic primitive
        cell's basis, which differs whenever the magnetic order lowers the
        lattice symmetry, so input and output matrices are kept distinct.
        """
        if not self.kpoints_data:
            raise ValueError("No custom KPOINTS path has been loaded.")

        b_standard = np.asarray(centroid_result['b_matrix'], dtype=float)
        rotation = np.asarray(
            centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
            dtype=float,
        )
        b_internal = b_standard @ rotation
        b_input = np.asarray(centroid_result['b_matrix_input'], dtype=float)
        b_output = np.asarray(
            centroid_result.get('b_matrix_output', b_input),
            dtype=float,
        )

        try:
            input_to_internal = b_input @ np.linalg.inv(b_internal)
            converted = []
            for point in self.kpoints_data:
                frac = np.asarray(point[:3], dtype=float) @ input_to_internal
                converted.append([float(frac[0]), float(frac[1]), float(frac[2]), point[3]])
        except Exception as exc:
            raise RuntimeError(f"Custom KPOINTS basis conversion failed: {exc}") from exc

        self.kpoints_data = converted
        self.kpoints_basis_matrix = b_standard
        self.kpoints_basis_rotation = rotation
        self.output_basis_matrix = b_output
        print("[Basis] Converted custom KPOINTS path from input-cell to standardized basis.")

    def _kpoint_for_output_basis(self, point: List) -> List:
        """Convert an internal k-point to the POSCAR reciprocal basis for VASP."""
        if self.kpoints_basis_matrix is None or self.output_basis_matrix is None:
            return point

        try:
            k_frac = np.array(point[:3], dtype=float)
            b_kpoints = np.array(self.kpoints_basis_matrix, dtype=float)
            if self.kpoints_basis_rotation is not None:
                b_kpoints = b_kpoints @ np.array(self.kpoints_basis_rotation, dtype=float)
            b_output = np.array(self.output_basis_matrix, dtype=float)
            k_out = k_frac @ b_kpoints @ np.linalg.inv(b_output)
        except Exception as exc:
            # Writing unconverted coordinates would silently reproduce the
            # oI3/221-P-d conventional-cell mismatch; stop instead.
            raise RuntimeError(
                f"Output-basis conversion failed for k-point '{point[3]}': {exc}. "
                "Refusing to write unconverted coordinates into KPOINTS."
            ) from exc
        k_out = [k_out[0], k_out[1], k_out[2]]
        if self.mode_2d:
            # The physical 2D path lies in the input vacuum k=0 plane; clear any
            # numerical residue along that axis so the written KPOINTS are exactly
            # in-plane.
            k_out[self.input_vacuum_axis] = 0.0
        return [k_out[0], k_out[1], k_out[2], point[3]]

    def load_flip_operations(self, filename: str = None) -> List[np.ndarray]:
        """Reads pre-calculated rotation matrices from file"""
        if filename is None:
            filename = os.path.join(OUTPUT_DIR, "spin_flip_operations.txt")
        matrices = []
        unique = []
        current_matrix = []
        if not os.path.exists(filename):
            legacy = os.path.join(
                os.path.dirname(filename),
                "preserve_spin_operations.txt"
                if os.path.basename(filename) == "spin_preserve_operations.txt"
                else "flip_spin_operations.txt",
            )
            if filename != legacy and os.path.exists(legacy):
                filename = legacy
        try:
            with open(filename, 'r', encoding='utf-8-sig') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Operation"):
                    continue
                if "|" in line:
                    line = line.split("|", 1)[0].strip()
                
                parts = line.split()
                if len(parts) == 3:
                    current_matrix.append([float(x) for x in parts])
                
                if len(current_matrix) == 3:
                    mat = np.array(current_matrix)
                    if not any(np.allclose(mat, ex, atol=1e-8) for ex in unique):
                        unique.append(mat)
                        matrices.append(mat)
                    current_matrix = []
            return matrices
        except FileNotFoundError:
            return []

    def load_preserve_operations(self, filename: str = None) -> List[np.ndarray]:
        """Reads pre-calculated spin-preserve rotation matrices from file."""
        if filename is None:
            filename = os.path.join(OUTPUT_DIR, "spin_preserve_operations.txt")
        return self.load_flip_operations(filename)

    def transform_kpoint(self, kpoint: List[float], transformation_matrix: np.ndarray) -> List[float]:
        """Apply transformation matrix to k-point: R^(-1)^T * k"""
        k = np.array(kpoint[:3])  # Take only x, y, z coordinates
        # Calculate R^(-1)^T
        R_inv_T = np.linalg.inv(transformation_matrix).T
        k_transformed = R_inv_T @ k
        return k_transformed.tolist()
    
    def insert_general_kpoints(self,
                               general_kpoint: List[float],
                               transformation_matrix: np.ndarray,
                               extra_general_points: Optional[List[List]] = None) -> List[List]:
        """
        Insert general k-points into every segment of the high symmetry path.
        Per-segment butterfly: for each (A,B) segment:
          - If A not yet opened: A-k | k'-A'
          - Connect A' to B'
          - Close: B'-k' | k-B
        Points already butterflied in earlier chains appear as plain segments.
        """
        if not self.kpoints_data:
            print("Error: No k-points data loaded. Please read KPOINTS file first.")
            return []

        k_prime = self.transform_kpoint(general_kpoint, transformation_matrix)
        kpt = general_kpoint
        kp  = k_prime

        def coords_eq(p, q, tol=1e-6):
            return abs(p[0]-q[0]) < tol and abs(p[1]-q[1]) < tol and abs(p[2]-q[2]) < tol

        def pt_key(p):
            return (round(p[0], 6), round(p[1], 6), round(p[2], 6))

        def is_gamma(p):
            label = str(p[3])
            return label.strip().upper() == 'GAMMA' or label == '\u0393'

        def get_prime(p):
            """Return primed version of p (Gamma stays unprimed)."""
            if is_gamma(p):
                return p.copy()
            tc = self.transform_kpoint(p, transformation_matrix)
            return [tc[0], tc[1], tc[2], f"{p[3]}'"]

        # --- Step 1: group flat kpoints_data into segment pairs ---
        raw = self.kpoints_data
        seg_pairs = [(raw[i], raw[i+1]) for i in range(0, len(raw) - 1, 2)]

        # --- Step 2: build connected chains ---
        chains = []
        current_chain = [seg_pairs[0][0], seg_pairs[0][1]]
        for sp_start, sp_end in seg_pairs[1:]:
            if coords_eq(current_chain[-1], sp_start):
                current_chain.append(sp_end)
            else:
                chains.append(current_chain)
                current_chain = [sp_start, sp_end]
        chains.append(current_chain)

        # --- Step 3: alternating plain / butterfly segments ---
        #
        # Pattern:
        #   Even-indexed segment --plain  (emit A, B)
        #   Odd-indexed segment  --butterfly (emit A, k, k', A', B', k', k, B)
        #
        # First chain: parity 0 --segment 0 is plain, segment 1 is butterfly, ...
        # Other chains: parity 1 --segment 0 is butterfly, segment 1 is plain, ...
        #
        # Consecutive segments share an endpoint; the duplicate is suppressed by
        # the same-label skip in write_kpoints_file and the display builder.
        # GAMMA is self-conjugate: get_prime(GAMMA) = GAMMA.
        path_sequence = []

        def emit_plain(A, B):
            path_sequence.append(A.copy())
            path_sequence.append(B.copy())

        def emit_butterfly(A, B):
            path_sequence.append(A.copy())
            path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
            path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
            path_sequence.append(get_prime(A))
            path_sequence.append(get_prime(B))
            # Universal rule: skip close if B already received butterfly treatment.
            # Applies to all points equally (GAMMA, X, W, etc.).
            if pt_key(B) not in butterflied:
                path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                path_sequence.append(B.copy())

        butterflied = set()  # pt_key of points that have received butterfly treatment

        for ci, chain in enumerate(chains):
            # Dedup within chain (by coordinates, not labels)
            unique = []
            for pt in chain:
                if not unique or not coords_eq(unique[-1], pt):
                    unique.append(pt)

            # Degenerate chain: all points share the same coordinates
            # (e.g. U_0 --T for certain lattice parameters). Skip silently.
            if len(unique) < 2:
                labels = [p[3] for p in chain]
                print(f"  [Note] Part {ci+1} ({' - '.join(labels)}) skipped: "
                      f"endpoints coincide in coordinates")
                continue

            if ci > 0:
                path_sequence.append(None)

            start_parity = 0 if ci == 0 else 1

            parity = start_parity  # tracks alternation independently of s
            for s, (A, B) in enumerate(zip(unique, unique[1:])):
                A_key = pt_key(A)
                B_key = pt_key(B)
                is_last = (s == len(unique) - 2)
                # If A was already butterflied, continue from the side that is
                # already active. If the previous butterfly ended at A', the next
                # plain segment must be A'->B, not A'->A->B.
                if A_key in butterflied:
                    A_start = A
                    if path_sequence:
                        prev_pt = path_sequence[-1]
                        A_prime = get_prime(A)
                        if (prev_pt is not None and prev_pt[3] == A_prime[3]
                                and A_prime[3] != A[3]):
                            A_start = A_prime
                    emit_plain(A_start, B)
                    # B is a new endpoint --partial butterfly B-k|k'-B'
                    if B_key not in butterflied and is_last:
                        path_sequence.append(B.copy())
                        path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                        path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                        path_sequence.append(get_prime(B))
                        butterflied.add(B_key)
                    # Only hold parity when butterflied-A overrides a butterfly slot;
                    # if parity is already even (plain slot), advance normally.
                    if parity % 2 == 0:
                        parity += 1
                elif parity % 2 == 1:
                    emit_butterfly(A, B)
                    butterflied.add(A_key)
                    butterflied.add(B_key)
                    parity += 1
                else:
                    emit_plain(A, B)
                    parity += 1

            # After processing all segments in this chain, check if the last
            # point still lacks butterfly treatment (happens when chain 0 has
            # only 1 segment --the single plain pair leaves B un-butterflied).
            last_key = pt_key(unique[-1])
            if last_key not in butterflied:
                last_pt = unique[-1]
                path_sequence.append(last_pt.copy())
                path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
                path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
                path_sequence.append(get_prime(last_pt))
                butterflied.add(last_key)

        # Remove trailing sentinel
        if path_sequence and path_sequence[-1] is None:
            path_sequence.pop()

        # Doubled-IBZ append-only anchors. These are project-only copied
        # vertices that should be sampled through the general point without
        # adding duplicated high-symmetry edges.
        extra_general_points = extra_general_points or []
        for pt in extra_general_points:
            if path_sequence:
                path_sequence.append(None)
            path_sequence.append(pt.copy())
            path_sequence.append([kpt[0], kpt[1], kpt[2], "k"])
            path_sequence.append([kp[0],  kp[1],  kp[2],  "k'"])
            path_sequence.append(get_prime(pt))

        # Print generated path as label string
        tokens = []  # list of label strings or '|' for breaks
        prev = None
        i = 0
        while i < len(path_sequence) - 1:
            cur = path_sequence[i]
            nxt = path_sequence[i + 1]
            # Chain boundary sentinel --insert break
            if cur is None or nxt is None:
                if tokens and tokens[-1] != '|':
                    tokens.append('|')
                i += 1
                continue
            # k↔k' connection --insert break
            if (cur[3] == "k" and nxt[3] == "k'") or \
               (cur[3] == "k'" and nxt[3] == "k"):
                if tokens and tokens[-1] != '|':
                    tokens.append('|')
                i += 1
                continue
            # same label --skip
            if cur[3] == nxt[3]:
                i += 1
                continue
            if prev != cur[3]:
                tokens.append(cur[3])
            tokens.append(nxt[3])
            prev = nxt[3]
            i += 1
        # Build display string: labels joined by '-', breaks become '|'.
        display_segments = []
        current_segment = []
        for t in tokens:
            if t == '|':
                if current_segment:
                    display_segments.append('-'.join(self._display_label(x) for x in current_segment))
                    current_segment = []
            else:
                current_segment.append(t)
        if current_segment:
            display_segments.append('-'.join(self._display_label(x) for x in current_segment))

        if len(display_segments) > 6:
            display = " | ".join(display_segments[:3] + ["..."] + display_segments[-3:])
        else:
            display = " | ".join(display_segments)
        print(f"Generated path: {display}")

        generated_segments = self._count_written_segments(path_sequence)
        generated_points = sum(1 for pt in path_sequence if pt is not None)
        print(f"Full path: {len(seg_pairs)} original segments -> "
              f"{generated_segments} generated segments, {generated_points} k-points")
        if extra_general_points:
            labels = ", ".join(str(pt[3]) for pt in extra_general_points)
            print(f"Added doubled-IBZ general anchors: {labels}")

        return path_sequence

    def insert_general_kpoint_anchors(self, kpoint: List[float],
                                      extra_general_points: Optional[List[List]] = None) -> List[List]:
        """Keep the ordinary path, then append compact high-symmetry/k comparisons."""
        kpt = [kpoint[0], kpoint[1], kpoint[2], "k"]
        raw = self.kpoints_data
        seg_pairs = [(raw[i], raw[i + 1]) for i in range(0, len(raw) - 1, 2)]
        path_sequence = []
        for idx, (start, end) in enumerate(seg_pairs):
            if idx:
                path_sequence.append(None)
            path_sequence.append(start.copy())
            path_sequence.append(end.copy())

        seen = set()
        anchors = []
        for pt in self.kpoints_data + (extra_general_points or []):
            if pt is None:
                continue
            key = (
                round(float(pt[0]), 8),
                round(float(pt[1]), 8),
                round(float(pt[2]), 8),
                str(pt[3]),
            )
            if key in seen:
                continue
            seen.add(key)
            anchors.append(pt)

        for idx in range(0, len(anchors), 2):
            if path_sequence:
                path_sequence.append(None)
            path_sequence.append(anchors[idx].copy())
            path_sequence.append(kpt.copy())
            if idx + 1 < len(anchors):
                path_sequence.append(anchors[idx + 1].copy())

        pair_count = len(anchors) // 2
        leftover = len(anchors) % 2
        print(
            f"Kept ordinary path and added {pair_count} A-k-B comparison segments"
            f"{' plus 1 A-k tail' if leftover else ''}."
        )
        return path_sequence
    
    def _valid_segment_pairs(self, new_kpoints: List[List]):
        """Collect the writable path segments shared by the VASP and QE writers.

        Skips chain-boundary sentinels (``None``), direct k<->k' connections,
        and zero-length (same-label) pairs. Returns a list of
        ``(start_out, end_out, break_before, index, end_raw_label)`` tuples:
        endpoints converted to the output basis, ``break_before`` marking a
        discontinuity with the previously kept segment, ``index`` the position
        of the start point in *new_kpoints*, and the raw end-point label.
        """
        pairs = []
        forced_break = False
        i = 0
        while i < len(new_kpoints) - 1:
            sp = new_kpoints[i]
            ep = new_kpoints[i + 1]
            if sp is None or ep is None:
                i += 1
                continue
            if (sp[3] == "k" and ep[3] == "k'") or (sp[3] == "k'" and ep[3] == "k"):
                forced_break = True
                i += 1
                continue
            if sp[3] == ep[3]:
                i += 1
                continue
            sp_out = self._kpoint_for_output_basis(sp)
            ep_out = self._kpoint_for_output_basis(ep)
            break_before = forced_break
            if pairs:
                previous_end = pairs[-1][1]
                break_before = break_before or not np.allclose(
                    previous_end[:3], sp_out[:3], atol=1e-10, rtol=0.0
                )
            pairs.append((sp_out, ep_out, break_before, i, ep[3]))
            forced_break = False
            i += 1
        return pairs

    def _general_kpoint_output_basis(self, general_kpoint) -> Optional[List[float]]:
        """Return the general point k in the output (input-cell) basis, or
        None when the conversion is unavailable. When the standardized and
        input bases coincide the returned value equals the input."""
        if general_kpoint is None:
            return None
        if self.kpoints_basis_matrix is None or self.output_basis_matrix is None:
            return None
        try:
            out = self._kpoint_for_output_basis([*general_kpoint[:3], "k"])
        except Exception:
            return None
        return [float(out[0]), float(out[1]), float(out[2])]

    def write_kpoints_file(self, new_kpoints: List[List], output_file: str = "KPOINTS_alter",
                           transformation_matrix: Optional[np.ndarray] = None,
                           transformation_label: Optional[str] = None):
        """Write modified KPOINTS file with proper Line-Mode format and discontinuity"""
        try:
            if transformation_matrix is not None:
                flat_matrix = np.array(transformation_matrix).flatten()
                matrix_str = " ".join(f"{x:.8f}" for x in flat_matrix)
                label = f" ({transformation_label})" if transformation_label else ""
                first_line = (
                    f"Selected spin-flip operation{label} in input-cell "
                    f"fractional basis: {matrix_str}\n"
                )
            else:
                first_line = f"{self.header_lines[0]}\n"

            lines = [
                first_line,
                "   30\n",
                f"{self.header_lines[2]}\n",
                f"{self.header_lines[3]}\n",
            ]

            segments = self._valid_segment_pairs(new_kpoints)
            for start_point_out, end_point_out, _break_before, i, end_raw_label in segments:
                start_label = self._kpoints_label(start_point_out[3])
                end_label = self._kpoints_label(end_point_out[3])
                lines.append(
                    f"   {_fmt_coord(start_point_out[0])}   {_fmt_coord(start_point_out[1])}   "
                    f"{_fmt_coord(start_point_out[2])}     {start_label}\n"
                )
                lines.append(
                    f"   {_fmt_coord(end_point_out[0])}   {_fmt_coord(end_point_out[1])}   "
                    f"{_fmt_coord(end_point_out[2])}     {end_label}\n"
                )
                if end_raw_label == "k" or i < len(new_kpoints) - 2:
                    lines.append("\n")

            if not segments:
                print("Error writing file: path contains no writable segments.")
                return False
            _atomic_write_text(output_file, "".join(lines))
            
            print(f"Modified KPOINTS file written to: {output_file}")
            return True
            
        except Exception as e:
            print(f"Error writing file: {e}")
            return False

    def write_kpoints_file_qe(self, new_kpoints: List[List],
                              output_file: str = "KPOINTS_alter_qe",
                              transformation_matrix: Optional[np.ndarray] = None,
                              transformation_label: Optional[str] = None,
                              ninterp: int = 30):
        """Write KPOINTS in QE K_POINTS crystal_b format."""
        try:
            valid_pairs = self._valid_segment_pairs(new_kpoints)
        except Exception as exc:
            print(f"Error writing QE KPOINTS: {exc}")
            return False

        if not valid_pairs:
            print("Error writing QE KPOINTS: path contains no writable segments.")
            return False

        # Build sequential waypoint list, deduplicating shared endpoints.
        # ninterp rule: dead-end k or last waypoint → 1; everything else → ninterp.
        waypoints = []
        for idx, (sp_out, ep_out, break_before, _i, _end_raw_label) in enumerate(valid_pairs):
            if idx == 0:
                waypoints.append((sp_out, ninterp))
            elif break_before:
                previous_point, _ = waypoints[-1]
                waypoints[-1] = (previous_point, 1)
                waypoints.append((sp_out, ninterp))
            waypoints.append((ep_out, ninterp))
        final_point, _ = waypoints[-1]
        waypoints[-1] = (final_point, 1)

        try:
            lines = ["K_POINTS crystal_b\n", f"  {len(waypoints)}\n"]
            for pt_out, ni in waypoints:
                lbl = self._kpoints_label(pt_out[3])
                lines.append(
                    f"  {_fmt_coord(pt_out[0])}  {_fmt_coord(pt_out[1])}  "
                    f"{_fmt_coord(pt_out[2])}  {ni:3d}  ! {lbl}\n"
                )
            if transformation_matrix is not None:
                flat = np.array(transformation_matrix).flatten()
                mat_str = " ".join(f"{x:.8f}" for x in flat)
                lbl = f" ({transformation_label})" if transformation_label else ""
                lines.append(
                    f"! Spin-flip operation{lbl} in input-cell fractional basis: "
                    f"{mat_str}\n"
                )
            _atomic_write_text(output_file, "".join(lines))
            print(f"QE KPOINTS written to: {output_file}")
            return True
        except Exception as e:
            print(f"Error writing QE KPOINTS: {e}")
            return False

    def _generate_spin_figures(self, centroid_result, struct_file, general_kpoint,
                               R_for_kpts, R_cart_for_plot, flip_ops_for_plot,
                               preserve_ops_for_plot, new_kpoints, display_figures,
                               save_pdf=False):
        """Generate Figures 2-4 (spin-flip / spin-BZ / kz=0 top view) for the
        selected spin-flip operation; append any created figures to
        display_figures.  Extracted from interactive_modify (phase 5)."""
        # Generate Figures 2-4 (spin-flip, spin-BZ, kz=0 top view).
        # One shared call scaffold; per-figure kwargs hold the
        # differences between the three plots.
        if centroid_result is not None and self.mode_2d:
            # 2D slab mode: render dedicated top-down 2D figures instead
            # of the tilted 3D BZ plate.
            basename = (os.path.splitext(os.path.basename(struct_file))[0]
                        if struct_file else 'output')
            if plot_2d_figures is not None and 'b_matrix' in centroid_result:
                try:
                    plot_2d_figures(
                        centroid_result, general_kpoint, R_for_kpts,
                        basename, output_dir=OUTPUT_DIR,
                        flip_ops_for_plot=(flip_ops_for_plot
                                           if flip_ops_for_plot else None),
                        save_pdf=save_pdf,
                    )
                except Exception as _e:
                    print(f"[Warning] Could not generate 2D figures: {_e}")
        elif centroid_result is not None:
            basename = (os.path.splitext(os.path.basename(struct_file))[0]
                        if struct_file else 'output')
            sc_type = centroid_result.get('sc_type', 'BZ')
            flip_kwargs = dict(
                flip_ops_frac=flip_ops_for_plot if flip_ops_for_plot else None,
                preserve_ops_frac=preserve_ops_for_plot if preserve_ops_for_plot else None,
            )
            view_kwargs = dict(
                bz_center=centroid_result.get('bz_center'),
                bz_span=centroid_result.get('bz_span'),
                elev=centroid_result.get('elev', 14),
                azim=centroid_result.get('azim', 20),
            )
            figure_specs = []
            if plot_spin_flip_figure is not None and 'b_matrix' in centroid_result:
                figure_specs.append((
                    plot_spin_flip_figure, 'spin-flip',
                    os.path.join(OUTPUT_DIR, f'{basename}_spinflip_{sc_type}.png'),
                    dict(
                        kpoints_data=self.kpoints_data,
                        ibz_kpoints_frac=centroid_result.get('ibz_kpoints_frac', {}),
                        centroid_frac=general_kpoint,
                        R_cart=R_cart_for_plot,
                        block=False,
                        path_sequence=new_kpoints,
                        unique_ops=centroid_result.get('unique_ops'),
                        **view_kwargs,
                    ),
                ))
            if 'unique_ops' in centroid_result:
                spin_bz_kwargs = dict(
                    unique_ops=centroid_result['unique_ops'],
                    centroid_cart=centroid_result.get('centroid_cart'),
                    hull_labels=centroid_result.get('hull_labels'),
                    **flip_kwargs,
                )
                top_view_z0 = 0.0
                # b_matrix[2, 2] (the z-component of the standardized
                # reciprocal frame's 3rd vector) is not a reliable cut
                # height: depending on how seekpath orients the primitive
                # cell, no single reciprocal vector may point along
                # Cartesian z at all, making that element exactly zero
                # (silently defeating any 0.5 * b_mat[2, 2] cut) even
                # though the actual 3D BZ solid still spans real Cartesian
                # kz range. Use the true kz half-extent of the BZ
                # boundary itself instead, which is robust to orientation.
                if (sc_type in ('hR1', 'hR2') or sc_type.startswith('c')) \
                        and 'bz_loops' in centroid_result:
                    # Rhombohedral (hR1/hR2) and cubic AM point groups
                    # (432/-43m/m-3m) reaching this branch have no vertical
                    # mirror, so a kz=0 (Gamma-centered) cut sits on a
                    # degenerate/low-symmetry height and shows a spurious
                    # or collapsed domain split instead of the true bulk
                    # pattern. Cut away from kz=0 instead: 0.5 of the BZ's
                    # kz extent for hR1/hR2, 0.25 for cubic (cubic's kz=0.5
                    # cut can coincide with a boundary/high-symmetry plane
                    # that makes the pattern look identical to the
                    # tetragonal case; 0.25 avoids that).
                    bz_pts = np.vstack(centroid_result['bz_loops'])
                    z_max = float(np.abs(bz_pts[:, 2]).max())
                    z_frac = 0.5 if sc_type in ('hR1', 'hR2') else 0.25
                    top_view_z0 = z_frac * z_max
                if plot_spin_bz_figure is not None:
                    figure_specs.append((
                        plot_spin_bz_figure, 'spin-BZ',
                        os.path.join(OUTPUT_DIR, f'{basename}_spinbz_{sc_type}.png'),
                        dict(z0=top_view_z0, **spin_bz_kwargs, **view_kwargs),
                    ))
                if plot_spin_bz_top_view_figure is not None:
                    figure_specs.append((
                        plot_spin_bz_top_view_figure, 'spin-BZ top-view',
                        os.path.join(OUTPUT_DIR, f'{basename}_spinbz_top_{sc_type}.png'),
                        dict(z0=top_view_z0, **spin_bz_kwargs),
                    ))
            for plot_fn, fig_name, fig_path, extra_kwargs in figure_specs:
                try:
                    fig = plot_fn(
                        b_matrix=centroid_result['b_matrix'],
                        bz_loops=centroid_result['bz_loops'],
                        hull_pts=centroid_result.get('hull_pts'),
                        hull_simplices=centroid_result.get('hull_simplices'),
                        R=R_for_kpts,
                        output_path=fig_path,
                        show_plot=True,
                        defer_show=True,
                        save_pdf=save_pdf,
                        **extra_kwargs,
                    )
                    if fig is not None:
                        display_figures.append(fig)
                except Exception as _e:
                    print(f"[Warning] Could not generate {fig_name} figure: {_e}")

    def _select_spin_flip_operation(
        self,
        flip_ops,
        centroid_result,
        preset_choice=None,
        operation_basis_label="operation-source structure",
    ):
        """Step 3: choose the spin-flip operation R (default Option 1 /
        numbered / 'list' / 'manual').  Returns (R, selected_transformation_label).
        Extracted from interactive_modify (phase 5). `preset_choice`, if given
        (from alterseek_input.toml's `flip_option`), supplies the answer for
        the first prompt instead of reading stdin; run() validates it (a
        positive integer within range) before calling, so a preset always
        selects its numbered option directly."""
        def _op_name(op_input):
            if describe_spinflip_op is None:
                return ""
            try:
                if centroid_result is not None and 'b_matrix' in centroid_result:
                    b_mat = np.array(centroid_result['b_matrix'], dtype=float)
                    b_in = np.array(centroid_result.get(
                        'b_matrix_input', centroid_result.get('b_matrix_conv', b_mat)),
                        dtype=float)
                    b_prim = b_mat @ np.array(
                        centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                        dtype=float)
                    R_arr = np.array(op_input, dtype=float)
                    R_cart_in = b_in.T @ np.linalg.inv(R_arr).T @ np.linalg.inv(b_in.T)
                    R_prim_inv_T = np.linalg.inv(b_prim.T) @ R_cart_in @ b_prim.T
                    R_prim = np.linalg.inv(R_prim_inv_T.T)
                    R_cart = b_mat.T @ np.linalg.inv(R_prim).T @ np.linalg.inv(b_mat.T)
                    return describe_spinflip_op(R_cart, b_mat)
                # No standardized basis available: type/order are still correct
                # from the invariants; omit the (basis-dependent) axis.
                return describe_spinflip_op(np.array(op_input, dtype=float), None)
            except Exception:
                return ""

        R = None
        selected_transformation_label = None
        if flip_ops:
            print(f"Found {len(flip_ops)} spin-flip operations R.")
            _names_available = (describe_spinflip_op is not None
                                and centroid_result is not None
                                and 'b_matrix' in centroid_result)
            if _names_available:
                print(
                    f"  Note: R is in the {operation_basis_label} "
                    "fractional basis;"
                )
                print("  rotation axis/mirror plane indices are in the reciprocal (b1,b2,b3) basis.")
            print("Default R: Option 1")
            _preset_pending = preset_choice is not None
            while R is None:
                print("Press [Enter] to use default, type a number, 'list' to show matrices, or 'manual': ", end='', flush=True)
                if _preset_pending:
                    choice = str(preset_choice).strip().lower()
                    _preset_pending = False
                    print(choice)
                else:
                    choice = input().strip().lower()

                if not choice:
                    R = flip_ops[0]
                    selected_transformation_label = "Option 1"
                    _nm = _op_name(flip_ops[0])
                    print("Selected: Option 1" + (f"  ({_nm})" if _nm else ""))
                elif choice == 'list':
                    for i, op in enumerate(flip_ops):
                        _nm = _op_name(op)
                        print(f"\n  Option {i+1}:" + (f"  {_nm}" if _nm else ""))
                        print(self._format_matrix(op))
                    print()
                elif choice == 'manual':
                    break
                else:
                    try:
                        idx = int(choice) - 1
                        if 0 <= idx < len(flip_ops):
                            R = flip_ops[idx]
                            selected_transformation_label = f"Option {idx+1}"
                            _nm = _op_name(flip_ops[idx])
                            print(f"Selected: Option {idx+1}"
                                  + (f"  ({_nm})" if _nm else ""))
                        else:
                            print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")
                    except ValueError:
                        print(f"Please choose 1-{len(flip_ops)}, 'list', or 'manual'.")

        if R is None:
            print("Enter custom transformation matrix R.")
            print("Enter row by row (3 numbers per row, space-separated):")
            transformation_matrix = []
            for i in range(3):
                while True:
                    try:
                        print(f"Row {i+1}: ", end='', flush=True)
                        row_input = input().strip().split()
                        if len(row_input) == 3:
                            transformation_matrix.append([float(x) for x in row_input])
                            break
                        print("Please enter exactly 3 numbers.")
                    except ValueError: pass
            R = np.array(transformation_matrix)
            selected_transformation_label = "manual"
        return R, selected_transformation_label

    def _convert_operation_to_primitive_basis(
        self,
        R,
        flip_ops,
        preserve_ops,
        centroid_result,
        operation_basis_label,
    ):
        """Convert operations from their source-cell fractional basis to the
        SeeK-path primitive basis and annotate the files with both bases."""
        # The ordinary route writes FindSpinGroup rotations in the submitted
        # structure's fractional basis. The magnetic-cell route replaces them
        # with rotations in the magnetic primitive fractional basis. IBZ
        # coordinates are in the SeeK-path primitive reciprocal basis, which
        # can differ for centered lattices (e.g. BCT, RHL).
        # Convert through Cartesian k-space so k', Figure 2, and the path all use
        # the same physical spin-flip operation:
        #   R_cart_k       = b_input.T @ inv(R_input).T @ inv(b_input.T)
        #   R_prim^{-T}    = inv(b_prim.T) @ R_cart_k @ b_prim.T
        R_cart_for_plot = None
        flip_ops_for_plot = flip_ops
        preserve_ops_for_plot = preserve_ops
        if centroid_result is not None and 'b_matrix' in centroid_result:
            _b_input = np.array(centroid_result.get('b_matrix_input',
                                                    centroid_result['b_matrix_conv']), dtype=float)
            _b_prim = np.array(centroid_result['b_matrix'], dtype=float) @ np.array(
                centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                dtype=float,
            )
            def _convert_input_frac_R_to_prim(_R):
                _R_arr = np.array(_R, dtype=float)
                _R_cart = _b_input.T @ np.linalg.inv(_R_arr).T @ np.linalg.inv(_b_input.T)
                _R_prim_inv_T = np.linalg.inv(_b_prim.T) @ _R_cart @ _b_prim.T
                return np.linalg.inv(_R_prim_inv_T.T), _R_cart

            R_for_kpts, _ = _convert_input_frac_R_to_prim(R)
            # Figure 2 draws the HPKOT hull in seekpath's standardized
            # Cartesian frame.  Let it reconstruct R_cart from R_for_kpts in
            # that same frame; the Cartesian matrix above remains in the
            # orientation of the input structure (notably different for MCIF).
            R_cart_for_plot = None
            if flip_ops:
                flip_ops_for_plot = [
                    _convert_input_frac_R_to_prim(op)[0] for op in flip_ops
                ]
            if preserve_ops:
                preserve_ops_for_plot = [
                    _convert_input_frac_R_to_prim(op)[0] for op in preserve_ops
                ]
            def _annotate_ops_with_standardized_basis(filename, input_ops, standardized_ops, label):
                try:
                    with open(filename, "w", encoding="utf-8", newline="\n") as f:
                        f.write(
                            f"# Left basis: {operation_basis_label} real-space "
                            "fractional basis (a1, a2, a3).\n"
                        )
                        f.write(
                            "# Right basis: SeeK-path standardized primitive "
                            "real-space fractional basis (a1, a2, a3).\n"
                        )
                        f.write(
                            "# k mapping: k' = R^(-T) k (mod G) in each "
                            "corresponding reciprocal basis (b1, b2, b3).\n"
                        )
                        f.write(f"# Found {len(input_ops)} inversion-extended {label} point operations\n")
                        for i, (input_op, std_op) in enumerate(zip(input_ops, standardized_ops), 1):
                            f.write(f"Operation_{i}\n")
                            input_int = np.rint(np.array(input_op, dtype=float)).astype(int)
                            std_int = np.rint(np.array(std_op, dtype=float)).astype(int)
                            for input_row, std_row in zip(input_int, std_int):
                                left = " ".join(str(int(x)) for x in input_row)
                                right = " ".join(str(int(x)) for x in std_row)
                                f.write(f"{left}    |    {right}\n")
                            f.write("\n")
                except Exception as exc:
                    print(f"[Warning] Could not write {filename}: {exc}")

            _annotate_ops_with_standardized_basis(
                os.path.join(OUTPUT_DIR, "spin_flip_operations.txt"),
                flip_ops,
                flip_ops_for_plot,
                "spin-flipping",
            )
            _annotate_ops_with_standardized_basis(
                os.path.join(OUTPUT_DIR, "spin_preserve_operations.txt"),
                preserve_ops,
                preserve_ops_for_plot,
                "spin-preserving",
            )
            # R_for_kpts (the selected operation in the standardized primitive
            # basis) is intentionally not printed: both bases are always saved to
            # spin_flip_operations.txt for reference (left = input, right =
            # standardized). To echo it on screen -- it differs from the Step-3
            # input-basis matrix only for centered cells (BCT/RHL/base-centered) --
            # re-enable:
            #     if not np.allclose(np.asarray(R, dtype=float), R_for_kpts, atol=1e-6):
            #         print("[Basis] Spin-flip operation R in the standardized "
            #               "primitive basis (converted from the input-cell basis; "
            #               "used for KPOINTS and figures):")
            #         print(self._format_matrix(R_for_kpts))
        else:
            R_for_kpts = R
        return R_for_kpts, R_cart_for_plot, flip_ops_for_plot, preserve_ops_for_plot

    def interactive_modify(self):
        """Interactive modification of KPOINTS file"""
        BOLD  = "\033[1m"
        RESET = "\033[0m"
        print("=== Altermagnetic K-Path Generator ===")

        try:
            input_config = _read_input_config()
        except ValueError as exc:
            print(f"[Error] {exc}")
            return False
        view_elev = float(input_config['view_elev']) if 'view_elev' in input_config else None
        view_azim = float(input_config['view_azim']) if 'view_azim' in input_config else None
        symprec = float(input_config['symprec']) if 'symprec' in input_config else None
        save_pdf = bool(input_config.get('save_pdf', False))

        def _ask(prompt_text, key):
            print(prompt_text, end='', flush=True)
            if key in input_config:
                val = str(input_config[key]).strip()
                print(val)
                return val
            return input().strip()

        def _choose_output_code():
            while True:
                choice = _ask("Output code ([vasp]/qe): ", "output_code").lower() or "vasp"
                if choice in {"vasp", "qe"}:
                    return choice
                print("Invalid output code. Enter 'vasp', 'qe', or press Enter for VASP.")

        # Step 0: Compute spin-flip operations from structure
        print(f"\n{BOLD}>>> Step 0: Spin symmetry{RESET}")
        struct_file = _ask(
            "Enter structure file (default: POSCAR, supports .vasp/.cif/.mcif): ",
            "structure",
        )
        if not struct_file: struct_file = "POSCAR"

        # None = Step 0 not run; True = file freshly written; False = ran but no flip ops found
        _step0_wrote_flip_file = None
        standard_path_reason = None
        standard_path_reason_reported = False
        centroid_result = None
        centroid_error = None
        centroid_struct_file = struct_file
        centroid_seekpath_type_numbers = None
        operation_basis_label = (
            f"submitted structure '{os.path.basename(struct_file)}'"
        )
        display_figures = []
        self.extra_general_points = []

        def _load_custom_path(custom_filename):
            """Read a custom KPATH.in/KPOINTS path and convert it from the
            input-cell basis to the standardized internal basis."""
            if not self.read_kpoints_file(custom_filename):
                return False
            if centroid_result is not None:
                try:
                    self.convert_custom_path_from_input_basis(centroid_result)
                except Exception as exc:
                    print(f"[Error] {exc}")
                    return False
            return True

        def _save_and_finish(path_points, R_matrix, R_label):
            """Step 5: choose the output code, write the KPOINTS file (plus
            its band-plot config), and show any deferred figures."""
            print(f"\n{BOLD}>>> Step 5: Save{RESET}")
            code_choice = _choose_output_code()
            if code_choice == "qe":
                write_ok = self.write_kpoints_file_qe(
                    path_points, "KPOINTS_alter_qe", R_matrix, R_label
                )
                if write_ok:
                    write_qe_bandplot_config()
            else:
                write_ok = self.write_kpoints_file(
                    path_points, "KPOINTS_alter", R_matrix, R_label
                )
                if write_ok and centroid_result is not None:
                    write_bandplot_lattice_config(
                        centroid_result.get('lattice_key', centroid_result.get('sc_type'))
                    )
            if not write_ok:
                print("[Error] KPOINTS output was not written.")
                return False
            print("\nDone.")
            if display_figures and plt is not None:
                print('Displaying generated figure(s)...')
                plt.show()
                for fig in display_figures:
                    save_after_show = getattr(fig, '_alterseek_save_after_show', None)
                    if save_after_show is not None:
                        save_after_show()
            return True

        if not os.path.exists(struct_file):
            print(f"[Error] Structure file '{struct_file}' not found. Aborting "
                  "(not falling back to a possibly stale spin_flip_operations.txt).")
            return False
        elif not FIND_SF_AVAILABLE:
            print("[Note] find_sf_operations.py not found. Skipping Step 0.")
            struct_file = None
            centroid_struct_file = None
        else:
            is_mcif = struct_file.lower().endswith('.mcif')
            spin_axis_cart = None
            sf_result = None
            if is_mcif:
                print("Detected .mcif file --magnetic moments will be read from file.")
                moments_str = ""
            else:
                spin_axis_cart = _ask(
                    "Spin axis in Cartesian coordinates (default: 0 0 1): ", "spin_axis"
                )
                moments_str = _ask(
                    "Magnetic moments along this axis (atom order, trailing atoms auto-fill to 0): ",
                    "moments",
                )
            if not is_mcif and not moments_str:
                standard_path_reason = "No magnetic moments entered."
                standard_path_reason_reported = True
                _step0_wrote_flip_file = False
                print(f"{BOLD}[Note] {standard_path_reason}{RESET} Ordinary structural path will be written.")
            else:
                try:
                    sf_result = find_sf_run(
                        struct_file,
                        moments_str,
                        verbose=False,
                        spin_axis_cart=spin_axis_cart,
                        symprec=symprec,
                        output_dir=OUTPUT_DIR,
                    )
                except Exception as e:
                    print(f"[Error] Spin-symmetry analysis failed: {e} Aborting.")
                    return False
                if not isinstance(sf_result, dict):
                    print("[Error] Spin-symmetry analysis failed. Aborting.")
                    return False
            if isinstance(sf_result, dict):
                magnetic_setting_counts = None
                magnetic_setting_outputs = None
                # Symmetry of the cell the path is actually built in. Stays None
                # in the ordinary path, where sf_result already describes the
                # submitted cell.
                working_cell_symmetry = None
                # The magnetic primitive cell is adopted only when the submitted
                # cell genuinely cannot carry the path: a changed space group or
                # an enlarged cell (see _magnetic_cell_needed). Otherwise the
                # ordinary route is kept, because re-deriving the geometry from
                # FindSpinGroup's cell only re-orients it -- its axes can come
                # back permuted and sign-flipped, which silently rewrites the
                # path as a symmetry-equivalent but different set of
                # coordinates while the user's own cell stops matching their
                # KPOINTS.
                magnetic_cell_needed = _magnetic_cell_needed(sf_result)
                if self.magnetic_setting and magnetic_cell_needed:
                    try:
                        mag_setting = prepare_magnetic_setting_files(
                            struct_file,
                            moments_str=moments_str,
                            spin_axis_cart=spin_axis_cart,
                            output_dir=OUTPUT_DIR,
                        )
                        centroid_struct_file = mag_setting["helper_path"]
                        centroid_seekpath_type_numbers = mag_setting["seekpath_type_numbers"]
                        magnetic_setting_counts = mag_setting
                        operation_basis_label = mag_setting["operation_basis_label"]
                        # The altermagnetism gate below must judge the cell the
                        # path is built in, and must judge it by its *magnetic*
                        # symmetry. Both matter for MnSe2 (MAGNDATA 1.0.47): it
                        # is submitted in a cubic Pa-3 parent whose Laue group
                        # m-3 forbids altermagnetism, and the magnetic cell that
                        # replaces it is still the cubic Pa-3 *crystal* once the
                        # moments are stripped -- it is orthorhombic only because
                        # the moments make it so. Re-detecting symmetry from its
                        # coordinates therefore answers the wrong question, and
                        # answers it tolerance-dependently. G0, the spatial part
                        # of the spin space group, is the symmetry that actually
                        # holds, and FindSpinGroup reports it directly.
                        working_cell_symmetry = _g0_symmetry(
                            sf_result, sites=mag_setting.get('magnetic_cell_sites'))
                    except Exception as e:
                        print(f"[Error] Magnetic primitive cell construction failed: {e}")
                        print(
                            "[Error] The default magnetic-state path cannot be "
                            "generated. Fix the input/dependency problem, or "
                            "explicitly rerun with --parent-setting to request "
                            "the nonmagnetic reference path. Aborting."
                        )
                        return False
                if CENTROID_AVAILABLE:
                    try:
                        centroid_result = compute_centroid(
                            centroid_struct_file, output_dir=OUTPUT_DIR, show_plot=True,
                            defer_show=True, verbose=False,
                            seekpath_type_numbers=centroid_seekpath_type_numbers,
                            mode_2d=self.mode_2d,
                            input_vacuum_axis=self.input_vacuum_axis,
                            view_elev=view_elev, view_azim=view_azim, symprec=symprec,
                            figure_basename=_figure_basename(struct_file),
                            save_pdf=save_pdf,
                        )
                        if self.magnetic_setting and magnetic_setting_counts is not None:
                            magnetic_setting_outputs = finalize_magnetic_setting_outputs(
                                magnetic_setting_counts,
                                centroid_result,
                                output_dir=OUTPUT_DIR,
                                verbose_output=self.output_verbose,
                            )
                            if magnetic_setting_outputs:
                                centroid_result["b_matrix_output"] = magnetic_setting_outputs[
                                    "b_matrix_output"
                                ]
                                # Say something only when the user has to act.
                                # The analysis running in the magnetic cell is
                                # not itself news; a changed calculation cell
                                # is, because their own POSCAR no longer matches
                                # the path.
                                if magnetic_setting_outputs.get("cell_changed"):
                                    calc_cell = magnetic_setting_outputs.get(
                                        "calculation_cell_path")
                                    calc_magmom = magnetic_setting_outputs.get(
                                        "calculation_magmom_path")
                                    species_order = magnetic_setting_outputs.get(
                                        "calculation_species_order")
                                    print(
                                        "[Cell] The magnetic order changes the cell; the path "
                                        "is written in the magnetic primitive cell's basis."
                                    )
                                    if calc_cell:
                                        print(
                                            f"[Cell] Run the band calculation with {calc_cell}"
                                        )
                                    if calc_magmom:
                                        print(
                                            f"[Cell] Matching magnetic moments: {calc_magmom}"
                                        )
                                    if species_order:
                                        print(
                                            "[Cell] Species order: "
                                            f"{' '.join(species_order)} "
                                            "(match POTCAR and species-indexed settings)."
                                        )
                                if self.output_verbose:
                                    print(
                                        "[SSG setting] Kept intermediates in "
                                        f"{magnetic_setting_outputs.get('intermediate_dir')}"
                                    )
                        if magnetic_setting_counts is None:
                            # The submitted cell already carries the magnetic
                            # order, but SeeK-path may still permute/rotate it
                            # into its standard setting. Write a spin-bearing
                            # MCIF beside the structural standardized VASP so
                            # that setting can be inspected directly.
                            try:
                                standard_vasp = centroid_result.get(
                                    "standardized_structure_path"
                                )
                                if standard_vasp:
                                    (magnetic_lattice,
                                     magnetic_positions,
                                     magnetic_elements,
                                     magnetic_moments,
                                     _) = _load_magnetic_input_data(
                                        struct_file,
                                        moments_str,
                                        spin_axis_cart,
                                    )
                                    standard_mcif = os.path.splitext(
                                        standard_vasp
                                    )[0] + ".mcif"
                                    _write_seekpath_standard_mcif(
                                        standard_vasp,
                                        standard_mcif,
                                        f"{_figure_basename(struct_file)}_seekpath_standard",
                                        magnetic_lattice,
                                        magnetic_positions,
                                        magnetic_elements,
                                        magnetic_moments,
                                        symprec=centroid_result.get("symprec"),
                                    )
                                    centroid_result[
                                        "standardized_magnetic_structure_path"
                                    ] = standard_mcif
                            except Exception as exc:
                                print(
                                    "[Warning] Could not write the SeeK-path-"
                                    f"standardized magnetic MCIF: {exc}"
                                )
                        display_figures.extend(centroid_result.get('display_figures', []))
                    except Exception as e:
                        centroid_error = e
                if magnetic_setting_counts is not None:
                    _step0_wrote_flip_file = (
                        magnetic_setting_counts.get('spin_flip_operations', 0) > 0
                    )
                else:
                    _step0_wrote_flip_file = sf_result.get('spin_flip_operations', 0) > 0
                gate_laue_group = (working_cell_symmetry or sf_result).get('laue_group')
                laue_no_altermag = _altermagnetism_gate(sf_result, working_cell_symmetry)
                spin_split_diagnostic = sf_result.get('spin_split_diagnostic', '')

                parent_recovery = None
                lattice_tag = None
                if centroid_result is not None:
                    parent_recovery = centroid_result.get("mcif_parent_recovery")
                    lattice_tag = centroid_result.get(
                        'sc_type', centroid_result.get('seekpath_bravais', 'unknown'))
                # The two cells are printed adjacently, in the same field order,
                # because the point of the block is whether they differ. The
                # lattice tag and site count ride on whichever cell the k-path
                # is actually built in.
                print(f"\nInput structure: {sf_result['structure_file']}, "
                      f"{sf_result['num_atoms']} atoms")
                cell_rows = [(
                    'Nonmagnetic primitive cell:',
                    sf_result['space_group'],
                    sf_result['point_group'],
                    sf_result['laue_group'],
                    _cell_suffix(sf_result.get('nonmagnetic_sites'),
                                 sf_result.get('nonmagnetic_lattice')),
                )]
                reported_g0_symmetry = working_cell_symmetry
                reported_g0_lattice = lattice_tag
                if reported_g0_symmetry is None:
                    # G0 is still a result when no separate magnetic-cell
                    # construction is needed. In that case the submitted cell
                    # is already primitive for the magnetic order, so its site
                    # count and lattice class describe both rows. If a distinct
                    # magnetic cell exists but was not constructed (for example
                    # under --parent-setting), report G0 without borrowing the
                    # parent cell's size or lattice label.
                    reported_g0_symmetry = _g0_symmetry(
                        sf_result,
                        sites=(None if magnetic_cell_needed else
                               sf_result.get('nonmagnetic_sites')),
                    )
                    if magnetic_cell_needed:
                        reported_g0_lattice = None
                if reported_g0_symmetry is not None:
                    cell_rows.append((
                        'Magnetic primitive cell (G0):',
                        reported_g0_symmetry['label'],
                        reported_g0_symmetry['point_group'],
                        reported_g0_symmetry['laue_group'],
                        _cell_suffix(
                            reported_g0_symmetry.get('sites'),
                            reported_g0_lattice,
                        ),
                    ))
                recovery_note = None
                if parent_recovery:
                    recovery_note = (
                        "recovered from the input cell at symprec="
                        f"{parent_recovery['symprec']:g} (index {parent_recovery['index']})")
                _print_cell_rows(cell_rows, note_after_first=recovery_note)
                print(f"Phase: {sf_result['magnetic_phase']}")
                print(f"Oriented SSG: {sf_result['ssg_index']}")
                print(f"SSG Symbol (Chen-Liu): {sf_result['ssg_symbol']}")
                print(f"MSG without SOC: {sf_result['magnetic_space_group_without_soc']}")

                if laue_no_altermag:
                    laue = laue_no_altermag.get('laue_group', gate_laue_group)
                    standard_path_reason = f"Laue group {laue}: no altermagnetism."
                    print(f"{BOLD}[Note] {standard_path_reason}{RESET} Default path will be written.")
                    standard_path_reason_reported = True
                else:
                    print("Spin operations: "
                          f"{sf_result['actual_spin_flip_point_operations']} flip, "
                          f"{sf_result['actual_spin_preserve_point_operations']} preserve")
                    if spin_split_diagnostic:
                        standard_path_reason = spin_split_diagnostic
                        print(f"{BOLD}[Note] {standard_path_reason}{RESET} Default path will be written.")
                        standard_path_reason_reported = True

                if STEP0_VERBOSE_SUMMARY:
                    print(f"Magnetic SG: {sf_result['magnetic_space_group']}")
                    print(f"G0: {sf_result['g0_symbol']} ({sf_result['g0_number']}), "
                          f"L0: {sf_result['l0_symbol']} ({sf_result['l0_number']}), "
                          f"EMPG: {sf_result['empg']}")
                    if sf_result.get('findspingroup_warning'):
                        print(f"FindSpinGroup warning: {sf_result['findspingroup_warning']}")
                    print(f"Spin axis: {sf_result['spin_group']}")
                    unique_ops = sf_result.get('unique_point_operations',
                                               sf_result['total_operations'])
                    print(f"Space-group operations: {sf_result['total_operations']} total")
                    print(f"Point operations: {unique_ops} unique")
                    print("Inversion-extended k operations: "
                          f"{sf_result['extended_spin_flip_point_operations']} spin-flip, "
                          f"{sf_result['extended_spin_preserve_point_operations']} spin-preserving "
                          f"({sf_result['extended_spin_flip_operations']} + "
                          f"{sf_result['extended_spin_preserve_operations']} with translations)")
                    print(f"Saved: {', '.join(sf_result['saved_files'])}")

        if centroid_struct_file and CENTROID_AVAILABLE:
            try:
                if centroid_result is None:
                    if centroid_error is not None:
                        raise centroid_error
                    centroid_result = compute_centroid(
                        centroid_struct_file, output_dir=OUTPUT_DIR, show_plot=True,
                        defer_show=True, verbose=False,
                        seekpath_type_numbers=centroid_seekpath_type_numbers,
                        mode_2d=self.mode_2d,
                        input_vacuum_axis=self.input_vacuum_axis,
                        view_elev=view_elev, view_azim=view_azim, symprec=symprec,
                        figure_basename=_figure_basename(struct_file),
                        save_pdf=save_pdf,
                    )
                    display_figures.extend(centroid_result.get('display_figures', []))
                    print(
                        "Lattice type: "
                        f"{centroid_result.get('sc_type', centroid_result.get('seekpath_bravais', 'unknown'))}"
                    )
                print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
                sp_path   = centroid_result['sp_path']
                sp_coords = centroid_result['sp_point_coords']
                displayed_path = centroid_result.get(
                    'band_kpath',
                    centroid_result.get('ibz_kpath', sp_path)
                )
                print(f"Path: {self._format_path(displayed_path)}")
                path_choice = _ask(
                    "Press [Enter] to use this path, or type a filename to load your own: ",
                    "path",
                )
                if not path_choice:
                    # Build kpoints_data in the same HPKOT/SeeK-path convention
                    # as Figure 1.  lattice_kpoints.py may include curated
                    # closure vertices for the hull, but ibz_kpath contains only
                    # the public band-path labels.  path_kpoints_frac keeps
                    # optional path-only labels such as H_2 available without
                    # adding them to the centroid hull.
                    self.kpoints_data = []
                    sc_type_auto = centroid_result.get('sc_type', '')
                    if (
                        ('band_kpath' in centroid_result and 'band_kpoints_frac' in centroid_result)
                        or ('ibz_kpath' in centroid_result and 'ibz_kpoints_frac' in centroid_result)
                    ):
                        self.header_lines = [f'K-Path generated by AlterSeeK-Path (HPKOT {sc_type_auto})',
                                             '20', 'Line-Mode', 'Reciprocal']
                        self.kpoints_basis_matrix = np.array(centroid_result['b_matrix'], dtype=float)
                        self.kpoints_basis_rotation = np.array(
                            centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                            dtype=float,
                        )
                        self.output_basis_matrix = np.array(
                            centroid_result.get(
                                'b_matrix_output',
                                centroid_result.get('b_matrix_input', centroid_result['b_matrix']),
                            ),
                            dtype=float,
                        )
                        # Prefer the selected band path when present. This
                        # keeps the prompt, Figure 1 path overlay, and KPOINTS
                        # path consistent.
                        auto_path = centroid_result.get(
                            'band_kpath',
                            centroid_result['ibz_kpath']
                        )
                        ibz_coords = centroid_result.get(
                            'band_kpoints_frac',
                            centroid_result.get(
                                'path_kpoints_frac',
                                centroid_result['ibz_kpoints_frac']
                            )
                        )
                        for seg_start, seg_end in auto_path:
                            for label in (seg_start, seg_end):
                                coords = ibz_coords[label]
                                self.kpoints_data.append([coords[0], coords[1], coords[2], label])
                        extra_vertices = centroid_result.get('extra_general_vertices', [])
                        self.extra_general_points = []
                        for label in extra_vertices:
                            if label in ibz_coords:
                                coords = ibz_coords[label]
                                self.extra_general_points.append([coords[0], coords[1], coords[2], label])
                        print(f"Using HPKOT {sc_type_auto} path ({len(auto_path)} segments, {len(self.kpoints_data)} k-points)")
                        if self.extra_general_points:
                            labels = ", ".join(str(pt[3]) for pt in self.extra_general_points)
                            print(f"Extra doubled-IBZ general anchors: {labels}")
                    else:
                        self.header_lines = ['K-Path generated by AlterSeeK-Path (seekpath)', '20', 'Line-Mode', 'Reciprocal']
                        self.kpoints_basis_matrix = np.array(centroid_result['b_matrix'], dtype=float)
                        self.kpoints_basis_rotation = np.array(
                            centroid_result.get('seekpath_rotation_matrix', np.eye(3)),
                            dtype=float,
                        )
                        self.output_basis_matrix = np.array(
                            centroid_result.get(
                                'b_matrix_output',
                                centroid_result.get('b_matrix_input', centroid_result['b_matrix']),
                            ),
                            dtype=float,
                        )
                        for seg_start, seg_end in sp_path:
                            for label in (seg_start, seg_end):
                                coords = sp_coords[label]
                                self.kpoints_data.append([coords[0], coords[1], coords[2], label])
                        print(f"Using auto-generated path ({len(sp_path)} segments, {len(self.kpoints_data)} k-points)")
                else:
                    if not _load_custom_path(path_choice):
                        return False
            except Exception as e:
                print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
                print(f"[Warning] Auto path generation failed: {e}")
                print("Falling back to manual file input.")
                print("Enter KPOINTS file name (default: KPATH.in): ", end='', flush=True)
                filename = input().strip()
                if not filename: filename = "KPATH.in"
                if not _load_custom_path(filename):
                    return False
        else:
            print(f"\n{BOLD}>>> Step 1: High-symmetry k-path{RESET}")
            print("Enter KPOINTS file name (default: KPATH.in): ", end='', flush=True)
            filename = input().strip()
            if not filename: filename = "KPATH.in"
            if not _load_custom_path(filename):
                return False

        # Laue groups -1, -3, and m-3 do not have a one-dimensional,
        # nonidentical inversion-even irrep, so no altermagnetic splitting is
        # possible. Write the ordinary IBZ path without butterfly insertion.
        no_altermag = None
        if centroid_result is not None:
            no_altermag = centroid_result.get('no_altermagnetism')
            if no_altermag is None and no_altermagnetism_reason is not None:
                no_altermag = no_altermagnetism_reason(
                    centroid_result.get('point_group'),
                    centroid_result.get('spacegroup'),
                )
        if standard_path_reason is None and no_altermag:
            laue = no_altermag.get('laue_group', 'unknown')
            standard_path_reason = f"Laue group {laue}: no altermagnetism."

        # Step 2: Auto-compute or manually enter general k-point
        print(f"\n{BOLD}>>> Step 2: General k-point{RESET}")
        general_kpoint = None

        if centroid_result is not None:
            try:
                c = centroid_result['centroid_frac']
                general_kpoint = [c[0], c[1], c[2]]
                print(f"IBZ centroid (standardized basis): [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
            except Exception as e:
                print(f"[Warning] Centroid retrieval failed: {e}")
        elif struct_file and CENTROID_AVAILABLE:
            try:
                result = compute_centroid(centroid_struct_file, output_dir=OUTPUT_DIR, show_plot=True,
                                          defer_show=True, verbose=False,
                                          seekpath_type_numbers=centroid_seekpath_type_numbers,
                                          mode_2d=self.mode_2d,
                                          input_vacuum_axis=self.input_vacuum_axis,
                                          view_elev=view_elev, view_azim=view_azim,
                                          symprec=symprec,
                                          figure_basename=_figure_basename(struct_file),
                                          save_pdf=save_pdf)
                display_figures.extend(result.get('display_figures', []))
                c = result['centroid_frac']
                general_kpoint = [c[0], c[1], c[2]]
                print(f"IBZ centroid (standardized basis): [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
            except Exception as e:
                print(f"[Warning] Centroid computation failed: {e}")

        # Append centroid to spin_operations.txt for reference
        if general_kpoint is not None:
            out_k = self._general_kpoint_output_basis(general_kpoint)
            if out_k is not None:
                print(
                    "IBZ centroid (KPOINTS output basis): "
                    f"[{out_k[0]:.6f}, {out_k[1]:.6f}, {out_k[2]:.6f}]"
                )
            try:
                with open(os.path.join(OUTPUT_DIR, "spin_operations.txt"),
                          "a", encoding="utf-8", newline="\n") as f:
                    f.write(f"\nGeneral k-point (IBZ centroid, standardized primitive basis): "
                            f"[{general_kpoint[0]:.6f}, {general_kpoint[1]:.6f}, {general_kpoint[2]:.6f}]\n")
                    if out_k is not None:
                        f.write(f"General k-point (IBZ centroid, KPOINTS output basis): "
                                f"[{out_k[0]:.6f}, {out_k[1]:.6f}, {out_k[2]:.6f}]\n")
            except Exception:
                pass

        if general_kpoint is None:
            print("Format: kx ky kz (space-separated)")
            while True:
                try:
                    print("Enter k-point: ", end='', flush=True)
                    k_input = input().strip().split()
                    if len(k_input) == 3:
                        general_kpoint = [float(x) for x in k_input]
                        break
                    else:
                        print("Please enter exactly 3 coordinates.")
                except ValueError:
                    print("Invalid input. Please enter three numbers.")

        def _write_ordinary_path_and_stop(reason, reason_reported):
            if not reason_reported:
                print(f"\n{BOLD}[Note] {reason}{RESET}")
            print("Writing ordinary k-path with non-spin-flip general-k anchors.")
            standard_general_path = self.insert_general_kpoint_anchors(
                general_kpoint,
                self.extra_general_points,
            )
            return _save_and_finish(standard_general_path, None, None)

        if standard_path_reason:
            return _write_ordinary_path_and_stop(
                standard_path_reason, standard_path_reason_reported
            )

        # Step 3: Input transformation matrix
        print(f"\n{BOLD}>>> Step 3: Spin-flip operation{RESET}")
        if _step0_wrote_flip_file is False:
            # Step 0 ran but found no flip ops for this structure.
            # Don't load a stale file from a previous run on a different structure.
            print("[Note] Step 0 found no spin-flip operations --skipping any existing spin_flip_operations.txt.")
            flip_ops = []
        else:
            flip_ops = self.load_flip_operations()
        preserve_ops = self.load_preserve_operations()

        # Always include inversion-extended spatial partners. The spin-flip
        # classification comes from the spin rotation in find_sf_operations.py;
        # multiplying the spatial operation by inversion does not change that
        # spin-flip status. Deduplicate after extension.
        def _inversion_extended(ops):
            expanded = list(ops)
            for op in ops:
                neg_op = -np.array(op, dtype=float)
                if not any(np.allclose(neg_op, ex, atol=1e-8) for ex in expanded):
                    expanded.append(neg_op)
            return expanded

        flip_ops = _inversion_extended(flip_ops)
        preserve_ops = _inversion_extended(preserve_ops)

        # 2D mode: keep only genuine in-plane spin-flip operations. A flip
        # operation whose in-plane 2x2 block is +-I (C2z/mz/inversion type)
        # produces no in-plane splitting, so it must not be chosen for the
        # 2D path. If none survive, the slab is not a 2D altermagnet.
        _flip_ops_emptied_2d = False
        if self.mode_2d and flip_ops:
            vax = self.input_vacuum_axis
            valid_flip_ops = [op for op in flip_ops
                              if is_valid_2d_spin_flip(op, vax)]
            n_excluded = len(flip_ops) - len(valid_flip_ops)
            if valid_flip_ops:
                print(f"[2D mode] In-plane spin splitting: YES "
                      f"({len(valid_flip_ops)} valid 2D spin-flip ops"
                      + (f", {n_excluded} trivial C2/m ops excluded" if n_excluded else "")
                      + ").")
                flip_ops = valid_flip_ops
            else:
                print("[2D mode] In-plane spin splitting: NO -- every spin-flip "
                      "operation acts as +-identity in-plane (C2z/mz/inversion "
                      "type). This slab is not a 2D altermagnet; writing the "
                      "ordinary in-plane path without a k' partner.")
                flip_ops = []
                _flip_ops_emptied_2d = True

        if not flip_ops:
            # No candidate spin-flip point operation survived (structural
            # symmetry search found none, or the 2D in-plane filter emptied
            # the list) -- this means "not altermagnetic" for this
            # configuration, same as the Laue-group/no-moments cases handled
            # above. Do not fall through to manual 3x3 matrix entry.
            return _write_ordinary_path_and_stop(
                "No spin-flip point operation available: not altermagnetic.",
                _flip_ops_emptied_2d,
            )

        preset_flip_choice = input_config.get('flip_option')
        if preset_flip_choice is not None and preset_flip_choice > len(flip_ops):
            print(
                f"[Error] flip_option={preset_flip_choice} is out of range; "
                f"available choices are 1-{len(flip_ops)}."
            )
            return False

        # Operation naming: the listed matrices are in the recorded operation-
        # source fractional basis, so the axis cannot be read from the integers
        # directly. Convert each through Cartesian (same path Figure 2 uses) and
        # classify there, so the name matches the figure labels.
        R, selected_transformation_label = self._select_spin_flip_operation(
            flip_ops,
            centroid_result,
            preset_choice=preset_flip_choice,
            operation_basis_label=operation_basis_label,
        )
        # Step 4: Process k-points
        print(f"\n{BOLD}>>> Step 4: Build altermagnetic path{RESET}")

        (R_for_kpts, R_cart_for_plot, flip_ops_for_plot,
         preserve_ops_for_plot) = self._convert_operation_to_primitive_basis(
            R,
            flip_ops,
            preserve_ops,
            centroid_result,
            operation_basis_label,
        )
        # Calculate and show k'
        k_prime = self.transform_kpoint(general_kpoint, R_for_kpts)
        print(f"k' = [{k_prime[0]:.4f}, {k_prime[1]:.4f}, {k_prime[2]:.4f}]")

        try:
            new_kpoints = self.insert_general_kpoints(
                general_kpoint, R_for_kpts, self.extra_general_points
            )

            if new_kpoints:
                self._generate_spin_figures(
                    centroid_result, struct_file, general_kpoint, R_for_kpts,
                    R_cart_for_plot, flip_ops_for_plot, preserve_ops_for_plot,
                    new_kpoints, display_figures, save_pdf=save_pdf)
                # Step 5: Save modified file
                return _save_and_finish(new_kpoints, R, selected_transformation_label)
            else:
                print("Error: Failed to process k-points.")
                return False
                
        except Exception as e:
            print(f"Error processing k-points: {e}")
            return False
