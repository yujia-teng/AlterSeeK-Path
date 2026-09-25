"""Read AlterSeeK-Path results into Python instead of parsing its output files.

Runs the scripted path workflow on several structures in one process and
collects the results in a table.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # figures are saved to files; no windows open

from alterseek import general_k_path

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
OUTPUT = HERE / "output"

# name, structure file, general_k_path settings
STRUCTURES = [
    ("MnF2", EXAMPLES / "MnF2" / "VASP" / "POSCAR",
     dict(moments="5 -5 4*0", flip_option=2)),
    ("GdAuGe", EXAMPLES / "GdAuGe" / "VASP" / "POSCAR",
     dict(moments="8 -8 4*0", flip_option=11)),
    # MCIF: moments and spin axis are read from the file
    ("La2NiO4", HERE / "structures" / "La2NiO4.mcif", dict()),
    # 2D slab with vacuum along c
    ("square-4m", HERE / "structures" / "square_4m_POSCAR",
     dict(moments="1 -1 -1 1 6*0", mode_2d=True, vacuum_axis="c")),
]

rows = []
for name, path, settings in STRUCTURES:
    kpath = general_k_path(
        str(path), output_dir=str(OUTPUT / name), **settings
    )
    if kpath["no_splitting_reason"]:
        print(f"{name}: ordinary path -- {kpath['no_splitting_reason']}")
        continue

    rows.append({
        "name": name,
        "phase": kpath["magnetic_phase"].replace("\n", " "),
        "msg": kpath["magnetic_space_group_without_soc"],
        "flips": kpath["flip_option_count"],
        "lattice": kpath["lattice"],
        "kpath": kpath,
    })

header = f"{'structure':10s} {'phase':22s} {'MSG without SOC':36s} {'#R':>3s} {'lat':6s} general k"
print(header)
print("-" * len(header))
for r in rows:
    k = r["kpath"]["k"]
    print(f"{r['name']:10s} {r['phase']:22s} {r['msg']:36s} "
          f"{r['flips']:3d} {r['lattice']:6s} "
          + " ".join(f"{x:8.5f}" for x in k))

for r in rows:
    kpath = r["kpath"]
    k, kp = kpath["k"], kpath["k_prime"]
    print()
    print(f"{r['name']}: {kpath['path']}")
    print(f"  k  = [{k[0]:.4f}, {k[1]:.4f}, {k[2]:.4f}]")
    print(f"  k' = [{kp[0]:.4f}, {kp[1]:.4f}, {kp[2]:.4f}]")
    print(f"  {len(kpath['segments'])} segments")
