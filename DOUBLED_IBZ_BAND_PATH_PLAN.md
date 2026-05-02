# Doubled-IBZ Band Path Plan

Planning note started 2026-04-30. This file records the intended band-path
policy for doubled-IBZ cases before code changes are made.

## General Policy

- Keep BZ/IBZ plotting and centroid geometry separate from the band path.
- For band calculations, start from the standard SeeK-path/HPKOT path whenever
  possible so the usual high-symmetry comparison path is preserved.
- For doubled-IBZ cases, append only the copied-vertex general-point anchors
  needed to sample the enlarged wedge, rather than replacing the usual
  SeeK-path path with a fully custom path or duplicating equivalent
  high-symmetry edges.
- Do not remove the existing project-curated path implementation until the new
  policy is tested case by case.

## hP Doubled 60-Degree Wedge: Trigonal 32/3m/-3m And Hexagonal 6/-6/6/m

Current issue:

- The existing project-curated doubled hP path covers the doubled wedge, but it
  changes the usual hexagonal path by replacing the standard `K-Gamma` part
  with a route through the doubled-wedge vertex family, e.g. `K-M_A-Gamma`.
- This loses the direct standard `Gamma-K`/`K-Gamma` comparison segment.

Implemented direction:

- Use the usual SeeK-path hexagonal path as the base:

```text
Gamma-M-K-Gamma-A-L-H-A | L-M | H-K
```

- Append only isolated general-point anchors for the added doubled-wedge
  vertices. Do not append the copied edge `L_A-M_A`, because it duplicates the
  usual-path edge `L-M`.

```text
M_A-k | k'-M_A'
L_A-k | k'-L_A'
```

Implementation decision:

- `compute_centroid_hybrid.py` returns `extra_general_vertices`.
- `alterseek_path.py` appends each extra vertex as `B_A-k | k'-B_A'`.
- These anchors are not ordinary path segments and must not create copied
  high-symmetry edges such as `L_A'-M_A'`.

Implementation scope:

- Applies to hP doubled 60-degree wedges:
  - trigonal `32`, `3m`, `-3m`: SG `149-167`
  - hexagonal `6`, `-6`, `6/m`: SG `168-176`
- Does not yet change the separate `3`/`-3` quadrupled hP case SG `143-148`.

## 4/m Tetragonal Doubled Wedge

- Apply the same anchor-only policy:
  standard SeeK-path tetragonal path first, then append isolated
  `B_A-k | k'-B_A'` anchors for the project-only copied vertices.
- Implemented copied vertices:
  - `tP1`: `X_A`, `R_A`
  - `tI1`: `M_p`, `N_p`, `Z_0_p`
  - `tI2`: `N_p`, `S_0_p`, `S_p`
- Do not append ordinary high-symmetry segments for these vertices unless a
  separate paper-path decision explicitly needs the duplicated edge.
