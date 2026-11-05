# Steel plate response cases

Six standalone linear-static analyses exercise the working-tree
`python_scripts/helper/responses.py`. No test-suite registration is required.

| Case | Analysis | Structured mesh | Element |
| --- | --- | --- | --- |
| `plane_stress_q4` | Plane stress | 10 × 4 | SmallDisplacementElement2D4N |
| `plane_stress_q8` | Plane stress | 10 × 4 | SmallDisplacementElement2D8N |
| `plane_strain_q4` | Plane strain | 10 × 4 | SmallDisplacementElement2D4N |
| `plane_strain_q8` | Plane strain | 10 × 4 | SmallDisplacementElement2D8N |
| `solid_h8` | 3D solid | 10 × 4 × 2 | SmallDisplacementElement3D8N |
| `solid_h20` | 3D solid | 10 × 4 × 2 | SmallDisplacementElement3D20N |

Each case contains `ProjectParameters.json`, `Materials.json`, and `plate.mdpa`.
The shared `Main.py` resolves their paths, so it can be launched from any folder.
Use the Python environment containing your Kratos build; NumPy must be available.
The solver is the core `skyline_lu_factorization` direct solver.

## Geometry, steel, and loading

Units are mm, N, MPa, and tonne. The plate spans x = 0–100 mm and y = 0–20 mm.
2D cases have 2 mm thickness; 3D cases span z = 0–2 mm.
Steel uses E = 210000 MPa, nu = 0.3, and density = 7.85e-9 tonne/mm³.
The small-displacement linear-elastic model has no initial stress or strain.

The x = 0 edge/face is fully fixed in the active displacement directions.
The x = 100 edge/face carries uniform tension in +x; other boundaries are free.
2D line load is 200 N/mm (100 MPa × 2 mm thickness). 3D surface load is 100 MPa.
Both produce a 4000 N total force. Quadratic boundary conditions use matching
quadratic line/surface load conditions, including their midside nodes.

All cases use one load step, `analysis_type: linear`, and `move_mesh_flag: false`.
The 3D model can contract through its thickness away from the fixed face; it is
therefore closer to the plane-stress case than the plane-strain case.

## Run

From this folder:

```bash
python3 Main.py
python3 Main.py --case plane_stress_q4
python3 Main.py --case plane_strain_q8
python3 Main.py --case solid_h20 --mode center_gp
```

`Main.py` loads the local source helper explicitly, so rebuilding Kratos is not
needed after editing `responses.py`. It runs all six cases by default and uses
`average`, meaning the arithmetic mean of the Gauss-point von Mises values.

Default first-order integration has 2 × 2 points in 2D and 2 × 2 × 2 in 3D,
with no center Gauss point. Default second-order integration has 3 × 3 points
or 3 × 3 × 3 points and includes the center. Thus `--mode center_gp` is expected
to warn and raise `NotImplementedError` for the first-order cases, because its
`center_element` fallback is intentionally unfinished. Use an individual
second-order case to exercise the successful center-GP path. Selecting
`--mode center_element` also raises `NotImplementedError` as expected.

## Inspect results

After a successful run, each case has a `results` folder containing:

- `von_mises_<mode>.csv`: element IDs, centers, GP counts, and extracted stresses.
- `nodal_results.csv`: coordinates, displacement, and reaction components.
- `summary_<mode>.json`: stress range, downstream stress mean, loaded-edge mean
  displacement, and total support reaction.

The runner checks that the support reaction balances -4000 N in x and raises
an error otherwise. Away from the clamp, ideal uniform-tension von Mises stress
is 100 MPa for plane stress and a freely contracting 3D solid. For plane strain,
sigma_yy approaches zero and sigma_zz approaches nu × sigma_xx, giving
100 × sqrt(1 - nu + nu²) = 88.8819 MPa. These are far-field references, not
exact values for every element: the fully fixed boundary creates transverse
stresses and corner concentrations. Compare the downstream mean and first- vs
second-order results rather than requiring the maximum to equal the reference.

The analyses are intended for manual execution. Input parsing, mesh connectivity,
and analysis initialization were checked during setup; the load steps were not run.
