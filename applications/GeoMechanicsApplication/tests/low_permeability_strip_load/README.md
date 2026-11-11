# Strip load on low-permeability soil

This benchmark checks that coupled U-Pw analyses with a very low permeability are solved
accurately. It is run by `test_low_permeability_strip_load.py`.

## Model

- A 2D plane-strain domain, 40 m deep, consisting of three linear elastic layers.
- The soil is modelled with fully coupled U-Pw elements, with a permeability of 1e-16 m².
- **Stage `1_Initial`:** gravity loading with the K0 procedure.
- **Stage `2_Load`:** a strip load of 50 kN/m². Displacements are reset at the start of this
  stage, so its results are caused by the strip load alone.
- Both stages use the static scheme (`"solution_type": "static"`).
- Hydraulic boundaries: water pressure is 0 at the top surface, and all other boundaries are
  impermeable.

## Why it is a useful benchmark

The force equations have entries of about 1e7, set by the soil stiffness. The flow equations
scale with the permeability and have entries of about 1e-13. An unscaled sparse LU solve loses the
flow equations to round-off. The residual-based convergence check is dominated by the force
residuals, so it still reports convergence, but the water pressures and displacements are wrong.

| Solve | Max settlement in stage `2_Load` |
|---|---|
| Unscaled (incorrect) | 0.1664 m |
| Scaled (correct) | 0.1749 m |

## Expected results

The static flow problem has an exact solution, so the benchmark doesn't need an external reference:

- **Water pressure** must be hydrostatic, `p = 10000 · y` [Pa], in both stages.
- **Displacements** must match the drained reference, in which the hydrostatic water pressure is
  prescribed at every node. With the pressures prescribed, only the force equations are solved,
  so the result doesn't depend on scaling.

The `<stage>__expected_results.csv` files contain the nodal total displacements of that drained
reference.

## Variants

| Test | Strategy | Scaling |
|---|---|---|
| `test_drained_reference` | Newton-Raphson, water pressure prescribed | off |
| `test_newton_raphson_with_scaling` | Newton-Raphson | `"scaling": true` in `linear_solver_settings` |
| `test_quasi_newton_lbfgs_with_scaling` | quasi-Newton (L-BFGS) | `"quasi_newton_scaling": true` |
| `test_quasi_newton_broyden_with_scaling` | quasi-Newton (Broyden) | `"quasi_newton_scaling": true` |
