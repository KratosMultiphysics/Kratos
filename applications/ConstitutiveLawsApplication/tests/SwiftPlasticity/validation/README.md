# Swift plasticity: independent one-element FE validation

This standalone benchmark compares an actual Kratos finite-element solution
using `SmallStrainSwiftJ2Plasticity3DLaw` with an independent analytical
uniaxial Swift solution. It exercises assembly, nonlinear equilibrium,
constitutive integration, history commitment and integration-point extraction.
It requires no external solver or reference data.

## FE model and material

Units are **N, mm, MPa**; strains are dimensionless. The model is one
1 mm cube with `SmallDisplacementElement3D8N`, using `Hexahedra3D8` and
`GI_GAUSS_2` (2 x 2 x 2 Gauss points). Connectivity is `[1,2,3,4,5,6,7,8]`:

```text
node:  1       2       3       4       5       6       7       8
xyz: (0,0,0) (1,0,0) (1,1,0) (0,1,0) (0,0,1) (1,0,1) (1,1,1) (0,1,1)
```

Symmetry conditions fix `ux` on `x=0`, `uy` on `y=0`, and `uz` on `z=0`.
The `x=1` face receives uniform prescribed `ux` from 0 to 0.02 mm in
**100 equal increments**. Transverse contraction remains free, producing
uniaxial stress. Axial strain is `ux / 1 mm`, ending at 0.02; the FE load
levels are used directly without resampling. Pseudo-time introduces no rate
effect, and the reference mesh stays fixed.

Both the FE runner and analytical comparison read
[`kratos/Materials.json`](kratos/Materials.json):

| Parameter | Value |
|---|---:|
| `YOUNG_MODULUS` | 210000 MPa |
| `POISSON_RATIO` | 0.3 |
| `SWIFT_COEFFICIENT` | 1000 MPa |
| `SWIFT_INITIAL_STRAIN` | 0.01 |
| `SWIFT_HARDENING_EXPONENT` | 0.2 |

The material uses infinitesimal strain, isotropic elasticity and associated
J2 plasticity with continuous Swift hardening `sigma_y = K*(epsilon_0+p)^n`.
Density is zero because there are no body or inertial loads. See the
[detailed law documentation](../README.md) for its formulation and API.

The solver settings are unchanged from the development benchmark:

| Setting | Value |
|---|---|
| Strategy | `ResidualBasedNewtonRaphsonStrategy` |
| Scheme | `ResidualBasedIncrementalUpdateStaticScheme` |
| Linear solver | `SkylineLUFactorizationSolver` |
| Builder and solver | `ResidualBasedBlockBuilderAndSolver` |
| Criterion | `ResidualCriteria(1e-10, 1e-12)` |
| Maximum iterations | 30 per increment |
| Reactions | Calculated |
| Reform DOFs / move mesh | False / False |

The residual criterion accepts a relative norm ratio <= 1e-10 **or** an
absolute norm < 1e-12 N. Its absolute norm is the residual L2 norm divided
by the number of considered DOFs. Each increment must converge before history
is finalized. The runner extracts `CAUCHY_STRESS_VECTOR` and committed
`ACCUMULATED_PLASTIC_STRAIN` at all eight integration points. Equivalent
stress is computed from the full 3D von Mises invariant. Before averaging,
point spreads must be <= 1e-9 MPa for stresses and <= 1e-14 for `p`.

## Independent analytical reference

For monotonic uniaxial tension, associated J2 flow gives `epsilon_p_xx = p`:

```text
sigma_y0 = K*epsilon_0^n
elastic: sigma = E*epsilon, p = 0       (E*epsilon <= sigma_y0)
plastic: F(p) = E*(epsilon-p) - K*(epsilon_0+p)^n = 0
         sigma = E*(epsilon-p), seqv = abs(sigma)
```

`compare_swift_validation.py` imports no Kratos modules. It bisects the unique
plastic root on `[0, epsilon]` because `F(0)>0`, `F(epsilon)<0` and
`F'(p)=-E-H(p)<0`. Up to 100 bisections, or floating-point bracket stagnation,
are followed by an independent yield-equation check within 1e-9 MPa.
This does not reproduce the production incremental return mapping.

## Execution and outputs

With a configured Kratos Python environment, run from the repository root:

```bash
OMP_NUM_THREADS=1 python3 -B applications/ConstitutiveLawsApplication/tests/SwiftPlasticity/validation/run_swift_validation.py --output-dir /tmp/kratos_swift_validation
```

This solves the FE problem and runs the analytical comparison. Omit
`--output-dir` to use `kratos_swift_validation` under Python's system temporary
directory (`/tmp` on this Linux environment). Relative output paths are resolved
from the working directory. Reruns replace generated outputs in that directory;
use separate directories for concurrent runs. All generated files are disposable
and are not version-controlled. Default execution writes no source-tree outputs.

For this repository's existing Release build, the exact invocation is:

```bash
OMP_NUM_THREADS=1 \
PYTHONPATH=build/Release/kratos:build/Release/applications/StructuralMechanicsApplication:build/Release/applications/ConstitutiveLawsApplication:bin/Release \
LD_LIBRARY_PATH=build/Release/kratos:build/Release/applications/StructuralMechanicsApplication:build/Release/applications/ConstitutiveLawsApplication:bin/Release/libs \
/usr/bin/python3 -B applications/ConstitutiveLawsApplication/tests/SwiftPlasticity/validation/run_swift_validation.py --output-dir /tmp/kratos_swift_validation
```

Use the Python interpreter matching the build. Native build-library paths above
precede the installed package tree to avoid loading stale installed libraries.
To repeat only the comparison (no Kratos Python environment needed):

```bash
python3 -B applications/ConstitutiveLawsApplication/tests/SwiftPlasticity/validation/compare_swift_validation.py --output-dir /tmp/kratos_swift_validation
```

| Generated file in the output directory | Contents |
|---|---|
| `swift_validation_kratos.csv` | 100 FE integration-point averages |
| `swift_validation_integration_points.csv` | All eight point states, 800 rows |
| `swift_validation_diagnostics.csv` | Iterations, residuals and point spreads |
| `swift_validation_comparison.csv` | Analytical values and per-step FE errors |
| `swift_validation_summary.json` | Aggregate errors, spreads and final values |

The raw FE columns are
`substep,eps_xx,sigma_xx,sigma_yy,sigma_zz,seqv,p_eq,ux`.
The comparison preserves the raw result file, reports maximum and RMS absolute
errors for axial/equivalent stress and `p`, and recomputes point spreads from
all 800 rows. Relative errors use analytical magnitudes as denominators, only
when >= 1 MPa for stress or >= 1e-8 for `p`; otherwise entries are blank/null.
Relative errors are fractions, not percentages.

If matplotlib is already installed, two overlays are also written:
`swift_validation_stress_strain.png` and `swift_validation_plastic_strain.png`.
No plotting dependency is installed by these scripts.

## Interpretation and limitations

Expect a homogeneous uniaxial-stress state and agreement with continuous Swift
to constitutive-integration and FE-equilibrium precision. The generated summary
reports actual errors; agreement thresholds are not fitted to observed results.
The cleanup verification run converged in all 100 increments with at most two
Newton iterations. Maximum/RMS axial-stress errors were 3.8432e-10/9.9861e-11
MPa; maximum `p` error was 2.0047e-15. Maximum transverse stress over all points
was 1.4106e-10 MPa, and maximum stress/`p` spreads were 1.9327e-12 MPa/1.0408e-17.
These are observed floating-point results, not portable acceptance limits.

The benchmark covers one monotonic homogeneous loading path, not general
multiaxial loading, unloading, mesh convergence, finite strain or damage.
The focused C++ material-point tests provide complementary coverage.
