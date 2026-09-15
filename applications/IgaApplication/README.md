## IGA Application
### Isogeometric beam mass formulations

`IsogeometricBeamElement::CalculateMassMatrix` requires an explicit property:

```python
properties.SetValue(IGA.BEAM_MASS_FORMULATION, "simplified")
# Or: properties.SetValue(IGA.BEAM_MASS_FORMULATION, "full")
```

In a materials JSON file, set `"BEAM_MASS_FORMULATION": "simplified"` or
`"BEAM_MASS_FORMULATION": "full"` under `Variables`. There is no implicit default.
Static calculations do not require this property.

Both options use consistent NURBS integration over the reference geometry, with
four DOFs per control point ordered as `[u_x, u_y, u_z, phi]`. They require positive,
finite `DENSITY`, `CROSS_AREA`, `I_N`, and `I_V`. The section is assumed to have its
mass centroid on the centerline and principal axes along its reference directors.
`I_V` weights offsets along N, and `I_N` weights offsets along V. Torsional mass
uses the polar area moment `I_N + I_V`, **not** the stiffness torsion constant `I_T`.

* `simplified`: integrates `rho*N_i*N_j*diag(A,A,A,I_N+I_V)` over reference arc
  length. It is independent of displacement and scalar twist, and neglects
  bending rotary inertia and director-motion coupling.
* `full`: integrates `rho*(A*B_r^T*B_r + I_V*B_n^T*B_n + I_N*B_v^T*B_v)`, where
  the B matrices differentiate the current centerline and directors with respect
  to all control-point DOFs. It includes bending and torsional rotary inertia and
  translation/twist coupling. The director derivatives are analytic. Initial
  coordinates plus historical displacement determine the current state; mesh
  coordinates are not used. Opposite reference/current tangents are singular and
  rejected, as in the displacement mapper.

Each element represents a single quadrature-point geometry. Assemble contributions
from all quadrature-point elements of the curve. Choose quadrature appropriate for
the mass integrand: degree-p polynomial shape-function products require at least
p+1 Gauss points per span; rational and configuration-dependent integrands require
quadrature convergence checks.

These are **mass-matrix implementations, not a complete transient formulation**.
In particular, using `full` in an ordinary `M*a` scheme omits the velocity-dependent
inertial residual and its consistent tangent. The scalar-twist time integration still needs validation before running transient
beam or FSI analyses.

`GetFirstDerivativesVector` returns `[VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
ANGULAR_VELOCITY_X]` per control point; `GetSecondDerivativesVector` uses the
corresponding acceleration variables. Allocate historical `VELOCITY`,
`ACCELERATION`, `ANGULAR_VELOCITY`, and `ANGULAR_ACCELERATION`. The angular X values
represent scalar twist derivatives, not global physical angular velocity or
acceleration. Both accessors honor the requested buffer step and reject missing
variables or invalid steps.

`CalculateDampingMatrix` returns Rayleigh damping `C = alpha*M + beta*K` of size
`4*n` by `4*n`, overwriting previous contents. Set `RAYLEIGH_ALPHA` (units 1/time)
and `RAYLEIGH_BETA` (units time) in properties, or in `ProcessInfo` as a fallback.
Properties take precedence, including explicit zero values; absent coefficients
default to zero. Coefficients must be finite and non-negative.

```python
properties.SetValue(IGA.RAYLEIGH_ALPHA, alpha)
properties.SetValue(IGA.RAYLEIGH_BETA, beta)
```

Only nonzero contributions are calculated. Nonzero alpha requires the mass
formulation and its properties; nonzero beta requires initialized material data.
`K` is the current tangent stiffness, not a cached initial stiffness. In a nonlinear
state with indefinite tangent stiffness, beta*K is not guaranteed dissipative.
This damping matrix does not supply the full formulation's additional inertial
terms or derivatives of state-dependent damping forces.

Run the matrix tests after building:

```bash
OMP_NUM_THREADS=1 build/Release/applications/IgaApplication/KratosIgaCoreTest \
  --gtest_filter='*IsogeometricBeam*Mass*'
```

Tests check exact quadratic consistent mass, reference-coordinate invariance,
rigid translation and pure-twist kinetic energy, bending rotary inertia, a finite
quarter-turn of an anisotropic section, the full matrix against finite differences
of deformed rational-curve kinematics, and invalid properties/states.

### Simplified-mass implicit vibration benchmark

`tests/test_isogeometric_beam_dynamics.py` runs the actual Newton strategy with
`ResidualBasedBossakDisplacementScheme(0.0)` (average-acceleration Newmark), no
external loads, and zero Rayleigh coefficients. A straight quadratic beam has all
DOFs fixed except the axial displacement of the final control point. Thus its
admissible displacement is `u(s,t)=(s/L)^2*q(t)`; this is a constrained semidiscrete
oscillator, not the continuum cantilever's fundamental axial mode.

In the small-amplitude limit, `m=rho*A*L/5`, `k=4*E*A/(3*L)`, and
`q(t)=q0*cos(sqrt(k/m)*t)` for initially zero velocity. Initial acceleration is
`-k/m*q0`. The test uses `L=2`, `E=1000`, `rho=3`, `A=2`, and `q0=1e-7` to make
geometric nonlinearity negligible. It runs two periods at 40 and 80 steps per
period, verifies displacement error and second-order time convergence, and checks
energy conservation. It validates translational dynamics only; twist integration
and the full inertia formulation remain separate work.

```bash
OMP_NUM_THREADS=1 python applications/IgaApplication/tests/test_isogeometric_beam_dynamics.py -v
```

The same file includes `test_simplified_mass_distributed_axial_vibration`, a
fixed-free rod with multiple free axial control-point DOFs and two excited modes:

```
u(s,t) = q0*[sin(pi*s/(2*L))*cos(omega*t)
             + 0.3*sin(3*pi*s/(2*L))*cos(3*omega*t)]
omega = pi/(2*L)*sqrt(E/rho)
```

This solves the continuum linear rod equation with `u(0,t)=0` and zero axial
traction at `s=L`. The initial displacement and acceleration fields are
interpolated at Greville points into a quadratic spline space; initial velocity
is zero. Small amplitude makes nonlinear effects negligible. The actual implicit
Kratos solve is checked against this continuum solution at 41 physical positions
over one fundamental period (three periods of the second mode). The test compares
8 spans/200 steps with 16 spans/400 steps and checks reduced error under combined
spatial and temporal refinement. This does not isolate the spatial convergence
order. As in the first benchmark, rotations and transverse motion are constrained.

### Constrained bending and torsion vibration

The dynamics test file also contains `test_simplified_mass_bending_free_vibration`
and `test_simplified_mass_torsion_free_vibration`. All coefficients except the
selected tip coefficient are fixed, so the admissible field is `(s/L)^2*q(t)`.
These are analytical semidiscrete oscillators, not continuum cantilever modes.
For initially zero velocity, `q(t)=q0*cos(sqrt(k/m)*t)` in the small-amplitude limit:

| Motion | Generalized mass m | Generalized stiffness k |
|---|---|---|
| Bending along Y | rho*A*L/5 | 4*E*I_N/L^3 |
| Scalar twist | rho*(I_N+I_V)*L/5 | 4*G*I_T/(3*L) |

Both check energy and time convergence over two periods, with 40 and 80 steps per
period. Bending uses the actual Kratos displacement Bossak/Newmark scheme. Torsion
uses a **test-only scalar Newmark/Newton driver** that assembles the actual beam
element residual, stiffness, and mass at every iteration and explicitly updates
twist velocity and acceleration. The core displacement Bossak scheme does not
update those angular derivatives. Thus this test validates torsional element
contributions, but does not establish production scheme support for scalar twist.

```bash
OMP_NUM_THREADS=1 python applications/IgaApplication/tests/test_isogeometric_beam_dynamics.py -v -k bending
OMP_NUM_THREADS=1 python applications/IgaApplication/tests/test_isogeometric_beam_dynamics.py -v -k torsion
```
