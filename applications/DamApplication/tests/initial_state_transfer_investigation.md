# Elastic initial-state transfer investigation

This is a controlled prototype, not a general transfer feature. It does not
change DamAnalysis, its solvers, or the meaning of existing nodal stress input.
Run `test_dam_elastic_initial_state_transfer.py` with the built Dam, Structural
Mechanics, Constitutive Laws, and Poromechanics applications.

## Model and measurements

A 1 x 2 plane-strain body of unit thickness is discretized as four linear
triangles (one point per element), or two bilinear quadrilaterals (four points
per element). The two bottom nodes are fixed in both directions. The SMA
`SmallDisplacementElement` is used, with E = 10000, nu = 0.2, density = 2,
and vertical volume acceleration = -10. Gravity produces a total vertical
load of -40; its support reactions total +40. The increment is horizontal
volume acceleration = +2, while gravity and constraints remain in place.
The benchmark is small strain and does not move the mesh.

The reference solves gravity and then the increment. A fresh model is
initialized with the same material, load, and constraints. Displacements use
Dam's existing mapper; constitutive states are assigned to initialized laws by
calling the core `SetInitialStateProcess2D.ExecuteInitializeSolutionStep`
once per element submodelpart, before the initial Gauss-point query. No process
is registered for automatic reapplication on subsequent steps. The assembled
element RHS is external minus internal force; free entries measure equilibrium,
and negative constrained entries give reactions. The test also obtains external
forces as the difference between loaded and zero-acceleration RHS. Absolute
comparison tolerance is 1e-10 in the chosen units.

For the triangles, all three methods match the reference Gauss-point stresses,
assembled internal/external balance, reactions and incremental displacement and
stress to 1e-10. Maximum actual discrepancies were:

| Model and method | Initial stress error | Initial free residual | First solve displacement | Increment displacement error | Increment stress error |
|---|---:|---:|---:|---:|---:|
| Triangles, displacement | 0 | 7.11e-15 | 1.95e-18 | 1.73e-18 | 1.42e-14 |
| Triangles, stress only | 0 | 7.11e-15 | 2.09e-18 | 1.73e-18 | 7.11e-15 |
| Triangles, combined | 0 | 7.11e-15 | 1.95e-18 | 1.73e-18 | 1.42e-14 |
| Thermoelastic triangles, displacement | 0 | 6.22e-15 | 1.73e-18 | 8.67e-19 | 8.88e-15 |
| Thermoelastic triangles, stress only | 8.88e-16 | 6.22e-15 | 1.63e-18 | 8.67e-19 | 7.11e-15 |
| Thermoelastic triangles, combined | 0 | 6.22e-15 | 1.73e-18 | 1.73e-18 | 1.07e-14 |
| Quadrilaterals, displacement | 0 | 8.88e-15 | 6.50e-18 | 3.47e-18 | 1.42e-14 |
| Quadrilaterals, combined | 1.78e-15 | 8.88e-15 | 7.07e-18 | 5.20e-18 | 2.18e-14 |

The direct stress-only formulation intentionally resets the displacement datum:
its first displacement is zero, and its subsequent displacement is an increment,
not the total displacement of the reference calculation. Constrained values and
loads must use that datum consistently.

For the quadrilateral stress-only experiment, applying the first point's stress
to every point of each element gives a 4.018 maximum initial stress error, a
2.319 maximum free residual and a 3.257e-4 artificial first-step displacement.
This demonstrates that one constant state per element cannot represent this
nonuniform pointwise stress field. The combined quadrilateral result is a special
elastic case: the selected point's stress and strain cancel as an offset under
the same linear constitutive tensor, while mapped displacement reconstructs the
varying field. It is not a demonstration of pointwise stress transfer.

The thermoelastic test uses the Dam plane-strain adapter, uniform temperature 30,
reference temperature 20 and expansion coefficient 1e-5. Both temperatures are
mapped. The law uses the plane-strain thermal vector
`1.2e-4 * [1, 1, 0]`. With sigma0 equal to the original stress, the direct
stress-only test sets initial strain to minus that thermal vector; the combined
test uses the original elastic strain as initial strain. In this specific SMA
output path `GREEN_LAGRANGE_STRAIN_VECTOR` returned by the thermoelastic law is
already thermally reduced. Do not subtract the thermal term again. This test
exercises the thermoelastic constitutive law and mechanical solve; it does not
solve the thermal diffusion equation or validate the Dam thermal solver lifecycle.

## Scope and potential reuse

The core state contains stress, strain and deformation gradient, not arbitrary
history variables. SMA's elastic law computes `C * (strain - initial_strain -
thermal_strain) + initial_stress`. Results depend on the concrete law, not just
on setting the state pointer. The standard core process can use element geometry
values but assigns one state to all integration points of an element. GeoMechanics'
`ApplyFinalStressesOfPreviousStageToInitialState` instead captures each point's
stress and assigns it to each law, but is designed for staged displacements, calls
`InitializeMaterial`, and has not been tested here with Dam or separate meshes.
It is a candidate for investigation, not a drop-in replacement.

Dam maps nodal displacement, temperature and reference temperature. Its mapped
nodal Cauchy stress is output data; `TransferInitialStressToMainModelPartUtility`
adds it to smoothed nodal output, not to solid constitutive forces. This output
compatibility path, the associated mapping and the scheme behavior should remain
until replacements are verified. Cohesive interfaces separately use nodal
`INITIAL_STRESS_TENSOR`; UP/P behavior is distinct and their nodal stress output
is carried over rather than reconstructed by these schemes. These formulations
were not changed by this prototype.

The ordinary-process `ExecuteBeforeSolutionLoop` hook in DamAnalysis runs after
solver initialization and before output hooks. The prototype assigns the core
state after strategy initialization and queries stress before solving; it does not
exercise an actual DamAnalysis output-process callback. A small one-shot adapter
could call the existing process at that hook for supported inputs, but real data
import, per-point state fidelity, interaction with output and first-step process
ordering must be tested first. No production adapter is included here.

Same-mesh restart is a candidate for history-dependent continuation: SMA solid
elements serialize their constitutive-law vectors and some Dam damage laws
serialize history components. Completeness depends on each law, element, solver,
and nonlocal data; exact nonlinear continuation was not tested. Mapping nodal
fields to a different mesh does not guarantee a statically admissible stress
field, and neither this test nor a nodal stress interpolation establishes one.
Pointwise remapping, equilibration, nonlinear history transfer, and changing
material properties require separate verification.

## Existing regression results

On the existing Release build, the Dam small suite passes 136 tests (including
these five). The Dam C++ fast suite passes 48 of 49 tests; the known
`UPCouplingPreservesExpectedSign` failure remains. A temporary copy of
`/home/cimne/Simulaciones/00_DEBBUGING_TESTS/DamApplication/TestTranferUtility2D.gid/`
completed successfully with its supplied input and zero-RHS warnings. This run
checks execution only; no independent numerical comparison against a baseline
was performed. The original simulation files were not modified.
