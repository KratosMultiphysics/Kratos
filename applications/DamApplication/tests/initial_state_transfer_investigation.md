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

## Lifecycle and pointwise assignment follow-up

A test-only Python process now invokes `SetInitialStateProcess2D` once from
`ExecuteBeforeSolutionLoop`, after the Dam thermomechanical solver initializes.
A real DamAnalysis writes the configured nonzero state into its ASCII GiD t=0
frame; the same stress is present in the next frame after the thermal and
mechanical solvers run. This is a fully constrained, zero-temperature-increment
case. It confirms the hook and output ordering, but it does not validate the
requested gravity equilibrium, load increment, or thermal diffusion response in
DamAnalysis. The adapter is test-only, not a supported production input format.
It does not execute the core process at each solution step.

A C++ Q4 test obtains the individual law pointers through
`CalculateOnIntegrationPoints(CONSTITUTIVE_LAW)` and assigns one core
`InitialState` per integration point. A nonuniform reference displacement field
produces different Gauss-point stresses; the test reproduces all four stresses,
the element internal-force residual and tangent, and the point stresses and
internal force after an identical prescribed displacement increment. Tolerance
is 1e-9 absolute. This is an element-level comparison, not a solved gravity
continuation or a test of reactions. No element or law changes are required for
C++ pointwise assignment. The Python `CalculateOnIntegrationPoints` overloads do
not expose the law-pointer vector, and the Python `ConstitutiveLaw` binding does
not expose `SetInitialState` or the `InitialState` class.

Assignment and transport are separate: output can record Gauss-point stresses,
but these output files are not a validated general initial-state import path.
Same-mesh transfer needs identification by element ID, integration rule and
point order, compatible materials and stress measures, and validation of state
sizes and equilibrium. `SetValuesOnIntegrationPoints` delegates to individual
law `Has`/`SetValue` implementations and cannot be assumed to populate
`InitialState`. GeoMechanics' staged process directly captures point stresses
and sets law states in C++, but it reinitializes material and is not established
as a Dam import mechanism. Nodal mapping and Gauss-to-node extrapolation do not
provide a guaranteed equilibrated Gauss-point state on a new mesh. Restart is a
separate same-mesh option whose history completeness still requires testing.

A rerun of the supplied `TestTranferUtility2D.gid` case from a temporary copy
produced a GiD result file byte-identical to the preceding temporary run. An
instrumented copy isolated all ten zero-RHS warnings to the thermal solver. The
builder emits that warning when its assembled RHS has exactly zero norm; the
output temperatures are zero initially and remain zero for all 11 frames, and
this case has no applied thermal temperature or heat-flux load in its process
lists. These observations are consistent with a stationary zero thermal state,
but they do not establish mechanical equilibrium or validate every inactive or
constrained degree of freedom. The original case files were not modified.

The updated Dam small suite passes 137 tests. The rebuilt Dam C++ fast suite
passes 49 of 50, with the existing `UPCouplingPreservesExpectedSign` failure.
Production integration still requires a validated input path for per-element or
per-point values, explicit supported element/law checks, and end-to-end gravity,
reaction, incremental-load and thermomechanical lifecycle tests. The current
benchmark's absolute tolerances are documented above; relative error acceptance
criteria remain to be added to the end-to-end continuation tests.

## Solved quadrilateral gravity continuation and restart

The C++ `QuadrilateralGravityPerGaussPointInitialStateContinuation` regression
now extends the previous element-level pointwise experiment to a complete
assembled solve. It uses the same 1 x 2 body, two Q4 elements, material,
constraints, gravity, and horizontal acceleration increment described above.
Stage A stores each four-point stress vector under its element ID and point
index. Stage B initializes its constitutive laws, starts with zero
nodal displacement, and calls core `ConstitutiveLaw::SetInitialState` separately
for every point; gravity and constraints are retained.

The test checks pointwise stress, separately reconstructed internal and
external force vectors, the assembled RHS, support reactions, first-solve
stress and displacement, and the response to the horizontal load increment.
It uses `abs_error <= 1e-10 + 1e-10 * abs(reference)` for stresses and forces,
and `1e-12 + 1e-10 * abs(reference)` for displacements. Stage A has a maximum
free-DOF residual of 8.88e-15 and vertical reactions summing to 40. Stage B's
maximum first-solve displacement is 6.50e-18, and the maximum error in the
incremental displacement compared with Stage A is 1.74e-18. The displacement
comparison is between increments, because Stage B uses a new datum. The
maximum assembled RHS difference after the increment is 1.42e-14, the
continuation free-DOF residual is 5.55e-15, and the maximum Gauss-point stress
difference is 1.92e-14. All pointwise
stresses and reactions pass their stated tolerances. This demonstrates a valid
same-mesh elastic case, not arbitrary input-field equilibration.

A separate Python regression serializes the complete Model using Kratos
`StreamSerializer` after quadrilateral gravity equilibrium, reloads it with
`IS_RESTARTED` set, and compares displacement, stress and residual before and
after applying the horizontal increment. It passes the existing 1e-10 absolute
checks. This is an in-memory serializer experiment using an SMA linear strategy,
not an end-to-end DamAnalysis restart. SMA's base solid element serializes both
its integration method and constitutive-law vector and skips material
reinitialization when `IS_RESTARTED` is set. Core's `RestartUtility` supplies a
file-based save/load mechanism and marks the loaded model as restarted. The
current Dam mechanical and thermomechanical solvers, however, accept only mdpa
input in `ImportModelPart`, and their `SaveRestart` methods are empty; a normal
DamAnalysis restart has not been established here.

The existing `RestartLocalDamage3DContinuation` and
`RestartNonlocalSimoJuContinuation` C++ tests pass for the examined thermal
Simo-Ju damage laws, as do the multi-properties rebinding and frozen legacy
damage-element tests. These are element/law serialization tests with prescribed
history and, for nonlocal damage, explicitly supplied nonlocal driving; they do
not prove a complete DamAnalysis restart or preservation of every law and
nonlocal neighborhood state. Exact continuation with the original displacement
reference and complete history is the restart use case. An import of prescribed
stresses with zero displacement serves a different use case: it creates a new
mechanical datum and requires an equilibrated, compatible stress field.

For a prospective import, the minimum evidence is a matching mesh and reference
configuration, element identity and connectivity, element formulation and
integration rule, point ordering, stress measure and component convention,
material/law and parameter compatibility, consistent loads and constraints, and
explicit rejection of unsupported history-dependent or thermomechanical cases.
Assignment can reuse the law pointers and `InitialState` demonstrated in this
test. Data transport and compatibility validation are not supplied by the
per-element core process, the current Python bindings, or the Dam nodal mapper.
A generic same-mesh importer could be considered after validating a durable
export format and the real DamAnalysis lifecycle; it should not perform
interpolation. For exact same-mesh continuation, investigate wiring the existing
restart mechanism into DamAnalysis before creating a separate stress importer.
No production integration is proposed as validated by these tests alone.

After this follow-up, the Dam small suite passes 138 tests. The rebuilt Dam C++
fast suite passes 50 of 51; `UPCouplingPreservesExpectedSign` remains the known
pre-existing failure.
