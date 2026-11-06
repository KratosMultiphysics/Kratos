### AI Disclosure
“This content was partially created using Copilot365. Copilot365 uses a Microsoft business license. Data processed within Copilot365 remains within Deltares’ secure Microsoft 365 environment and is not used to train external AI models. The author has reviewed the content, made adjustments where necessary, and remains responsible for its quality, accuracy, and final version.”

# Functional Test Matrix

The existing execution-time suites group tests by how they are run. In addition, use functional keywords to describe what each Python test case verifies. Apply keywords at the individual test-case level, rather than to an entire test file, so cases in the same file can have different coverage.

Execution-time suites and functional keywords are separate classifications: each test case belongs to one execution-time suite, but may have multiple functional keywords. Functional suites select cases by keyword and may overlap with one another.

The inventory below organizes the detailed test evidence by subject. The dimensions below define cross-cutting keyword categories; a test case can be tagged in several dimensions. Use a consistent `dimension:variant` format for keywords (for example, `analysis:transient`, `element:upw`, or `spatial:3d`).

Assess coverage against the applicable variants in the matrix. Mark unsupported or inapplicable variants as such rather than treating them as missing coverage. A variant is covered only when a test meaningfully verifies the associated behaviour; merely using an element, boundary condition, or output format does not count. Require evidence for each applicable variant, but do not require tests for every possible combination of variants.

| Keyword dimension | Variants to distinguish at minimum |
|---|---|
| Element family | `UPw`, `Pw`, interface, piping, thermal, structural elements |
| Spatial formulation (stress state) | plane strain, axisymmetric, 3D, line interface, surface interface |
| Analysis | steady-state, transient, staged, restart/state serialization |
| Constitutive behaviour | linear elastic, Mohr-Coulomb, other supported models, UMAT/UDSM integration |
| Saturation | saturated, partially saturated, Van Genuchten |
| Numerical robustness | mesh refinement, time-step refinement, solver variants, platform/compiler |
| Error handling | missing parameter, wrong type, physically invalid value, incompatible settings (mainly covered by C++ unit tests)|
| Boundary conditions | Dirichlet, Neumann, hydrostatic, seepage/switching, water loads (user-facing behaviour) |
| Output | JSON, GiD, VTK (user-facing behaviour) |

# Inventory by Subject

The subject inventories are complementary, not exclusive: a test case may provide evidence for more than one subject. **In progress** means that a detailed inventory has been started; **Planned** means that the subject has not yet been inventoried in detail. These statuses describe the inventory, not whether tests exist. Files listed below are candidates to review, not a claim that every test in them provides sufficient evidence for that subject.

- [Elements](#elements)
  - [UPw Elements](#upw-elements) — **In progress**
  - [Pw Elements](#pw-elements) — **Planned**
  - [Interface Elements](#interface-elements) — **Planned**
  - [Piping Elements](#piping-elements) — **Planned**
  - [Thermal Elements](#thermal-elements) — **Planned**
  - [Structural Elements](#structural-elements) — **Planned**
- [Constitutive Behaviour](#constitutive-behaviour) — **Planned**
- [Groundwater Flow and Saturation](#groundwater-flow-and-saturation) — **Planned**
- [Boundary Conditions and Loads](#boundary-conditions-and-loads) — **Planned**
- [Processes and Initial Conditions](#processes-and-initial-conditions) — **Planned**
- [Solvers, Analysis and Time Integration](#solvers-analysis-and-time-integration) — **Planned**
- [End-to-End Workflows and Benchmarks](#end-to-end-workflows-and-benchmarks) — **Planned**
- [Output and Post-processing](#output-and-post-processing) — **Planned**

## Elements
### UPw Elements
#### Intended behaviour
The UPw elements couple the displacement and the water pressure field within soil.
Relevant behaviour to assess includes:
- mechanical response,
- groundwater flow,
- hydro-mechanical coupling,
- effective stress,
- pore pressure response,
- drained behaviour,
- undrained behaviour,
- transient behaviour,
- body forces,
- boundary contributions,
- behaviour for supported geometries and dimensions.

#### Risk level
**Very high** impact on model outcomes. Incorrect behaviour can lead to incorrect displacements, stresses or pore pressures while the analysis may still converge successfully.

#### Testing requirement
According to the Deltares software development canvas, a product with TRL 7-8 should have 100% unit test coverage of logic.

#### Existing tests

This is an initial, Python-only inventory and is still being reviewed. It is
not a complete list; classifications should be confirmed against each test's
assertions and verification oracle.

| Test file | Focus | Test level | Evidence type | What is assessed? |
|---|---|---|---|---|
| `test_elements.py` | Linear-elastic response across element geometries | Python integration | Analytical | Checks UPw element response against analytical displacements, total/effective stresses and strains for triangles, quadrilaterals and tetrahedra, including different node orders and 2D/3D geometries. Nodal stress extrapolation/smoothing is not considered covered merely because nodal stress output is requested. |
| `test_compressibility.py` | Compressibility contributions in UPw formulations | Python workflow integration | Reference / expected numerical behaviour | Checks the pore-pressure response for small-strain, updated-Lagrangian, FIC and different-order formulations using expected pore pressures in a one-dimensional model. |
| `test_strain_measures.py` | Small and large deformation using linear/Hencky strain measures | Python integration | Analytical | Compares small-strain and updated-Lagrangian UPw elements, including different-order formulations, against analytical vertical displacement. |
| `test_single_element_with_Mohr_Coulomb.py` | Mohr-Coulomb response of UPw elements | Python integration | Reference / expected numerical behaviour | Checks mechanical response and shear capacity for small-strain and different-order elements using expected displacements, Cauchy stress and shear capacity. |
| `test_element_lab.py` | Geotechnical element response: triaxial, oedometer, DSS and CRS | Python integration | Reference / expected numerical behaviour | Checks drained/undrained and time-dependent response of small-strain different-order and updated-Lagrangian UPw elements using expected displacements, stresses, strains and pore pressure. |
| `test_apply_initial_uniform_stress_field.py` | Initial uniform stress applied to a UPw element | Python workflow integration (indirect) | Property / invariant | Checks that an imposed initial stress is correctly represented by the element through the Cauchy stress at the integration points. The primary focus is the initialisation process. |

Tests that merely use UPw elements while primarily assessing a boundary
condition, process, solver strategy or larger project scenario are not included
as direct UPw element evidence.

UPw interface elements belong under the Interface Elements subject below,
rather than in this bulk UPw element inventory.

### Pw Elements

**Inventory status: Planned.**

Candidate tests to review include `test_steady_state_groundwater_flow.py`,
`test_transient_groundwater_flow.py`, `test_water_pressure.py`. Separate
direct element evidence from tests whose primary purpose is groundwater-flow
or boundary-condition behaviour.

### Interface Elements

**Inventory status: Planned.**

Candidate tests to review include `test_interface.py`,
`test_line_interface_elements.py`, `test_surface_interface_elements.py`,
`test_upw_interface.py`, and `test_soil_structure_interactions.py`. Distinguish
interface formulations and geometries, and separate direct element evidence
from broader soil-structure interaction workflows.

### Piping Elements

**Inventory status: Planned.**

Candidate tests to review include `test_consecutive_pipe_lines.py`,
`test_three_dimensional_piping_validation.py`, `test_sellmeijers_rule.py`,
and `test_sellmeijers_rule_validation.py`. Classify element-level behaviour
separately from piping validation and complete-scenario tests.

### Thermal Elements

**Inventory status: Planned.**

Candidate tests to review include `test_transient_thermal.py` and
`test_transient_thermal_validation.py`. Identify which thermal behaviours and
formulations each test actually verifies.

### Structural Elements

**Inventory status: Planned.**

Candidate tests to review include `test_curved_beam_elements.py`,
`test_truss_backbone_mat.py`, and
`test_deactivation_with_structural_element.py`. Record element response
separately from structural-element use in a larger workflow.

## Constitutive Behaviour

**Inventory status: Planned.**

Candidate tests to review include `test_single_element_with_Mohr_Coulomb.py`,
`test_mohr_coulomb_with_tension_cutoff.py`,
`test_piecewise_linear_moment_capacity.py`, `test_truss_backbone_mat.py`, and
`test_udsm_serialization.py`. Organize coverage by supported law and behaviour
(for example, elastic, Mohr-Coulomb, tension cut-off, and external-law
integration/serialization).

## Groundwater Flow and Saturation

**Inventory status: Planned.**

Candidate tests to review include `test_elementary_groundwater_flow.py`,
`test_steady_state_groundwater_flow.py`,
`test_transient_groundwater_flow.py`, `test_partial_saturation.py`,
`test_nodal_hydraulic_head.py`, `test_pressure_line_element.py`,
`test_pressure_point_flux.py`, `test_seepage.py`, `test_Muskat.py`,
`test_sellmeijers_rule.py`, and `test_sellmeijers_rule_validation.py`.
Distinguish flow regime, saturation/retention behaviour, and piping-related
validation.

## Boundary Conditions and Loads

**Inventory status: Planned.**

Candidate tests to review include `test_conditions.py`, `test_dirichlet_u.py`,
`test_dirichlet_u_constant.py`, `test_dirichlet_release.py`,
`test_line_loads.py`, `test_normal_load_on_1d_element.py`,
`test_normal_load_on_hexa_element.py`, `test_water_pressure.py`,
`test_pressure_line_element.py`, `test_pressure_point_flux.py`, and
`test_seepage.py`. Tag the behaviour actually checked, rather than every
boundary condition or load merely present in a model.

## Processes and Initial Conditions

**Inventory status: Planned.**

Candidate tests to review include `test_apply_initial_uniform_stress_field.py`,
`test_k0_procedure_process.py`, `test_soil_weight.py`,
`test_column_changing_waterlevel.py`, `test_set_moving_load_process.py`,
`test_set_multiple_moving_load_process.py`, `test_reset_totals.py`, and
`test_parameter_field.py`. Identify initialization, state update, and process
behaviour independently where a test covers more than one.

## Solvers, Analysis and Time Integration

**Inventory status: Planned.**

Candidate tests to review include `test_geomechanics_solver.py`,
`test_dynamics.py`, `test_dynamics_long.py`, `test_avoid_small_end_step.py`,
`test_master_slave_constraints.py`, and `test_prescribed_derivatives.py`.
Inventory supported analysis modes and time-integration behaviour, including
transient, staged, and restart/serialization cases where present.

## End-to-End Workflows and Benchmarks

**Inventory status: Planned.**

Candidate tests to review include `test_building_pit.py`, `test_excavation.py`,
`test_settlement.py`, `test_dsettlement_validation.py`,
`test_soil_structure_interactions.py`, `test_CROW_case.py`,
`test_benchmark_set_1.py`, and `test_benchmark_set_2.py`. Record the scenario
and verification oracle for each case; using a feature in a scenario alone
does not establish its direct functional coverage.

## Output and Post-processing

**Inventory status: Planned.**

Candidate tests to review include `test_integration_node_extrapolation.py`,
`test_elements.py`, and `test_water_pressure.py`. Determine explicitly which
output formats and post-processing operations are asserted; requesting or
writing output without checking its values does not count as verification.
`