# Nodal Water Flow Visualization — Design

**Date:** 2026-08-26
**Application:** GeoMechanicsApplication
**Status:** Approved for planning

## Problem

`Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows` assembles, per node, the
WATER_PRESSURE right-hand-side contributions of every element. This nodal flow drives the
seepage boundary switching in `GeoSeepageNewtonRaphsonStrategy` (via `SwitchOneSeepageNode`),
where the sign convention is the key modelling assumption: a positive nodal flow means outflow
(`ShouldReleaseToNeumann` returns `NodalFlow > 0.0`).

Today this value is transient — computed each non-linear iteration and discarded. There is no way
to look at it. To verify (or invert) the sign convention against a validation case such as Muskat,
we need to expose the nodal flow as a nodal variable that can be written to output and visualized.

## Goal

Add a permanent, first-class nodal `Variable<double>` — **`NODAL_WATER_FLOW`** — that holds the
assembled WATER_PRESSURE right-hand-side at each node, computed once per converged time step, stored
in the historical (solution-step) database, and written through the normal `nodal_results` GiD path.

## Non-goals

- No per-iteration history of the flow. Only the converged-step value is stored.
- No separate output process or scheme-level machinery.
- No computation in non-seepage runs (see Scope below).

## Design

### Variable

A new nodal `Variable<double>` named `NODAL_WATER_FLOW`. Its value at a node is exactly the entry
that `CalculateNodalWaterFlows` produces for that node (sum of the element WATER_PRESSURE RHS
entries). Sign convention matches `ShouldReleaseToNeumann`: positive = outflow.

### Where the value is produced

Override `FinalizeSolutionStep()` in `GeoSeepageNewtonRaphsonStrategy`
(`custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`):

```cpp
void FinalizeSolutionStep() override
{
    MotherType::FinalizeSolutionStep();               // keep base behaviour
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(
        BaseType::GetModelPart(), BaseType::GetModelPart().GetProcessInfo());
    Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(BaseType::GetModelPart(), nodal_flows);
}
```

This runs after the step converges, so the stored value is the flow of the converged solution — the
exact map `SwitchOneSeepageNode` consumes. That makes the visualized picture consistent with the
switching decision.

### New helper: `AssignNodalWaterFlows`

Add to `seepage_boundary_utilities` (`.h` / `.cpp`):

```cpp
void AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows);
```

Behaviour:

1. Set `NODAL_WATER_FLOW` to `0.0` on every node in the model part.
2. For each entry in `rNodalFlows`, write the value onto the matching node via
   `FastGetSolutionStepValue(NODAL_WATER_FLOW)`.

Zeroing first ensures nodes not present in the map (e.g. nodes without a WATER_PRESSURE dof) hold a
defined value rather than stale data. Splitting this out of the strategy keeps it unit-testable
without constructing a solver.

### Registration wiring (mirrors `HYDRAULIC_DISCHARGE`)

| File | Change |
| --- | --- |
| `geo_mechanics_application_variables.h` | `KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_WATER_FLOW)` |
| `geo_mechanics_application_variables.cpp` | `KRATOS_CREATE_VARIABLE(double, NODAL_WATER_FLOW)` |
| `geo_mechanics_application.cpp` | `KRATOS_REGISTER_VARIABLE(NODAL_WATER_FLOW)` |
| `custom_python/geo_mechanics_python_application.cpp` | `KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_WATER_FLOW)` |
| `GeoMechanicsApplication.py` | `NODAL_WATER_FLOW = KratosGeo.NODAL_WATER_FLOW` |
| `tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.cpp` | `KRATOS_REGISTER_VARIABLE(NODAL_WATER_FLOW)` |

### Making it addable & outputtable

| File | Change |
| --- | --- |
| `python_scripts/geomechanics_solver.py` (`_add_water_variables`) | `AddNodalSolutionStepVariable(NODAL_WATER_FLOW)` |
| `custom_workflows/dgeoflow.cpp` | `rModelPart.AddNodalSolutionStepVariable(NODAL_WATER_FLOW)` |
| `custom_workflows/dgeosettlement.cpp` | `rModelPart.AddNodalSolutionStepVariable(NODAL_WATER_FLOW)` |
| `custom_workflows/geo_output_writer.cpp` (`WriteNodalOutput` map) | `{"NODAL_WATER_FLOW", MakeNodalResultWriterFor(NODAL_WATER_FLOW)}` |

Because the value already lives in the historical database, the output writer needs no special-case
computation (unlike `HYDRAULIC_HEAD`, which is a C++-only recomputation). The standard
`WriteNodalResults` path reads the historical value directly.

## Data flow

```
converged step
  -> GeoSeepageNewtonRaphsonStrategy::FinalizeSolutionStep
       -> CalculateNodalWaterFlows(model_part)      // NodalFlowMap
       -> AssignNodalWaterFlows(model_part, map)    // writes NODAL_WATER_FLOW (historical)
  -> GeoOutputWriter::WriteNodalOutput
       -> WriteNodalResults(NODAL_WATER_FLOW, ...)   // GiD, when listed in nodal_results
```

## Scope & trade-offs

- Only `newton_raphson_with_seepage` runs populate `NODAL_WATER_FLOW`. Other strategies leave it at
  its default (0.0). This keeps the extra RHS assembly where it is meaningful and guarantees the
  stored value equals the switching input.
- One extra full RHS assembly per converged step in seepage runs. Acceptable: it happens once per
  step, not once per iteration.

## Testing

- **C++ unit test** for `AssignNodalWaterFlows`: build a small model part with `NODAL_WATER_FLOW`
  as a solution-step variable, feed a `NodalFlowMap` covering some but not all nodes, and assert:
  - listed nodes hold the mapped value;
  - unlisted nodes are zeroed.
- **Integration (optional)**: extend an existing seepage integration test's
  `ProjectParameters.json` `nodal_results` to include `NODAL_WATER_FLOW`, and assert a known sign at
  a draining node — this doubles as verification of the sign convention.

## Open questions

None.

