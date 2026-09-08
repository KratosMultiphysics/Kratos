# Seepage Boundary Prototype — Step 3: Switching Criterion — Design

**Issue:** #14672 (prototype for a seepage boundary), bullet 3 of 4.
**Depends on:** step 1 `GeoSeepageCondition`, step 2 `GeoSeepageNewtonRaphsonStrategy`.
**Date:** 2026-08-25

## Goal

Implement `GeoSeepageNewtonRaphsonStrategy::UpdateSeepageBoundaryConditions()`, the seam left by
step 2. Each non-linear iteration it decides whether a seepage node must change between a Dirichlet
boundary (`WATER_PRESSURE` fixed at 0) and a zero-flux Neumann boundary, applies that change, and
reports whether anything changed.

## Sign convention (verified)

The criterion depends on the meaning of a positive `WATER_PRESSURE`. Verified in the retention laws:

- `custom_retention/saturated_below_phreatic_level_law.cpp:31` —
  `return (rParameters.GetFluidPressure() < 0.0) ? SATURATED_SATURATION : RESIDUAL_SATURATION;`
- `custom_retention/van_genuchten_law.cpp:36` — `if (p > 0.0) { ...unsaturated formula... } else { return SATURATED_SATURATION; }`

So **`WATER_PRESSURE > 0` means suction, i.e. an unsaturated (dry) node**, and `WATER_PRESSURE < 0`
means saturated. This makes the criterion physically coherent: a dry node has nothing to drain and
should carry no flux, while a node with water arriving must be allowed to let it out.

`PORE_PRESSURE_SIGN_FACTOR` is `1.0` (`geo_mechanics_application_constants.h:34`) and plays no part
in this decision.

## The algorithm

Executed once per non-linear iteration, from the seam established in step 2. **At most one node
changes state per iteration, across the whole model.**

1. **Compute the nodal flow** for every node (see next section).
2. **Take the cached set of seepage nodes**, collected from all `GeoSeepageCondition`s.
3. **Neumann to Dirichlet — this direction takes precedence.** Among seepage nodes that are currently
   *free* and have `WATER_PRESSURE > 0`, select the one with the **highest** `WATER_PRESSURE`. If one
   exists: set its `WATER_PRESSURE` to `0.0`, `Fix` it, and return `true`.
4. **Dirichlet to Neumann.** Otherwise, among seepage nodes that are currently *fixed* and have
   nodal flow `> 0`, select the one with the **largest** nodal flow. `Free` it and return `true`.
5. Otherwise return `false`.

Returning `true` drives the seams built in step 2: `RebuildSystem()` runs, and `is_converged` is
forced to `false` so the solver takes at least one more iteration under the new configuration.

Ties are broken by the lowest node id, so the algorithm is deterministic and reproducible.

### Initial state

Every seepage node starts **fixed** at `WATER_PRESSURE = 0`, applied in
`GeoSeepageCondition::InitializeSolutionStep`. The seepage face therefore starts over-prescribed and
shrinks as nodes are released, which is what the "largest outflow is released" rule expects.

## State representation: nodal fixity is the single source of truth

A node's `WATER_PRESSURE` fixity **is** its mode: fixed means Dirichlet, free means Neumann. There is
no separate mode variable to fall out of step with it.

This is a deliberate change from step 1, which stored a mode string per condition. It was chosen
because pressure and flow are nodal quantities and the true exit point falls *between* nodes, so a
per-condition mode cannot represent the boundary faithfully. It also removes the parallel data race
recorded in the step 2 spec: with state held per node there is no shared node with two disagreeing
owners.

## Computing the nodal flow

The nodal flow is the assembled element right-hand side: permeability flow plus compressibility flow
plus fluid body flow. It is obtained directly rather than through any derived quantity.

Per iteration:

1. Zero a nodal accumulator (a `std::unordered_map<std::size_t, double>` keyed by node id, or a nodal
   variable) for the seepage nodes.
2. For each element in the model part, call
   `Element::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo)`.
3. Align the RHS with the element's degrees of freedom and accumulate.

### Aligning RHS entries with nodes

An element's RHS is a flat vector, and for a U-Pw element displacement and pressure entries are
interleaved, so entry *i* does not correspond to node *i*. Assuming a pressure-only ordering would
silently read displacement residuals as flows.

Instead, call `Element::GetDofList(dofs, rCurrentProcessInfo)`, which returns the dofs in the same
order as the RHS entries, and accumulate only those entries whose dof variable is `WATER_PRESSURE`:

```cpp
auto dofs = Element::DofsVectorType{};
r_element.GetDofList(dofs, rCurrentProcessInfo);
r_element.CalculateRightHandSide(rhs, rCurrentProcessInfo);

for (std::size_t i = 0; i < dofs.size(); ++i) {
    if (dofs[i]->GetVariable() == WATER_PRESSURE) {
        nodal_flow[dofs[i]->Id()] += rhs[i];
    }
}
```

This is correct for both `PwElement` and `UPwSmallStrainElement` without special-casing either.

For the prototype every element in the model part is visited. If that proves too slow, it can later
be narrowed to elements adjacent to seepage conditions using the existing
`custom_processes/find_neighbour_elements_of_conditions_process.h`. That optimisation is explicitly
**not** part of this step.

## Isolating the sign assumption

The rule "the largest **positive** nodal flow is released to Neumann" is the least certain part of
this design and the one most likely to be wrong. It is therefore isolated in one named function so
that inverting it is a one-line change:

```cpp
// Returns true if a Dirichlet seepage node carrying this nodal flow should be released to a
// zero-flux Neumann boundary. This sign is the key modelling assumption of the prototype and is
// what the Muskat validation is meant to confirm.
[[nodiscard]] bool ShouldReleaseToNeumann(double NodalFlow) { return NodalFlow > 0.0; }
```

## Changes to step 1

Making nodal fixity the source of truth retires part of `GeoSeepageCondition`. These are **deleted**
rather than left unused, so there is only one source of truth:

- `InitializeNonLinearIteration` — the strategy now sets fixity directly.
- `Initialize` and its `Properties` clone, plus the `mHasOwnProperties` guard and its serialisation.
- `GetBoundaryType` and `SetBoundaryType`.
- The `GEO_SEEPAGE_BOUNDARY_TYPE` variable and its registration, together with
  `Geo::SeepageDirichletType` and `Geo::SeepageNeumannType`.

`GeoSeepageCondition` keeps a clear and smaller responsibility:

- identify which nodes are seepage nodes, through its geometry;
- set the initial state in `InitializeSolutionStep`, fixing every node at `WATER_PRESSURE = 0`;
- `Check`, reduced to validating the geometry and the presence of `WATER_PRESSURE` dofs;
- contribute nothing to the linear system.

## Collecting the seepage nodes

The strategy caches the seepage node pointers once, in `Initialize`, by iterating the model part's
conditions and selecting those that `dynamic_cast` to `GeoSeepageCondition`. Duplicate nodes shared
by adjacent conditions are removed, so each node is considered exactly once. Conditions do not change
during a solve, so a single collection pass is sufficient.

## Testing

Unit tests for the decision logic, which is pure and easy to test in isolation:

- a free node with `WATER_PRESSURE > 0` is fixed and its pressure set to 0;
- among several such nodes, only the highest-pressure one switches;
- a fixed node with positive nodal flow is freed;
- among several such nodes, only the largest-flow one switches;
- when both kinds of candidate exist, the Neumann to Dirichlet switch wins and only it happens;
- when no candidate exists, nothing changes and the result is `false`;
- ties are broken by lowest node id;
- the RHS-to-dof mapping picks `WATER_PRESSURE` entries from an interleaved U-Pw ordering.

To keep these tests free of a full solver, the decision step is factored into a method taking the
seepage nodes and their nodal flows, so it can be exercised without assembling anything.

## Risks

- **The outflow sign** is the central modelling assumption, isolated as described above.
- **Migration rate.** Only one node switches per iteration and every switch resets convergence, so a
  face needing *n* switches costs at least *n* extra iterations. With a typical `max_iterations` of
  about 30, a fine mesh could exhaust the iteration budget before the boundary settles.
- **Oscillation.** A node could in principle be fixed and released on alternating iterations. Nothing
  in this design prevents that; the `mMaxIterationNumber` guard only bounds it. If the Muskat run
  oscillates, the natural remedy is to allow each node at most one switch per solution step.

## Open questions for #14637

- Whether one switch per iteration is fast enough in practice, or whether several non-adjacent nodes
  should be allowed to switch together.
- Whether recomputing every element's RHS each iteration is acceptable, or whether the neighbour
  restriction is needed.
- Whether `RebuildSystem()` (step 4) is genuinely required under a block builder and solver, since
  fixity changes are re-applied there on every build.

