# Seepage Boundary Prototype — Step 1: The Mixed Seepage Condition

**Date:** 2026-08-25
**Application:** GeoMechanicsApplication
**Related issues:** seepage boundary prototype; follow-up #14637

## Context

A seepage boundary is a boundary whose type is not known in advance. Where the
soil is saturated and water wants to leave the domain, the boundary behaves as a
prescribed-pressure (Dirichlet) boundary. Where it is unsaturated, no water
crosses and it behaves as a zero-flux (Neumann) boundary. The switch point
migrates during the solve, so the boundary type must be decided per non-linear
iteration rather than up front.

The prototype exists to find out whether the intended implementation is suited,
to validate it against the Muskat case created by John van Esch, and to let us
detail out follow-up issue #14637 from real findings.

The overall prototype has four parts:

1. A mixed seepage condition that can act as Dirichlet or Neumann. **(this spec)**
2. A Newton-Raphson strategy deriving from the one in core.
3. Switching logic in that strategy, driven by pressure and flux.
4. Rebuilding the system through the Builder & Solver when a condition switched.

This document specifies part 1 only. Parts 2 to 4 get their own spec, plan and
implementation cycle.

## Goal

A condition that owns no boundary-condition physics of its own beyond fixity. It
is told, from outside, which of two modes to be in, and it applies that mode at
the start of each non-linear iteration. Zero-flux Neumann needs no
implementation, so the condition contributes nothing to the linear system in
either mode.

## Non-goals

- Non-zero flux on the Neumann branch.
- Deciding *when* to switch. That is part 3.
- Rebuilding the system after a switch. That is part 4.
- 3D geometries, axisymmetric variants, or displacement degrees of freedom.

## Design

### Scope and class

A new standalone class `GeoSeepageCondition`, deriving directly from
`Kratos::Condition`, managing only `WATER_PRESSURE` degrees of freedom. It is
not templated on dimension or node count; the node count is taken from the
geometry at runtime.

It is registered for two geometries:

- `GeoSeepageCondition2D2N` — `Line2D2`
- `GeoSeepageCondition2D3N` — `Line2D3`

2D only, because the Muskat validation case is 2D. Further geometries are
deferred until the prototype has proven the approach.

### The mode property

A new application variable:

```cpp
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, GEO_SEEPAGE_BOUNDARY_TYPE)
```

This follows the existing `GEO_DRAINAGE_TYPE` precedent in this application, and
keeps the value readable in `.mdpa` and JSON input.

The only legal values are `"Dirichlet"` and `"Neumann"`. Any other value is a
hard error raised from `Check`, so a typo in an input file fails at startup
rather than silently degrading to one of the two behaviours.

The value lives in `GetProperties()` and is read fresh on every non-linear
iteration. It is never cached in a member variable. That is what makes the mode
switchable per iteration.

### Per-condition Properties ownership

In Kratos a `Properties` object is normally shared by many conditions. If the
mode were written to a shared `Properties`, flipping one condition would flip the
entire seepage face at once. The Muskat case is precisely the situation where
the switch point migrates gradually along the exit face, so per-condition
switching is required.

Therefore `Initialize()` clones the shared `Properties` into an instance owned
solely by this condition, and calls `SetProperties()` on itself:

```cpp
void GeoSeepageCondition::Initialize(const ProcessInfo&)
{
    if (mHasOwnProperties) return;
    SetProperties(Kratos::make_shared<PropertiesType>(GetProperties()));
    mHasOwnProperties = true;
}
```

The `mHasOwnProperties` guard makes this idempotent, so a second `Initialize()`
cannot silently reset a mode that has already been switched from outside.

The clone is a detached `Properties` that is not added to the `ModelPart`'s
properties container. The condition holds the only reference to it and
serializes it as part of itself, which is sufficient for the prototype.

`mHasOwnProperties` is saved and loaded by the serializer alongside the base
class. A restored condition has therefore already cloned its `Properties` and
will not clone again, which keeps the guard's meaning intact across a
save/restore cycle.

### How the mode is changed from outside, per non-linear iteration

This is the central mechanism of the prototype and the part most likely to
produce findings for #14637.

**The mechanism.** `Condition::GetProperties()` has a non-const overload
returning `PropertiesType&`. Because each condition owns its own `Properties`,
any external code holding a reference to a condition can write the mode
directly:

```cpp
r_condition.GetProperties().SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, "Neumann");
```

To keep the string literals in one place and out of the future strategy, the
condition exposes two thin accessors that are pure pass-throughs to
`GetProperties()` and hold no state of their own:

```cpp
[[nodiscard]] const std::string& GetBoundaryType() const;
void SetBoundaryType(const std::string& rType);
```

`SetBoundaryType` validates its argument and raises an error on anything other
than the two legal values, so a bug in the switching logic surfaces at the point
of the mistake rather than one iteration later.

**The timing.** One iteration of
`ResidualBasedNewtonRaphsonStrategy::SolveSolutionStep` runs in this order:

1. `mpScheme->InitializeNonLinearIteration(...)`, forwarded to every condition.
   This is where the condition reads the property and applies fixity.
2. `BuildAndSolve(...)`, where the builder & solver consults `IsFixed()`.
3. `mpScheme->FinalizeNonLinearIteration(...)`.
4. Convergence check, then the next iteration.

An external switch therefore takes effect if it is written after point 1 of
iteration *n* and before point 1 of iteration *n+1*. The natural slot is
immediately after point 3, where the pressures and fluxes of the just-computed
iterate are available. That is where the part-3 strategy override will sit.

**The ordering guarantee.** `Initialize()` runs during the strategy's
`Initialize`, strictly before the first `InitializeNonLinearIteration`. By the
time any external code switches a mode, every condition already owns its own
`Properties`, so a write cannot leak to its neighbours.

**Concurrency.** Because no two seepage conditions share a `Properties` after
`Initialize`, the future switching loop can safely be a `block_for_each` over the
conditions; each iteration writes to a distinct object.

**In this step there is no strategy yet.** The contract is verified by unit
tests that play the role of the strategy, as described under Testing.

### Behaviour per mode

`InitializeNonLinearIteration` reads the mode and applies it:

- **Dirichlet** — for each node of the geometry, set `WATER_PRESSURE` to `0.0`
  and `Fix` the `WATER_PRESSURE` degree of freedom.
- **Neumann** — for each node, `Free` the `WATER_PRESSURE` degree of freedom.
  Nothing else. Zero flux means zero contribution, so there is deliberately no
  integration code on this branch.

In both modes:

- `GetDofList` and `EquationIdVector` return exactly one `WATER_PRESSURE` entry
  per node.
- `CalculateLocalSystem`, `CalculateLeftHandSide` and `CalculateRightHandSide`
  produce correctly-sized zero output. The condition never contributes to the
  linear system; it only manipulates fixity.

### Check

`Check` validates, and raises a descriptive error otherwise:

- the geometry has at least two nodes;
- every node has `WATER_PRESSURE` as a solution step variable and as a degree of
  freedom;
- the `Properties` has `GEO_SEEPAGE_BOUNDARY_TYPE`;
- its value is one of the two legal strings.

It returns 0 on a valid setup.

### Known interaction to carry into parts 2 to 4

Fixity is changed in `InitializeNonLinearIteration`, which runs before
`BuildAndSolve`. A block builder & solver reads `IsFixed()` in
`ApplyDirichletConditions` and therefore honours a mid-solve change. An
elimination builder & solver bakes fixity into the equation IDs during
`SetUpSystem` and therefore needs an explicit rebuild after a switch.

This is exactly the "rebuild the system via the Builder & Solver" bullet of the
prototype. It is recorded here as a finding and addressed in part 4, not in this
step.

## Testing

C++ unit tests in
`tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`, in the
existing `KratosGeoMechanicsFastSuiteWithoutKernel` suite:

- `Create` returns a condition of the expected type, and `Info()` is correct.
- Two conditions constructed with the same `Properties` hold different
  `Properties` after `Initialize`, and mutating one does not affect the other.
- A second `Initialize` after a mode switch does not reset the mode.
- `"Dirichlet"` fixes every node and sets `WATER_PRESSURE` to 0.
- `"Neumann"` frees every node.
- Flipping the property between two `InitializeNonLinearIteration` calls
  switches the fixity, in both directions. This is the test that stands in for
  the strategy.
- `GetDofList` and `EquationIdVector` return one entry per node.
- `CalculateLocalSystem` yields a zero left-hand side and right-hand side sized
  to the node count.
- `SetBoundaryType` errors on an invalid value.
- `Check` returns 0 on a valid setup, and errors on a missing property, an
  invalid value, and a missing `WATER_PRESSURE` degree of freedom.

## Files touched

| File | Change |
| --- | --- |
| `custom_conditions/geo_seepage_condition.h` | new condition, declaration |
| `custom_conditions/geo_seepage_condition.cpp` | new condition, definition |
| `geo_mechanics_application_variables.h` | define `GEO_SEEPAGE_BOUNDARY_TYPE` |
| `geo_mechanics_application_variables.cpp` | create `GEO_SEEPAGE_BOUNDARY_TYPE` |
| `geo_mechanics_application.h` | condition member prototypes for both geometries |
| `geo_mechanics_application.cpp` | register the variable and both conditions |
| `tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp` | new tests |

## Decisions and rationale

| Decision | Rationale |
| --- | --- |
| Derive from `Condition`, `WATER_PRESSURE` only | Self-contained and easiest to reason about for a prototype. Displacement coupling is not needed to run Muskat. |
| The condition fixes its own degrees of freedom | Keeps "how to be Dirichlet" inside the condition. The strategy only flips a flag. |
| Dirichlet prescribes `WATER_PRESSURE = 0` | A classic seepage face exits at atmospheric pressure. No extra configuration. |
| Mode stored as `Variable<std::string>` | Matches the existing `GEO_DRAINAGE_TYPE` precedent and stays readable in input files. Validation in `Check` covers the typo risk. |
| Each condition clones its own `Properties` | Without it, one switch flips the whole face. Per-condition switching is required by the physics of a migrating exit point. |
| The condition clones the `Properties` itself | Automatic, with no burden on the input files, and no separate setup step to forget. |
| One class with a branch, not two behaviour classes | The Neumann branch is empty. Polymorphism here would add three classes to carry no behaviour. Revisit if non-zero flux is needed. |
| 2D lines only | Enough to run the Muskat case. The condition performs no integration, so extending the geometry coverage later is cheap. |


