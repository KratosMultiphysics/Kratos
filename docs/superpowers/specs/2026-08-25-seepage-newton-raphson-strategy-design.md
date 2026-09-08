# Seepage Boundary Prototype — Step 2: Seepage Newton-Raphson Strategy — Design

**Issue:** #14672 (prototype for a seepage boundary), bullet 2 of 4.
**Depends on:** step 1, `GeoSeepageCondition` (`docs/superpowers/specs/2026-08-25-seepage-boundary-condition-design.md`).
**Date:** 2026-08-25

## Goal

Create a Newton-Raphson strategy that owns the *decision* to switch a seepage boundary between
Dirichlet and Neumann, and that can force further iterations when it does so. This step delivers the
strategy skeleton and the seams; the switching criterion (step 3) and the system rebuild (step 4)
plug into those seams.

## Background: what step 1 already gives us

Two facts, verified in the codebase, define the boundary of this step.

**Applying a switch needs no strategy code.** The Newton loop calls
`p_scheme->InitializeNonLinIteration(...)` on every iteration
(`kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h`, lines 942 and 996),
*before* `BuildAndSolve`. Both the core scheme (`kratos/solving_strategies/schemes/scheme.h:386`, via
`EntitiesUtilities::InitializeNonLinearIterationAllEntities`) and the GeoMechanics scheme
(`applications/GeoMechanicsApplication/custom_strategies/schemes/geomechanics_time_integration_scheme.hpp:102`,
via `BlockForEachActiveCondition`) forward that call to every active condition. So
`GeoSeepageCondition::InitializeNonLinearIteration` already runs once per iteration and already
applies the configured fixity. The strategy therefore only decides *what* the mode should be; the
condition applies it.

**The builder and solver determines whether a rebuild is required.** The elimination builder bakes
fixity into equation IDs (`residualbased_elimination_builder_and_solver.h:683`):

```cpp
if (dof_iterator->IsFixed()) dof_iterator->SetEquationId(--fix_id);
else                         dof_iterator->SetEquationId(free_id++);
BaseType::mEquationSystemSize = fix_id;
```

`SetUpSystem` runs once per *solution step* from `InitializeSolutionStep`, not per iteration. A
mid-iteration switch therefore leaves stale equation IDs and a wrong system size under the
elimination builder. The block builder assigns IDs to all dofs and re-applies Dirichlet conditions
on every build, so switches are honoured there without a rebuild. GeoMechanics defaults to
`block_builder: true` (`python_scripts/geomechanics_solver.py`), so the block path is the one the
prototype will exercise first. This is the finding recorded as open at the end of step 1.

## Architecture

A new header-only class template:

```
applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
```

```cpp
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoSeepageNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
```

Header-only `.hpp` matches every existing strategy in that folder. It derives from the **core**
strategy, as the issue specifies, rather than from `GeoMechanicsNewtonRaphsonStrategy`. The
consequence is that it takes no GeoMechanics `Parameters` block, so its constructor mirrors the core
one:

```cpp
GeoSeepageNewtonRaphsonStrategy(ModelPart&                                 rModelPart,
                                typename TSchemeType::Pointer              pScheme,
                                typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                typename TBuilderAndSolverType::Pointer    pNewBuilderAndSolver,
                                int  MaxIterations          = 30,
                                bool CalculateReactions     = false,
                                bool ReformDofSetAtEachStep = false,
                                bool MoveMeshFlag           = false);
```

## `SolveSolutionStep`: a copied loop with three seams

`SolveSolutionStep` is copied from the core strategy and modified only at named seam points. It is
copied rather than hooked because the prototype is expected to need further changes inside the loop,
and because one of the required behaviours — suppressing convergence — cannot be expressed from any
existing hook (see below).

All base members the loop touches (`mpA`, `mpDx`, `mpb`, `mpConvergenceCriteria`,
`mMaxIterationNumber`, `mCalculateReactionsFlag`, `mUseOldStiffnessInFirstIteration`,
`mStoreNonconvergedSolutionsFlag`, `mNonconvergedSolutionsMatrix`) are `protected` in the core
class, so the copy compiles. `GeoMechanicsNewtonRaphsonStrategy` already relies on this for
`mpA`/`mpb`/`mpDx`.

### Seam 1 and 2 — decide, then rebuild

A single `bool any_switched` is declared once, alongside the loop's other state variables (next to
`residual_is_updated`, core line 941), and **re-assigned** at each seam rather than redeclared. This
avoids one variable shadowing another between the first-iteration block and the iteration cycle.

The seams are inserted immediately after `UpdateDatabase(...)`, in both the first-iteration block and
the iteration cycle:

```cpp
UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

// Seam 1 (step 3): decide whether any seepage condition must flip.
any_switched = UpdateSeepageBoundaryConditions();

// Seam 2 (step 4): the dof set changed, so the system must be rebuilt.
if (any_switched) RebuildSystem();
```

This is the correct window. The solution has just been updated, so pressures and fluxes are current,
and the *next* iteration's `InitializeNonLinIteration` is what applies the new fixity to the nodes.

### Seam 3 — a switch must suppress convergence

When a condition has just switched, the current residual describes the *old* boundary configuration.
Declaring convergence at that point would end the step on a solution that does not satisfy the
boundary conditions now in force. So after each `PostCriteria` call the strategy forces another
iteration:

```cpp
if (is_converged) {
    // ... optional RHS refresh ...
    is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
}

// Seam 3: a switch invalidates this iteration's convergence verdict.
if (any_switched) is_converged = false;
```

This must be applied in **both** places the core loop evaluates `PostCriteria`: the first-iteration
block (core line 986) and the iteration cycle (core line 1064).

This requirement is why the loop is copied rather than hooked through the virtual `UpdateDatabase`:
`UpdateDatabase` is called *before* `PostCriteria`, so an override of it cannot suppress the
convergence verdict that follows.

Termination remains bounded by the existing `iteration_number++ < mMaxIterationNumber` guard, so a
pair of conditions that oscillate cannot loop forever; the step will exit through
`MaxIterationsExceeded()`.

## Step 2 scope

In this step both new methods are stubs, so the strategy is behaviourally identical to its base:

```cpp
protected:
    // Returns true if any seepage condition changed its boundary type this iteration.
    // Step 3 implements the criterion; for now nothing ever switches.
    virtual bool UpdateSeepageBoundaryConditions() { return false; }

    // Re-establishes the dof set and system after a switch changed the number of free dofs.
    // Step 4 implements this; under a block builder and solver it may prove unnecessary.
    virtual void RebuildSystem() {}
```

Explicitly **not** in this step: collecting the `GeoSeepageCondition`s from the model part and the
switching criterion (both step 3), and the rebuild itself (step 4).

## Python registration

Register in `applications/GeoMechanicsApplication/custom_python/add_custom_strategies_to_python.cpp`
following the existing pattern, minus the `Parameters` argument, exposing the class as
`GeoSeepageNewtonRaphsonStrategy`. Wiring it into `geomechanics_solver.py` is deferred until the
strategy does something observable (step 3).

## Testing

**By explicit decision, this step ships without tests.** It is a prototype, and in this step the
strategy is behaviourally identical to the base class, so a test would assert only that a copy is
faithful.

The risk this accepts is stated plainly: a transcription error in ~150 lines of intricate copied
code would surface later as a confusing convergence bug, and the copy will silently drift from core
on future Kratos merges. Mitigation is to copy verbatim and confine the diff to the three seams, so
that diffing the method against core stays a quick and meaningful check. Steps 3 and 4 will exercise
the loop for real, and that is where tests should return.

## Known issue deferred to step 3

The scheme applies conditions through `BlockForEachActiveCondition`, i.e. in parallel. Adjacent line
conditions share end nodes, so once neighbouring conditions can disagree about their mode, the
`Fix`/`Free` and the `WATER_PRESSURE` write on a shared node become a data race with order-dependent
results. This is latent from step 1 and harmless while every condition holds the same mode, but it
must be resolved in step 3, when the criterion first makes neighbours disagree along the seepage
face. Likely resolutions: decide all modes in the strategy and apply nodal fixity there in a single
threaded pass, or define a deterministic precedence rule for a shared node (for example, Dirichlet
wins).

## Open questions for #14637

- Whether suppressing convergence on every switch converges in practice, or whether the exit point
  oscillates between two adjacent conditions and exhausts the iteration budget. If it oscillates, a
  damping rule — such as allowing a condition to switch at most once per solution step — may be
  needed.
- Whether the block builder and solver alone is sufficient, making `RebuildSystem()` a no-op in
  practice and step 4 much smaller than the issue anticipates.


