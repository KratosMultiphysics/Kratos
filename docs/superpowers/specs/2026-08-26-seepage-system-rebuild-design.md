# Seepage Boundary Prototype — Step 4: System Rebuild After a Switch — Design

**Issue:** #14672 (prototype for a seepage boundary), bullet 4 of 4.
**Depends on:** step 2, `GeoSeepageNewtonRaphsonStrategy`
(`docs/superpowers/specs/2026-08-25-seepage-newton-raphson-strategy-design.md`), and step 3, the
switching criterion (`docs/superpowers/specs/2026-08-25-seepage-switching-criterion-design.md`).
**Date:** 2026-08-26

## Goal

Give `GeoSeepageNewtonRaphsonStrategy::RebuildSystem()` a real body. When a seepage node has just
switched between a Dirichlet boundary (`WATER_PRESSURE` fixed at zero) and a zero-flux Neumann
boundary (`WATER_PRESSURE` free), the set of *free* degrees of freedom has changed. The strategy must
re-shape the linear system to match, and re-solve so the current iteration's database is consistent
with the boundary configuration now in force.

This is the final stubbed piece of the prototype. After this step the strategy switches boundaries,
suppresses convergence on the iteration that switches (step 2, seam 3), and rebuilds the system here.

## Scope

- **In scope:** implement `RebuildSystem()`; add a focused C++ integration test that drives a fixity
  flip through the rebuild and asserts the system was re-shaped.
- **Out of scope (deferred to #14637):** wiring the strategy into `geomechanics_solver.py`, building a
  full seepage-boundary Muskat validation case, and judging the switching sign against a reference
  solution. Those need a running end-to-end solve and belong to the follow-up issue.

## Background: what a switch changes

A seepage switch never adds or removes a degree of freedom. Every seepage node always owns exactly
one `WATER_PRESSURE` dof. The switch only flips that dof between **fixed** (Dirichlet) and **free**
(Neumann). What that flip means for the linear system depends entirely on which builder and solver is
in use, and the two behave very differently. This distinction is the crux of step 4, so it is spelled
out in full below.

## The two builders and solvers, contrasted

Kratos offers two builder-and-solver implementations that matter here. They differ in how they turn
the dof set into a system of equations, and therefore in whether a mid-iteration fixity flip needs a
rebuild at all.

### Elimination builder and solver (`ResidualBasedEliminationBuilderAndSolver`)

- **Only free dofs enter the system.** Fixed (Dirichlet) dofs are *eliminated*: they are given no
  equation in the linear system.
- **Equation ids depend on fixity.** `SetUpSystem` walks the dof set and hands out equation ids like
  this (`residualbased_elimination_builder_and_solver.h`, around line 683):

  ```cpp
  if (dof_iterator->IsFixed()) dof_iterator->SetEquationId(--fix_id);
  else                         dof_iterator->SetEquationId(free_id++);
  BaseType::mEquationSystemSize = fix_id;   // number of FREE dofs
  ```

  So `mEquationSystemSize` equals the number of **free** dofs.
- **Consequence for a switch.** Flipping one dof from fixed to free (or back) changes the number of
  free dofs, hence changes `mEquationSystemSize`, hence changes the size of `A`, `Dx` and `b`, and
  renumbers the equation ids of *other* dofs too. The system genuinely must be rebuilt. If it is not,
  the next build writes into a matrix of the wrong size using stale equation ids — a silent
  corruption or a crash.
- **`SetUpSystem` runs once per solution step**, from `InitializeSolutionStep`, not per iteration. So
  a mid-iteration switch is exactly the case the base strategy never anticipated, and the case this
  step exists to handle.
- **The reshape-matrix flag must be on.** `ResizeAndInitializeVectors` only tolerates a changed
  system size when `GetReshapeMatrixFlag()` is `true` (`residualbased_elimination_builder_and_solver.h`
  line 736). If it is `false` and the size has changed, the builder throws
  `"The equation system size has changed during the simulation. This is not permitted."`. That flag is
  set from the strategy's `ReformDofSetAtEachStep` constructor argument (the base class forwards it via
  `SetReshapeMatrixFlag`). The Muskat project parameters already set `reform_dofs_at_each_step: true`,
  so this is satisfied in practice, but it is a hard prerequisite for `RebuildSystem()` under the
  elimination builder and the integration test must set it explicitly.

### Block builder and solver (`ResidualBasedBlockBuilderAndSolver`)

- **Every dof enters the system.** Both free and fixed dofs get an equation id, and
  `mEquationSystemSize` equals the **total** number of dofs.
- **Dirichlet is applied per build.** On every `BuildAndSolve` the builder re-applies the Dirichlet
  conditions to the assembled system (it zeroes the row/column and puts the prescribed value on the
  diagonal) rather than eliminating fixed dofs beforehand.
- **Consequence for a switch.** Because the size counts *all* dofs, a fixity flip does **not** change
  `mEquationSystemSize`. Equation ids do not shift, and `A`/`Dx`/`b` keep their size. The next
  `BuildAndSolve` simply applies the new fixity while assembling. In other words, under the block
  builder the switch is honoured **without any rebuild at all**, and `RebuildSystem()` reduces to a
  near no-op: `SetUpSystem` re-assigns the same ids, and `ResizeAndInitializeVectors` sees the same
  size and only zeroes the vectors.

### Why the implementation still does the full rebuild

The GeoMechanics solver defaults to the block builder (`block_builder: true` in
`geomechanics_solver.py`), so the prototype exercises the block path first, where a rebuild is not
strictly needed. But `RebuildSystem()` is written to be **correct for both** builders: the elimination
path *requires* the re-shape, and the block path *tolerates* it cheaply. A general, always-correct
rebuild is preferable to one that quietly assumes the block builder and breaks the moment someone
selects the elimination builder. The test therefore uses the **elimination** builder, because that is
the only builder under which the rebuild has an observable effect (a changed system size) to assert
against.

## A verified safety point: the system pointers are not swapped

`SolveSolutionStep` binds local references at the top of the method:

```cpp
TSystemMatrixType& rA  = *mpA;
TSystemVectorType& rDx = *mpDx;
TSystemVectorType& rb  = *mpb;
```

`RebuildSystem()` calls `ResizeAndInitializeVectors`, which could in principle reallocate and swap
`mpA`/`mpDx`/`mpb`, which would leave those references dangling. It does not, in this path: both
builders swap the pointer **only when it is null** (see
`residualbased_block_builder_and_solver.h` `ResizeAndInitializeVectors`, and the elimination
equivalent). By the time `RebuildSystem()` runs mid-`SolveSolutionStep`, the pointers are already
allocated, so the vectors are resized *in place* and the references stay valid. The rebuild is
therefore safe to call from inside the copied loop without re-binding `rA`/`rDx`/`rb`.

## Architecture

`RebuildSystem()` replaces the empty stub in
`applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`.
It mirrors the system-construction block of the base `InitializeSolutionStep`, then re-solves and
updates the database, mirroring the first-iteration build branch of `SolveSolutionStep`:

```cpp
void RebuildSystem()
{
    KRATOS_TRY

    auto  p_scheme             = this->GetScheme();
    auto  p_builder_and_solver = this->GetBuilderAndSolver();
    auto& r_model_part         = BaseType::GetModelPart();

    // A switch changed which WATER_PRESSURE dofs are fixed. Re-shape the system so the dof set,
    // the equation ids and the vector sizes match the new fixity. Under the elimination builder
    // this genuinely changes the system size; under the block builder it is a cheap no-op.
    p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
    p_builder_and_solver->SetUpSystem(r_model_part);
    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb, r_model_part);

    // Re-solve on the re-shaped system, mirroring the first-iteration build branch, so this
    // iteration's database already reflects the new boundary configuration rather than leaving
    // zeroed vectors for FinalizeNonLinIteration and the convergence criteria to see.
    TSystemMatrixType& rA  = *mpA;
    TSystemVectorType& rDx = *mpDx;
    TSystemVectorType& rb  = *mpb;
    TSparseSpace::SetToZero(rA);
    TSparseSpace::SetToZero(rDx);
    TSparseSpace::SetToZero(rb);
    p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
    this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

    KRATOS_CATCH("")
}
```

Notes:

- **It keeps `virtual`, not `override`.** The stub `virtual void RebuildSystem() {}` is declared in
  this same class, so there is no base declaration to override. The body is replaced in place, exactly
  as step 3 did for `UpdateSeepageBoundaryConditions`.
- **`SetUpDofSet` is included** even though a switch does not change the dof *list*. It is cheap, it
  makes the rebuild robust against future changes that do alter the dof set, and it matches the base
  `InitializeSolutionStep` sequence one-to-one, which keeps the method easy to read against its model.

## Seam interaction

`RebuildSystem()` is called from the two seams already present in `SolveSolutionStep` (step 2),
immediately after a switch is detected:

```cpp
any_switched = UpdateSeepageBoundaryConditions();
if (any_switched) RebuildSystem();
```

Seam 3 still forces `is_converged = false` on any iteration that switched, so a full further iteration
always follows. The re-solve inside `RebuildSystem()` is what makes the *current* iteration's database
consistent; the forced next iteration then re-evaluates convergence on the settled configuration.
Termination stays bounded by the existing `iteration_number++ < mMaxIterationNumber` guard.

## Testing

A focused C++ integration test isolates the rebuild from the switching criterion (which step 3 already
covers). It drives the fixity change directly rather than through `SwitchOneSeepageNode`, so the test
asserts one thing: that a fixity flip followed by `RebuildSystem()` re-shapes the system.

- **Builder:** `ResidualBasedEliminationBuilderAndSolver`, because it is the only builder under which
  a fixity flip changes the observable system size. Under the block builder the size would be constant
  and the assertion would be vacuous.
- **Reshape flag:** the builder's reshape-matrix flag must be set to `true` (via the strategy's
  `ReformDofSetAtEachStep` constructor argument, or directly with `SetReshapeMatrixFlag(true)`).
  Without it, `ResizeAndInitializeVectors` throws when the system size changes, and the rebuild cannot
  work. The test sets it deliberately and, ideally, also covers the negative case (a switch under a
  `false` flag surfaces the builder's error) if it is cheap to do so.
- **Model:** a minimal model part with a small number of `WATER_PRESSURE` dofs and enough elements to
  assemble a well-defined system, a GeoMechanics scheme, and the elimination builder and solver.
- **Access:** `RebuildSystem()` is protected. A minimal test-only subclass exposes it (a thin public
  `PublicRebuildSystem()` that forwards), so the test can call it without going through a full solve.
- **Assertions:** after establishing the system (dof set + system + vectors), record the equation
  system size; free one previously-fixed `WATER_PRESSURE` dof; call the rebuild; assert the equation
  system size and the sizes of `A`, `Dx` and `b` grew by one, and that the freed dof now has a valid
  free-range equation id.

The full suite must show zero failures and nothing previously passing now failing.

## Findings to record for #14637

- Whether the immediate re-solve inside `RebuildSystem()` is redundant with the forced next iteration
  (seam 3), which would let step 4 shrink to just the re-shape.
- Confirmation, from a real end-to-end run, that under the block builder the rebuild is indeed a
  no-op, so the block path never actually needed step 4.
- Whether re-shaping and re-solving on every switch is a noticeable cost when many nodes switch over a
  solution step.

## Build and test commands

```powershell
kp config                                    # only after adding a new .cpp test file
kp build
kp test -- --gtest_filter="*Seepage*"
kp test                                      # full suite
```



