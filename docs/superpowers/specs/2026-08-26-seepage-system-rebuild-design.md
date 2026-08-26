# Seepage Boundary Prototype — Step 4: Forcing a Rebuild After a Switch — Design

**Issue:** #14672 (prototype for a seepage boundary), bullet 4 of 4.
**Depends on:** step 2, `GeoSeepageNewtonRaphsonStrategy`
(`docs/superpowers/specs/2026-08-25-seepage-newton-raphson-strategy-design.md`), and step 3, the
switching criterion (`docs/superpowers/specs/2026-08-25-seepage-switching-criterion-design.md`).
**Date:** 2026-08-26

## Goal

Make `GeoSeepageNewtonRaphsonStrategy` re-solve with the correct boundary configuration when a seepage
node has just switched between a Dirichlet boundary (`WATER_PRESSURE` fixed at zero) and a zero-flux
Neumann boundary (`WATER_PRESSURE` free).

This is the final stubbed piece of the prototype. It turns out to need far less than the issue
anticipated: **no explicit system rebuild at all.** Under the block builder and solver, forcing the
existing build-and-solve step to reassemble the stiffness matrix is enough, because that step already
applies the Dirichlet conditions for the current fixity. So step 4 removes the `RebuildSystem()` stub
entirely and replaces it with a one-line steer of an existing flag.

## Scope

- **In scope:** move the switching seam to just before the build-and-solve branch, and on a switch set
  `mStiffnessMatrixIsBuilt = false` so that branch reassembles the stiffness matrix with the new
  fixity. Remove the `RebuildSystem()` method. Add a focused C++ test.
- **Block builder and solver only.** The GeoMechanics solver defaults to the block builder
  (`block_builder: true` in `geomechanics_solver.py`), and the prototype targets exactly that path.
- **Out of scope (deferred to #14637):** any support for the elimination builder and solver (which
  would need a genuine system resize, see below); wiring the strategy into `geomechanics_solver.py`;
  building a full seepage-boundary Muskat validation case; and judging the switching sign against a
  reference solution.

## Why no rebuild is needed: the two builders and solvers, contrasted

A seepage switch never adds or removes a degree of freedom. Every seepage node always owns exactly one
`WATER_PRESSURE` dof; the switch only flips that dof between **fixed** (Dirichlet) and **free**
(Neumann). What that flip costs depends entirely on the builder and solver, and the two differ
sharply.

### Block builder and solver (`ResidualBasedBlockBuilderAndSolver`) — the target

- **Every dof enters the system.** Free and fixed dofs alike get an equation id, and the system size
  `mEquationSystemSize` equals the **total** number of dofs.
- **Dirichlet is applied on every build.** Each `BuildAndSolve` assembles the full system and then
  re-applies the Dirichlet conditions for the *current* fixity (it zeroes the fixed dof's row and
  column and puts the prescribed value on the diagonal).
- **Consequence for a switch.** Because the size counts all dofs, a fixity flip does **not** change
  `mEquationSystemSize`; equation ids do not shift; `A`, `Dx` and `b` keep their sizes. The only thing
  required to honour the switch is to make the next `BuildAndSolve` actually reassemble the stiffness
  matrix, so the new fixity is baked in. That is what this step does — nothing more.

### Elimination builder and solver (`ResidualBasedEliminationBuilderAndSolver`) — deferred

- **Only free dofs enter the system.** Fixed dofs are eliminated and given no equation. `SetUpSystem`
  numbers free dofs first and sets `mEquationSystemSize` to the count of **free** dofs
  (`residualbased_elimination_builder_and_solver.h`, around line 683).
- **Consequence for a switch.** Flipping one dof changes the number of free dofs, hence the system
  size, hence the sizes of `A`, `Dx`, `b`, and renumbers other dofs' equation ids. Honouring a
  mid-iteration switch there requires a real `SetUpSystem` + `ResizeAndInitializeVectors` rebuild, and
  `ResizeAndInitializeVectors` even throws `"The equation system size has changed during the
  simulation. This is not permitted."` unless the reshape-matrix flag is on. None of that is needed
  for the block builder, so it is left to #14637.

The one-line block-builder steer is correct for the prototype and keeps the code minimal; supporting
the elimination builder later is an additive change, not a rework.

## Design: move the seam and set the flag

The change lives entirely in
`applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`,
inside the copied `SolveSolutionStep`. The guiding constraint, carried over from steps 2 and 3, is
that `SolveSolutionStep` must stay a faithful copy of the core method: a diff against core should show
**only** the seams. Step 4 keeps that property. The deviations from core become exactly:

1. **Seam 0** (unchanged): `bool any_switched = false;` declared once, before the first-iteration
   block.
2. **Seam 1** (moved): the switch decision moves from *after* `UpdateDatabase` to *just before* the
   build-and-solve branch, in both the first-iteration block and the iteration cycle. Its action
   changes from `RebuildSystem()` to setting the flag:

   ```cpp
   // ---- SEEPAGE SEAM 1: decide the switch, then force a stiffness rebuild -------------
   any_switched = UpdateSeepageBoundaryConditions();
   if (any_switched) this->SetStiffnessMatrixIsBuilt(false);
   // ------------------------------------------------------------------------------------

   // Function to perform the building and the solving phase.
   if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
       // ... unchanged core BuildAndSolve branch, which applies the new Dirichlet fixity ...
   }
   ```

3. **Seam 3** (unchanged): after each `PostCriteria`, `if (any_switched) is_converged = false;`.

The old "seam 2" (`if (any_switched) RebuildSystem();`) and the `RebuildSystem()` method are deleted.

### Why the seam moves before the build

With the seam *before* the build, the switch is applied in the **same** iteration that decides it:
the build reassembles with the new fixity and the solve already reflects the new boundary
configuration. The decision reads the current database — the previous iteration's solution (or the
initial condition on the very first pass) — which is the correct, latest available information.

### Why the flag, given the default rebuild level

The Newton-Raphson constructor calls `SetRebuildLevel(2)` (rebuild the stiffness matrix every
iteration), so under default settings `BuildAndSolve` already runs each iteration and the flag is
belt-and-suspenders. Setting `mStiffnessMatrixIsBuilt = false` makes the behaviour correct also at
lower rebuild levels, where an iteration would otherwise skip the LHS rebuild and fail to re-apply the
new Dirichlet fixity. `SetStiffnessMatrixIsBuilt(bool)` is a public method of the base
`ImplicitSolvingStrategy`, so no new access plumbing is needed.

### Why seam 3 is still required

Even though the switch is now applied within the same iteration, convergence must not be declared on
an iteration that switched: the freshly settled solution might itself warrant a further switch, and
only an iteration that makes **no** switch and is numerically converged represents a true fixed point.
Seam 3 forces one more iteration after any switch; the following iteration's seam re-checks and, if
nothing switches, convergence stands. Termination stays bounded by the existing
`iteration_number++ < mMaxIterationNumber` guard, so oscillating neighbours cannot loop forever.

## What this step explicitly does *not* touch

- It does not call `SetUpDofSet`, `SetUpSystem` or `ResizeAndInitializeVectors`. The system vectors and
  their pointers are left alone, which sidesteps any question of stale `rA`/`rDx`/`rb` references in
  the copied loop.
- It does not add any member variable. `mSeepageNodes` (step 3) is the only state the strategy holds.

## Testing

A focused C++ test drives a minimal block-builder Newton-Raphson solve through a small test subclass
that forces a single, deterministic switch, and checks the switch is honoured. Driving a real solve is
the honest way to exercise the change, because the seam lives inline in the copied `SolveSolutionStep`
and is not reachable through a hook.

- **Builder:** `ResidualBasedBlockBuilderAndSolver`, matching the target and the constant-size
  behaviour the design relies on.
- **Model:** a minimal model part with a small number of `WATER_PRESSURE` dofs and enough element(s)
  to assemble a solvable steady-state system, a GeoMechanics scheme, and a residual convergence
  criterion.
- **Forced switch:** a test-only subclass overrides `UpdateSeepageBoundaryConditions()` to free a
  chosen node and return `true` on its first call, then `false` afterwards. This makes exactly one
  switch happen deterministically, independent of the step-3 criterion (already tested on its own).
- **Assertions:**
  - the solve performs at least one extra iteration after the switch (seam 3 suppressed premature
    convergence); and
  - the freed node ends the solve as a free dof carrying its free-solve pressure rather than the fixed
    zero, proving the block builder reassembled with the new fixity.
- **Optional lower-rebuild-level check:** to show the flag itself matters (not just the default level-2
  rebuild), a variant may set the rebuild level to 1 and assert the freed node is still honoured — this
  is the case the `SetStiffnessMatrixIsBuilt(false)` line exists for.

The full suite must show zero failures, with nothing previously passing now failing.

## Findings to record for #14637

- Confirmation from a real end-to-end run that the block-builder steer honours switches correctly
  across a full solution step.
- Whether the elimination builder is ever wanted; if so, it needs the genuine
  `SetUpSystem` + `ResizeAndInitializeVectors` rebuild described above, plus `reform_dofs_at_each_step`.
- Whether reassembling the stiffness matrix on every switch is a noticeable cost when many nodes
  switch across a step.

## Build and test commands

```powershell
kp config                                    # only after adding a new .cpp test file
kp build
kp test -- --gtest_filter="*Seepage*"
kp test                                      # full suite
```

