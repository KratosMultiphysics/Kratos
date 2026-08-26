# Seepage System Rebuild — Step 4 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `GeoSeepageNewtonRaphsonStrategy` honour a seepage switch by forcing the existing build-and-solve step to reassemble the stiffness matrix with the new fixity, instead of doing a separate system rebuild.

**Architecture:** A small, inline edit to the copied `SolveSolutionStep` in the header-only strategy. The switching seam moves from *after* `UpdateDatabase` to *just before* the build-and-solve branch, and its action changes from calling `RebuildSystem()` to `this->SetStiffnessMatrixIsBuilt(false)`. The `RebuildSystem()` method is deleted. Under the block builder and solver (the GeoMechanics default) this is sufficient, because that builder keeps a constant system size and re-applies the Dirichlet conditions on every build.

**Tech Stack:** C++20, Kratos Multiphysics, CMake + Ninja.

**Spec:** `docs/superpowers/specs/2026-08-26-seepage-system-rebuild-design.md`

## Global Constraints

- Target application: `applications/GeoMechanicsApplication`. All paths are relative to `C:/checkouts/KratosProjects/dev`.
- C++20. The only file changed is the existing header `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`.
- **Block builder and solver only.** The elimination builder (which would need a genuine `SetUpSystem` + `ResizeAndInitializeVectors` resize) is out of scope and deferred to #14637.
- **`SolveSolutionStep` must stay a faithful copy of core.** A diff against `kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h` must show only the marked `SEEPAGE SEAM` blocks. Do not reformat, rename, or otherwise touch the copied core lines.
- **No new tests, by explicit decision.** Verification for this step is: it builds, a diff against core shows only the seams, and the full existing suite still passes with no regressions. The switching behaviour is proven later on the Muskat case (#14637).
- No new `.cpp` file, so no `kp config` is needed. Build with `kp build`.

## Build and test commands

```powershell
kp build     # header-only change; no new .cpp, so no kp config
kp test      # full C++ suite; must stay green with no regressions
```

---

### Task 1: Rewire the seepage seams to force a stiffness rebuild

Deliverable: the strategy decides the switch just before each build, forces a stiffness rebuild on a switch, and no longer has a `RebuildSystem()` method. Behaviour under the block builder: a switched node's new fixity is assembled into the very next solve.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`

**Interfaces:**
- Consumes: `UpdateSeepageBoundaryConditions()` (step 3, unchanged) and `ImplicitSolvingStrategy::SetStiffnessMatrixIsBuilt(bool)` (base class, public).
- Produces: no new public interface. Removes the `protected virtual void RebuildSystem()` method.

There are two structurally identical regions to edit: the **first-iteration block** (before the `while` loop) and the **iteration cycle** (inside the `while` loop). Each has one seam to move and one old seam to delete. Do both.

- [ ] **Step 1: Move seam 1 in the first-iteration block (insert before the build)**

In the first-iteration block, find:

```cpp
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
```

Replace it with (insert the seam between `PreCriteria` and the build comment):

```cpp
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // ---- SEEPAGE SEAM 1: decide the switch, then force a stiffness rebuild ----------------
        // Deciding before the build means this iteration's BuildAndSolve reassembles with the new
        // fixity, so the block builder applies the switched node's Dirichlet/Neumann state at once.
        any_switched = UpdateSeepageBoundaryConditions();
        if (any_switched) this->SetStiffnessMatrixIsBuilt(false);
        // --------------------------------------------------------------------------------------

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
```

- [ ] **Step 2: Delete the old post-`UpdateDatabase` seam in the first-iteration block**

In the first-iteration block, find:

```cpp
        // Updating the results stored in the database
        this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        // ---- SEEPAGE SEAMS 1 AND 2 -----------------------------------------------------------
        // The solution has just been updated, so pressures and fluxes are current. Decide now; the
        // next iteration's InitializeNonLinIteration is what applies the new fixity to the nodes.
        any_switched = UpdateSeepageBoundaryConditions();
        if (any_switched) RebuildSystem();
        // --------------------------------------------------------------------------------------

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
```

Replace it with (seam removed):

```cpp
        // Updating the results stored in the database
        this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
```

- [ ] **Step 3: Update the seam-3 comment in the first-iteration block**

The switch is now applied *before* the build, so the old rationale ("residual describes the OLD boundary configuration") is no longer accurate. Find:

```cpp
        // ---- SEEPAGE SEAM 3 ------------------------------------------------------------------
        // A switch means this iteration's residual describes the OLD boundary configuration, so its
        // convergence verdict is meaningless. Force at least one more iteration.
        if (any_switched) is_converged = false;
        // --------------------------------------------------------------------------------------
```

Replace it with:

```cpp
        // ---- SEEPAGE SEAM 3 ------------------------------------------------------------------
        // A switch was just applied and solved for. The settled solution might itself warrant a
        // further switch, so convergence may only be declared on an iteration that switches
        // nothing. Force at least one more iteration here.
        if (any_switched) is_converged = false;
        // --------------------------------------------------------------------------------------
```

- [ ] **Step 4: Move seam 1 in the iteration cycle (insert before the build)**

Inside the `while` loop, find:

```cpp
            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            // call the linear system solver to find the correction mDx for the
            // it is not called if there is no system to solve
            if (TSparseSpace::Size(rDx) != 0) {
```

Replace it with (insert the seam between `PreCriteria` and the solver comment):

```cpp
            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            // ---- SEEPAGE SEAM 1: decide the switch, then force a stiffness rebuild ------------
            any_switched = UpdateSeepageBoundaryConditions();
            if (any_switched) this->SetStiffnessMatrixIsBuilt(false);
            // ----------------------------------------------------------------------------------

            // call the linear system solver to find the correction mDx for the
            // it is not called if there is no system to solve
            if (TSparseSpace::Size(rDx) != 0) {
```

- [ ] **Step 5: Delete the old post-`UpdateDatabase` seam in the iteration cycle**

Inside the `while` loop, find:

```cpp
            // Updating the results stored in the database
            this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            // ---- SEEPAGE SEAMS 1 AND 2 -------------------------------------------------------
            any_switched = UpdateSeepageBoundaryConditions();
            if (any_switched) RebuildSystem();
            // ----------------------------------------------------------------------------------

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
```

Replace it with (seam removed):

```cpp
            // Updating the results stored in the database
            this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
```

- [ ] **Step 6: Delete the `RebuildSystem()` method**

In the `protected:` section, find and delete this entire method (it is no longer called):

```cpp
    // Re-establishes the degree of freedom set and the system vectors after a switch changed which
    // degrees of freedom are free.
    //
    // Step 4 of the prototype implements this. Note that under a block builder and solver this may
    // turn out to be unnecessary, because that builder gives every degree of freedom an equation id
    // and re-applies the Dirichlet conditions on every build.
    virtual void RebuildSystem() {}
```

Leave `UpdateSeepageBoundaryConditions()` and the `private:` `mSeepageNodes` member exactly as they are.

- [ ] **Step 7: Build**

```powershell
kp build
```

Expected: builds cleanly. The header is instantiated through the Python registration (from step 2), so the template body is fully type-checked. If you see an error about `SetStiffnessMatrixIsBuilt`, confirm it is called as `this->SetStiffnessMatrixIsBuilt(false)` — it is a public method of the dependent base `ImplicitSolvingStrategy`, so the `this->` qualification is required.

- [ ] **Step 8: Verify the diff against core shows only the seams**

Open `kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h` (the `SolveSolutionStep` starting around line 919) side by side with the edited method. Confirm the only differences are:
- seam 0: the one-time `bool any_switched = false;` declaration;
- seam 1: the `UpdateSeepageBoundaryConditions()` + `SetStiffnessMatrixIsBuilt(false)` block, now placed immediately before each build branch (first-iteration block and iteration cycle);
- seam 3: the `if (any_switched) is_converged = false;` block after each `PostCriteria`;
- plus the two harmless deviations already documented in step 2 (`TSparseSpace::Size` instead of `SparseSpaceType::Size`, and `this->` on inherited calls).

There must be **no** remaining reference to `RebuildSystem` anywhere in the file.

- [ ] **Step 9: Run the full C++ suite to check for regressions**

```powershell
kp test
```

Expected: zero failures, and no test that passed before this change now fails. Nothing wired the new seam into a running solve yet (that is #14637), so the suite total should be unchanged from the step-3 baseline.

- [ ] **Step 10: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
git commit -m "Force stiffness rebuild on a seepage switch instead of rebuilding the system"
```

---

## Findings to Record for the Follow-up Issue (#14637)

Fill this in during implementation and, later, during the Muskat validation:

- Confirmation from a real end-to-end run that the block-builder steer honours switches correctly across a full solution step.
- Whether the elimination builder is ever wanted; if so, it needs a genuine `SetUpSystem` + `ResizeAndInitializeVectors` rebuild plus `reform_dofs_at_each_step`.
- Whether reassembling the stiffness matrix on every switch is a noticeable cost when many nodes switch across a step.
- Whether deciding the switch before the first solve (on the initial condition) ever produces a worse starting point than deciding after the first solve.

