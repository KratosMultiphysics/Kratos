# Seepage Newton-Raphson Strategy — Step 2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `GeoSeepageNewtonRaphsonStrategy` to the GeoMechanicsApplication that reimplements the core Newton-Raphson `SolveSolutionStep` with three seams, so that steps 3 and 4 can decide seepage boundary switches, rebuild the system, and suppress convergence when a switch has just happened.

**Architecture:** A header-only class template deriving from the core `ResidualBasedNewtonRaphsonStrategy`. `SolveSolutionStep` is copied verbatim from core and modified only at three marked seams. In this step the two new hook methods are stubs that do nothing, so the strategy behaves exactly like its base class.

**Tech Stack:** C++20, Kratos Multiphysics, pybind11, CMake + Ninja.

**Spec:** `docs/superpowers/specs/2026-08-25-seepage-newton-raphson-strategy-design.md`

## Global Constraints

- Target application: `applications/GeoMechanicsApplication`. All paths are relative to `C:/checkouts/KratosProjects/dev`.
- C++20. New files carry the standard Kratos GeoMechanics header comment block.
- Derive from the **core** `ResidualBasedNewtonRaphsonStrategy`, not from `GeoMechanicsNewtonRaphsonStrategy`. The constructor therefore takes **no** `Parameters` argument.
- **This step ships without tests, by explicit decision.** The strategy is behaviourally identical to its base here. Verification is: it compiles, the full existing suite still passes, and the class is importable from Python.
- `SolveSolutionStep` must be a verbatim copy of core except the three seams. Do not reformat, rename variables, or "improve" the copied code. A diff against core must show only the seams.
- Do not modify `applications/GeoMechanicsApplication/CMakeLists.txt`. Strategies are header-only and `custom_python/*.cpp` is already globbed.
- Out of scope: collecting the seepage conditions and the switching criterion (step 3), the actual rebuild (step 4), and wiring into `geomechanics_solver.py` (step 3).

## Critical C++ Note: Dependent Base Classes

`GeoSeepageNewtonRaphsonStrategy` is a **class template** whose base is also a template. C++ does **not** look up unqualified names in dependent base classes. Every inherited member used in the copied loop must therefore be either:

- declared with a `using MotherType::member;` line (done for data members in Task 1), or
- called through `this->` (done for member functions).

If you skip this, you get errors like `'mpA': undeclared identifier`. The existing `GeoMechanicsNewtonRaphsonStrategy` uses the same `using MotherType::mpA;` technique.

## Build and Test Commands

Build (no `kp config` needed — no new `.cpp` in Task 1 or 2; Task 3 modifies an existing `.cpp`):

```powershell
kp build
```

Full C++ regression suite (baseline before this plan: **1128 passing**):

```powershell
kp test
```

---

### Task 1: Create the strategy skeleton with its two stub hooks

Deliverable: a compiling, registered-nowhere strategy class that inherits the core `SolveSolutionStep` unchanged. No loop copy yet, so a reviewer can judge the class shape in isolation.

**Files:**
- Create: `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`

**Interfaces:**
- Consumes: `Kratos::ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>` from `solving_strategies/strategies/residualbased_newton_raphson_strategy.h`.
- Produces:
  - `template <class TSparseSpace, class TDenseSpace, class TLinearSolver> class Kratos::GeoSeepageNewtonRaphsonStrategy`
  - Constructor `(ModelPart&, typename TSchemeType::Pointer, typename TConvergenceCriteriaType::Pointer, typename TBuilderAndSolverType::Pointer, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false)`
  - `protected virtual bool UpdateSeepageBoundaryConditions()` — returns `false` in this step. Step 3 overrides/implements it. Returns `true` when at least one `GeoSeepageCondition` changed its boundary type.
  - `protected virtual void RebuildSystem()` — empty in this step. Step 4 implements it.
  - `[[nodiscard]] std::string Info() const override` returning `"GeoSeepageNewtonRaphsonStrategy"`

- [ ] **Step 1: Create the header file**

Create `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`:

```cpp
// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

// A Newton-Raphson strategy that can switch seepage boundary conditions between Dirichlet and
// Neumann while iterating.
//
// Applying a switch is the condition's job: the scheme calls
// Condition::InitializeNonLinearIteration on every non-linear iteration, and GeoSeepageCondition
// fixes or frees its nodes there. This strategy owns the two things the condition cannot do:
// deciding when to switch, and making sure the solver does not declare convergence on an iteration
// whose boundary configuration has just changed underneath it.
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoSeepageNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoSeepageNewtonRaphsonStrategy);

    using BaseType   = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;

    // The base class is a dependent template base, so every inherited data member used by the
    // copied SolveSolutionStep must be pulled into scope explicitly.
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mMaxIterationNumber;
    using MotherType::mNonconvergedSolutionsMatrix;
    using MotherType::mpA;  // Tangent matrix
    using MotherType::mpb;  // Residual vector of iteration i
    using MotherType::mpConvergenceCriteria;
    using MotherType::mpDx; // Delta x of iteration i
    using MotherType::mStoreNonconvergedSolutionsFlag;
    using MotherType::mUseOldStiffnessInFirstIteration;

    GeoSeepageNewtonRaphsonStrategy(ModelPart&                                 rModelPart,
                                    typename TSchemeType::Pointer              pScheme,
                                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                    typename TBuilderAndSolverType::Pointer    pNewBuilderAndSolver,
                                    int  MaxIterations          = 30,
                                    bool CalculateReactions     = false,
                                    bool ReformDofSetAtEachStep = false,
                                    bool MoveMeshFlag           = false)
        : MotherType(rModelPart,
                     pScheme,
                     pNewConvergenceCriteria,
                     pNewBuilderAndSolver,
                     MaxIterations,
                     CalculateReactions,
                     ReformDofSetAtEachStep,
                     MoveMeshFlag)
    {
    }

    [[nodiscard]] std::string Info() const override { return "GeoSeepageNewtonRaphsonStrategy"; }

protected:
    // Decides whether any seepage condition must change its boundary type, and applies that
    // decision to the conditions' Properties. Returns true if at least one condition changed.
    //
    // Step 3 of the prototype implements this. Until then nothing ever switches, which makes this
    // strategy behave exactly like its base class.
    virtual bool UpdateSeepageBoundaryConditions() { return false; }

    // Re-establishes the degree of freedom set and the system vectors after a switch changed which
    // degrees of freedom are free.
    //
    // Step 4 of the prototype implements this. Note that under a block builder and solver this may
    // turn out to be unnecessary, because that builder gives every degree of freedom an equation id
    // and re-applies the Dirichlet conditions on every build.
    virtual void RebuildSystem() {}
}; // Class GeoSeepageNewtonRaphsonStrategy

} // namespace Kratos
```

- [ ] **Step 2: Build**

```powershell
kp build
```

Expected: builds cleanly. Nothing includes the new header yet, so this only proves the file is syntactically valid once it is first included — which happens in Task 3. To get compile feedback now, temporarily add `#include "custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp"` to `applications/GeoMechanicsApplication/custom_python/add_custom_strategies_to_python.cpp`, build, then remove it again. Task 3 adds it permanently.

- [ ] **Step 3: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
git commit -m "Add GeoSeepageNewtonRaphsonStrategy skeleton with seepage hooks"
```

---

### Task 2: Copy `SolveSolutionStep` from core and insert the three seams

Deliverable: the reimplemented loop. Behaviour is still identical to core, because `UpdateSeepageBoundaryConditions()` always returns `false`.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`

**Interfaces:**
- Consumes: `UpdateSeepageBoundaryConditions()` and `RebuildSystem()` from Task 1.
- Produces: `bool SolveSolutionStep() override`.

**Source of the copy:** `kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h`, lines 919–1105.

**The three seams:**

| Seam | Location | Purpose |
| --- | --- | --- |
| 0 | next to `residual_is_updated` | declare `bool any_switched = false;` once, so it is not redeclared and shadowed in the loop body |
| 1 & 2 | right after each `UpdateDatabase(...)` | decide switches, then rebuild if anything switched |
| 3 | right after each `PostCriteria(...)` block | force another iteration when a switch just happened |

Seams 1, 2 and 3 each appear **twice**: once in the first-iteration block and once in the iteration cycle. Missing the first-iteration copy would let a switch made on iteration 1 be reported as converged.

- [ ] **Step 1: Add the `SolveSolutionStep` override**

In `geo_seepage_newton_raphson_strategy.hpp`, insert this method into the `public:` section, directly below the `Info()` method:

```cpp
    // NOTE: This is a deliberate copy of
    // ResidualBasedNewtonRaphsonStrategy::SolveSolutionStep
    // (kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h, lines
    // 919-1105). The ONLY intended differences are the blocks marked "SEEPAGE SEAM". Please keep it
    // that way: a diff against the core method should show nothing else. It is copied rather than
    // hooked because seam 3 has to run after PostCriteria, and no existing virtual method of the
    // base class sits at that point.
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart&                             r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer          p_scheme     = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        auto&                                  r_dof_set    = p_builder_and_solver->GetDofSet();
        std::vector<Vector>                    NonconvergedSolutions;

        if (mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            this->GetCurrentSolution(r_dof_set, initial);
            NonconvergedSolutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        // initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number                      = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated                           = false;

        // ---- SEEPAGE SEAM 0 ------------------------------------------------------------------
        // Declared once here, and only re-assigned at the seams below. Declaring it inside both
        // the first-iteration block and the loop body would shadow one with the other.
        bool any_switched = false;
        // --------------------------------------------------------------------------------------

        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (mUseOldStiffnessInFirstIteration) {
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(
                    p_scheme, r_model_part, rA, rDx, rb, BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx); // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        this->EchoInfo(iteration_number);

        // Updating the results stored in the database
        this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        // ---- SEEPAGE SEAMS 1 AND 2 -----------------------------------------------------------
        // The solution has just been updated, so pressures and fluxes are current. Decide now; the
        // next iteration's InitializeNonLinIteration is what applies the new fixity to the nodes.
        any_switched = UpdateSeepageBoundaryConditions();
        if (any_switched) RebuildSystem();
        // --------------------------------------------------------------------------------------

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            Vector first;
            this->GetCurrentSolution(r_dof_set, first);
            NonconvergedSolutions.push_back(first);
        }

        if (is_converged) {
            if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        // ---- SEEPAGE SEAM 3 ------------------------------------------------------------------
        // A switch means this iteration's residual describes the OLD boundary configuration, so its
        // convergence verdict is meaningless. Force at least one more iteration.
        if (any_switched) is_converged = false;
        // --------------------------------------------------------------------------------------

        // Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < mMaxIterationNumber) {
            // setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            // call the linear system solver to find the correction mDx for the
            // it is not called if there is no system to solve
            if (TSparseSpace::Size(rDx) != 0) {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false) {
                    if (this->GetKeepSystemConstantDuringIterations() == false) {
                        // A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    } else {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                } else {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            this->EchoInfo(iteration_number);

            // Updating the results stored in the database
            this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            // ---- SEEPAGE SEAMS 1 AND 2 -------------------------------------------------------
            any_switched = UpdateSeepageBoundaryConditions();
            if (any_switched) RebuildSystem();
            // ----------------------------------------------------------------------------------

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (mStoreNonconvergedSolutionsFlag == true) {
                Vector ith;
                this->GetCurrentSolution(r_dof_set, ith);
                NonconvergedSolutions.push_back(ith);
            }

            residual_is_updated = false;

            if (is_converged == true) {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }

            // ---- SEEPAGE SEAM 3 --------------------------------------------------------------
            if (any_switched) is_converged = false;
            // ----------------------------------------------------------------------------------
        }

        // plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            this->MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("GeoSeepageNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIterationNumber << " iterations" << std::endl;
        }

        // calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            mNonconvergedSolutionsMatrix = Matrix(r_dof_set.size(), NonconvergedSolutions.size());
            for (std::size_t i = 0; i < NonconvergedSolutions.size(); ++i) {
                block_for_each(r_dof_set, [&](const auto& r_dof) {
                    mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) =
                        NonconvergedSolutions[i](r_dof.EquationId());
                });
            }
        }

        return is_converged;
    }
```

Note two intentional, harmless deviations from the literal core text, both forced by the language rather than by design:

1. Core writes `SparseSpaceType::Size(rDx)`; the copy writes `TSparseSpace::Size(rDx)`. Same type, but the template parameter name is what is in scope here.
2. Core calls `GetScheme()`, `EchoInfo(...)`, `UpdateDatabase(...)` etc. unqualified. The copy calls them through `this->`, which is required for dependent base classes.

Core's dead `if (residual_is_updated == false) { /* commented-out body */ }` block is omitted, since its body is entirely commented out in core. `residual_is_updated` is still assigned, matching core.

- [ ] **Step 2: Build**

Temporarily add the include to `add_custom_strategies_to_python.cpp` (Task 3 makes it permanent), then:

```powershell
kp build
```

Expected: builds cleanly. If you see `'mpA': undeclared identifier` or similar, a `using MotherType::...;` line from Task 1 is missing — re-read the "Critical C++ Note" above.

- [ ] **Step 3: Verify the copy against core**

Open `kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h` at line 919 side by side with the new method. Confirm the only differences are the seam blocks and the two deviations listed in Step 1. This manual diff is the substitute for the test we are not writing.

- [ ] **Step 4: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
git commit -m "Reimplement SolveSolutionStep with seepage switching seams"
```

---

### Task 3: Register the strategy with Python

Deliverable: `GeoSeepageNewtonRaphsonStrategy` is constructible from Python, which is how steps 3 and 4 will drive the Muskat case.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_python/add_custom_strategies_to_python.cpp` (include near line 28; alias near line 82; registration near line 138)

**Interfaces:**
- Consumes: `Kratos::GeoSeepageNewtonRaphsonStrategy` from Tasks 1 and 2.
- Produces: Python class `KratosMultiphysics.GeoMechanicsApplication.GeoSeepageNewtonRaphsonStrategy`, constructible with `(ModelPart, Scheme, ConvergenceCriteria, BuilderAndSolver, int, bool, bool, bool)`.

Note the deliberate difference from its siblings: **no `Parameters&` argument**, because this strategy derives from the core class rather than from `GeoMechanicsNewtonRaphsonStrategy`.

- [ ] **Step 1: Add the include**

In `add_custom_strategies_to_python.cpp`, in the strategies include block (the one containing `geo_mechanics_newton_raphson_strategy.hpp`), add:

```cpp
#include "custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp"
```

- [ ] **Step 2: Add the type alias**

Directly below the `GeoMechanicsNewtonRaphsonErosionProcessStrategyType` alias (around line 82):

```cpp
    using GeoSeepageNewtonRaphsonStrategyType =
        GeoSeepageNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
```

- [ ] **Step 3: Add the registration**

Directly below the `GeoMechanicsNewtonRaphsonErosionProcessStrategy` registration (which ends around line 138):

```cpp
    py::class_<GeoSeepageNewtonRaphsonStrategyType, typename GeoSeepageNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType>(
        m, "GeoSeepageNewtonRaphsonStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer,
                      BuilderAndSolverType::Pointer, int, bool, bool, bool>());
```

- [ ] **Step 4: Build**

```powershell
kp build
```

Expected: builds cleanly.

- [ ] **Step 5: Verify the class is importable from Python**

```powershell
python -c "import KratosMultiphysics; import KratosMultiphysics.GeoMechanicsApplication as Geo; print(Geo.GeoSeepageNewtonRaphsonStrategy)"
```

Expected: prints something like `<class 'KratosMultiphysics.GeoMechanicsApplication.GeoSeepageNewtonRaphsonStrategy'>`.

If Python cannot find `KratosMultiphysics`, the installed package directory is not on `PYTHONPATH`. Locate it under `bin/FullDebug` and prepend it, for example:

```powershell
$env:PYTHONPATH = "C:\checkouts\KratosProjects\dev\bin\FullDebug;$env:PYTHONPATH"
```

- [ ] **Step 6: Run the full C++ suite to check for regressions**

```powershell
kp test
```

Expected: `[  PASSED  ] 1128 tests.` — the same count as before this plan, since nothing existing changed behaviour.

- [ ] **Step 7: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_python/add_custom_strategies_to_python.cpp
git commit -m "Register GeoSeepageNewtonRaphsonStrategy with Python"
```

---

## Findings to Record for the Follow-up Issue (#14637)

Fill this in during implementation:

- Whether the copied loop needed any deviation from core beyond the two listed in Task 2. Each extra deviation is a maintenance cost to report.
- Whether `mStoreNonconvergedSolutionsFlag` and the other protected members were all reachable as expected, or whether any needed a workaround.
- Confirm the strategy is a drop-in replacement: with both hooks stubbed, a case run with `GeoSeepageNewtonRaphsonStrategy` should produce results identical to `GeoMechanicsNewtonRaphsonStrategy`. Worth spot-checking once step 3 wires it into `geomechanics_solver.py`.

