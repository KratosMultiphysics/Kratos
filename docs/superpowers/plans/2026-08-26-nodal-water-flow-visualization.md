# NODAL_WATER_FLOW Visualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose the per-node assembled WATER_PRESSURE right-hand-side as a registered nodal variable `NODAL_WATER_FLOW`, computed once per converged step in the seepage strategy and writable to GiD output, so the seepage sign convention can be visualized and verified.

**Architecture:** Register a new nodal `Variable<double>` (`NODAL_WATER_FLOW`) mirroring `HYDRAULIC_DISCHARGE`. Add a small, unit-tested helper `AssignNodalWaterFlows` to `SeepageBoundaryUtilities` that stores a `NodalFlowMap` onto the historical database (zeroing all nodes first). Call it from a new `FinalizeSolutionStep` override in `GeoSeepageNewtonRaphsonStrategy`. Add the variable to the solver/workflow variable lists and to the GiD output writer's nodal map.

**Tech Stack:** C++20, Kratos Multiphysics, GeoMechanicsApplication, GTest (`KratosGeoMechanicsCoreTest`), pybind11, Python.

## Global Constraints

- Variable name is exactly `NODAL_WATER_FLOW`, type `double`, nodal (historical/solution-step).
- Sign convention is unchanged: positive value = outflow (matches `ShouldReleaseToNeumann`).
- Only `newton_raphson_with_seepage` runs populate the variable; other strategies leave it at default `0.0`.
- Registration must mirror the existing `HYDRAULIC_DISCHARGE` pattern in every location.
- All new C++ files carry the standard GeoMechanics license header block.
- C++ standard is C++20; namespace for the helper is `Kratos::Geo::SeepageBoundaryUtilities`.

---

### Task 1: Register the `NODAL_WATER_FLOW` variable

Declares, creates, and registers the new nodal variable everywhere the application and its C++ test kernel need it. No behavior yet — this task's deliverable is that the variable exists and the code still builds. The unit test in Task 2 depends on the test-suite registration done here.

**Files:**
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.h` (near line 31)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp` (near line 27)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application.cpp` (near line 386)
- Modify: `applications/GeoMechanicsApplication/custom_python/geo_mechanics_python_application.cpp` (near line 57)
- Modify: `applications/GeoMechanicsApplication/GeoMechanicsApplication.py` (near line 93)
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.cpp` (near line 41)

**Interfaces:**
- Consumes: nothing.
- Produces: `Kratos::NODAL_WATER_FLOW` — a `Variable<double>` available in C++ (`geo_mechanics_application_variables.h`), in the Python module as `KratosGeo.NODAL_WATER_FLOW` and `GeoMechanicsApplication.NODAL_WATER_FLOW`, and registered in the `KratosGeoMechanicsFastSuiteWithoutKernel` C++ test kernel.

- [ ] **Step 1: Declare the variable**

In `geo_mechanics_application_variables.h`, add after the `HYDRAULIC_DISCHARGE` line (line 31):

```cpp
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, HYDRAULIC_DISCHARGE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, double, NODAL_WATER_FLOW)
```

- [ ] **Step 2: Create the variable**

In `geo_mechanics_application_variables.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 27):

```cpp
KRATOS_CREATE_VARIABLE(double, HYDRAULIC_DISCHARGE)
KRATOS_CREATE_VARIABLE(double, NODAL_WATER_FLOW)
```

- [ ] **Step 3: Register the variable in the application**

In `geo_mechanics_application.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 386):

```cpp
    KRATOS_REGISTER_VARIABLE(HYDRAULIC_DISCHARGE)
    KRATOS_REGISTER_VARIABLE(NODAL_WATER_FLOW)
```

- [ ] **Step 4: Register the variable in Python**

In `custom_python/geo_mechanics_python_application.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 57):

```cpp
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HYDRAULIC_DISCHARGE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_WATER_FLOW)
```

- [ ] **Step 5: Expose the Python alias**

In `GeoMechanicsApplication.py`, add after the `HYDRAULIC_HEAD` alias (line 94), keeping the alphabetical grouping:

```python
HYDRAULIC_DISCHARGE = KratosGeo.HYDRAULIC_DISCHARGE
HYDRAULIC_HEAD = KratosGeo.HYDRAULIC_HEAD
NODAL_WATER_FLOW = KratosGeo.NODAL_WATER_FLOW
```

- [ ] **Step 6: Register the variable in the C++ test kernel**

In `tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 41):

```cpp
    KRATOS_REGISTER_VARIABLE(HYDRAULIC_DISCHARGE)
    KRATOS_REGISTER_VARIABLE(NODAL_WATER_FLOW)
```

- [ ] **Step 7: Build to verify the registration compiles**

Run (from repo root, adjust build dir if different):

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCore KratosGeoMechanicsCoreTest
```

Expected: build succeeds, no undefined-symbol errors for `NODAL_WATER_FLOW`.

- [ ] **Step 8: Commit**

```powershell
git add applications/GeoMechanicsApplication/geo_mechanics_application_variables.h `
        applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp `
        applications/GeoMechanicsApplication/geo_mechanics_application.cpp `
        applications/GeoMechanicsApplication/custom_python/geo_mechanics_python_application.cpp `
        applications/GeoMechanicsApplication/GeoMechanicsApplication.py `
        applications/GeoMechanicsApplication/tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.cpp
git commit -m "feat(geo): register NODAL_WATER_FLOW nodal variable"
```

---

### Task 2: Add and unit-test `AssignNodalWaterFlows`

Adds the pure helper that writes a `NodalFlowMap` onto the historical `NODAL_WATER_FLOW` database, zeroing every node first. This is the unit-tested core.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h`
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp`
- Test: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp`

**Interfaces:**
- Consumes: `Kratos::NODAL_WATER_FLOW` (Task 1); `Geo::SeepageBoundaryUtilities::NodalFlowMap` (existing, `std::unordered_map<std::size_t, double>`).
- Produces: `void Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows)` — sets `NODAL_WATER_FLOW` to `0.0` on every node of `rModelPart`, then writes `rNodalFlows` values onto the nodes whose ids appear as keys.

- [ ] **Step 1: Write the failing test**

Append inside the `namespace Kratos::Testing` block in `test_seepage_boundary_utilities.cpp`, before its closing brace (line 253). The helper `CreateModelPartWithNodes` (line 24) already adds `WATER_PRESSURE`; add `NODAL_WATER_FLOW` to the nodes in the test directly:

```cpp
KRATOS_TEST_CASE_IN_SUITE(AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(NODAL_WATER_FLOW);
    for (std::size_t i = 1; i <= 3; ++i) {
        r_model_part.CreateNewNode(static_cast<int>(i), static_cast<double>(i), 0.0, 0.0);
    }
    // Pre-seed a stale value on node 3 to prove it gets zeroed.
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW) = 99.0;

    // Node 2 is deliberately absent from the map and must end up at 0.0.
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 4.0}, {3, -2.0}};

    Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(r_model_part, nodal_flows);

    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW), -2.0);
}
```

- [ ] **Step 2: Run the test to verify it fails**

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCoreTest
./build/FullDebug/applications/GeoMechanicsApplication/tests/KratosGeoMechanicsCoreTest --gtest_filter=*AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest*
```

Expected: compile error (`AssignNodalWaterFlows` not declared) — this counts as a failing test.

- [ ] **Step 3: Declare the helper**

In `seepage_boundary_utilities.h`, add after the `CalculateNodalWaterFlows` declaration (line 42):

```cpp
// Writes the nodal water flows onto the NODAL_WATER_FLOW solution-step variable of the model part.
// Every node is set to zero first, so nodes absent from rNodalFlows (e.g. nodes without a
// WATER_PRESSURE degree of freedom) hold a defined value rather than stale data. This makes the
// assembled flow visualisable through the normal nodal output path.
void KRATOS_API(GEO_MECHANICS_APPLICATION)
    AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows);
```

- [ ] **Step 4: Implement the helper**

In `seepage_boundary_utilities.cpp`, first add the include for the application variables near the existing includes (after line 15, `#include "includes/variables.h"`):

```cpp
#include "geo_mechanics_application_variables.h"
```

Then add the implementation after `CalculateNodalWaterFlows` (after line 51, the closing brace of that function):

```cpp
void AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_WATER_FLOW) = 0.0;
    }

    for (const auto& [node_id, flow] : rNodalFlows) {
        rModelPart.GetNode(node_id).FastGetSolutionStepValue(NODAL_WATER_FLOW) = flow;
    }
}
```

- [ ] **Step 5: Run the test to verify it passes**

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCoreTest
./build/FullDebug/applications/GeoMechanicsApplication/tests/KratosGeoMechanicsCoreTest --gtest_filter=*AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest*
```

Expected: 1 test, PASS.

- [ ] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h `
        applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp `
        applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp
git commit -m "feat(geo): add AssignNodalWaterFlows helper for NODAL_WATER_FLOW"
```

---

### Task 3: Store the flow each converged step in the seepage strategy

Overrides `FinalizeSolutionStep` in the seepage strategy to assemble the nodal flows once per converged step and store them via the Task 2 helper.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`

**Interfaces:**
- Consumes: `Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(ModelPart&, const ProcessInfo&)` (existing) and `Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(ModelPart&, const NodalFlowMap&)` (Task 2).
- Produces: `NODAL_WATER_FLOW` populated on all nodes at the end of every converged solution step, only for `GeoSeepageNewtonRaphsonStrategy`.

- [ ] **Step 1: Add the `FinalizeSolutionStep` override**

In `geo_seepage_newton_raphson_strategy.hpp`, in the `protected:` section, add above the existing `UpdateSeepageBoundaryConditions` method (before line 306):

```cpp
    // After the step converges, store the assembled nodal water flow on the nodes so it can be
    // visualised. This is exactly the map that drives the boundary switching, which is what makes
    // it useful for verifying the sign convention in ShouldReleaseToNeumann.
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        MotherType::FinalizeSolutionStep();

        auto&       r_model_part  = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const auto  nodal_flows =
            Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(r_model_part, r_process_info);
        Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(r_model_part, nodal_flows);

        KRATOS_CATCH("")
    }

```

- [ ] **Step 2: Build to verify it compiles**

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCore
```

Expected: build succeeds.

- [ ] **Step 3: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
git commit -m "feat(geo): store NODAL_WATER_FLOW each converged step in seepage strategy"
```

---

### Task 4: Make `NODAL_WATER_FLOW` addable and outputtable

Adds the variable to the nodal solution-step variable lists (Python solver + both C++ workflows) and to the GiD output writer's nodal map, so a project can list `NODAL_WATER_FLOW` in `nodal_results`.

**Files:**
- Modify: `applications/GeoMechanicsApplication/python_scripts/geomechanics_solver.py` (near line 364)
- Modify: `applications/GeoMechanicsApplication/custom_workflows/dgeoflow.cpp` (near line 433)
- Modify: `applications/GeoMechanicsApplication/custom_workflows/dgeosettlement.cpp` (near line 287)
- Modify: `applications/GeoMechanicsApplication/custom_workflows/geo_output_writer.cpp` (near line 120)

**Interfaces:**
- Consumes: `Kratos::NODAL_WATER_FLOW` (Task 1); the variable is populated by Task 3.
- Produces: `NODAL_WATER_FLOW` present in the historical database for standard GeoMechanics and workflow runs, and accepted as a `nodal_results` output item that writes via `WriteNodalResults`.

- [ ] **Step 1: Add the variable in the Python solver**

In `geomechanics_solver.py`, `_add_water_variables`, add after the `HYDRAULIC_DISCHARGE` line (line 364):

```python
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.HYDRAULIC_DISCHARGE)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_WATER_FLOW)
```

- [ ] **Step 2: Add the variable in the dgeoflow workflow**

In `dgeoflow.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 433):

```cpp
    rModelPart.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    rModelPart.AddNodalSolutionStepVariable(NODAL_WATER_FLOW);
```

- [ ] **Step 3: Add the variable in the dgeosettlement workflow**

In `dgeosettlement.cpp`, add after the `HYDRAULIC_DISCHARGE` line (line 287):

```cpp
    rModelPart.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    rModelPart.AddNodalSolutionStepVariable(NODAL_WATER_FLOW);
```

- [ ] **Step 4: Add the variable to the nodal output writer map**

In `geo_output_writer.cpp`, `WriteNodalOutput`, add after the `HYDRAULIC_DISCHARGE` entry (line 120):

```cpp
        {"HYDRAULIC_DISCHARGE", MakeNodalResultWriterFor(HYDRAULIC_DISCHARGE)},
        {"NODAL_WATER_FLOW", MakeNodalResultWriterFor(NODAL_WATER_FLOW)},
```

- [ ] **Step 5: Build to verify it compiles**

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCore KratosGeoMechanicsApplication
```

Expected: build succeeds.

- [ ] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/python_scripts/geomechanics_solver.py `
        applications/GeoMechanicsApplication/custom_workflows/dgeoflow.cpp `
        applications/GeoMechanicsApplication/custom_workflows/dgeosettlement.cpp `
        applications/GeoMechanicsApplication/custom_workflows/geo_output_writer.cpp
git commit -m "feat(geo): expose NODAL_WATER_FLOW as a nodal output result"
```

---

### Task 5: Verify the full C++ suite still passes

Guards against regressions in the seepage utilities and output writer.

**Files:**
- No source changes.

**Interfaces:**
- Consumes: everything from Tasks 1–4.
- Produces: nothing.

- [ ] **Step 1: Build the test target**

```powershell
cmake --build build/FullDebug --target KratosGeoMechanicsCoreTest
```

Expected: build succeeds.

- [ ] **Step 2: Run the seepage utility tests**

```powershell
./build/FullDebug/applications/GeoMechanicsApplication/tests/KratosGeoMechanicsCoreTest --gtest_filter=*SeepageBoundary*:*NodalWaterFlow*:*AssignNodalWaterFlows*
```

Expected: all matched tests PASS, including `AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest` from Task 2.

- [ ] **Step 3: (No commit — verification only)**

If any test fails, fix the offending task before proceeding.

---

## Notes

- The GiD writer path needs no special-case computation (unlike `HYDRAULIC_HEAD`): because Task 3 stores `NODAL_WATER_FLOW` in the historical database each step, `MakeNodalResultWriterFor(NODAL_WATER_FLOW)` → `WriteNodalResults` reads it directly.
- Adjust `build/FullDebug` in commands to your actual build directory (e.g. `build/Release-unity`) and the test runner path to wherever `KratosGeoMechanicsCoreTest` is emitted, if different.

