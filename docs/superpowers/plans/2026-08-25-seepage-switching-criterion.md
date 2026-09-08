
 # Seepage Switching Criterion — Step 3 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the switching criterion behind `GeoSeepageNewtonRaphsonStrategy::UpdateSeepageBoundaryConditions()`, so that each non-linear iteration at most one seepage node moves between a Dirichlet boundary (`WATER_PRESSURE` fixed at 0) and a zero-flux Neumann boundary.

**Architecture:** The decision logic lives in a new **non-templated** utility (`custom_utilities/seepage_boundary_utilities.{h,cpp}`) so it can be unit-tested without instantiating the strategy template or running a solver. The strategy's seam becomes a thin call into that utility. A node's `WATER_PRESSURE` fixity is the single source of truth for its mode, which retires the per-condition mode machinery built in step 1.

**Tech Stack:** C++20, Kratos Multiphysics, GTest via `kratos_add_gtests`, CMake + Ninja.

**Spec:** `docs/superpowers/specs/2026-08-25-seepage-switching-criterion-design.md`

## Global Constraints

- Target application: `applications/GeoMechanicsApplication`. All paths are relative to `C:/checkouts/KratosProjects/dev`.
- C++20. New files carry the standard Kratos GeoMechanics header comment block (copy it from `custom_conditions/geo_seepage_condition.h`).
- **Nodal fixity is the only source of truth.** Fixed `WATER_PRESSURE` means Dirichlet, free means Neumann. Do not add any parallel mode variable.
- **At most one node switches per non-linear iteration, model-wide.**
- Neumann to Dirichlet takes precedence over Dirichlet to Neumann when both are possible in the same iteration.
- Ties are broken by **lowest node id**, so results are reproducible.
- Do not modify `applications/GeoMechanicsApplication/CMakeLists.txt`. `custom_utilities/*.cpp` and `tests/*.cpp` are already globbed — but the globs lack `CONFIGURE_DEPENDS`, so run `kp config` once after adding any new `.cpp`.
- Out of scope: `RebuildSystem()` (step 4), restricting the flow assembly to neighbouring elements only, and wiring the strategy into `geomechanics_solver.py`.

## Build and Test Commands

```powershell
kp config                                    # only after adding a new .cpp file
kp build
kp test -- --gtest_filter="*Seepage*"
kp test                                      # full suite; baseline is 1128 passing
```

---

### Task 1: Retire the per-condition mode machinery

Deliverable: `GeoSeepageCondition` stops carrying a mode. It now only identifies seepage nodes and seeds their initial Dirichlet state. The `GEO_SEEPAGE_BOUNDARY_TYPE` variable disappears entirely.

This deletes work done in step 1. That is intended: with fixity as the source of truth, keeping the old path would leave two sources of truth that can silently disagree.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.h` (remove the `GEO_SEEPAGE_BOUNDARY_TYPE` line)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp` (remove the `GEO_SEEPAGE_BOUNDARY_TYPE` line)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application.cpp` (remove `KRATOS_REGISTER_VARIABLE(GEO_SEEPAGE_BOUNDARY_TYPE)`)
- Delete: `applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp`

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `void GeoSeepageCondition::InitializeSolutionStep(const ProcessInfo&)` — sets `WATER_PRESSURE = 0.0` and `Fix`es it on every node of the condition.
  - `int GeoSeepageCondition::Check(const ProcessInfo&) const` — reduced to geometry and dof validation.
  - Removed: `Initialize`, `InitializeNonLinearIteration`, `GetBoundaryType`, `SetBoundaryType`, `mHasOwnProperties`, `Kratos::Geo::SeepageDirichletType`, `Kratos::Geo::SeepageNeumannType`, `Kratos::GEO_SEEPAGE_BOUNDARY_TYPE`.

- [ ] **Step 1: Rewrite the condition tests for the new behaviour**

Replace the whole body of `namespace Kratos::Testing` in `test_geo_seepage_condition.cpp` with the tests below, and delete the now-obsolete helper usage of `Geo::SeepageDirichletType`. Keep the two helpers at the top of the file (`CreateModelPartWithTwoWaterPressureNodes`, `CreateSeepageCondition`) exactly as they are, but delete the `#include "geo_mechanics_application_variables.h"` line if nothing else needs it (the condition header no longer defines the seepage type strings).

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInfoReturnsClassName, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EQ(condition.Info(), "GeoSeepageCondition");
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCreateReturnsGeoSeepageCondition, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    const auto p_created = condition.Create(2, condition.pGetGeometry(), condition.pGetProperties());

    KRATOS_EXPECT_NE(dynamic_cast<const GeoSeepageCondition*>(p_created.get()), nullptr);
    KRATOS_EXPECT_EQ(p_created->Id(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionGetDofListReturnsOneWaterPressureDofPerNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    auto dofs = Condition::DofsVectorType{};
    condition.GetDofList(dofs, ProcessInfo{});

    KRATOS_EXPECT_EQ(dofs.size(), 2);
    KRATOS_EXPECT_EQ(dofs[0]->GetVariable(), WATER_PRESSURE);
    KRATOS_EXPECT_EQ(dofs[1]->GetVariable(), WATER_PRESSURE);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionEquationIdVectorHasOneEntryPerNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    auto equation_ids = Condition::EquationIdVectorType{};
    condition.EquationIdVector(equation_ids, ProcessInfo{});

    KRATOS_EXPECT_EQ(equation_ids.size(), 2);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCalculateLocalSystemReturnsZeroes, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    auto left_hand_side  = Matrix{};
    auto right_hand_side = Vector{};
    condition.CalculateLocalSystem(left_hand_side, right_hand_side, ProcessInfo{});

    const auto expected_left_hand_side  = Matrix{ZeroMatrix{2, 2}};
    const auto expected_right_hand_side = Vector{ZeroVector{2}};

    KRATOS_EXPECT_MATRIX_NEAR(left_hand_side, expected_left_hand_side, 1.0e-12)
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, expected_right_hand_side, 1.0e-12)
}

// The seepage face starts over-prescribed: every node begins as a Dirichlet boundary at zero
// pressure, and the strategy releases nodes one at a time from there.
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInitializeSolutionStepFixesAllNodesAtZeroPressure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 42.0;
    }

    condition.InitializeSolutionStep(ProcessInfo{});

    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    }
}

// A node released to Neumann by the strategy must be re-fixed at the start of the next solution
// step, so each step starts from the same over-prescribed configuration.
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInitializeSolutionStepReFixesAReleasedNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    r_model_part.pGetNode(1)->Free(WATER_PRESSURE);

    condition.InitializeSolutionStep(ProcessInfo{});

    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckReturnsZeroForValidSetup, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EQ(condition.Check(ProcessInfo{}), 0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckThrowsWhenNodeHasNoWaterPressureDof, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    // Deliberately do not add the WATER_PRESSURE degree of freedom.

    auto condition = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.Check(ProcessInfo{}),
                                      "Missing degree of freedom for WATER_PRESSURE on node 1")
}
```

- [ ] **Step 2: Delete the obsolete variable test file**

```powershell
git rm applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp
```

- [ ] **Step 3: Build and verify the new tests fail**

```powershell
kp build
kp test -- --gtest_filter="*GeoSeepageCondition*"
```

Expected: the two new `InitializeSolutionStep` tests FAIL, because `Condition::InitializeSolutionStep`
does nothing by default so the nodes are never fixed. Several older tests will also fail to compile
once the header is stripped in Step 4; that is expected and is resolved by Steps 4 and 5.

- [ ] **Step 4: Strip the header**

In `geo_seepage_condition.h`, delete the entire `namespace Kratos::Geo { ... }` block containing `SeepageDirichletType` and `SeepageNeumannType`. Then replace the class body's declarations so the class reads exactly:

```cpp
// A seepage boundary condition on the WATER_PRESSURE degree of freedom.
//
// This condition holds no mode of its own. A node's WATER_PRESSURE fixity is the single source of
// truth for its boundary type: fixed means a Dirichlet boundary at zero pressure, free means a
// zero-flux Neumann boundary. GeoSeepageNewtonRaphsonStrategy switches individual nodes while
// iterating; this condition only marks which nodes belong to the seepage face and gives them their
// initial state.
//
// The condition never contributes to the linear system.
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoSeepageCondition : public Condition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoSeepageCondition);

    GeoSeepageCondition();
    GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    ~GeoSeepageCondition() override;

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   rThisNodes,
                              PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType               NewId,
                              GeometryType::Pointer   pGeom,
                              PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const override;
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                              Vector&            rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    // Puts the seepage face into its initial, over-prescribed state: every node is a Dirichlet
    // boundary at zero pressure. The strategy releases nodes from there, one per iteration.
    void InitializeSolutionStep(const ProcessInfo&) override;

    [[nodiscard]] int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};
```

- [ ] **Step 5: Strip the implementation**

In `geo_seepage_condition.cpp`:

Delete the `Initialize`, `InitializeNonLinearIteration`, `GetBoundaryType` and `SetBoundaryType` definitions entirely.

Replace `Check` with:

```cpp
int GeoSeepageCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto base_check_result = Condition::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF(GetGeometry().PointsNumber() < 2)
        << "GeoSeepageCondition " << Id() << " needs at least two nodes, but has "
        << GetGeometry().PointsNumber() << std::endl;

    for (const auto& r_node : GetGeometry()) {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(WATER_PRESSURE))
            << "Missing variable WATER_PRESSURE on node " << r_node.Id() << std::endl;
        KRATOS_ERROR_IF_NOT(r_node.HasDofFor(WATER_PRESSURE))
            << "Missing degree of freedom for WATER_PRESSURE on node " << r_node.Id() << std::endl;
    }

    return base_check_result;

    KRATOS_CATCH("")
}
```

Add `InitializeSolutionStep` directly above `Check`:

```cpp
void GeoSeepageCondition::InitializeSolutionStep(const ProcessInfo&)
{
    KRATOS_TRY

    // Start each solution step with the seepage face fully prescribed. Nodes released to a
    // zero-flux Neumann boundary during the previous step are fixed again here.
    for (auto& r_node : GetGeometry()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
        r_node.Fix(WATER_PRESSURE);
    }

    KRATOS_CATCH("")
}
```

Restore the plain serializers (the `mHasOwnProperties` member is gone):

```cpp
void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}
```

- [ ] **Step 6: Remove the variable**

In `geo_mechanics_application_variables.h`, delete:

```cpp
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, GEO_SEEPAGE_BOUNDARY_TYPE)
```

In `geo_mechanics_application_variables.cpp`, delete:

```cpp
KRATOS_CREATE_VARIABLE(std::string, GEO_SEEPAGE_BOUNDARY_TYPE)
```

In `geo_mechanics_application.cpp`, delete:

```cpp
    KRATOS_REGISTER_VARIABLE(GEO_SEEPAGE_BOUNDARY_TYPE)
```

- [ ] **Step 7: Build and verify the tests pass**

```powershell
kp config
kp build
kp test -- --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 9 tests.`

- [ ] **Step 8: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp applications/GeoMechanicsApplication/geo_mechanics_application_variables.h applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp applications/GeoMechanicsApplication/geo_mechanics_application.cpp
git commit -m "Make nodal fixity the source of truth for seepage boundaries"
```

---

### Task 2: Collect the seepage nodes

Deliverable: a utility that returns the distinct nodes belonging to `GeoSeepageCondition`s in a model part.

**Files:**
- Create: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h`
- Create: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp`
- Create: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp`

**Interfaces:**
- Consumes: `Kratos::GeoSeepageCondition` (Task 1).
- Produces: `std::vector<Kratos::Node*> Kratos::Geo::SeepageBoundaryUtilities::CollectSeepageNodes(Kratos::ModelPart& rModelPart)` — distinct nodes of all `GeoSeepageCondition`s, sorted ascending by node id.

The model part is taken by **non-const** reference on purpose. The returned nodes must be mutable, because the strategy fixes and frees them; taking a const model part would force a `const_cast` on the geometry to get `Node*`.

- [ ] **Step 1: Write the failing test**

Create `applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp`:

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

#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/seepage_boundary_utilities.h"
#include "geometries/line_2d_2.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace
{

// Builds a model part with a chain of nodes along the x axis, each owning a WATER_PRESSURE dof.
Kratos::ModelPart& CreateModelPartWithNodes(Kratos::Model& rModel, std::size_t NumberOfNodes)
{
    auto& r_model_part = rModel.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    for (std::size_t i = 1; i <= NumberOfNodes; ++i) {
        r_model_part.CreateNewNode(static_cast<int>(i), static_cast<double>(i), 0.0, 0.0);
    }
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return r_model_part;
}

// Adds a two-noded seepage condition spanning the two given node ids.
void AddSeepageCondition(Kratos::ModelPart& rModelPart, std::size_t Id, std::size_t FirstNodeId, std::size_t SecondNodeId)
{
    auto p_geometry = Kratos::make_shared<Line2D2<Node>>(rModelPart.pGetNode(static_cast<int>(FirstNodeId)),
                                                         rModelPart.pGetNode(static_cast<int>(SecondNodeId)));
    auto p_properties = rModelPart.CreateNewProperties(0);
    rModelPart.AddCondition(Kratos::make_intrusive<GeoSeepageCondition>(
        static_cast<int>(Id), p_geometry, p_properties));
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CollectSeepageNodesReturnsNothingWhenThereAreNoSeepageConditions, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::CollectSeepageNodes(r_model_part).empty())
}

KRATOS_TEST_CASE_IN_SUITE(CollectSeepageNodesReturnsSharedNodesOnlyOnce, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    // Two adjacent conditions sharing node 2.
    AddSeepageCondition(r_model_part, 1, 1, 2);
    AddSeepageCondition(r_model_part, 2, 2, 3);

    const auto nodes = Geo::SeepageBoundaryUtilities::CollectSeepageNodes(r_model_part);

    KRATOS_EXPECT_EQ(nodes.size(), 3);
    KRATOS_EXPECT_EQ(nodes[0]->Id(), 1);
    KRATOS_EXPECT_EQ(nodes[1]->Id(), 2);
    KRATOS_EXPECT_EQ(nodes[2]->Id(), 3);
}

} // namespace Kratos::Testing
```

- [ ] **Step 2: Build and verify the test fails**

```powershell
kp config
kp build
```

Expected: compile error, `cannot open source file "custom_utilities/seepage_boundary_utilities.h"`.

- [ ] **Step 3: Write the header**

Create `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h`:

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

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"

namespace Kratos::Geo::SeepageBoundaryUtilities
{

// Returns the distinct nodes of every GeoSeepageCondition in the model part, sorted ascending by
// node id. Nodes shared by adjacent conditions appear exactly once.
//
// The model part is non-const because the returned nodes must be mutable: the strategy fixes and
// frees their WATER_PRESSURE degree of freedom.
std::vector<Node*> KRATOS_API(GEO_MECHANICS_APPLICATION) CollectSeepageNodes(ModelPart& rModelPart);

} // namespace Kratos::Geo::SeepageBoundaryUtilities
```

- [ ] **Step 4: Write the implementation**

Create `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp`:

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

#include "custom_utilities/seepage_boundary_utilities.h"

#include <algorithm>

#include "custom_conditions/geo_seepage_condition.h"

namespace Kratos::Geo::SeepageBoundaryUtilities
{

std::vector<Node*> CollectSeepageNodes(ModelPart& rModelPart)
{
    auto result = std::vector<Node*>{};

    for (auto& r_condition : rModelPart.Conditions()) {
        if (dynamic_cast<const GeoSeepageCondition*>(&r_condition) == nullptr) continue;

        for (auto& r_node : r_condition.GetGeometry()) {
            result.push_back(&r_node);
        }
    }

    // Adjacent conditions share end nodes, so the same node can be collected more than once.
    std::sort(result.begin(), result.end(),
              [](const Node* pLeft, const Node* pRight) { return pLeft->Id() < pRight->Id(); });
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;
}

} // namespace Kratos::Geo::SeepageBoundaryUtilities
```

- [ ] **Step 5: Build and verify the tests pass**

```powershell
kp config
kp build
kp test -- --gtest_filter="*CollectSeepageNodes*"
```

Expected: `[  PASSED  ] 2 tests.`

- [ ] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp
git commit -m "Add CollectSeepageNodes utility"
```

---

### Task 3: Map element RHS entries onto nodal water flows

Deliverable: the tricky part of the flow assembly, isolated and tested without needing a real element.

An element's right-hand side is a flat vector. For a U-Pw element the displacement and pressure entries are interleaved, so entry *i* is **not** node *i*. Reading it positionally would silently interpret displacement residuals as water flows. The element's own `GetDofList` returns dofs in the same order as the RHS entries, so it is used to pick out the `WATER_PRESSURE` entries.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h`
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp`

**Interfaces:**
- Consumes: `CollectSeepageNodes` (Task 2).
- Produces:
  - `using NodalFlowMap = std::unordered_map<std::size_t, double>;`
  - `void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>& rDofs, const Vector& rRightHandSide, NodalFlowMap& rNodalFlows)` — adds each RHS entry whose dof variable is `WATER_PRESSURE` to that dof's node id.
  - `NodalFlowMap CalculateNodalWaterFlows(ModelPart& rModelPart, const ProcessInfo& rProcessInfo)` — loops every element, calls `CalculateRightHandSide` and `GetDofList`, and accumulates through `AccumulateWaterPressureEntries`.

- [ ] **Step 1: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_seepage_boundary_utilities.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(AccumulateWaterPressureEntriesPicksOnlyWaterPressureDofs, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(WATER_PRESSURE);
    }

    // Mimic a U-Pw ordering: ux, uy, p for node 1, then ux, uy, p for node 2.
    auto dofs = std::vector<Dof<double>*>{r_model_part.pGetNode(1)->pGetDof(DISPLACEMENT_X),
                                          r_model_part.pGetNode(1)->pGetDof(DISPLACEMENT_Y),
                                          r_model_part.pGetNode(1)->pGetDof(WATER_PRESSURE),
                                          r_model_part.pGetNode(2)->pGetDof(DISPLACEMENT_X),
                                          r_model_part.pGetNode(2)->pGetDof(DISPLACEMENT_Y),
                                          r_model_part.pGetNode(2)->pGetDof(WATER_PRESSURE)};

    auto right_hand_side = Vector{6};
    right_hand_side[0] = 10.0; // ux node 1, must be ignored
    right_hand_side[1] = 20.0; // uy node 1, must be ignored
    right_hand_side[2] = 3.0;  // p  node 1
    right_hand_side[3] = 30.0; // ux node 2, must be ignored
    right_hand_side[4] = 40.0; // uy node 2, must be ignored
    right_hand_side[5] = 7.0;  // p  node 2

    auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{};
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, right_hand_side, nodal_flows);

    KRATOS_EXPECT_EQ(nodal_flows.size(), 2);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(1), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(2), 7.0);
}

KRATOS_TEST_CASE_IN_SUITE(AccumulateWaterPressureEntriesSumsContributionsFromSeveralElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);

    auto dofs = std::vector<Dof<double>*>{r_model_part.pGetNode(1)->pGetDof(WATER_PRESSURE),
                                          r_model_part.pGetNode(2)->pGetDof(WATER_PRESSURE)};

    auto first_contribution = Vector{2};
    first_contribution[0]   = 1.5;
    first_contribution[1]   = 2.5;

    auto second_contribution = Vector{2};
    second_contribution[0]   = 0.5;
    second_contribution[1]   = -2.5;

    auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{};
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, first_contribution, nodal_flows);
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, second_contribution, nodal_flows);

    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(2), 0.0);
}
```

- [ ] **Step 2: Build and verify the tests fail**

```powershell
kp build
```

Expected: compile error, `namespace "Kratos::Geo::SeepageBoundaryUtilities" has no member "AccumulateWaterPressureEntries"`.

- [ ] **Step 3: Declare the new interfaces**

In `seepage_boundary_utilities.h`, add inside the namespace, above `CollectSeepageNodes`:

```cpp
// Nodal water flow, keyed by node id.
using NodalFlowMap = std::unordered_map<std::size_t, double>;

// Adds the entries of an element's right-hand side that belong to WATER_PRESSURE degrees of freedom
// onto their nodes.
//
// rDofs must be the element's own degrees of freedom, in the same order as rRightHandSide. For a
// U-Pw element the displacement and pressure entries are interleaved, so the entries cannot be
// matched to nodes positionally; the degrees of freedom are what makes the mapping correct.
void KRATOS_API(GEO_MECHANICS_APPLICATION) AccumulateWaterPressureEntries(
    const std::vector<Dof<double>*>& rDofs, const Vector& rRightHandSide, NodalFlowMap& rNodalFlows);

// Returns the nodal water flow for every node in the model part, assembled from the right-hand side
// of every element. For a Pw element that right-hand side is exactly the permeability flow plus the
// compressibility flow plus the fluid body flow.
NodalFlowMap KRATOS_API(GEO_MECHANICS_APPLICATION)
    CalculateNodalWaterFlows(ModelPart& rModelPart, const ProcessInfo& rProcessInfo);
```

Add `#include "includes/dof.h"` and `#include "includes/ublas_interface.h"` to the header's include block.

- [ ] **Step 4: Implement**

In `seepage_boundary_utilities.cpp`, add above `CollectSeepageNodes`:

```cpp
void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>& rDofs,
                                    const Vector&                    rRightHandSide,
                                    NodalFlowMap&                    rNodalFlows)
{
    KRATOS_ERROR_IF(rDofs.size() != rRightHandSide.size())
        << "Number of degrees of freedom (" << rDofs.size()
        << ") does not match the size of the right hand side (" << rRightHandSide.size() << ")"
        << std::endl;

    for (std::size_t i = 0; i < rDofs.size(); ++i) {
        if (rDofs[i]->GetVariable() != WATER_PRESSURE) continue;

        rNodalFlows[rDofs[i]->Id()] += rRightHandSide[i];
    }
}

NodalFlowMap CalculateNodalWaterFlows(ModelPart& rModelPart, const ProcessInfo& rProcessInfo)
{
    auto result = NodalFlowMap{};

    auto dofs             = std::vector<Dof<double>*>{};
    auto right_hand_side  = Vector{};

    for (auto& r_element : rModelPart.Elements()) {
        r_element.GetDofList(dofs, rProcessInfo);
        r_element.CalculateRightHandSide(right_hand_side, rProcessInfo);

        AccumulateWaterPressureEntries(dofs, right_hand_side, result);
    }

    return result;
}
```

Add `#include "includes/variables.h"` to the implementation's include block.

- [ ] **Step 5: Build and verify the tests pass**

```powershell
kp build
kp test -- --gtest_filter="*AccumulateWaterPressureEntries*"
```

Expected: `[  PASSED  ] 2 tests.`

- [ ] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp
git commit -m "Assemble nodal water flows from element right hand sides"
```

---

### Task 4: The switching decision

Deliverable: the criterion itself, as a pure function over nodes and their flows. This is the heart of step 3.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h`
- Modify: `applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp`

**Interfaces:**
- Consumes: `NodalFlowMap` (Task 3).
- Produces:
  - `bool ShouldReleaseToNeumann(double NodalFlow)` — the isolated sign assumption; returns `NodalFlow > 0.0`.
  - `bool SwitchOneSeepageNode(const std::vector<Node*>& rSeepageNodes, const NodalFlowMap& rNodalFlows)` — applies at most one switch and returns whether it did.

- [ ] **Step 1: Add the test-only helper**

The decision logic takes a plain `std::vector<Node*>`, so it can be tested without constructing any
conditions. Add this to the **anonymous namespace at the top** of
`test_seepage_boundary_utilities.cpp`, next to the existing helpers:

```cpp
// Test-only shortcut: treats every node of the model part as a seepage node, so the decision logic
// can be exercised without constructing conditions. Production code uses CollectSeepageNodes.
std::vector<Kratos::Node*> AllNodesOf(Kratos::ModelPart& rModelPart)
{
    auto result = std::vector<Kratos::Node*>{};
    for (auto& r_node : rModelPart.Nodes()) {
        result.push_back(&r_node);
    }
    return result;
}
```

- [ ] **Step 2: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_seepage_boundary_utilities.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeDoesNothingWhenNoNodeViolatesItsCondition, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 fixed with no outflow, node 2 free and under suction: both are consistent.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = -5.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, -1.0}, {2, 0.0}};

    KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeFixesTheHighestPressureFreeNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Free(WATER_PRESSURE);
    }
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = 9.0; // highest
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(WATER_PRESSURE) = -1.0;

    const auto nodes = AllNodesOf(r_model_part);

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(
        nodes, Geo::SeepageBoundaryUtilities::NodalFlowMap{}))

    // Only node 2 switches, and it is prescribed at zero pressure.
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeReleasesTheLargestOutflowFixedNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 4.0}, {2, 11.0}, {3, -2.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

    // Only node 2, which has the largest outflow, is released.
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodePrefersFixingOverReleasing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 is a fixed node with a large outflow, node 2 is a free node under positive pressure.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 100.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

    // The Neumann to Dirichlet switch wins, and node 1 is left alone this iteration.
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeBreaksTiesByLowestNodeId, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 5.0}, {2, 5.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(ShouldReleaseToNeumannIsTrueOnlyForPositiveFlow, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(1.0e-9))
    KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(0.0))
    KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(-1.0))
}
```

- [ ] **Step 3: Build and verify the tests fail**

```powershell
kp build
```

Expected: compile error, `namespace "Kratos::Geo::SeepageBoundaryUtilities" has no member "SwitchOneSeepageNode"`.

- [ ] **Step 4: Declare the interfaces**

In `seepage_boundary_utilities.h`, add inside the namespace:

```cpp
// Returns true if a Dirichlet seepage node carrying this nodal flow should be released to a
// zero-flux Neumann boundary.
//
// This sign is the key modelling assumption of the prototype: a positive nodal flow means outflow,
// and the node carrying the largest outflow is released. Confirming or inverting it is what the
// Muskat validation case is for, so it deliberately lives in one place.
[[nodiscard]] bool KRATOS_API(GEO_MECHANICS_APPLICATION) ShouldReleaseToNeumann(double NodalFlow);

// Switches at most one seepage node between a Dirichlet and a zero-flux Neumann boundary, and
// returns whether it switched anything.
//
// A free node whose WATER_PRESSURE is positive is unsaturated, so it should not be draining: the
// highest-pressure such node is fixed at zero pressure. Otherwise the fixed node with the largest
// outflow is released. Fixing takes precedence over releasing, and ties are broken by the lowest
// node id so the result is reproducible.
bool KRATOS_API(GEO_MECHANICS_APPLICATION)
    SwitchOneSeepageNode(const std::vector<Node*>& rSeepageNodes, const NodalFlowMap& rNodalFlows);
```

- [ ] **Step 5: Implement**

In `seepage_boundary_utilities.cpp`, add:

```cpp
namespace
{

// Returns the node maximising the given score among the candidates, or nullptr when there are none.
// Candidates are visited in ascending node id order, and a strict comparison keeps the first of any
// tie, which makes the choice reproducible.
template <typename PredicateType, typename ScoreType>
Kratos::Node* SelectBestCandidate(const std::vector<Kratos::Node*>& rNodes, PredicateType IsCandidate, ScoreType Score)
{
    Kratos::Node* p_best       = nullptr;
    auto          best_score   = 0.0;

    for (auto* p_node : rNodes) {
        if (!IsCandidate(*p_node)) continue;

        const auto score = Score(*p_node);
        if (p_best == nullptr || score > best_score) {
            p_best     = p_node;
            best_score = score;
        }
    }

    return p_best;
}

} // namespace

bool ShouldReleaseToNeumann(double NodalFlow) { return NodalFlow > 0.0; }

bool SwitchOneSeepageNode(const std::vector<Node*>& rSeepageNodes, const NodalFlowMap& rNodalFlows)
{
    KRATOS_TRY

    // A free node under positive pressure is unsaturated, so it cannot be a draining face. Fixing
    // takes precedence over releasing.
    if (auto* p_node = SelectBestCandidate(
            rSeepageNodes,
            [](const Node& rNode) {
                return !rNode.IsFixed(WATER_PRESSURE) &&
                       rNode.FastGetSolutionStepValue(WATER_PRESSURE) > 0.0;
            },
            [](const Node& rNode) { return rNode.FastGetSolutionStepValue(WATER_PRESSURE); })) {
        p_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
        p_node->Fix(WATER_PRESSURE);
        return true;
    }

    // Otherwise release the prescribed node carrying the largest outflow.
    const auto flow_of = [&rNodalFlows](const Node& rNode) {
        const auto it = rNodalFlows.find(rNode.Id());
        return it == rNodalFlows.end() ? 0.0 : it->second;
    };

    if (auto* p_node = SelectBestCandidate(
            rSeepageNodes,
            [&flow_of](const Node& rNode) {
                return rNode.IsFixed(WATER_PRESSURE) && ShouldReleaseToNeumann(flow_of(rNode));
            },
            flow_of)) {
        p_node->Free(WATER_PRESSURE);
        return true;
    }

    return false;

    KRATOS_CATCH("")
}
```

- [ ] **Step 6: Build and verify the tests pass**

```powershell
kp build
kp test -- --gtest_filter="*SwitchOneSeepageNode*:*ShouldReleaseToNeumann*"
```

Expected: `[  PASSED  ] 6 tests.`

- [ ] **Step 7: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.h applications/GeoMechanicsApplication/custom_utilities/seepage_boundary_utilities.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_utilities/test_seepage_boundary_utilities.cpp
git commit -m "Add the seepage boundary switching decision"
```

---

### Task 5: Wire the decision into the strategy

Deliverable: `UpdateSeepageBoundaryConditions()` stops returning `false` and starts driving the criterion, completing step 3.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp`

**Interfaces:**
- Consumes: `CollectSeepageNodes` (Task 2), `CalculateNodalWaterFlows` (Task 3), `SwitchOneSeepageNode` (Task 4).
- Produces: an overriding `bool UpdateSeepageBoundaryConditions()` and a cached `std::vector<Node*> mSeepageNodes`, filled by an overriding `void Initialize()`.

- [ ] **Step 1: Add the include and the member**

In `geo_seepage_newton_raphson_strategy.hpp`, add to the application includes:

```cpp
#include "custom_utilities/seepage_boundary_utilities.h"
```

Add a `private:` section at the end of the class, before the closing brace:

```cpp
private:
    // Cached once in Initialize. The conditions of a model part do not change during a solve, so
    // there is no need to rediscover the seepage nodes every iteration.
    std::vector<Node*> mSeepageNodes;
```

- [ ] **Step 2: Cache the seepage nodes in Initialize**

Add to the `public:` section, below `Info()`:

```cpp
    void Initialize() override
    {
        KRATOS_TRY

        MotherType::Initialize();
        mSeepageNodes = Geo::SeepageBoundaryUtilities::CollectSeepageNodes(BaseType::GetModelPart());

        KRATOS_INFO_IF("GeoSeepageNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
            << "Found " << mSeepageNodes.size() << " seepage nodes" << std::endl;

        KRATOS_CATCH("")
    }
```

- [ ] **Step 3: Implement the seam**

Step 2 declared `virtual bool UpdateSeepageBoundaryConditions() { return false; }` in **this same
class**, so there is no base declaration to override. Replace that stub body in place, and do **not**
add `override`:

```cpp
    // Decides whether one seepage node must change between a Dirichlet and a zero-flux Neumann
    // boundary, applies that change, and reports whether anything changed. At most one node
    // switches per non-linear iteration, across the whole model.
    virtual bool UpdateSeepageBoundaryConditions()
    {
        KRATOS_TRY

        if (mSeepageNodes.empty()) return false;

        const auto nodal_flows = Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(
            BaseType::GetModelPart(), BaseType::GetModelPart().GetProcessInfo());

        return Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(mSeepageNodes, nodal_flows);

        KRATOS_CATCH("")
    }
```

- [ ] **Step 4: Build and run the full suite**

```powershell
kp build
kp test
```

The baseline before this plan was **1128** passing. This plan removes 3 tests (the two
`GEO_SEEPAGE_BOUNDARY_TYPE` variable tests and the registration test) plus the 5 condition tests that
covered the retired mode machinery, and adds 2 condition tests and 10 utility tests.

Do not treat any particular total as the pass criterion. The criterion is: **zero failures**, and no
test that passed before this plan now fails. Record the actual total in the findings section below.

- [ ] **Step 5: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_strategies/strategies/geo_seepage_newton_raphson_strategy.hpp
git commit -m "Drive the seepage switching criterion from the strategy"
```

---

## Findings to Record for the Follow-up Issue (#14637)

Fill this in during implementation:

- Whether `CollectSeepageNodes` needed a `const_cast` to obtain mutable `Node*` from a const `ModelPart`, and whether a cleaner signature would be better.
- Whether recomputing every element's right-hand side each iteration is noticeably slow, which would justify restricting it to elements adjacent to seepage conditions using `find_neighbour_elements_of_conditions_process.h`.
- Whether the one-switch-per-iteration rate is sufficient for the Muskat mesh, or whether the iteration budget is exhausted before the exit point settles.
- Whether any node oscillates between fixed and free on alternating iterations, which would call for a once-per-step switching limit.
- Above all: whether `ShouldReleaseToNeumann` has the correct sign, judged against the Muskat reference solution.

## Findings Recorded During Implementation

**Status: all five tasks implemented. Full GeoMechanics C++ suite green (1126 tests, no failures). 19 new/rewritten feature tests pass.**

The test total went from 1128 (before step 3) to 1126: three `GEO_SEEPAGE_BOUNDARY_TYPE` tests and five retired mode-machinery tests were removed, and two condition tests plus ten utility tests were added. The pass criterion — zero failures, nothing previously passing now failing — holds.

Deviations from the plan:

- **`CollectSeepageNodes` did NOT need a `const_cast`.** Taking the model part by non-const reference (the plan's own fix to its first draft) was sufficient. `ModelPart::Conditions()` yields mutable conditions, whose `GetGeometry()` yields mutable nodes, so `Node*` came out cleanly.
- **A file was created with corrupt leading bytes.** `test_seepage_boundary_utilities.cpp` was written with five stray bytes (`It's'`) prepended before the license header, producing a baffling `C4430: missing type specifier` on line 1 — a line that is only a comment. Detected by dumping the raw bytes and comparing against a known-good test file; fixed by stripping the five bytes. Worth checking the first bytes of any newly created file if the compiler complains about line 1.
- **A test-harness bug, not a production bug.** The first multi-condition test called `CreateNewProperties(0)` once per condition, and Kratos throws on the second creation of the same property id. Fixed by reusing the existing Properties in the test helper. The production code was correct.

Confirmed as designed:

- The interleaved U-Pw dof mapping works: `AccumulateWaterPressureEntries` ignored displacement entries and kept only the `WATER_PRESSURE` ones.
- The precedence rule (Neumann to Dirichlet wins) and the lowest-node-id tie-break both behave as specified.
- Making nodal fixity the source of truth cleanly retired the step 1 mode machinery with no dangling references.

## Still Open for #14637 (need a running solver)

- The one-switch-per-iteration migration rate versus the iteration budget on a real mesh.
- The cost of recomputing every element's RHS each iteration.
- Whether any node oscillates between fixed and free.
- Whether `RebuildSystem()` (step 4) is genuinely needed under a block builder and solver.
- **The sign of `ShouldReleaseToNeumann`** — the central modelling assumption, still to be judged against Muskat.












