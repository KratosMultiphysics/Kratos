# Seepage Boundary Prototype — Step 1: Mixed Seepage Condition — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Build a `GeoSeepageCondition` for the GeoMechanicsApplication that acts as either a Dirichlet (fixed `WATER_PRESSURE = 0`) or a zero-flux Neumann boundary, switchable from outside the condition on every non-linear iteration.

**Architecture:** A single non-templated condition deriving directly from `Kratos::Condition`, managing only `WATER_PRESSURE` degrees of freedom. The mode is stored as a `std::string` in the condition's `Properties` and read fresh in `InitializeNonLinearIteration`, which fixes or frees the nodal degrees of freedom. Each condition clones its `Properties` in `Initialize` so that conditions can switch independently. The condition never contributes to the linear system.

**Tech Stack:** C++20, Kratos Multiphysics, GTest via `kratos_add_gtests`, CMake + Ninja.

**Spec:** `docs/superpowers/specs/2026-08-25-seepage-boundary-condition-design.md`

## Global Constraints

- Target application: `applications/GeoMechanicsApplication`. All paths below are relative to `C:/checkouts/KratosProjects/dev`.
- C++20. All new files carry the standard Kratos GeoMechanics header comment block (copy it from `custom_conditions/Pw_point_flux_condition.h`).
- The only legal values of the mode are the exact strings `"Dirichlet"` and `"Neumann"`. Any other value is a hard error.
- Dirichlet mode prescribes `WATER_PRESSURE = 0.0`. Neumann mode is zero flux and contributes nothing — do not write integration code for it.
- The condition must never contribute to the linear system. `CalculateLocalSystem` returns zeroes in both modes.
- 2D only: `Line2D2` and `Line2D3`. Do not add 3D geometries.
- The class is not templated on `TDim` / `TNumNodes`. Node count comes from `GetGeometry().PointsNumber()` at runtime.
- Export the class with `KRATOS_API(GEO_MECHANICS_APPLICATION)`.
- Do not modify `applications/GeoMechanicsApplication/CMakeLists.txt`. It already globs `custom_conditions/*.cpp` and `tests/*.cpp`.
- Out of scope for this plan: the Newton-Raphson strategy, the switching criterion, and rebuilding the system. Those are steps 2 to 4 of the prototype.

## Build and Test Commands

Because `CMakeLists.txt` uses `file(GLOB_RECURSE ...)` **without** `CONFIGURE_DEPENDS`, adding a new `.cpp` file requires re-running CMake once before it is picked up.

Re-run CMake (only needed after adding a new `.cpp` file), use the custom tooling:

```powershell
kp config
```

Build the test binary (custom tooling):

```powershell
kp b
```

Run only this feature's tests (custom tooling):

```powershell
kp t --gtest_filter="*GeoSeepageCondition*"
```

---

### Task 1: Register the `GEO_SEEPAGE_BOUNDARY_TYPE` variable

**Files:**
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.h` (near line 45, next to `GEO_DRAINAGE_TYPE`)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp` (near line 86, next to `GEO_DRAINAGE_TYPE`)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application.cpp` (near line 442, next to `KRATOS_REGISTER_VARIABLE(GEO_DRAINAGE_TYPE)`)
- Create: `applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp`

**Interfaces:**
- Consumes: nothing.
- Produces: `Kratos::GEO_SEEPAGE_BOUNDARY_TYPE`, a `Variable<std::string>`, registered under the component name `"GEO_SEEPAGE_BOUNDARY_TYPE"`.

- [x] **Step 1: Write the failing test**

Create `applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp`:

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

#include "geo_mechanics_application_variables.h"
#include "includes/kratos_components.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionBoundaryTypeVariableIsRegistered, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Variable<std::string>>::Has("GEO_SEEPAGE_BOUNDARY_TYPE"))
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionBoundaryTypeVariableHasExpectedName, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_EQ(GEO_SEEPAGE_BOUNDARY_TYPE.Name(), "GEO_SEEPAGE_BOUNDARY_TYPE");
}

} // namespace Kratos::Testing
```

- [x] **Step 2: Re-run CMake, build, and verify the test fails**

```powershell
cmake C:/checkouts/KratosProjects/dev/build/FullDebug
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
```

Expected: compile error, `'GEO_SEEPAGE_BOUNDARY_TYPE': undeclared identifier`.

- [x] **Step 3: Declare the variable**

In `geo_mechanics_application_variables.h`, directly below the existing `GEO_DRAINAGE_TYPE` line (line 45):

```cpp
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, GEO_DRAINAGE_TYPE)
KRATOS_DEFINE_APPLICATION_VARIABLE(GEO_MECHANICS_APPLICATION, std::string, GEO_SEEPAGE_BOUNDARY_TYPE)
```

- [x] **Step 4: Create the variable**

In `geo_mechanics_application_variables.cpp`, directly below the existing `GEO_DRAINAGE_TYPE` line (line 86):

```cpp
KRATOS_CREATE_VARIABLE(std::string, GEO_DRAINAGE_TYPE)
KRATOS_CREATE_VARIABLE(std::string, GEO_SEEPAGE_BOUNDARY_TYPE)
```

- [x] **Step 5: Register the variable**

In `geo_mechanics_application.cpp`, in `KratosGeoMechanicsApplication::Register()`, directly below the existing `KRATOS_REGISTER_VARIABLE(GEO_DRAINAGE_TYPE)` line (line 442):

```cpp
    KRATOS_REGISTER_VARIABLE(GEO_DRAINAGE_TYPE)
    KRATOS_REGISTER_VARIABLE(GEO_SEEPAGE_BOUNDARY_TYPE)
```

- [x] **Step 6: Build and verify the tests pass**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageConditionBoundaryTypeVariable*"
```

Expected: `[  PASSED  ] 2 tests.`

- [x] **Step 7: Commit**

```powershell
git add applications/GeoMechanicsApplication/geo_mechanics_application_variables.h applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp applications/GeoMechanicsApplication/geo_mechanics_application.cpp applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp
git commit -m "Add GEO_SEEPAGE_BOUNDARY_TYPE variable"
```

---

### Task 2: Create the `GeoSeepageCondition` skeleton and register it

Deliverable: a constructible, registered condition that reports its `WATER_PRESSURE` degrees of freedom and contributes nothing to the linear system. No mode handling yet.

**Files:**
- Create: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Create: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application.h` (include near line 32; members near line 803)
- Modify: `applications/GeoMechanicsApplication/geo_mechanics_application.cpp` (registration near line 265)
- Create: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`

**Interfaces:**
- Consumes: `Kratos::GEO_SEEPAGE_BOUNDARY_TYPE` from Task 1 (via the `geo_mechanics_application_variables.h` include only; not used yet). `Geo::DofUtilities::ExtractDofsFromNodes(const NodeRange&, const Variable<double>&)` and `Geo::DofUtilities::ExtractEquationIdsFrom(const std::vector<Dof<double>*>&)` from `custom_utilities/dof_utilities.hpp`.
- Produces:
  - `class Kratos::GeoSeepageCondition : public Kratos::Condition`
  - `GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)`
  - `Condition::Pointer Create(IndexType, NodesArrayType const&, PropertiesType::Pointer) const`
  - `Condition::Pointer Create(IndexType, GeometryType::Pointer, PropertiesType::Pointer) const`
  - `void GetDofList(DofsVectorType&, const ProcessInfo&) const`
  - `void EquationIdVector(EquationIdVectorType&, const ProcessInfo&) const`
  - `void CalculateLocalSystem(Matrix&, Vector&, const ProcessInfo&)`
  - `std::string Info() const` returning `"GeoSeepageCondition"`
  - `inline const std::string Kratos::Geo::SeepageDirichletType` = `"Dirichlet"`
  - `inline const std::string Kratos::Geo::SeepageNeumannType` = `"Neumann"`
  - Registered condition names `"GeoSeepageCondition2D2N"` and `"GeoSeepageCondition2D3N"`

- [x] **Step 1: Write the failing test**

Create `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`:

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

#include "custom_conditions/geo_seepage_condition.h"
#include "containers/model.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace
{

// Creates a model part holding two nodes that own a WATER_PRESSURE degree of freedom.
Kratos::ModelPart& CreateModelPartWithTwoWaterPressureNodes(Kratos::Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return r_model_part;
}

// Creates a two-noded seepage condition on the given model part, with its own Properties.
Kratos::GeoSeepageCondition CreateSeepageCondition(Kratos::ModelPart& rModelPart)
{
    auto p_geometry = std::make_shared<Line2D2<Node>>(rModelPart.pGetNode(1), rModelPart.pGetNode(2));
    auto p_properties = Kratos::make_shared<Properties>(1);
    return GeoSeepageCondition{1, p_geometry, p_properties};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInfoReturnsClassName, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model      = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition  = CreateSeepageCondition(r_model_part);

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

    KRATOS_EXPECT_MATRIX_NEAR(left_hand_side, ZeroMatrix(2, 2), 1.0e-12)
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, ZeroVector(2), 1.0e-12)
}

} // namespace Kratos::Testing
```

- [x] **Step 2: Re-run CMake, build, and verify the test fails**

```powershell
cmake C:/checkouts/KratosProjects/dev/build/FullDebug
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
```

Expected: compile error, `cannot open source file "custom_conditions/geo_seepage_condition.h"`.

- [x] **Step 3: Write the header**

Create `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`:

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

#include "includes/condition.h"
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos::Geo
{

// The only legal values of GEO_SEEPAGE_BOUNDARY_TYPE.
inline const std::string SeepageDirichletType = "Dirichlet";
inline const std::string SeepageNeumannType   = "Neumann";

} // namespace Kratos::Geo

namespace Kratos
{

// A mixed seepage boundary condition on the WATER_PRESSURE degree of freedom.
//
// The condition acts either as a Dirichlet boundary (WATER_PRESSURE fixed at zero) or as a
// zero-flux Neumann boundary (WATER_PRESSURE free). The mode is stored in the condition's
// Properties under GEO_SEEPAGE_BOUNDARY_TYPE and is re-read on every non-linear iteration, so
// it can be switched from outside the condition while solving.
//
// The condition never contributes to the linear system; it only changes the fixity of its nodes.
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

    [[nodiscard]] std::string Info() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
```

- [x] **Step 4: Write the implementation**

Create `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`:

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

#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/dof_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/variables.h"

namespace Kratos
{

GeoSeepageCondition::GeoSeepageCondition() : GeoSeepageCondition(0, nullptr, nullptr) {}

GeoSeepageCondition::GeoSeepageCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoSeepageCondition(NewId, pGeometry, nullptr)
{
}

GeoSeepageCondition::GeoSeepageCondition(IndexType               NewId,
                                         GeometryType::Pointer   pGeometry,
                                         PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

GeoSeepageCondition::~GeoSeepageCondition() = default;

Condition::Pointer GeoSeepageCondition::Create(IndexType               NewId,
                                               NodesArrayType const&   rThisNodes,
                                               PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

Condition::Pointer GeoSeepageCondition::Create(IndexType               NewId,
                                               GeometryType::Pointer   pGeom,
                                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoSeepageCondition>(NewId, pGeom, pProperties);
}

void GeoSeepageCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

void GeoSeepageCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void GeoSeepageCondition::CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
                                               Vector& rRightHandSideVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void GeoSeepageCondition::CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo&)
{
    KRATOS_TRY

    // A seepage condition never contributes to the system matrix: in Dirichlet mode the fixed
    // degrees of freedom are handled by the builder and solver, and in Neumann mode the flux is zero.
    const auto number_of_nodes = GetGeometry().PointsNumber();
    if (rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes, number_of_nodes);

    KRATOS_CATCH("")
}

void GeoSeepageCondition::CalculateRightHandSide(Vector& rRightHandSideVector, const ProcessInfo&)
{
    KRATOS_TRY

    const auto number_of_nodes = GetGeometry().PointsNumber();
    if (rRightHandSideVector.size() != number_of_nodes) {
        rRightHandSideVector.resize(number_of_nodes, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);

    KRATOS_CATCH("")
}

Condition::DofsVectorType GeoSeepageCondition::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
}

std::string GeoSeepageCondition::Info() const { return "GeoSeepageCondition"; }

void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

} // namespace Kratos
```

- [x] **Step 5: Register the condition for both geometries**

In `geo_mechanics_application.h`, add the include directly below the existing `#include "custom_conditions/Pw_point_flux_condition.h"` (line 32):

```cpp
#include "custom_conditions/Pw_point_flux_condition.h"
#include "custom_conditions/geo_seepage_condition.h"
```

In the same file, add the two members directly below `mPwPointFluxCondition3D1N` (around line 803):

```cpp
    const PwPointFluxCondition<3, 1> mPwPointFluxCondition3D1N{
        0, Kratos::make_shared<Point3D<NodeType>>(Condition::GeometryType::PointsArrayType(1))};

    const GeoSeepageCondition mGeoSeepageCondition2D2N{
        0, Kratos::make_shared<Line2D2<NodeType>>(Condition::GeometryType::PointsArrayType(2))};
    const GeoSeepageCondition mGeoSeepageCondition2D3N{
        0, Kratos::make_shared<Line2D3<NodeType>>(Condition::GeometryType::PointsArrayType(3))};
```

In `geo_mechanics_application.cpp`, in `Register()`, directly below the `PwPointFluxCondition` registrations (around line 265):

```cpp
    KRATOS_REGISTER_CONDITION("PwPointFluxCondition3D1N", mPwPointFluxCondition3D1N)

    KRATOS_REGISTER_CONDITION("GeoSeepageCondition2D2N", mGeoSeepageCondition2D2N)
    KRATOS_REGISTER_CONDITION("GeoSeepageCondition2D3N", mGeoSeepageCondition2D3N)
```

- [x] **Step 6: Add a registration test**

Append to `applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp`, inside `namespace Kratos::Testing`, and add `#include "includes/kratos_components.h"` and `#include "includes/condition.h"` at the top if not already present:

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionIsRegisteredForBothLineGeometries, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has("GeoSeepageCondition2D2N"))
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has("GeoSeepageCondition2D3N"))
}
```

- [x] **Step 7: Build and verify all tests pass**

```powershell
cmake C:/checkouts/KratosProjects/dev/build/FullDebug
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 8 tests.`

- [x] **Step 8: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/geo_mechanics_application.h applications/GeoMechanicsApplication/geo_mechanics_application.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/test_geo_seepage_boundary_type_variable.cpp
git commit -m "Add GeoSeepageCondition skeleton and register it for Line2D2 and Line2D3"
```

---

### Task 3: Add validated boundary type accessors

Deliverable: the mode can be read and written through the condition, with invalid values rejected at the point of the mistake.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`

**Interfaces:**
- Consumes: `Kratos::GEO_SEEPAGE_BOUNDARY_TYPE` (Task 1); `GeoSeepageCondition` (Task 2); `Kratos::Geo::SeepageDirichletType` and `Kratos::Geo::SeepageNeumannType` (Task 2).
- Produces:
  - `const std::string& GeoSeepageCondition::GetBoundaryType() const` — errors if the property is absent.
  - `void GeoSeepageCondition::SetBoundaryType(const std::string& rBoundaryType)` — errors on any value other than the two legal strings, otherwise writes `GEO_SEEPAGE_BOUNDARY_TYPE` into `GetProperties()`.

- [x] **Step 1: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_geo_seepage_condition.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionSetBoundaryTypeIsReadBackByGetBoundaryType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    condition.SetBoundaryType(Geo::SeepageDirichletType);
    KRATOS_EXPECT_EQ(condition.GetBoundaryType(), Geo::SeepageDirichletType);

    condition.SetBoundaryType(Geo::SeepageNeumannType);
    KRATOS_EXPECT_EQ(condition.GetBoundaryType(), Geo::SeepageNeumannType);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionSetBoundaryTypeThrowsOnUnknownValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.SetBoundaryType("Robin"),
                                      "Unknown seepage boundary type 'Robin' for GeoSeepageCondition 1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionGetBoundaryTypeThrowsWhenPropertyIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto& r_type = condition.GetBoundaryType(),
        "GEO_SEEPAGE_BOUNDARY_TYPE is not defined for GeoSeepageCondition 1")
}
```

- [x] **Step 2: Build and verify the tests fail**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
```

Expected: compile error, `class "Kratos::GeoSeepageCondition" has no member "SetBoundaryType"`.

- [x] **Step 3: Declare the accessors**

In `geo_seepage_condition.h`, add directly above `[[nodiscard]] std::string Info() const override;`:

```cpp
    // Returns the currently configured boundary type. Errors if GEO_SEEPAGE_BOUNDARY_TYPE is absent.
    [[nodiscard]] const std::string& GetBoundaryType() const;

    // Sets the boundary type on this condition's Properties. Errors on any value other than
    // Geo::SeepageDirichletType or Geo::SeepageNeumannType.
    void SetBoundaryType(const std::string& rBoundaryType);
```

- [x] **Step 4: Implement the accessors**

In `geo_seepage_condition.cpp`, add directly above `std::string GeoSeepageCondition::Info() const`:

```cpp
const std::string& GeoSeepageCondition::GetBoundaryType() const
{
    KRATOS_ERROR_IF_NOT(GetProperties().Has(GEO_SEEPAGE_BOUNDARY_TYPE))
        << "GEO_SEEPAGE_BOUNDARY_TYPE is not defined for GeoSeepageCondition " << Id() << std::endl;

    return GetProperties()[GEO_SEEPAGE_BOUNDARY_TYPE];
}

void GeoSeepageCondition::SetBoundaryType(const std::string& rBoundaryType)
{
    KRATOS_ERROR_IF(rBoundaryType != Geo::SeepageDirichletType && rBoundaryType != Geo::SeepageNeumannType)
        << "Unknown seepage boundary type '" << rBoundaryType << "' for GeoSeepageCondition " << Id()
        << ", expected '" << Geo::SeepageDirichletType << "' or '" << Geo::SeepageNeumannType << "'"
        << std::endl;

    GetProperties().SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, rBoundaryType);
}
```

- [x] **Step 5: Build and verify the tests pass**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 11 tests.`

- [x] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp
git commit -m "Add validated boundary type accessors to GeoSeepageCondition"
```

---

### Task 4: Give each condition its own Properties on Initialize

Deliverable: two conditions constructed from one shared `Properties` can hold different modes after `Initialize`. This is what makes per-condition switching possible.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`

**Interfaces:**
- Consumes: `GeoSeepageCondition::SetBoundaryType` and `GetBoundaryType` (Task 3).
- Produces: `void GeoSeepageCondition::Initialize(const ProcessInfo&)` — clones the `Properties` exactly once per condition instance, guarded by the serialized member `bool mHasOwnProperties`.

- [x] **Step 1: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_geo_seepage_condition.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInitializeGivesEachConditionItsOwnProperties, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);

    auto p_geometry   = std::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    auto p_properties = Kratos::make_shared<Properties>(1);
    p_properties->SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, Geo::SeepageDirichletType);

    auto first_condition  = GeoSeepageCondition{1, p_geometry, p_properties};
    auto second_condition = GeoSeepageCondition{2, p_geometry, p_properties};

    first_condition.Initialize(ProcessInfo{});
    second_condition.Initialize(ProcessInfo{});

    // Switching one condition must not affect the other.
    first_condition.SetBoundaryType(Geo::SeepageNeumannType);

    KRATOS_EXPECT_EQ(first_condition.GetBoundaryType(), Geo::SeepageNeumannType);
    KRATOS_EXPECT_EQ(second_condition.GetBoundaryType(), Geo::SeepageDirichletType);
    KRATOS_EXPECT_EQ((*p_properties)[GEO_SEEPAGE_BOUNDARY_TYPE], Geo::SeepageDirichletType);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionSecondInitializeDoesNotResetTheBoundaryType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);

    auto p_geometry   = std::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    auto p_properties = Kratos::make_shared<Properties>(1);
    p_properties->SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, Geo::SeepageDirichletType);

    auto condition = GeoSeepageCondition{1, p_geometry, p_properties};

    condition.Initialize(ProcessInfo{});
    condition.SetBoundaryType(Geo::SeepageNeumannType);
    condition.Initialize(ProcessInfo{});

    KRATOS_EXPECT_EQ(condition.GetBoundaryType(), Geo::SeepageNeumannType);
}
```

- [x] **Step 2: Build and verify the tests fail**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageConditionInitializeGivesEach*:*GeoSeepageConditionSecondInitialize*"
```

Expected: FAIL. Without the clone, both conditions share one `Properties`, so `second_condition.GetBoundaryType()` returns `"Neumann"` instead of `"Dirichlet"`.

- [x] **Step 3: Declare `Initialize` and the ownership guard**

In `geo_seepage_condition.h`, add directly above `[[nodiscard]] const std::string& GetBoundaryType() const;`:

```cpp
    // Gives this condition its own copy of its Properties, so that its boundary type can be
    // switched independently of the other conditions that originally shared those Properties.
    void Initialize(const ProcessInfo&) override;
```

And in the `private:` section, directly above `friend class Serializer;`:

```cpp
    bool mHasOwnProperties = false;
```

- [x] **Step 4: Implement `Initialize` and serialize the guard**

In `geo_seepage_condition.cpp`, add directly above `const std::string& GeoSeepageCondition::GetBoundaryType() const`:

```cpp
void GeoSeepageCondition::Initialize(const ProcessInfo&)
{
    KRATOS_TRY

    // Properties are normally shared by many conditions. Since the boundary type is switched per
    // condition, each condition needs its own copy. The guard keeps this idempotent, so that a
    // second Initialize cannot undo a switch that has already been made.
    if (mHasOwnProperties) return;

    SetProperties(Kratos::make_shared<PropertiesType>(GetProperties()));
    mHasOwnProperties = true;

    KRATOS_CATCH("")
}
```

Replace the existing `save` and `load` in the same file with:

```cpp
void GeoSeepageCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
    rSerializer.save("HasOwnProperties", mHasOwnProperties);
}

void GeoSeepageCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
    rSerializer.load("HasOwnProperties", mHasOwnProperties);
}
```

- [x] **Step 5: Build and verify the tests pass**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 13 tests.`

- [x] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp
git commit -m "Give each GeoSeepageCondition its own Properties on Initialize"
```

---

### Task 5: Apply the boundary type on every non-linear iteration

Deliverable: the behaviour the whole prototype exists for. Dirichlet fixes the nodes at zero pressure, Neumann frees them, and a switch made from outside takes effect on the next iteration.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`

**Interfaces:**
- Consumes: `GeoSeepageCondition::Initialize` (Task 4), `GetBoundaryType` / `SetBoundaryType` (Task 3).
- Produces: `void GeoSeepageCondition::InitializeNonLinearIteration(const ProcessInfo&)` — in Dirichlet mode sets `WATER_PRESSURE` to `0.0` and fixes the degree of freedom on every node; in Neumann mode frees it. Errors on an unknown boundary type.

- [x] **Step 1: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_geo_seepage_condition.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionDirichletFixesWaterPressureAtZero, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    // Give the nodes a non-zero pressure, to prove the condition overwrites it.
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 42.0;
    }

    condition.Initialize(ProcessInfo{});
    condition.SetBoundaryType(Geo::SeepageDirichletType);
    condition.InitializeNonLinearIteration(ProcessInfo{});

    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionNeumannFreesWaterPressure, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    condition.Initialize(ProcessInfo{});
    condition.SetBoundaryType(Geo::SeepageNeumannType);
    condition.InitializeNonLinearIteration(ProcessInfo{});

    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.IsFixed(WATER_PRESSURE))
    }
}

// This test plays the role of the Newton-Raphson strategy that will be built in step 2 of the
// prototype: it switches the boundary type between two non-linear iterations and checks that the
// condition picks the change up.
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionSwitchesBetweenNonLinearIterations, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    condition.Initialize(ProcessInfo{});

    condition.SetBoundaryType(Geo::SeepageDirichletType);
    condition.InitializeNonLinearIteration(ProcessInfo{});
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
    }

    condition.SetBoundaryType(Geo::SeepageNeumannType);
    condition.InitializeNonLinearIteration(ProcessInfo{});
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.IsFixed(WATER_PRESSURE))
    }

    condition.SetBoundaryType(Geo::SeepageDirichletType);
    condition.InitializeNonLinearIteration(ProcessInfo{});
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInitializeNonLinearIterationThrowsOnUnknownBoundaryType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    condition.Initialize(ProcessInfo{});
    // Bypass the validating setter, to mimic a value that came straight from an input file.
    condition.GetProperties().SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, "Robin");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.InitializeNonLinearIteration(ProcessInfo{}),
                                      "Unknown seepage boundary type 'Robin' for GeoSeepageCondition 1")
}
```

- [x] **Step 2: Build and verify the tests fail**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageConditionDirichletFixes*"
```

Expected: FAIL. `Condition::InitializeNonLinearIteration` does nothing by default, so `IsFixed(WATER_PRESSURE)` is false.

- [x] **Step 3: Declare `InitializeNonLinearIteration`**

In `geo_seepage_condition.h`, add directly below the `Initialize` declaration:

```cpp
    // Applies the currently configured boundary type to the nodes of this condition. Called once
    // per non-linear iteration, so that a boundary type switched from outside takes effect here.
    void InitializeNonLinearIteration(const ProcessInfo&) override;
```

- [x] **Step 4: Implement `InitializeNonLinearIteration`**

In `geo_seepage_condition.cpp`, add directly below the `Initialize` implementation:

```cpp
void GeoSeepageCondition::InitializeNonLinearIteration(const ProcessInfo&)
{
    KRATOS_TRY

    const auto& r_boundary_type = GetBoundaryType();
    KRATOS_ERROR_IF(r_boundary_type != Geo::SeepageDirichletType && r_boundary_type != Geo::SeepageNeumannType)
        << "Unknown seepage boundary type '" << r_boundary_type << "' for GeoSeepageCondition "
        << Id() << ", expected '" << Geo::SeepageDirichletType << "' or '"
        << Geo::SeepageNeumannType << "'" << std::endl;

    const auto is_dirichlet = r_boundary_type == Geo::SeepageDirichletType;
    for (auto& r_node : GetGeometry()) {
        if (is_dirichlet) {
            // A seepage face lets water out at atmospheric pressure.
            r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
            r_node.Fix(WATER_PRESSURE);
        } else {
            // Zero flux: freeing the degree of freedom is all that is needed, since a zero flux
            // makes no contribution to the system.
            r_node.Free(WATER_PRESSURE);
        }
    }

    KRATOS_CATCH("")
}
```

- [x] **Step 5: Build and verify the tests pass**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 17 tests.`

- [x] **Step 6: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp
git commit -m "Apply Dirichlet or Neumann boundary type per non-linear iteration"
```

---

### Task 6: Validate the setup in `Check`

Deliverable: a misconfigured seepage condition fails loudly at startup instead of behaving unexpectedly during the solve.

**Files:**
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h`
- Modify: `applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp`
- Modify: `applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp`

**Interfaces:**
- Consumes: everything from Tasks 1 to 5.
- Produces: `int GeoSeepageCondition::Check(const ProcessInfo&) const` — returns `0` on a valid setup, otherwise errors.

- [x] **Step 1: Write the failing tests**

Append inside `namespace Kratos::Testing` in `test_geo_seepage_condition.cpp`:

```cpp
KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckReturnsZeroForValidSetup, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    condition.Initialize(ProcessInfo{});
    condition.SetBoundaryType(Geo::SeepageDirichletType);

    KRATOS_EXPECT_EQ(condition.Check(ProcessInfo{}), 0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckThrowsWhenBoundaryTypeIsMissing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.Check(ProcessInfo{}),
                                      "GEO_SEEPAGE_BOUNDARY_TYPE is not defined for GeoSeepageCondition 1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckThrowsOnUnknownBoundaryType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition    = CreateSeepageCondition(r_model_part);

    condition.GetProperties().SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, "Robin");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.Check(ProcessInfo{}),
                                      "Unknown seepage boundary type 'Robin' for GeoSeepageCondition 1")
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
    condition.SetBoundaryType(Geo::SeepageDirichletType);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(condition.Check(ProcessInfo{}),
                                      "Missing degree of freedom for WATER_PRESSURE on node 1")
}
```

- [x] **Step 2: Build and verify the tests fail**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageConditionCheck*"
```

Expected: FAIL. The base `Condition::Check` returns 0 and raises none of the expected messages.

- [x] **Step 3: Declare `Check`**

In `geo_seepage_condition.h`, add directly below the `InitializeNonLinearIteration` declaration:

```cpp
    // Validates that this condition can be used: the geometry, the WATER_PRESSURE degrees of
    // freedom, and a legal GEO_SEEPAGE_BOUNDARY_TYPE.
    [[nodiscard]] int Check(const ProcessInfo& rCurrentProcessInfo) const override;
```

- [x] **Step 4: Implement `Check`**

In `geo_seepage_condition.cpp`, add directly below the `InitializeNonLinearIteration` implementation:

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

    const auto& r_boundary_type = GetBoundaryType();
    KRATOS_ERROR_IF(r_boundary_type != Geo::SeepageDirichletType && r_boundary_type != Geo::SeepageNeumannType)
        << "Unknown seepage boundary type '" << r_boundary_type << "' for GeoSeepageCondition "
        << Id() << ", expected '" << Geo::SeepageDirichletType << "' or '"
        << Geo::SeepageNeumannType << "'" << std::endl;

    return base_check_result;

    KRATOS_CATCH("")
}
```

- [x] **Step 5: Build and verify the whole feature's tests pass**

```powershell
cmake --build C:/checkouts/KratosProjects/dev/build/FullDebug --target KratosGeoMechanicsCoreTest
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe --gtest_filter="*GeoSeepageCondition*"
```

Expected: `[  PASSED  ] 21 tests.`

- [x] **Step 6: Run the full GeoMechanics C++ test suite to check for regressions**

```powershell
C:/checkouts/KratosProjects/dev/bin/FullDebug/test/KratosGeoMechanicsCoreTest.exe
```

Expected: no new failures compared to before this plan was started.

- [x] **Step 7: Commit**

```powershell
git add applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.h applications/GeoMechanicsApplication/custom_conditions/geo_seepage_condition.cpp applications/GeoMechanicsApplication/tests/cpp_tests/custom_conditions/test_geo_seepage_condition.cpp
git commit -m "Validate GeoSeepageCondition setup in Check"
```

---

## Findings Recorded During Implementation

**Status: all six tasks implemented, 21 feature tests green, full GeoMechanics C++ suite green (1128 tests, no regressions).**

Deviations from the plan as written:

- **`KRATOS_EXPECT_VECTOR_NEAR` cannot take a ublas zero-expression.** The plan's test used
  `KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, ZeroVector(2), 1.0e-12)`. That macro expands to
  `EXPECT_THAT(a, Pointwise(DoubleNear(tol), b))`, and gmock's `Pointwise` iterates its second
  argument. Iterating a `zero_vector` throws `"bad index"` from
  `boost/numeric/ublas/vector.hpp:1813`. Fixed by materializing the expected values into concrete
  `Matrix` / `Vector` objects first. Worth knowing for any future test in this application.
- Everything else went as planned. The per-condition `Properties` clone in `Initialize` behaves as
  designed: the failing test before the fix showed both conditions *and* the shared `Properties`
  object picking up a switch made on one condition, which is exactly the failure mode the clone
  exists to prevent.

Confirmed as designed:

- No `CMakeLists.txt` change was needed; the existing globs picked up both new `.cpp` files. But
  `kp config` had to be re-run each time a new `.cpp` was added, because the globs lack
  `CONFIGURE_DEPENDS`.
- `Geo::DofUtilities::ExtractDofsFromNodes` and `ExtractEquationIdsFrom` were reusable as-is.

## Open Questions to Answer in Steps 2 to 4 (feeding #14637)

These could not be settled at the level of a single condition and remain open:

- **Block vs. elimination builder and solver.** Still unverified. The condition changes dof fixity
  in `InitializeNonLinearIteration`; whether that is honoured without a rebuild depends on the
  builder and solver in use. This is the central question for step 4 and needs an integration-level
  test, not a unit test.
- **The detached `Properties` clone.** It is deliberately not registered in the `ModelPart`'s
  properties container. No trouble observed in the unit tests, but output processes and restart were
  not exercised.
- **Prescribing `WATER_PRESSURE = 0`.** Sufficient for the unit tests by construction, but whether
  the Muskat case needs an elevation-dependent value is unknown until the validation case is run.


