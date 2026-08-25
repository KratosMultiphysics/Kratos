// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include <string>

// Project includes
#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
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

    // Materialize the expected values: the gmock matcher behind KRATOS_EXPECT_VECTOR_NEAR
    // iterates its argument, which a ublas zero-expression does not support.
    const auto expected_left_hand_side  = Matrix{ZeroMatrix{2, 2}};
    const auto expected_right_hand_side = Vector{ZeroVector{2}};

    KRATOS_EXPECT_MATRIX_NEAR(left_hand_side, expected_left_hand_side, 1.0e-12)
    KRATOS_EXPECT_VECTOR_NEAR(right_hand_side, expected_right_hand_side, 1.0e-12)
}

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

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionInitializeGivesEachConditionItsOwnProperties, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithTwoWaterPressureNodes(model);

    auto p_geometry = std::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
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

    auto p_geometry = std::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(2));
    auto p_properties = Kratos::make_shared<Properties>(1);
    p_properties->SetValue(GEO_SEEPAGE_BOUNDARY_TYPE, Geo::SeepageDirichletType);

    auto condition = GeoSeepageCondition{1, p_geometry, p_properties};

    condition.Initialize(ProcessInfo{});
    condition.SetBoundaryType(Geo::SeepageNeumannType);
    condition.Initialize(ProcessInfo{});

    KRATOS_EXPECT_EQ(condition.GetBoundaryType(), Geo::SeepageNeumannType);
}

} // namespace Kratos::Testing




