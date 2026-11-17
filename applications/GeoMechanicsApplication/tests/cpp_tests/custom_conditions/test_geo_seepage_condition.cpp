// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Wijtze Pieter Kikstra

#include <string>

#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
#include "geometries/line_2d_2.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;
using namespace std::string_literals;

namespace
{

// Creates a model part with two nodes that have a WATER_PRESSURE degree of freedom.
ModelPart& CreateModelPartWithTwoWaterPressureNodes(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("Main"s);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return r_model_part;
}

// Creates a two-noded seepage condition on the given model part, with its own Properties.
GeoSeepageCondition CreateSeepageCondition(ModelPart& rModelPart)
{
    auto p_geometry = std::make_shared<Line2D2<Node>>(rModelPart.pGetNode(1), rModelPart.pGetNode(2));
    auto p_properties = std::make_shared<Properties>(1);
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
    auto& r_model_part = model.CreateModelPart("Main"s);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    // Deliberately do not add the WATER_PRESSURE degree of freedom.

    auto condition = CreateSeepageCondition(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN([[maybe_unused]] const auto result = condition.Check(ProcessInfo{}),
                                      "Missing degree of freedom for WATER_PRESSURE on node 1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionCheckThrowsWhenItHasMoreThanOneNeighbouringElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model              = Model{};
    auto& r_model_part       = CreateModelPartWithTwoWaterPressureNodes(model);
    auto  condition          = CreateSeepageCondition(r_model_part);
    auto  neighbour_elements = GlobalPointersVector<Element>{};
    // Add two null pointers, to pretend the condition has two neighbouring elements
    neighbour_elements.push_back({});
    neighbour_elements.push_back({});
    condition.SetValue(NEIGHBOUR_ELEMENTS, neighbour_elements);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN([[maybe_unused]] const auto result = condition.Check(ProcessInfo{}), "The seepage condition with ID 1 has more than one neighbouring element, which is not allowed")
}

} // namespace Kratos::Testing
