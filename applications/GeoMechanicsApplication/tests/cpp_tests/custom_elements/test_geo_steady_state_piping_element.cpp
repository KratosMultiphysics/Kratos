// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

// #include <boost/numeric/ublas/assignment.hpp>
// #include <cstddef>

namespace
{

using namespace Kratos;

PointerVector<Node> CreateNodes()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    return result;
}

PointerVector<Node> CreateIdenticalNodes()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    return result;
}

ModelPart& CreateModelPartWithWaterPressureVariable(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);

    return r_result;
}

GeoSteadyStatePwPipingElement<2, 2> CreateGeoSteadyStatePwPipingElementWithPWDofs(
    const Properties::Pointer& rProperties, const Geometry<Node>::Pointer& rGeometry)
{
    auto result = GeoSteadyStatePwPipingElement<2, 2>{1, rGeometry, rProperties};
    for (auto& node : result.GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }

    return result;
}

GeoSteadyStatePwPipingElement<2, 2> CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithWaterPressureVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    return CreateGeoSteadyStatePwPipingElementWithPWDofs(rProperties, p_geometry);
}

GeoSteadyStatePwPipingElement<2, 2> CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithoutPWDofs(
    Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithWaterPressureVariable(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    return GeoSteadyStatePwPipingElement<2, 2>{1, p_geometry, rProperties};
}

} // namespace

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const GeoSteadyStatePwPipingElement<2, 2> element;
    auto p_casted_element = dynamic_cast<const Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(p_casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementCanCreateInstanceWithGeometryInput,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const GeoSteadyStatePwPipingElement<2, 2> element;
    const auto p_geometry   = std::make_shared<Line2D2<Node>>(CreateNodes());
    const auto p_properties = std::make_shared<Properties>();

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementCanCreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodes();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<Line2D2<Node>>(CreateNodes());
    const GeoSteadyStatePwPipingElement<2, 2> element(0, p_geometry, p_properties);

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    const auto element = CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(model, p_properties);

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    element.GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 2);
    for (auto dof : degrees_of_freedom) {
        KRATOS_EXPECT_EQ(dof->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedEquationIdVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto element = CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(model, p_properties);

    unsigned int i = 0;
    for (const auto& node : element.GetGeometry()) {
        ++i;
        node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    element.EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const std::vector<int> expected_ids = {1, 2};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedIntegrationMethod,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const GeoSteadyStatePwPipingElement<2, 2> element;
    const auto p_geometry        = std::make_shared<Line2D2<Node>>(CreateNodes());
    const auto p_properties      = std::make_shared<Properties>();
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Act
    const auto p_integration_method = p_created_element->GetIntegrationMethod();

    // Assert
    const auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementCheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const GeoSteadyStatePwPipingElement<2, 2> element;
    auto       p_geometry        = std::make_shared<Line2D2<Node>>(CreateIdenticalNodes());
    const auto p_properties      = std::make_shared<Properties>();
    auto       p_created_element = element.Create(1, p_geometry, p_properties);

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_created_element->Check(dummy_process_info),
                                      "Error: DomainSize (0) is smaller than 1e-15 for element 1")

    p_geometry        = std::make_shared<Line2D2<Node>>(CreateNodes());
    p_created_element = element.Create(2, p_geometry, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_created_element->Check(dummy_process_info),
                                      "Error: Missing variable WATER_PRESSURE on node 1")

    Model model;
    auto new_element = CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithoutPWDofs(model, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(new_element.Check(dummy_process_info),
                                      "Error: Missing degree of freedom for WATER_PRESSURE on node 1")

    Model model1;
    auto element1 = CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(model1, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: DENSITY_WATER does not exist in the properties of element 1")
    element1.GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: DENSITY_WATER (-1000) has an invalid value at element 1")
    element1.GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: DYNAMIC_VISCOSITY does not exist in the properties of element 1")
    element1.GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: DYNAMIC_VISCOSITY (-0.01) has an invalid value at element 1")
    element1.GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: PIPE_HEIGHT does not exist in the properties of element 1")
    element1.GetProperties().SetValue(PIPE_HEIGHT, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: PIPE_HEIGHT (-1) has an invalid value at element 1")
    element1.GetProperties().SetValue(PIPE_HEIGHT, 1.0);

}

} // namespace Kratos::Testing