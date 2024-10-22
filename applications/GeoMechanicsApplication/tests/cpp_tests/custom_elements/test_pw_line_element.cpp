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
//

#include "custom_elements/transient_Pw_line_element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

PointerVector<Node> CreateNodes1()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 1.0, 1.0, 0.0));
    return result;
}

PointerVector<Node> CreateIdenticalNodes1()
{
    PointerVector<Node> result;
    result.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(Kratos::make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    return result;
}

ModelPart& CreateModelPartWithWaterPressureVariableAndVolumeAcceleration1(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

TransientPwLineElement<2, 2> TransientPwLineElementWithPWDofs(const Properties::Pointer& rProperties,
                                                              const Geometry<Node>::Pointer& rGeometry)
{
    auto result = TransientPwLineElement<2, 2>{1, rGeometry, rProperties};
    for (auto& node : result.GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
        node.AddDof(DT_WATER_PRESSURE);
    }

    return result;
}

TransientPwLineElement<2, 2> TransientPwLineElementWithPWDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration1(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(1, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    return TransientPwLineElementWithPWDofs(rProperties, p_geometry);
}

TransientPwLineElement<2, 2> TransientPwLineElementWithoutPWDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration1(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    return TransientPwLineElement<2, 2>{1, p_geometry, rProperties};
}

} // namespace

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementIsAnElement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const TransientPwLineElement<2, 2> element;
    auto                               p_casted_element = dynamic_cast<const Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(p_casted_element, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementCanCreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const TransientPwLineElement<2, 2> element;
    const auto                         p_geometry = std::make_shared<Line2D2<Node>>(CreateNodes1());
    const auto                         p_properties = std::make_shared<Properties>();

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementCanCreateInstanceWithNodeInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto nodes        = CreateNodes1();
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto                         p_geometry = std::make_shared<Line2D2<Node>>(CreateNodes1());
    const TransientPwLineElement<2, 2> element(0, p_geometry, p_properties);

    // Act
    const auto p_created_element = element.Create(1, nodes, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementReturnsTheExpectedDoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model      model;
    const auto element = TransientPwLineElementWithPWDofs(model, p_properties);

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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementReturnsTheExpectedEquationIdVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto  element = TransientPwLineElementWithPWDofs(model, p_properties);

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

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementReturnsTheExpectedIntegrationMethod, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const TransientPwLineElement<2, 2> element;
    const auto                         p_geometry = std::make_shared<Line2D2<Node>>(CreateNodes1());
    const auto                         p_properties = std::make_shared<Properties>();
    const auto p_created_element                    = element.Create(1, p_geometry, p_properties);

    // Act
    const auto p_integration_method = p_created_element->GetIntegrationMethod();

    // Assert
    const auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElementReturnsTheExpectedLeftHandSideAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       process_info = ProcessInfo{};
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0); // small to get a significant value for the compressibility term
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    element.GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_properties->SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(RESIDUAL_SATURATION, 0.000000e+00);
    p_properties->SetValue(CROSS_AREA, 1.0);
    p_properties->SetValue(IGNORE_UNDRAINED, false);

    process_info[DT_PRESSURE_COEFFICIENT] = 1.5;

    Model model;
    auto  element = TransientPwLineElementWithPWDofs(model, p_properties);
    element.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    element.GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    // Create a head gradient of -10.
    element.GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE)    = 10.0;
    element.GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE)    = 0.0;
    element.GetGeometry()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 4.0;
    element.GetGeometry()[1].FastGetSolutionStepValue(DT_WATER_PRESSURE) = 5.0;

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    element.Initialize(process_info);
    element.CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, process_info);

    // Assert
    auto expected_left_hand_side = Matrix{2, 2, 0.0006423358000298597};
    expected_left_hand_side <<= -0.00099588919125952972, 0.00046555910441502474,
        0.00046555910441502474, -0.00099588919125952972;
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{2};
    expected_right_hand_side <<= 6.43131, -6.42813;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}
} // namespace Kratos::Testing