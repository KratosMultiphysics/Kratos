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

namespace
{

using namespace Kratos;

PointerVector<Node> CreateTwoNodes()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    return result;
}

PointerVector<Node> CreateNodesOnModelPart(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 1.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 0.0, 0.0, 0.0));
    return result;
}

PointerVector<Node> CreateCoincidentNodes()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 0.0, 0.0, 0.0));
    return result;
}

ModelPart& CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

Element::IndexType NextElementNumber(const ModelPart& rModelPart)
{
    return rModelPart.NumberOfElements() + 1;
}

intrusive_ptr<GeoSteadyStatePwPipingElement<2, 2>> CreateGeoSteadyStatePwPipingElementWithPWDofs(
    const ModelPart& rModelPart, const Properties::Pointer& rProperties, const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        NextElementNumber(rModelPart), rGeometry, rProperties);
    for (auto& node : p_result->GetGeometry()) {
        node.AddDof(WATER_PRESSURE);
    }

    return p_result;
}

intrusive_ptr<GeoSteadyStatePwPipingElement<2, 2>> CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(
    ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Line2D2<Node>>(CreateNodesOnModelPart(rModelPart));
    auto p_element = CreateGeoSteadyStatePwPipingElementWithPWDofs(rModelPart, rProperties, p_geometry);

    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<GeoSteadyStatePwPipingElement<2, 2>> CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithoutPWDofs(
    ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Line2D2<Node>>(CreateNodesOnModelPart(rModelPart));
    auto       p_element  = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        NextElementNumber(rModelPart), p_geometry, rProperties);

    rModelPart.AddElement(p_element);
    return p_element;
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
    const auto p_geometry   = std::make_shared<Line2D2<Node>>(CreateTwoNodes());
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
    const auto p_properties = std::make_shared<Properties>();

    // The source element needs to have a geometry, otherwise the version of the
    // Create method with a node input will fail.
    const auto p_geometry = std::make_shared<Line2D2<Node>>(CreateTwoNodes());
    const GeoSteadyStatePwPipingElement<2, 2> element(0, p_geometry, p_properties);

    // Act
    const auto p_created_element = element.Create(1, CreateTwoNodes(), p_properties);

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

    Model      model;
    auto&      r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    const auto p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 2);
    for (auto p_dof : degrees_of_freedom) {
        KRATOS_EXPECT_EQ(p_dof->GetVariable(), WATER_PRESSURE);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedEquationIdVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    unsigned int i = 0;
    for (const auto& node : p_element->GetGeometry()) {
        ++i;
        node.pGetDof(WATER_PRESSURE)->SetEquationId(i);
    }

    // Act
    const auto                    dummy_process_info = ProcessInfo{};
    Element::EquationIdVectorType equation_id_vector;
    p_element->EquationIdVector(equation_id_vector, dummy_process_info);

    // Assert
    const Element::EquationIdVectorType expected_ids = {1, 2};
    KRATOS_EXPECT_VECTOR_EQ(equation_id_vector, expected_ids)
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedIntegrationMethod,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const GeoSteadyStatePwPipingElement<2, 2> element;

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    const auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementCheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry   = std::make_shared<Line2D2<Node>>(CreateCoincidentNodes());
    const auto p_properties = std::make_shared<Properties>();
    const GeoSteadyStatePwPipingElement<2, 2> element(1, p_geometry, p_properties);

    // Act and Assert
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element.Check(dummy_process_info),
                                      "Error: DomainSize (0) is smaller than 1e-15 for element 1")

    const auto p_geometry2 = std::make_shared<Line2D2<Node>>(CreateTwoNodes());
    const GeoSteadyStatePwPipingElement<2, 2> element1(1, p_geometry2, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(element1.Check(dummy_process_info),
                                      "Error: Missing variable WATER_PRESSURE on node 1")

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_new_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithoutPWDofs(r_model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_new_element->Check(dummy_process_info),
        "Error: Missing degree of freedom for WATER_PRESSURE on node 1")

    auto p_element2 =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: DENSITY_WATER does not exist in the properties of element 2")
    p_element2->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: DENSITY_WATER (-1000) is not in the range [0,-> at element 2")
    p_element2->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: DYNAMIC_VISCOSITY does not exist in the properties of element 2")
    p_element2->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: DYNAMIC_VISCOSITY (-0.01) is not in the range [0,-> at element 2")
    p_element2->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: PIPE_HEIGHT does not exist in the properties of element 2")
    p_element2->GetProperties().SetValue(PIPE_HEIGHT, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "Error: PIPE_HEIGHT (-1) is not in the range [0,-> at element 2")
    p_element2->GetProperties().SetValue(PIPE_HEIGHT, 1.0);

    p_element2->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element2->Check(dummy_process_info),
                                      "Error: Node with non-zero Z coordinate found. Id: 1")
    p_element2->GetGeometry().begin()->Z() = 0;

    // No exceptions on correct input
    KRATOS_EXPECT_EQ(p_element2->Check(dummy_process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsTheExpectedLeftHandSideAndRightHandSide,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);
    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_element->SetValue(PIPE_HEIGHT, 1.0E-1);
    // Set gravity perpendicular to the line ( so no fluid body flow vector from this )
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    // Create a head gradient of -10.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = 10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    auto expected_left_hand_side  = Matrix{ScalarMatrix{2, 2, 0.1 / 12.}};
    expected_left_hand_side(0, 0) = -expected_left_hand_side(0, 0);
    expected_left_hand_side(1, 1) = -expected_left_hand_side(1, 1);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{ScalarVector{2, 1. / 12.}};
    expected_right_hand_side(1)   = -expected_right_hand_side(1);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementHasValuesAfterInitialize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    // Act
    p_element->Initialize(dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(p_element->GetValue(PIPE_ELEMENT_LENGTH), 1.);
    KRATOS_EXPECT_EQ(p_element->GetValue(PIPE_EROSION), false);
    const double quite_small = 1.E-10;
    KRATOS_EXPECT_EQ(p_element->GetValue(PIPE_HEIGHT), quite_small);
    KRATOS_EXPECT_EQ(p_element->GetValue(PREV_PIPE_HEIGHT), quite_small);
    KRATOS_EXPECT_EQ(p_element->GetValue(DIFF_PIPE_HEIGHT), 0.);
    KRATOS_EXPECT_EQ(p_element->GetValue(PIPE_ACTIVE), false);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsEquilibriumHeightForHeadGradient,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);
    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    // Set gravity perpendicular to the line ( so no fluid body flow vector from this )
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    // Create a head gradient of 0.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;

    // Act
    auto pipe_height = p_element->CalculateEquilibriumPipeHeight(p_element->GetProperties(),
                                                                 p_element->GetGeometry(), 0.);

    // Assert
    const double infinite_pipe_height = 1.0e10;
    KRATOS_EXPECT_DOUBLE_EQ(pipe_height, infinite_pipe_height);

    // Create a head gradient of 1.E-3.
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    // Add other necessary material parameters
    // DENSITY_SOLID such that rho_s/rho_w - 1. = 1.
    p_element->GetProperties().SetValue(DENSITY_SOLID, 2.0E3);
    // PIPE_ETA such that multiplication with Pi/3 becomes 1.
    p_element->GetProperties().SetValue(PIPE_ETA, 3.0 / Globals::Pi);
    p_element->GetProperties().SetValue(PIPE_MODEL_FACTOR, 1.0);
    // slope = 0. PIPE_THETA = 45 deg. such that sin(theta + slope) / cos(theta) becomes 1.
    p_element->GetProperties().SetValue(PIPE_THETA, 45.0);
    // no PIPE_MODIFIED_D so the equilibrium height becomes PIPE_D_70 / |dh/dx| = 7.e-3 / |-1.e-3| = 7.
    p_element->GetProperties().SetValue(PIPE_D_70, 7.e-3);

    // Act
    pipe_height = p_element->CalculateEquilibriumPipeHeight(p_element->GetProperties(),
                                                            p_element->GetGeometry(), 0.);

    // Assert
    KRATOS_EXPECT_NEAR(pipe_height, 7.0, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsPipeActive, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    p_element->SetValue(PIPE_ACTIVE, false);

    std::vector<bool> pipe_active(
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod()));
    p_element->CalculateOnIntegrationPoints(PIPE_ACTIVE, pipe_active, dummy_process_info);
    // Assert
    KRATOS_EXPECT_FALSE(pipe_active[0])

    p_element->SetValue(PIPE_ACTIVE, true);
    p_element->CalculateOnIntegrationPoints(PIPE_ACTIVE, pipe_active, dummy_process_info);
    // Assert
    KRATOS_EXPECT_TRUE(pipe_active[0])
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsPipeHeight, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    const double pipe_height = 1.234E-5;
    p_element->SetValue(PIPE_HEIGHT, pipe_height);

    std::vector<double> pipe_heights(
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod()));
    p_element->CalculateOnIntegrationPoints(PIPE_HEIGHT, pipe_heights, dummy_process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(pipe_heights[0], pipe_height);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsPermeabilityMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    const double pipe_height = 1.E-1;
    p_element->SetValue(PIPE_HEIGHT, pipe_height);

    std::vector<Matrix> permeability_matrices(
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod()));
    p_element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, permeability_matrices, dummy_process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(permeability_matrices[0](0, 0), std::pow(pipe_height, 3) / 12.0);

    p_element->CalculateOnIntegrationPoints(LOCAL_PERMEABILITY_MATRIX, permeability_matrices, dummy_process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(permeability_matrices[0](0, 0), std::pow(pipe_height, 3) / 12.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePwPipingElementReturnsFluidFluxVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);
    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_element->GetProperties().SetValue(DENSITY_SOLID, 2.0E3);
    // PIPE_ETA such that multiplication with Pi/3 becomes 1.
    p_element->GetProperties().SetValue(PIPE_ETA, 3.0 / Globals::Pi);
    p_element->GetProperties().SetValue(PIPE_MODEL_FACTOR, 1.0);
    // slope = 0. PIPE_THETA = 45 deg. such that sin(theta + slope) / cos(theta) becomes 1.
    p_element->GetProperties().SetValue(PIPE_THETA, 45.0);
    // no PIPE_MODIFIED_D so the equilibrium height becomes PIPE_D_70 / |dh/dx| = 7.e-3 / |-1.e-3| = 7.
    p_element->GetProperties().SetValue(PIPE_D_70, 7.e-3);

    // Set gravity perpendicular to the line ( so no fluid body flow vector from this )
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) =
        array_1d<double, 3>{0.0, -10.0, 0.0};
    // Create a head gradient of -10/(density_water*|volume_acceleration|) = -1.E-3.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    p_element->SetValue(PIPE_HEIGHT, 1.0E-1);

    // Act
    std::vector<array_1d<double, 3>> fluid_fluxes(
        p_element->GetGeometry().IntegrationPointsNumber(p_element->GetIntegrationMethod()));
    p_element->CalculateOnIntegrationPoints(LOCAL_FLUID_FLUX_VECTOR, fluid_fluxes, dummy_process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(fluid_fluxes[0](0), 1.E-3 * 0.1 / 12.0);

    // element tangetial axis is opposite global X, so global flux has negative X component
    p_element->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, fluid_fluxes, dummy_process_info);
    KRATOS_EXPECT_DOUBLE_EQ(fluid_fluxes[0](0), -1.E-3 * 0.1 / 12.0);
}

} // namespace Kratos::Testing