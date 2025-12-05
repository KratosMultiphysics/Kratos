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

#include "containers/model.h"
#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "includes/expect.h"
#include "includes/model_part.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
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

intrusive_ptr<GeoSteadyStatePwPipingElement<3, 2>> CreateGeoSteadyStatePwPipingElement3D2NWithPWDofs(
    const ModelPart& rModelPart, const Properties::Pointer& rProperties, const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result = make_intrusive<GeoSteadyStatePwPipingElement<3, 2>>(
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

intrusive_ptr<GeoSteadyStatePwPipingElement<3, 2>> CreateHorizontalUnitLengthGeoSteadyStatePwPipingElement3D2NWithPWDofs(
    ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    const auto p_geometry = std::make_shared<Line3D2<Node>>(CreateNodesOnModelPart(rModelPart));
    auto p_element = CreateGeoSteadyStatePwPipingElement3D2NWithPWDofs(rModelPart, rProperties, p_geometry);

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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementIsAnElement)
{
    const GeoSteadyStatePwPipingElement<2, 2> element;
    auto p_casted_element = dynamic_cast<const Element*>(&element);
    KRATOS_CHECK_NOT_EQUAL(p_casted_element, nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementCanCreateInstanceWithGeometryInput)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementCanCreateInstanceWithNodeInput)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsTheExpectedDoFList)
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
    EXPECT_EQ(degrees_of_freedom.size(), 2);
    for (auto p_dof : degrees_of_freedom) {
        EXPECT_EQ(p_dof->GetVariable(), WATER_PRESSURE);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsTheExpectedEquationIdVector)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsTheExpectedIntegrationMethod)
{
    // Arrange
    const GeoSteadyStatePwPipingElement<2, 2> element;

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    const auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    EXPECT_EQ(p_integration_method, expected_integration_method);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementCheckThrowsOnFaultyInput)
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
                                      "Missing variable WATER_PRESSURE on nodes 1 2")

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_new_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithoutPWDofs(r_model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_new_element->Check(dummy_process_info),
        "Missing the DoF for the variable WATER_PRESSURE on nodes 1 2")

    auto p_element2 =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element2->Check(dummy_process_info),
                                      "DENSITY_WATER does not exist in the material properties at "
                                      "element with Id 0 at element with Id 2.")
    p_element2->GetProperties().SetValue(DENSITY_WATER, -1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "DENSITY_WATER in the material properties at element with Id 0 at element with Id 2 has an "
        "invalid value: -1000 is out of the range [0, -).")
    p_element2->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element2->Check(dummy_process_info),
                                      "DYNAMIC_VISCOSITY does not exist in the material properties "
                                      "at element with Id 0 at element with Id 2.")
    p_element2->GetProperties().SetValue(DYNAMIC_VISCOSITY, -1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "DYNAMIC_VISCOSITY in the material properties at element with Id 0 at element with Id 2 "
        "has an invalid value: -0.01 is out of the range [0, -).")
    p_element2->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element2->Check(dummy_process_info),
                                      "PIPE_HEIGHT does not exist in the material properties at "
                                      "element with Id 0 at element with Id 2.")
    p_element2->GetProperties().SetValue(PIPE_HEIGHT, -1.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element2->Check(dummy_process_info),
        "PIPE_HEIGHT in the material properties at element with Id 0 at element with Id 2 has an "
        "invalid value: -1 is out of the range [0, -).")
    p_element2->GetProperties().SetValue(PIPE_HEIGHT, 1.0);

    p_element2->GetGeometry().begin()->Z() += 1;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element2->Check(dummy_process_info),
                                      "Node with Id: 1 has non-zero Z coordinate.")
    p_element2->GetGeometry().begin()->Z() = 0;

    // No exceptions on correct input
    EXPECT_EQ(p_element2->Check(dummy_process_info), 0);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsTheExpectedLeftHandSideAndRightHandSide)
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
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       GeoSteadyStatePwPipingElement3D2NReturnsTheExpectedLeftHandSideAndRightHandSide)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElement3D2NWithPWDofs(r_model_part, p_properties);
    p_element->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    p_element->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    p_element->GetProperties().SetValue(PIPE_WIDTH_FACTOR, 1.5);
    p_element->SetValue(PIPE_HEIGHT, 1.0E-1);
    // Set gravity perpendicular to the line ( so no fluid body flow vector from this )
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    // Create a head gradient of -10.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = 10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    // Act
    Vector actual_right_hand_side;
    Matrix actual_left_hand_side;
    p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);

    // Assert
    auto expected_left_hand_side  = Matrix{ScalarMatrix{2, 2, 0.015 / 12.}};
    expected_left_hand_side(0, 0) = -expected_left_hand_side(0, 0);
    expected_left_hand_side(1, 1) = -expected_left_hand_side(1, 1);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance)

    auto expected_right_hand_side = Vector{ScalarVector{2, 0.15 / 12.}};
    expected_right_hand_side(1)   = -expected_right_hand_side(1);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementHasValuesAfterInitialize)
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
    EXPECT_EQ(p_element->GetValue(PIPE_ELEMENT_LENGTH), 1.);
    EXPECT_EQ(p_element->GetValue(PIPE_EROSION), false);
    const double quite_small = 1.E-10;
    EXPECT_EQ(p_element->GetValue(PIPE_HEIGHT), quite_small);
    EXPECT_EQ(p_element->GetValue(PREV_PIPE_HEIGHT), quite_small);
    EXPECT_EQ(p_element->GetValue(DIFF_PIPE_HEIGHT), 0.);
    EXPECT_EQ(p_element->GetValue(PIPE_ACTIVE), false);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsEquilibriumHeightForHeadGradient)
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
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    // Create a head gradient of 0.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;

    // Act
    auto pipe_height = p_element->CalculateEquilibriumPipeHeight(p_element->GetProperties(),
                                                                 p_element->GetGeometry(), 0.);

    // Assert
    const double infinite_pipe_height = 1.0e10;
    EXPECT_DOUBLE_EQ(pipe_height, infinite_pipe_height);

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
    EXPECT_NEAR(pipe_height, 7.0, 1e-10);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsPipeActive)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    p_element->SetValue(PIPE_ACTIVE, false);

    std::vector<bool> pipe_active_states;
    p_element->CalculateOnIntegrationPoints(PIPE_ACTIVE, pipe_active_states, dummy_process_info);
    // Assert
    EXPECT_EQ(pipe_active_states, (std::vector<bool>{false, false}));

    p_element->SetValue(PIPE_ACTIVE, true);
    pipe_active_states.clear();
    p_element->CalculateOnIntegrationPoints(PIPE_ACTIVE, pipe_active_states, dummy_process_info);
    // Assert
    EXPECT_EQ(pipe_active_states, (std::vector<bool>{true, true}));
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsPipeHeight)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    constexpr auto pipe_height = 1.234E-5;
    p_element->SetValue(PIPE_HEIGHT, pipe_height);

    std::vector<double> pipe_heights;
    p_element->CalculateOnIntegrationPoints(PIPE_HEIGHT, pipe_heights, dummy_process_info);

    // Assert
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
        pipe_heights, (std::vector<double>{pipe_height, pipe_height}), Defaults::relative_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsPermeabilityMatrix)
{
    // Arrange
    const auto dummy_process_info = ProcessInfo{};
    const auto p_properties       = std::make_shared<Properties>();

    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    auto  p_element =
        CreateHorizontalUnitLengthGeoSteadyStatePwPipingElementWithPWDofs(r_model_part, p_properties);

    constexpr auto pipe_height = 1.E-1;
    p_element->SetValue(PIPE_HEIGHT, pipe_height);

    std::vector<Matrix> permeability_matrices;
    p_element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, permeability_matrices, dummy_process_info);

    // Assert
    EXPECT_EQ(permeability_matrices.size(), 2);
    const auto expected_permeability_matrix = Matrix{1, 1, std::pow(pipe_height, 3) / 12.0};
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(permeability_matrices[0], expected_permeability_matrix,
                                       Defaults::relative_tolerance);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(permeability_matrices[1], expected_permeability_matrix,
                                       Defaults::relative_tolerance);

    permeability_matrices.clear();
    p_element->CalculateOnIntegrationPoints(LOCAL_PERMEABILITY_MATRIX, permeability_matrices, dummy_process_info);

    // Assert
    EXPECT_EQ(permeability_matrices.size(), 2);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(permeability_matrices[0], expected_permeability_matrix,
                                       Defaults::relative_tolerance);
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(permeability_matrices[1], expected_permeability_matrix,
                                       Defaults::relative_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeoSteadyStatePwPipingElementReturnsFluidFluxVector)
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
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    p_element->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
    // Create a head gradient of -10/(density_water*|volume_acceleration|) = -1.E-3.
    p_element->GetGeometry()[0].FastGetSolutionStepValue(WATER_PRESSURE) = -10.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    p_element->SetValue(PIPE_HEIGHT, 1.0E-1);

    // Act
    std::vector<array_1d<double, 3>> fluid_fluxes;
    p_element->CalculateOnIntegrationPoints(LOCAL_FLUID_FLUX_VECTOR, fluid_fluxes, dummy_process_info);

    // Assert
    auto expected_fluid_flux_array = array_1d<double, 3>{ZeroVector{3}};
    expected_fluid_flux_array[0]   = 1.E-3 * 0.1 / 12.0;
    EXPECT_EQ(fluid_fluxes.size(), 2);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(fluid_fluxes[0], expected_fluid_flux_array, Defaults::relative_tolerance);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(fluid_fluxes[1], expected_fluid_flux_array, Defaults::relative_tolerance);

    // element tangential axis is opposite global X, so global flux has negative X component
    fluid_fluxes.clear();
    p_element->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, fluid_fluxes, dummy_process_info);
    expected_fluid_flux_array[0] = -1.E-3 * 0.1 / 12.0;
    EXPECT_EQ(fluid_fluxes.size(), 2);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(fluid_fluxes[0], expected_fluid_flux_array, Defaults::relative_tolerance);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(fluid_fluxes[1], expected_fluid_flux_array, Defaults::relative_tolerance);
}

} // namespace Kratos::Testing