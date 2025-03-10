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

#include "containers/variable.h"
#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/plane_strain_stress_state.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

ModelPart& CreateModelPartWithUPwSolutionStepVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.AddNodalSolutionStepVariable(VELOCITY);
    r_result.AddNodalSolutionStepVariable(ACCELERATION);
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_result.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);

    return r_result;
}

std::shared_ptr<Properties> SetProperties()
{
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<GeoIncrementalLinearElasticLaw>(
                                                 std::make_unique<PlaneStrain>()));
    p_properties->SetValue(YOUNG_MODULUS, 1.000000e+07);
    p_properties->SetValue(POISSON_RATIO, 0.000000e+00);
    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    p_properties->SetValue(POROSITY, 1.000000e-01);
    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    p_properties->SetValue(BULK_MODULUS_FLUID, 200.0); // small to get a significant value for the compressibility term
    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_YY, 9.084000e-06);
    p_properties->SetValue(PERMEABILITY_XY, 0.000000e+00);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    // Biot alpha = 0, no coupling
    p_properties->SetValue(BIOT_COEFFICIENT, 0.000000e+00);
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);
    p_properties->SetValue(IGNORE_UNDRAINED, false);

    return p_properties;
}

void SetSolutionStepValues(Kratos::intrusive_ptr<Kratos::UPwSmallStrainElement<2, 3>>& rElement)
{
    const auto zero_array_1d = array_1d<double, 3>{0.0, 0.0, 0.0};
    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = zero_array_1d;
        r_node.FastGetSolutionStepValue(VELOCITY)     = zero_array_1d;
        // Zero acceleration -> no Fluid Body Flow
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = zero_array_1d;
        // Zero pressure gradient -> no permeability flow
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)    = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE) = 0.0;
    }
    rElement->GetGeometry()[2].FastGetSolutionStepValue(WATER_PRESSURE) = 2.0E4;
}

template <unsigned int TDim, unsigned int TNumNodes>
intrusive_ptr<UPwSmallStrainElement<TDim, TNumNodes>> UPwSmallStrainElementWithUPwDofs(
    const Properties::Pointer& rProperties, const Geometry<Node>::Pointer& rGeometry)
{
    auto p_result = make_intrusive<UPwSmallStrainElement<TDim, TNumNodes>>(
        1, rGeometry, rProperties, std::make_unique<PlaneStrainStressState>());
    for (auto& r_node : p_result->GetGeometry()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        if constexpr (TDim == 3) r_node.AddDof(DISPLACEMENT_Z);
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        if constexpr (TDim == 3) r_node.AddDof(VELOCITY_Z);
        r_node.AddDof(WATER_PRESSURE);
        r_node.AddDof(DT_WATER_PRESSURE);
    }

    return p_result;
}

auto UPwSmallStrainElementWithUPwDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithUPwSolutionStepVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(nodes);
    return UPwSmallStrainElementWithUPwDofs<2, 3>(rProperties, p_geometry);
}

} // namespace

namespace Kratos::Testing
{
using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CreateInstanceWithGeometryInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry =
        std::make_shared<Triangle2D3<Node>>(ModelSetupUtilities::Create2D3NTriangleGeometry());
    const auto                        p_properties = std::make_shared<Properties>();
    const UPwSmallStrainElement<2, 3> element(0, p_geometry, p_properties,
                                              std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_created_element = element.Create(1, p_geometry, p_properties);

    // Assert
    EXPECT_NE(p_created_element, nullptr);
    EXPECT_EQ(p_created_element->Id(), 1);
    EXPECT_NE(p_created_element->pGetGeometry(), nullptr);
    EXPECT_NE(p_created_element->pGetProperties(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_DoFList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model      model;
    const auto p_element = UPwSmallStrainElementWithUPwDofs(model, std::make_shared<Properties>());

    // Act
    const auto              dummy_process_info = ProcessInfo{};
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, dummy_process_info);

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 9);
    std::vector<std::string> variable_names;
    std::transform(degrees_of_freedom.cbegin(), degrees_of_freedom.cend(), std::back_inserter(variable_names),
                   [](const auto& rpDof) { return rpDof->GetVariable().Name(); });
    const std::vector<std::string> desired_variable_list{
        "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_X",
        "DISPLACEMENT_Y", "WATER_PRESSURE", "WATER_PRESSURE", "WATER_PRESSURE"};

    KRATOS_EXPECT_EQ(variable_names, desired_variable_list);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_IntegrationMethod, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry =
        std::make_shared<Triangle2D3<Node>>(ModelSetupUtilities::Create2D3NTriangleGeometry());
    const auto                        p_properties = std::make_shared<Properties>();
    const UPwSmallStrainElement<2, 3> element(0, p_geometry, p_properties,
                                              std::make_unique<PlaneStrainStressState>());

    // Act
    const auto p_integration_method = element.GetIntegrationMethod();

    // Assert
    constexpr auto expected_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    KRATOS_EXPECT_EQ(p_integration_method, expected_integration_method);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementCheckDoesNotThrowOnCorrectInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info                     = ProcessInfo{};
    process_info[DT_PRESSURE_COEFFICIENT] = 1.0;
    process_info[VELOCITY_COEFFICIENT]    = 1.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    SetSolutionStepValues(element);

    // Act, no exceptions on correct input
    KRATOS_EXPECT_EQ(element->Check(process_info), 0);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementCalculatesSteadyStateRightHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 0.0;
    process_info[VELOCITY_COEFFICIENT]    = 0.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    SetSolutionStepValues(element);

    // Act
    element->Initialize(process_info);
    Vector actual_right_hand_side;
    element->CalculateRightHandSide(actual_right_hand_side, process_info);

    // Assert
    Vector expected_right_hand_side = ZeroVector(9);
    expected_right_hand_side[7]     = -4.542;
    expected_right_hand_side[8]     = +4.542;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance);
    std::vector<array_1d<double, 3>> calculated_fluid_flux_at_integration_points;
    element->CalculateOnIntegrationPoints(
        FLUID_FLUX_VECTOR, calculated_fluid_flux_at_integration_points, process_info);
    for (auto i = std::size_t{0}; i < calculated_fluid_flux_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_fluid_flux_at_integration_points[i],
                                           (array_1d<double, 3>{0., 9.084, 0.}), Defaults::relative_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementCalculatesSteadyStateLeftHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 0.0;
    process_info[VELOCITY_COEFFICIENT]    = 0.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    SetSolutionStepValues(element);

    // Act
    element->Initialize(process_info);
    Matrix actual_left_hand_side;
    element->CalculateLeftHandSide(actual_left_hand_side, process_info);

    // Assert
    Matrix expected_left_hand_side(9, 9);
    // clang-format off
    expected_left_hand_side <<= 5000000,0,-5000000,0,0,0,0,0,0,
                                0,2500000,2500000,-2500000,-2500000,0,0,0,0,
                               -5000000,2500000,7500000,-2500000,-2500000,0,0,0,0,
                                0,-2500000,-2500000,7500000,2500000,-5000000,0,0,0,
                                0,-2500000,-2500000,2500000,2500000,0,0,0,0,
                                0,0,0,-5000000,0,5000000,0,0,0,
                                0,0,0,0,0,0,-0.00045419999999999998,0.00045419999999999998,0,
                                0,0,0,0,0,0,0.00045419999999999998,-0.00090839999999999996,0.00045419999999999998,
                                0,0,0,0,0,0,0,0.00045419999999999998,-0.00045419999999999998;
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementInitializeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 0.0;
    process_info[VELOCITY_COEFFICIENT]    = 0.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    SetSolutionStepValues(element);

    // Act, no exceptions on correct input
    element->Initialize(process_info);
    element->InitializeSolutionStep(process_info);

    // Assert
    for (const auto& r_node : element->GetGeometry()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    }
}

} // namespace Kratos::Testing