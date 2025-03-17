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

void SetSolutionStepValuesForFluidFluxCheck(const intrusive_ptr<UPwSmallStrainElement<2, 3>>& rElement)
{
    const auto zero_values = array_1d<double, 3>{0.0, 0.0, 0.0};
    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = zero_values;
        r_node.FastGetSolutionStepValue(VELOCITY)     = zero_values;
        // Zero acceleration -> no Fluid Body Flow
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = zero_values;
        // Zero pressure gradient -> no permeability flow
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)    = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE) = 0.0;
    }
    rElement->GetGeometry()[2].FastGetSolutionStepValue(WATER_PRESSURE) = 2.0E4;
}

void SetSolutionStepValuesForGeneralCheck(const intrusive_ptr<UPwSmallStrainElement<2, 3>>& rElement)
{
    const auto zero_values = array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto gravity     = array_1d<double, 3>{0.0, -10.0, 0.0};

    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VELOCITY)            = zero_values;
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = 0.0;
    }
    rElement->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{-0.015, 0.0, 0.0};
    rElement->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.015, 0.00, 0.0};
    rElement->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.015, 0.0};
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
    SetSolutionStepValuesForGeneralCheck(element);

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
    SetSolutionStepValuesForFluidFluxCheck(element);

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
    SetSolutionStepValuesForGeneralCheck(element);

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
    SetSolutionStepValuesForGeneralCheck(element);

    // Act, no exceptions on correct input
    element->Initialize(process_info);
    element->InitializeSolutionStep(process_info);

    // Assert
    for (const auto& r_node : element->GetGeometry()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementInitializeNonLinearIterationAndCalculateOnIntegrationPointsVectors,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 1.0;
    process_info[VELOCITY_COEFFICIENT]    = 1.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    element->Initialize(process_info);

    // Act and Assert
    element->InitializeNonLinearIteration(process_info);

    std::vector<Vector> calculated_values_at_integration_points;
    element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    Vector expected_values_at_integration_point(4);
    expected_values_at_integration_point <<= 300000, 150000, 0, -75000;
    ;
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points[i],
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 310000, 160000, 10000, -75000;
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points[i],
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.03, 0.015, 0, -0.015;
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points[i],
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points[i],
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    // getting a value from properties
    Vector initial_strain_vectop(4);
    initial_strain_vectop <<= 10000, -10000, 5000, -5000;
    ;
    element->GetProperties().SetValue(INITIAL_STRAIN_VECTOR, initial_strain_vectop);
    element->CalculateOnIntegrationPoints(INITIAL_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points[i],
                                           initial_strain_vectop, Defaults::relative_tolerance);
    }

    // getting values from a constitutive law
    calculated_values_at_integration_points.clear();
    element->CalculateOnIntegrationPoints(KIRCHHOFF_STRESS_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (auto i = std::size_t{0}; i < calculated_values_at_integration_points.size(); ++i) {
        KRATOS_EXPECT_EQ(calculated_values_at_integration_points[i].empty(), true);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementCalculateOnIntegrationPointsVariables, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 1.0;
    process_info[VELOCITY_COEFFICIENT]    = 1.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    element->Initialize(process_info);

    // Act and Assert
    element->InitializeNonLinearIteration(process_info);

    std::vector<double> calculated_values_at_integration_points;
    element->CalculateOnIntegrationPoints(VON_MISES_STRESS, calculated_values_at_integration_points, process_info);
    Vector expected_values_at_integration_point(3);
    expected_values_at_integration_point <<= 290474, 290474, 290474;
    ;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(MEAN_EFFECTIVE_STRESS,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 150000, 150000, 150000;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(MEAN_STRESS, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 160000, 160000, 160000;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(ENGINEERING_VON_MISES_STRAIN,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.0244949, 0.0244949, 0.0244949;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(ENGINEERING_VOLUMETRIC_STRAIN,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.045, 0.045, 0.045;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_VON_MISES_STRAIN,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.0244949, 0.0244949, 0.0244949;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_VOLUMETRIC_STRAIN,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.045, 0.045, 0.045;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(DEGREE_OF_SATURATION,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 1, 1, 1;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(HYDRAULIC_HEAD, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= -0.833333, -0.833333, -0.3333333;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(CONFINED_STIFFNESS, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 1e+07, 1e+07, 1e+07;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    element->CalculateOnIntegrationPoints(SHEAR_STIFFNESS, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 5e+06, 5e+06, 5e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    // get value from properties
    element->CalculateOnIntegrationPoints(BULK_MODULUS_FLUID, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 200, 200, 200;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);

    // get value from constitutive law
    calculated_values_at_integration_points.clear();
    element->CalculateOnIntegrationPoints(STRAIN_ENERGY, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementFinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 0.0;
    process_info[VELOCITY_COEFFICIENT]    = 0.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForFluidFluxCheck(element);
    element->Initialize(process_info);
    element->InitializeNonLinearIteration(process_info);

    // Act
    element->FinalizeSolutionStep(process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[0].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[1].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), -4.542);
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[2].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 4.542);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElementSetValuesOnIntegrationPointsMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto process_info = ProcessInfo{};
    // No storage, no dynamics, only statics and steady state
    process_info[DT_PRESSURE_COEFFICIENT] = 0.0;
    process_info[VELOCITY_COEFFICIENT]    = 0.0;

    Model model;
    auto  element = UPwSmallStrainElementWithUPwDofs(model, SetProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    element->Initialize(process_info);

    // Act
    Vector cauchy_stress_vector(4);
    cauchy_stress_vector <<= 1000.0, 2000.0, 3000.0, 4000.0;
    std::vector<Vector> cauchy_stress_vectors;
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);
    cauchy_stress_vector += Vector(4, 1000.0);
    cauchy_stress_vectors.push_back(cauchy_stress_vector);

    element->SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_stress_vectors, process_info);

    // Assert
    std::vector<Vector> calculated_cauchy_stress_vectors;
    element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, calculated_cauchy_stress_vectors, process_info);

    KRATOS_EXPECT_EQ(calculated_cauchy_stress_vectors.size(), cauchy_stress_vectors.size());
    for (auto i = std::size_t{0}; i < calculated_cauchy_stress_vectors.size(); ++i) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_cauchy_stress_vectors[i],
                                           cauchy_stress_vectors[i], Defaults::relative_tolerance);
    }

    std::vector<Matrix> calculated_tensors;
    element->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, calculated_tensors, process_info);
    Matrix expected_tensor(3, 3);
    expected_tensor <<= 1000, 4000, 0, 4000, 2000, 0, 0, 0, 3000;
    std::vector<Matrix> expected_tensors;
    expected_tensors.push_back(expected_tensor);
    expected_tensor <<= 2000, 5000, 0, 5000, 3000, 0, 0, 0, 4000;
    expected_tensors.push_back(expected_tensor);
    expected_tensor <<= 3000, 6000, 0, 6000, 4000, 0, 0, 0, 5000;
    expected_tensors.push_back(expected_tensor);
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_tensors[i], Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(TOTAL_STRESS_TENSOR, calculated_tensors, process_info);
    expected_tensors.clear();
    expected_tensor <<= 11000, 4000, 0, 4000, 12000, 0, 0, 0, 13000;
    expected_tensors.push_back(expected_tensor);
    expected_tensor <<= 12000, 5000, 0, 5000, 13000, 0, 0, 0, 14000;
    expected_tensors.push_back(expected_tensor);
    expected_tensor <<= 13000, 6000, 0, 6000, 14000, 0, 0, 0, 15000;
    expected_tensors.push_back(expected_tensor);
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_tensors[i], Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(ENGINEERING_STRAIN_TENSOR, calculated_tensors, process_info);
    expected_tensor <<= 0.03, -0.0075, 0, -0.0075, 0.015, 0, 0, 0, 0;
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_tensor, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, calculated_tensors, process_info);
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_tensor, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, calculated_tensors, process_info);
    Matrix expected_matrix(2, 2);
    expected_matrix <<= 9.084e-06, 0, 0, 9.084e-06;
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_matrix, Defaults::relative_tolerance);
    }

    // getting from constitutive law
    element->CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX, calculated_tensors, process_info);
    expected_matrix <<= 0, 0, 0, 0;
    for (auto i = std::size_t{0}; i < calculated_tensors.size(); ++i) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_tensors[i], expected_matrix, Defaults::relative_tolerance);
    }
}

} // namespace Kratos::Testing