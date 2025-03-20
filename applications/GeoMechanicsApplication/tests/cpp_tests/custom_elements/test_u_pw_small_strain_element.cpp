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
//                   Gennady Markelov
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_elements/plane_strain_stress_state.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
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

std::shared_ptr<Properties> CreateProperties()
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
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 1.0E4;
        r_node.FastGetSolutionStepValue(DT_WATER_PRESSURE)   = 0.0;
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
intrusive_ptr<UPwSmallStrainElement<TDim, TNumNodes>> CreateUPwSmallStrainElementWithUPwDofs(
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

auto CreateUPwSmallStrainElementWithUPwDofs(Model& rModel, const Properties::Pointer& rProperties)
{
    auto& r_model_part = CreateModelPartWithUPwSolutionStepVariables(rModel);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0));
    const auto p_geometry = std::make_shared<Triangle2D3<Node>>(nodes);
    return CreateUPwSmallStrainElementWithUPwDofs<2, 3>(rProperties, p_geometry);
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
    Model model;
    const auto p_element = CreateUPwSmallStrainElementWithUPwDofs(model, std::make_shared<Properties>());

    // Act
    Element::DofsVectorType degrees_of_freedom;
    p_element->GetDofList(degrees_of_freedom, ProcessInfo{});

    // Assert
    KRATOS_EXPECT_EQ(degrees_of_freedom.size(), 9);
    std::vector<std::string> variable_names;
    std::transform(degrees_of_freedom.cbegin(), degrees_of_freedom.cend(), std::back_inserter(variable_names),
                   [](const auto& rpDof) { return rpDof->GetVariable().Name(); });
    const std::vector<std::string> expected_variable_names{
        "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_X",
        "DISPLACEMENT_Y", "WATER_PRESSURE", "WATER_PRESSURE", "WATER_PRESSURE"};

    KRATOS_EXPECT_EQ(variable_names, expected_variable_names);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_IntegrationMethod, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_geometry =
        std::make_shared<Triangle2D3<Node>>(ModelSetupUtilities::Create2D3NTriangleGeometry());
    const auto                        p_properties = std::make_shared<Properties>();
    const UPwSmallStrainElement<2, 3> element(0, p_geometry, p_properties,
                                              std::make_unique<PlaneStrainStressState>());

    // Act and Assert
    KRATOS_EXPECT_EQ(element.GetIntegrationMethod(), GeometryData::IntegrationMethod::GI_GAUSS_2);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CheckDoesNotThrowOnCorrectInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    SetSolutionStepValuesForGeneralCheck(element);

    // Act, no exceptions on correct input
    KRATOS_EXPECT_EQ(element->Check(ProcessInfo{}), 0);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CalculatesSteadyStateRightHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    SetSolutionStepValuesForFluidFluxCheck(element);
    const auto process_info = ProcessInfo{};

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
    for (const auto& calculated_fluid_flux : calculated_fluid_flux_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(
            calculated_fluid_flux, (array_1d<double, 3>{0., 9.084, 0.}), Defaults::relative_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CalculatesSteadyStateLeftHandSide, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    SetSolutionStepValuesForGeneralCheck(element);
    const auto process_info = ProcessInfo{};

    // Act
    element->Initialize(process_info);
    Matrix actual_left_hand_side;
    element->CalculateLeftHandSide(actual_left_hand_side, process_info);

    // Assert
    Matrix expected_left_hand_side(9, 9);
    // clang-format off
    expected_left_hand_side <<= 5000000,       0, -5000000,        0,        0,         0,          0,          0,          0,
                                0,       2500000,  2500000, -2500000, -2500000,         0,          0,          0,          0,
                               -5000000, 2500000,  7500000, -2500000, -2500000,         0,          0,          0,          0,
                                0,      -2500000, -2500000,  7500000,  2500000, -5000000,           0,          0,          0,
                                0,      -2500000, -2500000,  2500000,  2500000,         0,          0,          0,          0,
                                0,             0,        0, -5000000,        0,  5000000,           0,          0,          0,
                                0,             0,        0,        0,        0,         0, -0.0004542,  0.0004542,          0,
                                0,             0,        0,        0,        0,         0,  0.0004542, -0.0009084,  0.0004542,
                                0,             0,        0,        0,        0,         0,          0,  0.0004542, -0.0004542;
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(actual_left_hand_side, expected_left_hand_side, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_InitializeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    SetSolutionStepValuesForGeneralCheck(element);
    const auto process_info = ProcessInfo{};

    // Act, no exceptions on correct input
    element->Initialize(process_info);
    EXPECT_NO_THROW(element->InitializeSolutionStep(process_info));

    // Assert
    for (const auto& r_node : element->GetGeometry()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_InitializeNonLinearIterationAndCalculateOnIntegrationPointsVectors,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    const auto process_info = ProcessInfo{};
    element->Initialize(process_info);

    // Act and Assert
    element->InitializeNonLinearIteration(process_info);

    std::vector<Vector> calculated_values_at_integration_points;
    element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    Vector expected_values_at_integration_point(4);
    expected_values_at_integration_point <<= 300000, 150000, 0, -75000;
    for (const auto& calculated_cauchy_stress_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_cauchy_stress_vector, expected_values_at_integration_point,
                                           Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(TOTAL_STRESS_VECTOR, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 310000, 160000, 10000, -75000;
    for (const auto& calculated_total_stress_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_total_stress_vector, expected_values_at_integration_point,
                                           Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(ENGINEERING_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0.03, 0.015, 0, -0.015;
    for (const auto& calculated_engineering_strain_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_engineering_strain_vector,
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (const auto& calculated_green_lagrange_strain_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_green_lagrange_strain_vector,
                                           expected_values_at_integration_point, Defaults::relative_tolerance);
    }

    // getting a value from properties
    Vector initial_strain_vector(4);
    initial_strain_vector <<= 10000, -10000, 5000, -5000;
    element->GetProperties().SetValue(INITIAL_STRAIN_VECTOR, initial_strain_vector);
    element->CalculateOnIntegrationPoints(INITIAL_STRAIN_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (const auto& calculated_initial_strain_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_initial_strain_vector, initial_strain_vector,
                                           Defaults::relative_tolerance);
    }

    // getting values from a constitutive law: does nothing right now
    calculated_values_at_integration_points.clear();
    element->CalculateOnIntegrationPoints(KIRCHHOFF_STRESS_VECTOR,
                                          calculated_values_at_integration_points, process_info);
    for (const auto& kirchhoff_stress_vector : calculated_values_at_integration_points) {
        KRATOS_EXPECT_TRUE(kirchhoff_stress_vector.empty())
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CalculateOnIntegrationPointsVariables, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    const auto process_info = ProcessInfo{};
    element->Initialize(process_info);

    // Act and Assert
    element->InitializeNonLinearIteration(process_info);

    std::vector<double> calculated_values_at_integration_points;
    element->CalculateOnIntegrationPoints(VON_MISES_STRESS, calculated_values_at_integration_points, process_info);
    Vector expected_values_at_integration_point(3);
    expected_values_at_integration_point <<= 290474, 290474, 290474;
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

    // get value from constitutive law: currently returns zero values
    calculated_values_at_integration_points.clear();
    element->CalculateOnIntegrationPoints(STRAIN_ENERGY, calculated_values_at_integration_points, process_info);
    expected_values_at_integration_point <<= 0, 0, 0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(calculated_values_at_integration_points,
                                       expected_values_at_integration_point, Defaults::relative_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_FinalizeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForFluidFluxCheck(element);
    const auto process_info = ProcessInfo{};
    element->Initialize(process_info);
    element->FinalizeNonLinearIteration(process_info);

    // Act
    element->FinalizeSolutionStep(process_info);

    // Assert
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[0].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[1].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), -4.542);
    KRATOS_EXPECT_DOUBLE_EQ(element->GetGeometry()[2].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 4.542);
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_SetValuesOnIntegrationPointsMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto  element = CreateUPwSmallStrainElementWithUPwDofs(model, CreateProperties());
    element->GetProperties().SetValue(BIOT_COEFFICIENT, 1.000000e+00);
    SetSolutionStepValuesForGeneralCheck(element);
    const auto process_info = ProcessInfo{};
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
    for (const auto& calculated_engineering_strain_tensor : calculated_tensors) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_engineering_strain_tensor, expected_tensor,
                                           Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, calculated_tensors, process_info);
    for (const auto& calculated_green_lagrange_strain_tensor : calculated_tensors) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_green_lagrange_strain_tensor, expected_tensor,
                                           Defaults::relative_tolerance);
    }

    element->CalculateOnIntegrationPoints(PERMEABILITY_MATRIX, calculated_tensors, process_info);
    Matrix expected_matrix(2, 2);
    expected_matrix <<= 9.084e-06, 0, 0, 9.084e-06;
    for (const auto& calculated_permeability_matrix : calculated_tensors) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_permeability_matrix, expected_matrix,
                                           Defaults::relative_tolerance);
    }

    // returns zero matrix for non-implemented outputs
    element->CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX, calculated_tensors, process_info);
    expected_matrix <<= 0, 0, 0, 0;
    for (const auto& calculated_constitutive_matrix : calculated_tensors) {
        KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(calculated_constitutive_matrix, expected_matrix,
                                           Defaults::relative_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(UPwSmallStrainElement_CalculateShearCapacityUsingUMatParameters,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto element_id   = std::size_t{1};
    auto           p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    auto umat_parameters = Vector{2};
    umat_parameters <<= 2.0, 0.0;
    p_properties->SetValue(UMAT_PARAMETERS, umat_parameters);
    p_properties->SetValue(INDEX_OF_UMAT_C_PARAMETER, 1);
    p_properties->SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);

    auto element = UPwSmallStrainElement<2, 3>{
        element_id, std::make_shared<Triangle2D3<Node>>(ModelSetupUtilities::Create2D3NTriangleGeometry()),
        p_properties, std::make_unique<PlaneStrainStressState>()};
    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    auto stress_vector = Vector{4};
    stress_vector <<= -2.0, 0.0, 2.0, 0.0;
    element.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR,
                                         std::vector<Vector>{3, stress_vector}, dummy_process_info);

    // Act
    auto actual_shear_capacity_values = std::vector<double>{};
    element.CalculateOnIntegrationPoints(GEO_SHEAR_CAPACITY, actual_shear_capacity_values, dummy_process_info);

    // Assert
    auto expected_shear_capacity_values = Vector{3};
    expected_shear_capacity_values <<= 1.0, 1.0, 1.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_shear_capacity_values, expected_shear_capacity_values,
                              Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing