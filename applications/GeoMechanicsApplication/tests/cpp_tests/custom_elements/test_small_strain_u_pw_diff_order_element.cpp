// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/small_strain_U_Pw_diff_order_element.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_retention/saturated_law.h"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "test_setup_utilities/element_setup_utilities.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
#include "tests/cpp_tests/test_utilities.h"

#include <string>
#include <tuple>
#include <vector>

namespace
{

using namespace Kratos;
using namespace std::string_literals;

std::shared_ptr<Properties> CreatePropertiesForUPwDiffOrderElementTest()
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

void SetSolutionStepValuesForGeneralCheck(const Element::Pointer& rElement)
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

auto CreateSmallStrainUPwDiffOrderElementWithUPwDofs(const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 0.0, -1.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.0, -0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(5, 0.5, -0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(6, 0.5, 0.05, 0.0));

    auto result = Testing::ElementSetupUtilities::Create2D6NDiffOrderElement(nodes, rProperties);
    const auto solution_step_variables = Geo::ConstVariableDataRefs{
        std::cref(WATER_PRESSURE),     std::cref(DT_WATER_PRESSURE), std::cref(DISPLACEMENT),
        std::cref(VELOCITY),           std::cref(ACCELERATION),      std::cref(VOLUME_ACCELERATION),
        std::cref(HYDRAULIC_DISCHARGE)};
    const auto degrees_of_freedom =
        Geo::ConstVariableRefs{std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT_X),
                               std::cref(DISPLACEMENT_Y), std::cref(DISPLACEMENT_Z)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(result, solution_step_variables, degrees_of_freedom);

    for (auto& r_node : nodes) {
        r_node.SetBufferSize(2);
    }
    return result;
}

Matrix ExpectedLeftHandSide()
{
    return UblasUtilities::CreateMatrix(
        {{8531959.3787335716, -2715844.4022770403, 932795.698924731, -52972.802024035424,
          1814814.8148148144, -833333.33333333302, -3649940.2628434869, 268817.20430107455,
          -586618.87694146018, 154965.21189120878, -7043010.7526881732, 3178368.1214421256, 0, 0, 0},
         {-2715844.4022770403, 7883152.3648886075, -886306.13535736862, 1687934.8513598982, 0,
          907407.40740740718, 3602150.537634409, -6663679.8088410972, 154965.21189120901,
          -388186.09881228721, -154965.21189120901, -3426628.7160025295, 0, 0, 0},
         {932795.698924731, -886306.13535736862, 2414874.5519713252, -790.63883617966258,
          111111.11111111108, 833333.33333333314, -3236559.1397849452, 3387096.7741935472,
          46594.982078852947, -3380771.6635041102, -268817.20430107508, 47438.330170777976, 0, 0, 0},
         {-52972.802024035424, 1687934.8513598982, -790.63883617966258, 4824609.9515074827, 0,
          55555.55555555554, 53763.440860215087, -6456989.2473118249, -47438.330170777976,
          118174.15138098202, 47438.330170777976, -229285.26249209302, 0, 0, 0},
         {1814814.8148148144, 0, 111111.11111111108, 0, 5629629.6296296297, 0, -296296.29629629571,
          0, -592592.59259259282, 0, -6666666.666666666, 0, 0, 0, 0},
         {-833333.33333333302, 907407.40740740718, 833333.33333333314, 55555.55555555554, 0,
          2814814.8148148148, -1.1334696144634405e-10, -148148.14814814785, -3333333.3333333326,
          -296296.29629629641, 3333333.3333333326, -3333333.333333333, 0, 0, 0},
         {-3649940.2628434869, 3602150.537634409, -3236559.1397849452, 53763.440860215087,
          -296296.29629629641, -1.1334696144634405e-10, 19923536.43966547, -3655913.9784946213,
          -13385902.03106332, 3225806.4516129009, 645161.29032257735, -3225806.4516129037, 0, 0, 0},
         {268817.20430107467, -6663679.8088410972, 3387096.7741935472, -6456989.2473118249, 0,
          -148148.1481481482, -3655913.9784946213, 19639187.57467144, 3225806.4516129023,
          -6692951.0155316582, -3225806.4516129023, 322580.64516128716, 0, 0, 0},
         {-586618.87694145949, 154965.21189120889, 46594.982078852947, -47438.330170777976,
          -592592.59259259235, -3333333.3333333326, -13385902.031063316, 3225806.4516129023,
          19464755.077658299, -2846299.810246679, -4946236.5591397816, 2846299.8102466771, 0, 0, 0},
         {154965.21189120878, -388186.09881228651, -3380771.6635041102, 118174.15138098202, 0,
          -296296.29629629618, 3225806.4516129009, -6692951.0155316582, -2846299.8102466785,
          18650783.610935409, 2846299.8102466785, -11391524.351676149, 0, 0, 0},
         {-7043010.7526881713, -154965.21189120889, -268817.20430107508, 47438.330170777976,
          -6666666.666666666, 3333333.3333333326, 645161.29032257712, -3225806.4516129023,
          -4946236.5591397816, 2846299.810246679, 18279569.892473117, -2846299.8102466771, 0, 0, 0},
         {3178368.1214421256, -3426628.7160025286, 47438.330170777976, -229285.26249209302, 0,
          -3333333.333333333, -3225806.4516129037, 322580.64516128669, 2846299.8102466771,
          -11391524.351676149, -2846299.8102466771, 18058191.018342815, 0, 0, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.00096896000000000005, 0.00048448000000000002, 0.00048448000000000002},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00048448000000000002, -0.00048448000000000002, 0},
         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00048448000000000002, 0, -0.00048448000000000002}});
}

Vector RightHandSideRegressionValues()
{
    return UblasUtilities::CreateVector({127876, -41008.2, -35286.7, 15096.6, -13055.6, -67314.3, 62688.2, -65544.4,
                                         -23942.7, 136350, -118280, 9166.82, -116.855, 97.39, 154.882});
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateShearCapacity, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_element = ElementSetupUtilities::Create2D6NDiffOrderElement();

    auto& p_properties = p_element->GetProperties();
    p_properties.SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    p_properties.SetValue(GEO_COHESION, 2.0);
    p_properties.SetValue(GEO_FRICTION_ANGLE, 0.0);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    auto stress_vector = UblasUtilities::CreateVector({-1.5, 0.0, 1.5, 0.0});
    p_element->SetValuesOnIntegrationPoints(
        CAUCHY_STRESS_VECTOR, std::vector<Vector>{3, stress_vector}, dummy_process_info);

    // Act
    auto actual_shear_capacity_values = std::vector<double>{};
    p_element->CalculateOnIntegrationPoints(GEO_SHEAR_CAPACITY, actual_shear_capacity_values, dummy_process_info);

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(actual_shear_capacity_values, (Vector{ScalarVector{3, 0.75}}),
                              Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateLHS, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_element =
        CreateSmallStrainUPwDiffOrderElementWithUPwDofs(CreatePropertiesForUPwDiffOrderElementTest());

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    // Act
    auto actual_lhs_values = Matrix{};
    p_element->CalculateLeftHandSide(actual_lhs_values, dummy_process_info);

    KRATOS_EXPECT_MATRIX_NEAR(ExpectedLeftHandSide(), actual_lhs_values, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateLHS_WithSaveAndLoad, KratosGeoMechanicsFastSuite)
{
    ScopedSerializerRegistration registration(
        std::make_pair("SaturatedLaw"s, SaturatedLaw{}), std::make_pair("PlaneStrain"s, PlaneStrain{}),
        std::make_pair("PlaneStrainStressState"s, PlaneStrainStressState{}));

    // Arrange
    auto p_element =
        CreateSmallStrainUPwDiffOrderElementWithUPwDofs(CreatePropertiesForUPwDiffOrderElementTest());

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    auto serializer = StreamSerializer{};
    serializer.save("test_tag", p_element);

    // Act
    auto p_loaded_element = make_intrusive<SmallStrainUPwDiffOrderElement<2, 6>>();
    serializer.load("test_tag", p_loaded_element);
    auto actual_lhs_values = Matrix{};
    p_loaded_element->CalculateLeftHandSide(actual_lhs_values, dummy_process_info);

    KRATOS_EXPECT_MATRIX_NEAR(ExpectedLeftHandSide(), actual_lhs_values, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateRHS, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = CreatePropertiesForUPwDiffOrderElementTest();
    p_properties->SetValue(BIOT_COEFFICIENT, 1.0); // to get RHS contributions of coupling
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    auto& r_geometry = p_element->GetGeometry();
    for (int counter = 0; auto& node : r_geometry) {
        node.FastGetSolutionStepValue(WATER_PRESSURE)    = counter * 1.0e5;
        node.FastGetSolutionStepValue(DT_WATER_PRESSURE) = counter * 5.0e5;
        ++counter;
    }

    // Act
    auto actual_rhs_values = Vector{};
    p_element->CalculateRightHandSide(actual_rhs_values, dummy_process_info);

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(RightHandSideRegressionValues(), actual_rhs_values, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_RHSEqualsUnbalanceVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = CreatePropertiesForUPwDiffOrderElementTest();
    p_properties->SetValue(BIOT_COEFFICIENT, 1.0); // to get RHS contributions of coupling
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);

    auto& r_geometry = p_element->GetGeometry();
    for (int counter = 0; auto& node : r_geometry) {
        node.FastGetSolutionStepValue(WATER_PRESSURE)    = counter * 1.0e5;
        node.FastGetSolutionStepValue(DT_WATER_PRESSURE) = counter * 5.0e5;
        ++counter;
    }

    // Act
    auto actual_rhs_values = Vector{};
    p_element->CalculateRightHandSide(actual_rhs_values, dummy_process_info);

    auto internal_forces = Vector{};
    p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, dummy_process_info);

    auto external_forces = Vector{};
    p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, dummy_process_info);

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_rhs_values, (external_forces - internal_forces), 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateThrowsDebugErrorForUnknownVectorVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto p_properties = CreatePropertiesForUPwDiffOrderElementTest();
    p_properties->SetValue(BIOT_COEFFICIENT, 1.0); // to get RHS contributions of coupling
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    Vector output;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Calculate(CAUCHY_STRAIN_VECTOR, output, ProcessInfo{}),
        "Variable CAUCHY_STRAIN_VECTOR is unknown for element with Id 1.");
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CheckThrowsOnFaultyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(RETENTION_LAW, "SaturatedLaw");
    p_properties->SetValue(SATURATED_SATURATION, 1.000000e+00);

    // Zero domain size (all nodes coincident)
    PointerVector<Node> coincident_nodes;
    for (int i = 0; i < 6; ++i)
        coincident_nodes.push_back(make_intrusive<Node>(i + 1, 0.0, 0.0, 0.0));
    auto p_element = Testing::ElementSetupUtilities::Create2D6NDiffOrderElement(coincident_nodes, p_properties);
    const auto dummy_process_info = ProcessInfo{};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Element 1 has non-positive size 0");

    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(3, 0.0, 1.0, 0.0));
    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(5, 0.5, 0.5, 0.0));
    nodes.push_back(make_intrusive<Node>(6, 0.0, 0.5, 0.0));
    p_element = Testing::ElementSetupUtilities::Create2D6NDiffOrderElement(nodes, p_properties);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable DISPLACEMENT on nodes 1 2 3 4 5 6");

    auto solution_step_variables = Geo::ConstVariableDataRefs{std::cref(DISPLACEMENT)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VELOCITY on nodes 1 2 3 4 5 6");

    solution_step_variables = Geo::ConstVariableDataRefs{std::cref(DISPLACEMENT), std::cref(VELOCITY)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable ACCELERATION on nodes 1 2 3 4 5 6");

    solution_step_variables =
        Geo::ConstVariableDataRefs{std::cref(DISPLACEMENT), std::cref(VELOCITY), std::cref(ACCELERATION)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable WATER_PRESSURE on nodes 1 2 3 4 5 6");

    solution_step_variables = Geo::ConstVariableDataRefs{
        std::cref(DISPLACEMENT), std::cref(VELOCITY), std::cref(ACCELERATION), std::cref(WATER_PRESSURE)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable DT_WATER_PRESSURE on nodes 1 2 3 4 5 6");

    solution_step_variables =
        Geo::ConstVariableDataRefs{std::cref(DISPLACEMENT), std::cref(VELOCITY), std::cref(ACCELERATION),
                                   std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "Missing variable VOLUME_ACCELERATION on nodes 1 2 3 4 5 6");

    solution_step_variables = Geo::ConstVariableDataRefs{
        std::cref(DISPLACEMENT),   std::cref(VELOCITY),          std::cref(ACCELERATION),
        std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "Missing the DoF for the variable DISPLACEMENT_X on nodes 1 2 3 4 5 6");

    const auto degrees_of_freedom =
        Geo::ConstVariableRefs{std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y),
                               std::cref(DISPLACEMENT_Z), std::cref(WATER_PRESSURE)};
    Testing::ElementSetupUtilities::AddVariablesToEntity(p_element, solution_step_variables, degrees_of_freedom);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_SOLID does not exist in the material properties with Id 0 at element with Id 1.");

    p_properties->SetValue(DENSITY_SOLID, 2.650000e+03);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "DENSITY_WATER does not exist in the material properties with Id 0 at element with Id 1.");

    p_properties->SetValue(DENSITY_WATER, 1.000000e+03);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(dummy_process_info),
                                      "BULK_MODULUS_SOLID does not exist in the material "
                                      "properties with Id 0 at element with Id 1.");

    p_properties->SetValue(BULK_MODULUS_SOLID, 1.000000e+12);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "POROSITY does not exist in the material properties with Id 0 at element with Id 1.");

    p_properties->SetValue(POROSITY, 1.000000e-01);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "IGNORE_UNDRAINED does not exist in the parameter list with Id 0 at element with Id 1.");

    p_properties->SetValue(IGNORE_UNDRAINED, false);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XX does not exist in the parameter list with Id 0 at element with Id 1.");

    p_properties->SetValue(PERMEABILITY_XX, 9.084000e-06);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_YY does not exist in the parameter list with Id 0 at element with Id 1.");

    p_properties->SetValue(PERMEABILITY_YY, 9.084000e-06);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "PERMEABILITY_XY does not exist in the parameter list with Id 0 at element with Id 1.");

    p_properties->SetValue(PERMEABILITY_XY, 0.000000e+00);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Check(dummy_process_info),
        "CONSTITUTIVE_LAW does not exist in the parameter list with Id 0 at element with Id 1.");

    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    p_element->Initialize(dummy_process_info);
    KRATOS_EXPECT_EQ(p_element->Check(dummy_process_info), 0);
}

} // namespace Kratos::Testing

using namespace Kratos;

namespace
{

struct DiffOrderElementTestParam {
    std::string                                                             name;
    std::function<Element::Pointer(ModelPart&, const Properties::Pointer&)> create_element;
    std::vector<std::size_t>                                                intermediate_indices;
    std::vector<double> expected_intermediate_pressures;
};

// Helper functions for each element type
Element::Pointer Create2D6(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    nodes.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(5, 0.5, 0.5, 0.0));
    nodes.push_back(rModelPart.CreateNewNode(6, 0.0, 0.5, 0.0));
    return Testing::ElementSetupUtilities::Create2D6NDiffOrderElement(nodes, rProperties);
}

Element::Pointer Create2D8(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 8; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<2, 8>>(
        1, Geometry<Node>::Pointer(new Quadrilateral2D8<Node>(nodes)), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer Create2D9(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 9; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<2, 9>>(
        1, Geometry<Node>::Pointer(new Quadrilateral2D9<Node>(nodes)), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer Create2D10(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 10; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<2, 10>>(
        1, Geometry<Node>::Pointer(new Triangle2D10<Node>(nodes)), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer Create2D15(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 15; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<2, 15>>(
        1, Geometry<Node>::Pointer(new Triangle2D15<Node>(nodes)), rProperties,
        std::make_unique<PlaneStrainStressState>());
}

Element::Pointer Create3D10(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 10; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<3, 10>>(
        1, Geometry<Node>::Pointer(new Tetrahedra3D10<Node>(nodes)), rProperties,
        std::make_unique<ThreeDimensionalStressState>());
}

Element::Pointer Create3D20(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 20; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<3, 20>>(
        1, Geometry<Node>::Pointer(new Hexahedra3D20<Node>(nodes)), rProperties,
        std::make_unique<ThreeDimensionalStressState>());
}

Element::Pointer Create3D27(ModelPart& rModelPart, const Properties::Pointer& rProperties)
{
    PointerVector<Node> nodes;
    for (int i = 0; i < 27; ++i)
        nodes.push_back(rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0));
    return make_intrusive<SmallStrainUPwDiffOrderElement<3, 27>>(
        1, Geometry<Node>::Pointer(new Hexahedra3D27<Node>(nodes)), rProperties,
        std::make_unique<ThreeDimensionalStressState>());
}

// Parameter set for each instantiation (for brevity, only 2D6 and 2D8 are fully filled; others can be extended similarly)
const std::vector<DiffOrderElementTestParam> diff_order_element_params = {
    {"2D6", Create2D6, {3, 4, 5}, {5.0, 15.0, 10.0}},
    {"2D8", Create2D8, {4, 5, 6, 7}, {5.0, 15.0, 25.0, 15.0}},
    {"2D9", Create2D9, {4, 5, 6, 7, 8}, {5.0, 15.0, 25.0, 15.0, 15.0}},
    {"2D10", Create2D10, {3, 4, 5, 6, 7, 8, 9}, {25.5556, 28.8889, 35.5556, 38.8889, 48.8889, 42.2222, 50.0}},
    {"2D15",
     Create2D15,
     {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
     {23.5938, 38.75, 37.0312, 42.0312, 60.0, 55.4688, 59.2969, 83.125, 70.3906, 84.8438, 78.2031, 88.4375}},
    {"3D10", Create3D10, {4, 5, 6, 7, 8, 9}, {5.0, 15.0, 10.0, 15.0, 20.0, 25.0}},
    {"3D20",
     Create3D20,
     {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
     {5.0, 15.0, 25.0, 15.0, 20.0, 30.0, 40.0, 50.0, 45.0, 55.0, 65.0, 55.0}},
    {"3D27",
     Create3D27,
     {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
     {5.0, 15.0, 25.0, 15.0, 20.0, 30.0, 40.0, 50.0, 45.0, 55.0, 65.0, 35.0, 15.0, 25.0, 35.0, 45.0,
      35.0, 55.0, 35.0}},
};

} // namespace

class FinalizeSolutionStep : public ::testing::TestWithParam<DiffOrderElementTestParam>
{
};

TEST_P(FinalizeSolutionStep, ReturnsIntermediateNodePressures)
{
    const auto& param = GetParam();
    Model       model;
    ModelPart&  r_model_part = model.CreateModelPart("Test");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto  properties = CreatePropertiesForUPwDiffOrderElementTest();
    auto  p_element  = param.create_element(r_model_part, properties);
    auto& r_geometry = p_element->GetGeometry();
    // Set input pressures: 10 * index
    for (std::size_t i = 0; i < r_geometry.size(); ++i) {
        r_geometry[i].SetBufferSize(2); // or higher if your test needs more steps
        r_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE)    = static_cast<double>(10 * i);
        r_geometry[i].FastGetSolutionStepValue(WATER_PRESSURE, 1) = 0.0;
    }
    // Act: use public API
    p_element->Initialize(ProcessInfo{});
    p_element->FinalizeSolutionStep(ProcessInfo{});
    // Assert
    constexpr auto tolerance = 0.0001;
    for (std::size_t i = 0; i < param.intermediate_indices.size(); ++i) {
        EXPECT_NEAR(r_geometry[param.intermediate_indices[i]].FastGetSolutionStepValue(WATER_PRESSURE),
                    param.expected_intermediate_pressures[i], tolerance)
            << "Failed for " << param.name << " at node " << param.intermediate_indices[i];
    }
}

INSTANTIATE_TEST_SUITE_P(AllDiffOrderElementTypes,
                         FinalizeSolutionStep,
                         ::testing::ValuesIn(diff_order_element_params),
                         [](const ::testing::TestParamInfo<FinalizeSolutionStep::ParamType>& info) {
                             return info.param.name;
                         });
