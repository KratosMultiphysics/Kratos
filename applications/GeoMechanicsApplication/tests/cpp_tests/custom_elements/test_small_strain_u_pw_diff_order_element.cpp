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
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_retention/saturated_law.h"
#include "custom_utilities/registration_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
#include "tests/cpp_tests/test_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <string>

namespace
{

using namespace Kratos;
using namespace std::string_literals;

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
    Matrix result(15, 15);
    result <<= 8531959.3787335716, -2715844.4022770403, 932795.698924731, -52972.802024035424,
        1814814.8148148144, -833333.33333333302, -3649940.2628434869, 268817.20430107455,
        -586618.87694146018, 154965.21189120878, -7043010.7526881732, 3178368.1214421256, 0, 0, 0,
        -2715844.4022770403, 7883152.3648886075, -886306.13535736862, 1687934.8513598982, 0,
        907407.40740740718, 3602150.537634409, -6663679.8088410972, 154965.21189120901,
        -388186.09881228721, -154965.21189120901, -3426628.7160025295, 0, 0, 0, 932795.698924731,
        -886306.13535736862, 2414874.5519713252, -790.63883617966258, 111111.11111111108,
        833333.33333333314, -3236559.1397849452, 3387096.7741935472, 46594.982078852947,
        -3380771.6635041102, -268817.20430107508, 47438.330170777976, 0, 0, 0, -52972.802024035424,
        1687934.8513598982, -790.63883617966258, 4824609.9515074827, 0, 55555.55555555554,
        53763.440860215087, -6456989.2473118249, -47438.330170777976, 118174.15138098202,
        47438.330170777976, -229285.26249209302, 0, 0, 0, 1814814.8148148144, 0, 111111.11111111108,
        0, 5629629.6296296297, 0, -296296.29629629571, 0, -592592.59259259282, 0, -6666666.666666666,
        0, 0, 0, 0, -833333.33333333302, 907407.40740740718, 833333.33333333314, 55555.55555555554,
        0, 2814814.8148148148, -1.1334696144634405e-10, -148148.14814814785, -3333333.3333333326,
        -296296.29629629641, 3333333.3333333326, -3333333.333333333, 0, 0, 0, -3649940.2628434869,
        3602150.537634409, -3236559.1397849452, 53763.440860215087, -296296.29629629641,
        -1.1334696144634405e-10, 19923536.43966547, -3655913.9784946213, -13385902.03106332,
        3225806.4516129009, 645161.29032257735, -3225806.4516129037, 0, 0, 0, 268817.20430107467,
        -6663679.8088410972, 3387096.7741935472, -6456989.2473118249, 0, -148148.1481481482,
        -3655913.9784946213, 19639187.57467144, 3225806.4516129023, -6692951.0155316582,
        -3225806.4516129023, 322580.64516128716, 0, 0, 0, -586618.87694145949, 154965.21189120889,
        46594.982078852947, -47438.330170777976, -592592.59259259235, -3333333.3333333326,
        -13385902.031063316, 3225806.4516129023, 19464755.077658299, -2846299.810246679,
        -4946236.5591397816, 2846299.8102466771, 0, 0, 0, 154965.21189120878, -388186.09881228651,
        -3380771.6635041102, 118174.15138098202, 0, -296296.29629629618, 3225806.4516129009,
        -6692951.0155316582, -2846299.8102466785, 18650783.610935409, 2846299.8102466785,
        -11391524.351676149, 0, 0, 0, -7043010.7526881713, -154965.21189120889, -268817.20430107508,
        47438.330170777976, -6666666.666666666, 3333333.3333333326, 645161.29032257712, -3225806.4516129023,
        -4946236.5591397816, 2846299.810246679, 18279569.892473117, -2846299.8102466771, 0, 0, 0,
        3178368.1214421256, -3426628.7160025286, 47438.330170777976, -229285.26249209302, 0,
        -3333333.333333333, -3225806.4516129037, 322580.64516128669, 2846299.8102466771,
        -11391524.351676149, -2846299.8102466771, 18058191.018342815, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, -0.00096896000000000005, 0.00048448000000000002, 0.00048448000000000002, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00048448000000000002, -0.00048448000000000002, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00048448000000000002, 0, -0.00048448000000000002;
    return result;
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

    auto stress_vector = Vector{4};
    stress_vector <<= -1.5, 0.0, 1.5, 0.0;
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
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(CreateProperties());

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
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(CreateProperties());

    SetSolutionStepValuesForGeneralCheck(p_element);

    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    auto serializer = StreamSerializer{};
    serializer.save("test_tag", p_element);

    // Act
    auto p_loaded_element = make_intrusive<SmallStrainUPwDiffOrderElement>();
    serializer.load("test_tag", p_loaded_element);
    auto actual_lhs_values = Matrix{};
    p_loaded_element->CalculateLeftHandSide(actual_lhs_values, dummy_process_info);

    KRATOS_EXPECT_MATRIX_NEAR(ExpectedLeftHandSide(), actual_lhs_values, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateRHS, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = CreateProperties();
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
    auto p_properties = CreateProperties();
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
    auto p_properties = CreateProperties();
    p_properties->SetValue(BIOT_COEFFICIENT, 1.0); // to get RHS contributions of coupling
    auto p_element = CreateSmallStrainUPwDiffOrderElementWithUPwDofs(p_properties);

    Vector output;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        p_element->Calculate(CAUCHY_STRAIN_VECTOR, output, ProcessInfo{}),
        "Variable CAUCHY_STRAIN_VECTOR is unknown for element with Id 1.");
}

} // namespace Kratos::Testing
