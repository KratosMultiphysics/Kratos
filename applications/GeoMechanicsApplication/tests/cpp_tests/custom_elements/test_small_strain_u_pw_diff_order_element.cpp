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

#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUPwDiffOrderElement_CalculateShearCapacity, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto element_id   = std::size_t{1};
    auto           p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CONSTITUTIVE_LAW, std::make_shared<StubConstitutiveLaw>());
    p_properties->SetValue(GEO_COHESION, 2.0);
    p_properties->SetValue(GEO_FRICTION_ANGLE, 0.0);

    auto element = SmallStrainUPwDiffOrderElement{
        element_id, std::make_shared<Triangle2D6<Node>>(ModelSetupUtilities::Create2D6NTriangleGeometry()),
        p_properties, std::make_unique<PlaneStrainStressState>()};
    const auto dummy_process_info = ProcessInfo{};
    element.Initialize(dummy_process_info);

    auto stress_vector = Vector{4};
    stress_vector <<= -1.5, 0.0, 1.5, 0.0;
    element.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR,
                                         std::vector<Vector>{3, stress_vector}, dummy_process_info);

    // Act
    auto actual_shear_capacity_values = std::vector<double>{};
    element.CalculateOnIntegrationPoints(GEO_SHEAR_CAPACITY, actual_shear_capacity_values, dummy_process_info);

    // Assert
    auto expected_shear_capacity_values = Vector{3};
    expected_shear_capacity_values <<= 0.75, 0.75, 0.75;
    KRATOS_EXPECT_VECTOR_NEAR(actual_shear_capacity_values, expected_shear_capacity_values,
                              Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
