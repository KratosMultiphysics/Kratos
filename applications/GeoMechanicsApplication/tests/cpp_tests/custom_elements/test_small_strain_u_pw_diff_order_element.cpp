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

#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_constitutive_law.h"
#include "tests/cpp_tests/test_utilities.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

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

} // namespace Kratos::Testing
