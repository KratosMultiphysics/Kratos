// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#include "custom_utilities/function_object_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_aliases.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FunctionObjectUtilities_MakeConstantFunction, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto value = 3.0;

    // Act
    auto kappa = 0.0;
    auto cf    = FunctionObjectUtilities::MakeConstantFunction(value);

    // Assert
    KRATOS_EXPECT_NEAR(cf(kappa), value, Defaults::absolute_tolerance);

    // Act
    kappa = 10.0;
    cf    = FunctionObjectUtilities::MakeConstantFunction(value);

    // Assert
    KRATOS_EXPECT_NEAR(cf(kappa), value, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FunctionObjectUtilities_MakeLinearFunction, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto value       = 3.0;
    const auto     coefficient = 30.0;
    const auto     cf          = FunctionObjectUtilities::MakeLinearFunction(value, coefficient);

    // Act & Assert
    auto kappa = 4.0;
    KRATOS_EXPECT_NEAR(cf(kappa), 123.0, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
