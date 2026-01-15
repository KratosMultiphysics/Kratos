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
    constexpr double value = 3.0;

    // Act
    auto       kappa = 0.0;
    const auto cf    = FunctionObjectUtilities::MakeConstantFunction(value);

    // Assert
    KRATOS_CHECK_NEAR(cf(kappa), value, Defaults::absolute_tolerance);

    // Act
    // kappa = 10.0;
    // cf = FunctionObjectUtilities::MakeConstantFunction(value);

    // Assert
    // KRATOS_CHECK_NEAR(cf(kappa), value, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FunctionObjectUtilities_MakeLinearFunction, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr double value = 3.0;
    const auto       array = UblasUtilities::CreateVector({30.0});

    // Act
    auto       kappa = 4.0;
    const auto cf    = FunctionObjectUtilities::MakeLinearFunction(value, array[0]);

    // Assert
    KRATOS_CHECK_NEAR(cf(kappa), 123.0, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
