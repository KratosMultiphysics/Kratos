// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_utilities/math_utilities.h"
#include "geo_mechanics_fast_suite.h"
#include "includes/checks.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateDeterminants_ReturnsCorrectResults, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<Matrix> matrices = {ScalarMatrix(1, 1, 3.0), ZeroMatrix(3, 3), IdentityMatrix(3) * 2.0};

    const std::vector<double> results = GeoMechanicsMathUtilities::CalculateDeterminants(matrices);

    const std::vector<double> expected_results = {3.0, 0.0, 8.0};
    KRATOS_CHECK_VECTOR_NEAR(results, expected_results, 1.0e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDeterminants_ReturnsEmptyVectorForEmptyInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<Matrix> matrices = {};

    const std::vector<double> results = GeoMechanicsMathUtilities::CalculateDeterminants(matrices);

    KRATOS_EXPECT_TRUE(results.empty())
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_ReturnsNormalizedVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector vector = ScalarVector(3, 2);
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::Normalized(vector),
                              Vector{ScalarVector(3, 1 / std::sqrt(3))}, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_Throws_WhenInputtingZeroVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector vector = ZeroVector(3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto normalized = GeoMechanicsMathUtilities::Normalized(vector),
        "A zero vector cannot be normalized")
}

} // namespace Kratos::Testing