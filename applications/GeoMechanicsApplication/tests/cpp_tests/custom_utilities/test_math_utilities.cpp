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
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"
#include <boost/numeric/ublas/assignment.hpp>

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
    const Vector vector = ScalarVector(3, 2.0);
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::Normalized(vector),
                              Vector{ScalarVector(3, 1 / std::sqrt(3))}, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_ReturnsNormalizedVector_ForAllNegativeVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector vector = ScalarVector(3, -2.0);
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::Normalized(vector),
                              Vector{ScalarVector(3, -1 / std::sqrt(3))}, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_Throws_WhenInputtingZeroVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const Vector vector = ZeroVector(3);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto normalized = GeoMechanicsMathUtilities::Normalized(vector),
        "A zero vector cannot be normalized")
}

KRATOS_TEST_CASE_IN_SUITE(CheckRotateTensor, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Matrix stress_tensor = ZeroMatrix(3, 3);
    // clang-format off
    stress_tensor <<= 10.0, 40.0, 0.0,
                      40.0, 50.0, 0.0,
                      0.0, 0.0, 20.0;
    // clang-format on

    const auto angle           = MathUtils<>::DegreesToRadians(30.0);
    Matrix     rotation_matrix = ZeroMatrix(3, 3);
    rotation_matrix <<= std::cos(angle), -std::sin(angle), 0.0, std::sin(angle), std::cos(angle),
        0.0, 0.0, 0.0, 1.0;
    const auto result = GeoMechanicsMathUtilities::RotateSecondOrderTensor(stress_tensor, rotation_matrix);

    Matrix expected_result = ZeroMatrix(3, 3);
    // clang-format off
    expected_result <<= -14.641016151377542, 2.6794919243112303, 0.0,
                         2.6794919243112303, 74.64101615137753, 0.0,
                         0.0, 0.0, 20.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing