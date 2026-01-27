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

#include "custom_utilities/math_utilities.hpp"
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
    const auto vector = Vector{ScalarVector(3, 2.0)};
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::Normalized(vector),
                              Vector{ScalarVector(3, 1 / std::sqrt(3))}, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_ReturnsNormalizedVector_ForAllNegativeVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto vector = Vector{ScalarVector(3, -2.0)};
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::Normalized(vector),
                              Vector{ScalarVector(3, -1 / std::sqrt(3))}, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(Normalized_Throws_WhenInputtingZeroVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto vector = Vector{ZeroVector(3)};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto normalized = GeoMechanicsMathUtilities::Normalized(vector),
        "A zero vector cannot be normalized")
}

KRATOS_TEST_CASE_IN_SUITE(CheckRotateTensor, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Matrix stress_tensor = ZeroMatrix(3, 3);
    // clang-format off
    stress_tensor <<= 10.0, 40.0,  0.0,
                      40.0, 50.0,  0.0,
                       0.0,  0.0, 20.0;
    // clang-format on

    const auto angle           = MathUtils<>::DegreesToRadians(30.0);
    auto       rotation_matrix = Matrix(3, 3);
    // clang-format off
    rotation_matrix <<= std::cos(angle), -std::sin(angle), 0.0,
                        std::sin(angle),  std::cos(angle), 0.0,
                        0.0,              0.0,             1.0;
    // clang-format on
    const auto result = GeoMechanicsMathUtilities::RotateSecondOrderTensor(stress_tensor, rotation_matrix);

    auto expected_result = Matrix(3, 3);
    // clang-format off
    expected_result <<= -14.641016151377542,   2.6794919243112303,  0.0,
                          2.6794919243112303, 74.64101615137753,    0.0,
                          0.0,                 0.0,                20.0;
    // clang-format on

    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CheckVectorToDiagonalMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto vector = Vector(4);
    vector <<= 3.0, 4.0, 5.0, 6.0;

    // Act & assert
    auto expected_result = Matrix(4, 4);
    // clang-format off
    expected_result <<= 3.0, 0.0, 0.0, 0.0,
                        0.0, 4.0, 0.0, 0.0,
                        0.0, 0.0, 5.0, 0.0,
                        0.0, 0.0, 0.0, 6.0;
    // clang-format on
    KRATOS_EXPECT_MATRIX_EQ(GeoMechanicsMathUtilities::VectorToDiagonalMatrix(vector), expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(CheckDiagonalMatrixToVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto matrix = Matrix(3, 3);
    // clang-format off
    matrix <<= 10.0, 40.0,  0.0,
               40.0, 50.0,  0.0,
                0.0,  0.0, 20.0;
    // clang-format on

    // Act & assert
    auto expected_result = Vector(3);
    expected_result <<= 10.0, 50.0, 20.0;
    KRATOS_EXPECT_VECTOR_EQ(GeoMechanicsMathUtilities::DiagonalOfMatrixToVector(matrix), expected_result);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = GeoMechanicsMathUtilities::DiagonalOfMatrixToVector(Matrix(3, 2)), "Error: Attempting to convert diagonal of matrix to vector, but matrix has fewer columns than rows");
}
} // namespace Kratos::Testing