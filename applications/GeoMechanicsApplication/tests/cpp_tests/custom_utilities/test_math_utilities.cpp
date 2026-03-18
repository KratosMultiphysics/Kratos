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
#include "custom_utilities/ublas_utilities.h"
#include "includes/checks.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

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
    const auto stress_tensor =
        UblasUtilities::CreateMatrix({{10.0, 40.0, 0.0}, {40.0, 50.0, 0.0}, {0.0, 0.0, 20.0}});

    const auto angle           = MathUtils<>::DegreesToRadians(30.0);
    const auto rotation_matrix = UblasUtilities::CreateMatrix(
        {{std::cos(angle), -std::sin(angle), 0.0}, {std::sin(angle), std::cos(angle), 0.0}, {0.0, 0.0, 1.0}});
    const auto result = GeoMechanicsMathUtilities::RotateSecondOrderTensor(stress_tensor, rotation_matrix);

    const auto expected_result =
        UblasUtilities::CreateMatrix({{-14.641016151377542, 2.6794919243112303, 0.0},
                                      {2.6794919243112303, 74.64101615137753, 0.0},
                                      {0.0, 0.0, 20.0}});

    KRATOS_EXPECT_MATRIX_NEAR(result, expected_result, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(CheckVectorToDiagonalMatrix, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto vector = UblasUtilities::CreateVector({3.0, 4.0, 5.0, 6.0});

    // Act & assert
    const auto expected_result = UblasUtilities::CreateMatrix(
        {{3.0, 0.0, 0.0, 0.0}, {0.0, 4.0, 0.0, 0.0}, {0.0, 0.0, 5.0, 0.0}, {0.0, 0.0, 0.0, 6.0}});
    KRATOS_EXPECT_MATRIX_EQ(GeoMechanicsMathUtilities::VectorToDiagonalMatrix(vector), expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(CheckDiagonalMatrixToVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto matrix = UblasUtilities::CreateMatrix({{10.0, 40.0, 0.0}, {40.0, 50.0, 0.0}, {0.0, 0.0, 20.0}});

    // Act & assert
    auto expected_result = UblasUtilities::CreateVector({10.0, 50.0, 20.0});
    KRATOS_EXPECT_VECTOR_EQ(GeoMechanicsMathUtilities::DiagonalOfMatrixToVector(matrix), expected_result);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = GeoMechanicsMathUtilities::DiagonalOfMatrixToVector(Matrix(3, 2)), "Error: Attempting to convert diagonal of matrix to vector, but matrix has fewer columns than rows");
}

KRATOS_TEST_CASE_IN_SUITE(CheckRootsOfSecondOrderEquation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto A = 1.0;
    auto B = -0.9;
    auto C = -3.22;

    // Act & assert
    auto expected_result = UblasUtilities::CreateVector({-1.4, 2.3});
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::RootsOfSecondOrderEquation(A, B, C), expected_result, 1.0e-12);

    // Arrange
    A = 1.0;
    B = -2.4;
    C = 1.44;

    // Act & assert
    expected_result = UblasUtilities::CreateVector({1.2});
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::RootsOfSecondOrderEquation(A, B, C), expected_result, 1.0e-12);

    // Arrange
    A = 1.0;
    B = 2.0;
    C = 3.0;

    // Act & assert
    expected_result = UblasUtilities::CreateVector({});
    KRATOS_EXPECT_VECTOR_NEAR(GeoMechanicsMathUtilities::RootsOfSecondOrderEquation(A, B, C), expected_result, 1.0e-12);

    // Arrange
    A = 0.0;

    // Act & assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        [[maybe_unused]] const auto unused = GeoMechanicsMathUtilities::RootsOfSecondOrderEquation(A, B, C), "A == 0 does not define a second order equation.");
}

} // namespace Kratos::Testing