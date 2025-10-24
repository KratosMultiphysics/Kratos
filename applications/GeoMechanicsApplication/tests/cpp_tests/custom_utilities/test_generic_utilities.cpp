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
//                   Wijtze Pieter Kikstra
//

#include "custom_utilities/generic_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CheckVectorPermutation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto vector = Vector(5);
    vector <<= 1.0, 2.0, 3.0, 4.0, 5.0;
    const auto indices = std::vector<int>{4, 3, 2, 1, 0};

    // Act & assert
    auto expected_result = Vector(5);
    expected_result <<= 5.0, 4.0, 3.0, 2.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(GenericUtilities::PermutedVector(vector, indices), expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(CheckMatrixPermutation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto matrix = Matrix(4, 4);
    // clang-format off
    matrix <<= 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 1.0;
    // clang-format on
    const auto indices = std::vector<int>{3, 2, 1, 0};

    // Act & assert
    auto expected_result = Matrix(4, 4);
    // clang-format off
    expected_result <<= 0.0, 0.0, 0.0, 1.0,
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 1.0, 0.0, 0.0,
                        1.0, 0.0, 0.0, 0.0;
    // clang-format on
    KRATOS_EXPECT_MATRIX_EQ(GenericUtilities::MatrixWithPermutedColumns(matrix, indices), expected_result);
}

} // namespace Kratos::Testing