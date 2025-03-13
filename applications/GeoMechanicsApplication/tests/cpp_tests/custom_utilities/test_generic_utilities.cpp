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
    Vector vector(5);
    vector <<= 1.0, 2.0, 3.0, 4.0, 5.0;
    std::vector<int> indices = {4, 3, 2, 1, 0};

    // Act
    Vector result = GenericUtilities::PermutedVector(vector, indices);

    // Assert
    Vector expected_result(5);
    expected_result <<= 5.0, 4.0, 3.0, 2.0, 1.0;
    KRATOS_EXPECT_VECTOR_EQ(result, expected_result);
}

KRATOS_TEST_CASE_IN_SUITE(CheckMatrixPermutation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Matrix matrix(4, 4);
    matrix <<= 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 1.0;
    std::vector<int> indices = {3, 2, 1, 0};

    // Act
    Matrix result = GenericUtilities::MatrixWithPermutedColumns(matrix, indices);

    // Assert
    Matrix expected_result(4, 4);
    expected_result <<=  0.0, 0.0, 0.0, 1.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 1.0, 0.0, 0.0,
                         1.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_MATRIX_EQ(result, expected_result);
}


} // namespace Kratos::Testing