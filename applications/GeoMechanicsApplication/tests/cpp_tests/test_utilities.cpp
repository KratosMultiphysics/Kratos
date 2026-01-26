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

#include "tests/cpp_tests/test_utilities.h"
#include "includes/expect.h"

namespace Kratos::Testing
{

void AssertLHSMatrixBlocksAreNear(const Matrix& rActualLHSMatrix,
                                  const Matrix& rExpectedUUBlockMatrix,
                                  const Matrix& rExpectedUPBlockMatrix,
                                  const Matrix& rExpectedPUBlockMatrix,
                                  const Matrix& rExpectedPPBlockMatrix,
                                  std::size_t   NumberOfUDofs,
                                  std::size_t   NumberOfPwDofs,
                                  double        AbsoluteTolerance)
{
    ASSERT_EQ(rActualLHSMatrix.size1(), NumberOfUDofs + NumberOfPwDofs);
    ASSERT_EQ(rActualLHSMatrix.size2(), NumberOfUDofs + NumberOfPwDofs);

    AssertUUBlockMatrixIsNear(rActualLHSMatrix, rExpectedUUBlockMatrix, NumberOfUDofs, AbsoluteTolerance);
    AssertUPBlockMatrixIsNear(rActualLHSMatrix, rExpectedUPBlockMatrix, NumberOfUDofs,
                              NumberOfPwDofs, AbsoluteTolerance);
    AssertPUBlockMatrixIsNear(rActualLHSMatrix, rExpectedPUBlockMatrix, NumberOfUDofs,
                              NumberOfPwDofs, AbsoluteTolerance);
    AssertPPBlockMatrixIsNear(rActualLHSMatrix, rExpectedPPBlockMatrix, NumberOfUDofs,
                              NumberOfPwDofs, AbsoluteTolerance);
}

void AssertUUBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedUUBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedUUBlockMatrix.size1(), NumberOfUDofs);
    ASSERT_EQ(rExpectedUUBlockMatrix.size2(), NumberOfUDofs);

    KRATOS_EXPECT_MATRIX_NEAR(subrange(rActualLHSMatrix, 0, 0 + NumberOfUDofs, 0, 0 + NumberOfUDofs),
                              rExpectedUUBlockMatrix, AbsoluteTolerance)
}

void AssertUPBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedUPBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedUPBlockMatrix.size1(), NumberOfUDofs);
    ASSERT_EQ(rExpectedUPBlockMatrix.size2(), NumberOfPwDofs);

    KRATOS_EXPECT_MATRIX_NEAR(
        subrange(rActualLHSMatrix, 0, 0 + NumberOfUDofs, NumberOfUDofs, NumberOfUDofs + NumberOfPwDofs),
        rExpectedUPBlockMatrix, AbsoluteTolerance)
}

void AssertPUBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedPUBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedPUBlockMatrix.size1(), NumberOfPwDofs);
    ASSERT_EQ(rExpectedPUBlockMatrix.size2(), NumberOfUDofs);

    KRATOS_EXPECT_MATRIX_NEAR(
        subrange(rActualLHSMatrix, NumberOfUDofs, NumberOfUDofs + NumberOfPwDofs, 0, 0 + NumberOfUDofs),
        rExpectedPUBlockMatrix, AbsoluteTolerance)
}

void AssertPPBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedPPBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedPPBlockMatrix.size1(), NumberOfPwDofs);
    ASSERT_EQ(rExpectedPPBlockMatrix.size2(), NumberOfPwDofs);

    KRATOS_EXPECT_MATRIX_NEAR(subrange(rActualLHSMatrix, NumberOfUDofs, NumberOfUDofs + NumberOfPwDofs,
                                       NumberOfUDofs, NumberOfUDofs + NumberOfPwDofs),
                              rExpectedPPBlockMatrix, AbsoluteTolerance)
}

void AssertRHSVectorBlocksAreNear(const Vector& rActualRHSVector,
                                  const Vector& rExpectedUBlockVector,
                                  const Vector& rExpectedPBlockVector,
                                  std::size_t   NumberOfUDofs,
                                  std::size_t   NumberOfPwDofs,
                                  double        AbsoluteTolerance)
{
    ASSERT_EQ(rActualRHSVector.size(), NumberOfUDofs + NumberOfPwDofs);

    AssertUBlockVectorIsNear(rActualRHSVector, rExpectedUBlockVector, NumberOfUDofs, AbsoluteTolerance);
    AssertPBlockVectorIsNear(rActualRHSVector, rExpectedPBlockVector, NumberOfUDofs, NumberOfPwDofs,
                             AbsoluteTolerance);
}

void AssertUBlockVectorIsNear(const Vector& rActualRHSVector,
                              const Vector& rExpectedUBlockVector,
                              std::size_t   NumberOfUDofs,
                              double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedUBlockVector.size(), NumberOfUDofs);

    KRATOS_EXPECT_VECTOR_NEAR(subrange(rActualRHSVector, 0, 0 + NumberOfUDofs), rExpectedUBlockVector, AbsoluteTolerance)
}

void AssertPBlockVectorIsNear(const Vector& rActualRHSVector,
                              const Vector& rExpectedPBlockVector,
                              std::size_t   NumberOfUDofs,
                              std::size_t   NumberOfPwDofs,
                              double        AbsoluteTolerance)
{
    ASSERT_EQ(rExpectedPBlockVector.size(), NumberOfPwDofs);

    KRATOS_EXPECT_VECTOR_NEAR(subrange(rActualRHSVector, NumberOfUDofs, NumberOfUDofs + NumberOfPwDofs),
                              rExpectedPBlockVector, AbsoluteTolerance)
}

} // namespace Kratos::Testing
