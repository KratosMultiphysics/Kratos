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

#pragma once

#include "includes/ublas_interface.h"

#include <gtest/gtest.h>

namespace Kratos::Testing::Defaults
{

constexpr auto absolute_tolerance = 1.0e-12;
constexpr auto relative_tolerance = 1.0e-6;

} // namespace Kratos::Testing::Defaults

namespace Kratos::Testing
{

template <typename PointContainerType1, typename PointContainerType2>
void ExpectPointsAreNear(const PointContainerType1& rPoints1,
                         const PointContainerType2& rPoints2,
                         double AbsoluteTolerance = Defaults::absolute_tolerance)
{
    ASSERT_EQ(rPoints1.size(), rPoints2.size());
    for (auto i = std::size_t{0}; i < rPoints1.size(); ++i) {
        EXPECT_NEAR(rPoints1[i].X(), rPoints2[i].X(), AbsoluteTolerance);
        EXPECT_NEAR(rPoints1[i].Y(), rPoints2[i].Y(), AbsoluteTolerance);
        EXPECT_NEAR(rPoints1[i].Z(), rPoints2[i].Z(), AbsoluteTolerance);
    }
}

void AssertLHSMatrixBlocksAreNear(const Matrix& rActualLHSMatrix,
                                  const Matrix& rExpectedUUBlockMatrix,
                                  const Matrix& rExpectedUPBlockMatrix,
                                  const Matrix& rExpectedPUBlockMatrix,
                                  const Matrix& rExpectedPPBlockMatrix,
                                  std::size_t   NumberOfUDofs,
                                  std::size_t   NumberOfPwDofs,
                                  double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertUUBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedUUBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertUPBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedUPBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertPUBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedPUBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertPPBlockMatrixIsNear(const Matrix& rActualLHSMatrix,
                               const Matrix& rExpectedPPBlockMatrix,
                               std::size_t   NumberOfUDofs,
                               std::size_t   NumberOfPwDofs,
                               double        AbsoluteTolerance = Defaults::absolute_tolerance);

void AssertRHSVectorBlocksAreNear(const Vector& rActualRHSVector,
                                  const Vector& rExpectedUBlockVector,
                                  const Vector& rExpectedPBlockVector,
                                  std::size_t   NumberOfUDofs,
                                  std::size_t   NumberOfPwDofs,
                                  double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertUBlockVectorIsNear(const Vector& rActualRHSVector,
                              const Vector& rExpectedUBlockVector,
                              std::size_t   NumberOfUDofs,
                              double        AbsoluteTolerance = Defaults::absolute_tolerance);
void AssertPBlockVectorIsNear(const Vector& rActualRHSVector,
                              const Vector& rExpectedPBlockVector,
                              std::size_t   NumberOfUDofs,
                              std::size_t   NumberOfPwDofs,
                              double        AbsoluteTolerance = Defaults::absolute_tolerance);

} // namespace Kratos::Testing
