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
//                   Wijtze Pieter Kikstra
//

#include "custom_constitutive/p_q_theta.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

#include <vector>

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_HasZeroesAsValuesWhenDefaultConstructed)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PQTheta().Values(), (std::vector{0.0, 0.0, 0.0}), Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::PQTheta{}.P(), 0.0, Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::PQTheta{}.Q(), 0.0, Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::PQTheta{}.Theta(), 0.0, Defaults::absolute_tolerance);
}

template <typename T>
class TestPQThetaFixture : public ::testing::Test
{
};

using VectorTypesForTypedPQThetaTest =
    ::testing::Types<Vector, BoundedVector<double, 3>, std::vector<double>>;
TYPED_TEST_SUITE(TestPQThetaFixture, VectorTypesForTypedPQThetaTest);

TYPED_TEST(TestPQThetaFixture, PQTheta_CanBeConstructedFromAnyVectorWithSizeOf3)
{
    // Arrange
    auto initialization_vector = TypeParam(3);
    initialization_vector[0]   = 0.1;
    initialization_vector[1]   = 0.2;
    initialization_vector[2]   = 0.3;

    // Act
    const auto stress_state = Geo::PQTheta{initialization_vector};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), initialization_vector, Defaults::absolute_tolerance);
    EXPECT_NEAR(stress_state.P(), initialization_vector[0], Defaults::absolute_tolerance);
    EXPECT_NEAR(stress_state.Q(), initialization_vector[1], Defaults::absolute_tolerance);
    EXPECT_NEAR(stress_state.Theta(), initialization_vector[2], Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       PQTheta_RaisesADebugErrorWhenAttemptingToConstructFromAVectorWithSizeUnequalTo3)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    // Arrange
    const auto too_short = UblasUtilities::CreateVector({1.0});
    const auto too_long  = UblasUtilities::CreateVector({2.0, 3.0, 4.0, 5.0});

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        Geo::PQTheta{too_short},
        "Cannot construct a PQTheta instance: expected 3 values, but got 1 value(s)");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        Geo::PQTheta{too_long},
        "Cannot construct a PQTheta instance: expected 3 values, but got 4 value(s)");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_RaisesADebugErrorWhenAttemptingToConstructFromAVectorWithWrongValues)
{
    // Arrange
    const auto vector_with_wrong_q     = UblasUtilities::CreateVector({1.0, -2.0, 0.3});
    const auto vector_with_wrong_theta = UblasUtilities::CreateVector({2.0, 3.0, -1.0});

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(Geo::PQTheta{vector_with_wrong_q},
                                      "Invalid PQTheta: Q = -2 is negative.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        Geo::PQTheta{vector_with_wrong_theta},
        "Invalid PQTheta: Theta = -1 is out of the range [-pi/6, pi/6].");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_CanBeConstructedFromAVectorExpression)
{
    // Arrange
    const auto some_matrix =
        UblasUtilities::CreateMatrix({{1.0, 2.0, 0.0}, {3.0, 4.0, 0.0}, {0.0, 0.0, 1.0}});
    const auto     some_vector = UblasUtilities::CreateVector({2.0, 3.0, 0.04});
    constexpr auto some_scalar = 1.5;

    // Act & Assert
    // The following code won't compile when the template constructor of class `PQTheta` (which receives
    // a vector-like object) uses `std::ranges` algorithms when a UBlas vector expression is given.
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PQTheta{some_scalar * prod(some_matrix, some_vector)}.Values(),
                              (std::vector{12.0, 27.0, 0.06}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_CanBeConstructedFromThreeValues)
{
    EXPECT_NEAR((Geo::PQTheta{0.1, 0.2, 0.3}.P()), 0.1, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::PQTheta{0.1, 0.2, 0.3}.Q()), 0.2, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::PQTheta{0.1, 0.2, 0.3}.Theta()), 0.3, Defaults::absolute_tolerance);
}

TYPED_TEST(TestPQThetaFixture, PQTheta_CanBeCopiedToAnyVectorTypeWithSizeOf3)
{
    // Arrange
    const auto stress_state = Geo::PQTheta{1.0, 2.0, 0.3};

    // Act & Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.CopyTo<TypeParam>(), (std::vector{1.0, 2.0, 0.3}),
                              Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_SupportsCompoundAssignment)
{
    // Arrange
    auto stress_state = Geo::PQTheta{1.0, 2.0, 0.1};

    // Act
    stress_state += Geo::PQTheta{3.0, 4.0, 0.2};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), (std::vector{4.0, 6.0, 0.3}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_SupportsAdditionOfTwoInstances)
{
    // Arrange & Act
    const auto summed_stress_state = Geo::PQTheta{1.0, 3.0, 0.05} + Geo::PQTheta{2.0, 4.0, 0.06};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(summed_stress_state.Values(), (std::vector{3.0, 7.0, 0.11}),
                              Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing