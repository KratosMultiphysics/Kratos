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

#include "custom_constitutive/p_q.hpp"
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
class TestPQFixture : public ::testing::Test
{
};

using VectorTypesForTypedPQTest = ::testing::Types<Vector, BoundedVector<double, 3>, std::vector<double>>;
TYPED_TEST_SUITE(TestPQFixture, VectorTypesForTypedPQTest);

TYPED_TEST(TestPQFixture, PQTheta_CanBeConstructedFromAnyVectorWithSizeOf3)
{
    // Arrange
    auto initialization_vector = TypeParam(3);
    initialization_vector[0]   = 1.0;
    initialization_vector[1]   = 2.0;
    initialization_vector[1]   = 3.0;

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
    EXPECT_THROW(Geo::PQTheta{too_short}, Exception);
    EXPECT_THROW(Geo::PQTheta{too_long}, Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_CanBeConstructedFromAVectorExpression)
{
    // Arrange
    const auto some_matrix =
        UblasUtilities::CreateMatrix({{1.0, 2.0, 0.0}, {3.0, 4.0, 0.0}, {0.0, 0.0, 1.0}});
    const auto     some_vector = UblasUtilities::CreateVector({2.0, 3.0, 4.0});
    constexpr auto some_scalar = 1.5;

    // Act & Assert
    // The following code won't compile when the template constructor of class `PQTheta` (which receives
    // a vector-like object) uses `std::ranges` algorithms when a UBlas vector expression is given.
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PQTheta{some_scalar * prod(some_matrix, some_vector)}.Values(),
                              (std::vector{12.0, 27.0, 6.0}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_CanBeConstructedFromFromThreeValues)
{
    EXPECT_NEAR((Geo::PQTheta{1.0, 2.0, 3.0}.P()), 1.0, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::PQTheta{1.0, 2.0, 3.0}.Q()), 2.0, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::PQTheta{1.0, 2.0, 3.0}.Theta()), 3.0, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_ComponentsCanBeModifiedDirectly)
{
    // Arrange
    auto stress_state = Geo::PQTheta{1.0, 2.0, 3.0};

    // Act
    stress_state.P() = 3.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), (std::vector{3.0, 2.0, 3.0}), Defaults::absolute_tolerance);

    // Act
    stress_state.Q() = 4.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), (std::vector{3.0, 4.0, 3.0}), Defaults::absolute_tolerance);

    // Act
    stress_state.Theta() = 5.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), (std::vector{3.0, 4.0, 5.0}), Defaults::absolute_tolerance);
}

TYPED_TEST(TestPQFixture, PQTheta_CanBeCopiedToAnyVectorTypeWithSizeOf3)
{
    // Arrange
    const auto stress_state = Geo::PQTheta{1.0, 2.0, 3.0};

    // Act & Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.CopyTo<TypeParam>(), (std::vector{1.0, 2.0, 3.0}),
                              Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_SupportsCompoundAssignment)
{
    // Arrange
    auto stress_state = Geo::PQTheta{1.0, 2.0, 3.0};

    // Act
    stress_state += Geo::PQTheta{3.0, 4.0, 5.0};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_state.Values(), (std::vector{4.0, 6.0, 8.0}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQTheta_SupportsAdditionOfTwoInstances)
{
    // Arrange & Act
    const auto summed_stress_state = Geo::PQTheta{1.0, 3.0, 5.0} + Geo::PQTheta{2.0, 4.0, 6.0};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(summed_stress_state.Values(), (std::vector{3.0, 7.0, 11.0}),
                              Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing