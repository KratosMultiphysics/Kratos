// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#include "custom_constitutive/p_q.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQ_HasZeroesAsValuesWhenDefaultConstructed)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PQ().Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::PQ{}.P(), 0.0, Defaults::absolute_tolerance);
    EXPECT_NEAR(Geo::PQ{}.Q(), 0.0, Defaults::absolute_tolerance);
}

template <typename T>
class TestPQFixture : public ::testing::Test
{
};

using TestVectorTypes = ::testing::Types<Vector, BoundedVector<double, 2>, std::vector<double>>;
TYPED_TEST_SUITE(TestPQFixture, TestVectorTypes);

TYPED_TEST(TestPQFixture, PQ_CanBeConstructedFromAnyVectorWithSizeOf2)
{
    // Arrange
    auto initialization_vector = TypeParam(2);
    initialization_vector[0]   = 1.0;
    initialization_vector[1]   = 2.0;

    // Act
    const auto stress_invariant = Geo::PQ{initialization_vector};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_invariant.Values(), initialization_vector, Defaults::absolute_tolerance);
    EXPECT_NEAR(stress_invariant.P(), initialization_vector[0], Defaults::absolute_tolerance);
    EXPECT_NEAR(stress_invariant.Q(), initialization_vector[1], Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQ_RaisesADebugErrorWhenAttemptingToConstructFromAVectorWithSizeUnequalTo2)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    // Arrange
    const auto too_short = UblasUtilities::CreateVector({1.0});
    const auto too_long  = UblasUtilities::CreateVector({2.0, 3.0, 4.0});

    // Act & Assert
    EXPECT_THROW(Geo::PQ{too_short}, Exception);
    EXPECT_THROW(Geo::PQ{too_long}, Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQ_CanBeConstructedFromAStdInitializerListWithSize2)
{
    EXPECT_NEAR((Geo::PQ{1.0, 2.0}.P()), 1.0, Defaults::absolute_tolerance);
    EXPECT_NEAR((Geo::PQ{1.0, 2.0}.Q()), 2.0, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       PQ_RaisesADebugErrorWhenAttemptingToConstructFromAnStdInitializerListWithSizeUnequalTo2)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    EXPECT_THROW(Geo::PQ{1.0}, Exception);
    EXPECT_THROW((Geo::PQ{2.0, 3.0, 4.0}), Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PQ_ComponentsCanBeModifiedDirectly)
{
    // Arrange
    auto stress_invariant = Geo::PQ{1.0, 2.0};

    // Act
    stress_invariant.P() = 3.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_invariant.Values(), (std::vector{3.0, 2.0}), Defaults::absolute_tolerance);

    // Act
    stress_invariant.Q() = 4.0;

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_invariant.Values(), (std::vector{3.0, 4.0}), Defaults::absolute_tolerance);
}

TYPED_TEST(TestPQFixture, PQ_CanBeCopiedToAnyVectorTypeWithSizeOf2)
{
    // Arrange
    const auto stress_invariant = Geo::PQ{1.0, 2.0};

    // Act & Assert
    KRATOS_EXPECT_VECTOR_NEAR(stress_invariant.CopyTo<TypeParam>(), (std::vector{1.0, 2.0}),
                              Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing