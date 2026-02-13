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

#include "custom_constitutive/principal_stresses.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_HasZeroesAsValues_WhenDefaultConstructed)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PrincipalStresses{}.Values(), Vector(3, 0.0), Defaults::absolute_tolerance);
}

template <typename T>
class TestPrincipalStressFixture : public ::testing::Test
{
};

using TestVectorTypesPrincipalStress =
    ::testing::Types<Vector, BoundedVector<double, 3>, std::vector<double>>;
TYPED_TEST_SUITE(TestPrincipalStressFixture, TestVectorTypesPrincipalStress);

TYPED_TEST(TestPrincipalStressFixture, PrincipalStresses_CanBeConstructedFromAnyVectorWithSizeOf3)
{
    TypeParam initialization_vector(3);
    initialization_vector[0] = 1.0;
    initialization_vector[1] = 2.0;
    initialization_vector[2] = 3.0;

    const auto principal_stresses = Geo::PrincipalStresses{initialization_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), initialization_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_ThrowsWhenSizeIsIncorrect)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    const auto too_short = {1.0, 2.0};
    const auto too_long  = {1.0, 2.0, 3.0, 4.0};
    EXPECT_THROW(Geo::PrincipalStresses{too_long}, Exception);
    EXPECT_THROW(Geo::PrincipalStresses{too_short}, Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedFromAStdInitializerListWithSize3)
{
    KRATOS_EXPECT_VECTOR_NEAR((Geo::PrincipalStresses{1.0, 2.0, 3.0}.Values()),
                              (std::vector{1.0, 2.0, 3.0}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       PrincipalStress_RaisesADebugErrorWhenAttemptingToConstructFromANonEmptyStdInitializerListWithSizeUnequalTo3)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    EXPECT_THROW((Geo::PrincipalStresses{1.0, 2.0}), Exception);
    EXPECT_THROW((Geo::PrincipalStresses{1.0, 2.0, 3.0, 4.0}), Exception);
}

TYPED_TEST(TestPrincipalStressFixture, PrincipalStressesCanBeCopiedToAnyVectorWithSize3)
{
    Geo::PrincipalStresses principal_stresses{std::vector{1.0, 2.0, 3.0}};

    const auto copied_vector = principal_stresses.CopyTo<TypeParam>();

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), copied_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeChangedDirectly)
{
    auto stresses = Geo::PrincipalStresses{1.0, 2.0, 3.0};

    stresses.Values()[0] = 4.0;

    KRATOS_EXPECT_VECTOR_NEAR(stresses.Values(), (std::vector{4.0, 2.0, 3.0}), Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_SupportsCompoundAssignment)
{
    // Arrange
    auto principal_stresses = Geo::PrincipalStresses{3.0, 2.0, 1.0};

    // Act
    principal_stresses += Geo::PrincipalStresses{5.0, 4.0, 3.0};

    // Assert
    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), (std::vector{8.0, 6.0, 4.0}), Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing