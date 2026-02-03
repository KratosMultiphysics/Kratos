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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedWithAVectorWithSize3)
{
    const auto ublas_vector       = UblasUtilities::CreateVector({1.0, 2.0, 3.0});
    const auto principal_stresses = Geo::PrincipalStresses{ublas_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), ublas_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedWithABoundedVectorWithSize3)
{
    BoundedVector<double, 3> bounded_vector;
    bounded_vector[0] = 1.0;
    bounded_vector[1] = 2.0;
    bounded_vector[2] = 3.0;

    const auto principal_stresses = Geo::PrincipalStresses{bounded_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), bounded_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedWithAStdVectorWithSize3)
{
    const std::vector std_vector{1.0, 2.0, 3.0};
    const auto        principal_stresses = Geo::PrincipalStresses{std_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), std_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedWithAInitializerListWithSize3)
{
    const std::initializer_list initializer_list{1.0, 2.0, 3.0};
    const auto                  principal_stresses = Geo::PrincipalStresses{initializer_list};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), initializer_list, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_ThrowsWhenSizeIsIncorrect)
{
#ifndef KRATOS_DEBUG
    GTEST_SKIP() << "This test requires a debug build";
#endif

    const auto too_short = {1.0, 2.0};
    const auto too_long  = {1.0, 2.0, 3.0, 4.0};
    EXPECT_THROW((Geo::PrincipalStresses(too_long)), Exception);
    EXPECT_THROW((Geo::PrincipalStresses(too_short)), Exception);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStressesCanBeCopiedToVector) {}

template <typename T>
class TestCopyTo : public ::testing::Test
{
};

using ListTypes = ::testing::Types<Vector, BoundedVector<double, 3>, std::vector<double>>;
TYPED_TEST_SUITE(TestCopyTo, ListTypes);

TYPED_TEST(TestCopyTo, PrincipalStressesCanBeCopiedToListType)
{
    Geo::PrincipalStresses principal_stresses{std::vector{1.0, 2.0, 3.0}};

    const auto copied_vector = principal_stresses.CopyTo<TypeParam>();

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), copied_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing