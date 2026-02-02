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

#include "custom_constitutive/principal_stresses.h"
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

    const auto  principal_stresses = Geo::PrincipalStresses{bounded_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), bounded_vector, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_CanBeConstructedWithAStdVectorWithSize3)
{
    std::vector std_vector{1.0, 2.0, 3.0};
    const auto  principal_stresses = Geo::PrincipalStresses{std_vector};

    KRATOS_EXPECT_VECTOR_NEAR(principal_stresses.Values(), std_vector, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing