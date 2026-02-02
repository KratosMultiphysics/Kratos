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
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, PrincipalStresses_HasZeroesAsValues_WhenDefaultConstructed)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::PrincipalStresses{}.Values(), Vector(3, 0.0), Defaults::absolute_tolerance);
}
} // namespace Kratos::Testing