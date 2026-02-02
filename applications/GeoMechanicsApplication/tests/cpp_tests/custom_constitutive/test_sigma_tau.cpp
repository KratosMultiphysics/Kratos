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

#include "custom_constitutive/geo_sigma_tau.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SigmaTau_HasZeroesAsValuesWhenDefaultConstructed, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau().Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(Geo::SigmaTau{}.Values(), Vector(2, 0.0), Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(Geo::SigmaTau{}.Sigma(), 0.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(Geo::SigmaTau{}.Tau(), 0.0, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing