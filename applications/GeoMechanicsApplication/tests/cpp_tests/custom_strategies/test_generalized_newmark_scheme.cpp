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

#include "custom_strategies/schemes/generalized_newmark_scheme.hpp"
#include "spaces/default_spaces.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;
using SparseSpaceType = TDefaultSparseSpace<double>;
using LocalSpaceType  = TDefaultDenseSpace<double>;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ForInvalidTheta_CheckNewmarkScheme_Throws, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    constexpr double invalid_theta = -2.0;
    using SchemeType               = GeneralizedNewmarkScheme<SparseSpaceType, LocalSpaceType>;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeType scheme({}, invalid_theta),
                                      "Theta must be larger than zero, but got -2")
}

} // namespace Kratos::Testing
