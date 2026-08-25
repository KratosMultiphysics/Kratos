// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "geo_mechanics_application_variables.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionBoundaryTypeVariableIsRegistered, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Variable<std::string>>::Has("GEO_SEEPAGE_BOUNDARY_TYPE"))
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionBoundaryTypeVariableHasExpectedName, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_EQ(GEO_SEEPAGE_BOUNDARY_TYPE.Name(), "GEO_SEEPAGE_BOUNDARY_TYPE");
}

KRATOS_TEST_CASE_IN_SUITE(GeoSeepageConditionIsRegisteredForBothLineGeometries, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has("GeoSeepageCondition2D2N"))
    KRATOS_EXPECT_TRUE(KratosComponents<Condition>::Has("GeoSeepageCondition2D3N"))
}

} // namespace Kratos::Testing
