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

#include "testing/testing.h"
#include "custom_workflows/dgeosettlement.h"

using namespace Kratos;

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(TestCreatingKratosGeoSettlementDoesNotThrow, WorkInProgress) {
    KratosGeoSettlement settlement;

//    settlement.RunStage("",
//                        "",
//                        [](const char*){},
//                        [](const double){},
//                        [](const char*){},
//                        [](){return true;});

    KRATOS_EXPECT_TRUE(true); // No other way to check that the constructor does not throw
}
}