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
#include "custom_utilities/input_utility_stub.h"

using namespace Kratos;

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(CreatingKratosGeoSettlementDoesNotThrow, WorkInProgress) {
    InputUtilityStub inputUtility;
    KratosGeoSettlement settlement(inputUtility);
    KRATOS_EXPECT_TRUE(true); // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(RunStageDoesNotThrow, WorkInProgress) {
    InputUtilityStub inputUtility;
    KratosGeoSettlement settlement(inputUtility);

    settlement.RunStage("",
                        "",
                        [](const char*){},
                        [](const double){},
                        [](const char*){},
                        [](){return true;});

    KRATOS_EXPECT_TRUE(true); // No other way to check that there is no throw
}


}