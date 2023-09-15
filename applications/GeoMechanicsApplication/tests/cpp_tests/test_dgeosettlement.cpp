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

void RunStage(KratosGeoSettlement &settlement) {
    settlement.RunStage("",
                        "",
                        [](const char *) {},
                        [](const double) {},
                        [](const char *) {},
                        []() { return true; });
}

void ExpectNumberOfReadCallsIsEqualToOne(const KratosGeoSettlement &settlement) {
    const auto inputUtilityFromSettlement = dynamic_cast<const InputUtilityStub *>(settlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(inputUtilityFromSettlement, nullptr);
    KRATOS_EXPECT_EQ(inputUtilityFromSettlement->numberOfReadCalls(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(CreatingKratosGeoSettlementDoesNotThrow, WorkInProgress) {
    auto inputUtility = std::make_unique<InputUtilityStub>();
    KratosGeoSettlement settlement(std::move(inputUtility));
    KRATOS_EXPECT_TRUE(true); // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(RunStageMakesRelevantCallsOnce, WorkInProgress) {
    auto inputUtility = std::make_unique<InputUtilityStub>(); // invalid after the constructor
    KratosGeoSettlement settlement(std::move(inputUtility));

    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    const auto inputUtilityFromSettlement = dynamic_cast<const InputUtilityStub *>(settlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(inputUtilityFromSettlement, nullptr);
    KRATOS_EXPECT_EQ(inputUtilityFromSettlement->numberOfMaterialCalls(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageDoesPerformMaterialCallWhenNotSpecified, WorkInProgress) {
    auto inputUtility = std::make_unique<InputUtilityStub>(); // invalid after the constructor
    inputUtility->mParameterJsonString = "{"
                                         "\"solver_settings\":"
                                         "{"
                                         "\"model_part_name\":\"test\","
                                         "\"model_import_settings\":"
                                         "{"
                                         "\"input_type\": \"mdpa\","
                                         "\"input_filename\": \"mesh_stage1\""
                                         "}"
                                         "}"
                                         "}"; // no material import settings

    KratosGeoSettlement settlement(std::move(inputUtility));

    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    const auto inputUtilityFromSettlement = dynamic_cast<const InputUtilityStub *>(settlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(inputUtilityFromSettlement, nullptr);
    KRATOS_EXPECT_EQ(inputUtilityFromSettlement->numberOfMaterialCalls(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageTwiceOnlyCallsReadFromModelOnce, WorkInProgress) {
    auto inputUtility = std::make_unique<InputUtilityStub>();
    KratosGeoSettlement settlement(std::move(inputUtility));

    RunStage(settlement);
    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
}

}