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
#include "input_utility_stub.h"

using namespace Kratos;

namespace Kratos::Testing {

const std::string parameter_json_settings = "{"
                                            "\"solver_settings\":"
                                            "{"
                                            "\"model_part_name\":\"test\","
                                            "\"model_import_settings\":"
                                            "{"
                                            "\"input_type\": \"mdpa\","
                                            "\"input_filename\": \"mesh_stage1\""
                                            "},"
                                            "\"material_import_settings\": "
                                            "{"
                                            "\"materials_filename\": \"MaterialParameters1.json\""
                                            "}"
                                            "}"
                                            "}"; // these have material_import settings


void RunStage(KratosGeoSettlement& rSettlement) {
    rSettlement.RunStage("",
                         "",
                         [](const char *) {}, // kept empty as a stub method
                        [](const double) {}, // kept empty as a stub method
                        [](const char *) {}, // kept empty as a stub method
                        []() { return true; });
}

void ExpectNumberOfReadCallsIsEqualToOne(const KratosGeoSettlement &rSettlement) {
    const auto input_utility_from_settlement = dynamic_cast<const InputUtilityStub*>(rSettlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(input_utility_from_settlement, nullptr);
    KRATOS_EXPECT_EQ(input_utility_from_settlement->NumberOfReadCalls(), 1);
}

void ExpectNumberOfMaterialCallsEqualTo(const int expectedNumberOfMaterialCalls, const KratosGeoSettlement& rSettlement) {
    const auto input_utility_from_settlement = dynamic_cast<const InputUtilityStub*>(rSettlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(input_utility_from_settlement, nullptr);
    KRATOS_EXPECT_EQ(input_utility_from_settlement->NumberOfMaterialCalls(), expectedNumberOfMaterialCalls);
}

KRATOS_TEST_CASE_IN_SUITE(CreatingKratosGeoSettlementDoesNotThrow, KratosGeoMechanicsFastSuite) {
    auto input_utility = std::make_unique<InputUtilityStub>(parameter_json_settings);

    bool has_thrown = false;
    try
    {
        KratosGeoSettlement settlement(std::move(input_utility));
    }
    catch (...)
    {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown); // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(RunStageMakesRelevantCallsOnce, KratosGeoMechanicsFastSuite) {
    auto input_utility = std::make_unique<InputUtilityStub>(parameter_json_settings); // invalid after the constructor
    KratosGeoSettlement settlement(std::move(input_utility));

    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    ExpectNumberOfMaterialCallsEqualTo(1, settlement);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageDoesNotPerformMaterialCallWhenNotSpecified, KratosGeoMechanicsFastSuite) {
    const std::string parameter_json_string_without_material_import_settings = "{"
                                                                               "\"solver_settings\":"
                                                                               "{"
                                                                               "\"model_part_name\":\"test\","
                                                                               "\"model_import_settings\":"
                                                                               "{"
                                                                               "\"input_type\": \"mdpa\","
                                                                               "\"input_filename\": \"mesh_stage1\""
                                                                               "}"
                                                                               "}"
                                                                               "}";

    auto input_utility = std::make_unique<InputUtilityStub>(parameter_json_string_without_material_import_settings); // invalid after the constructor

    KratosGeoSettlement settlement(std::move(input_utility));

    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    ExpectNumberOfMaterialCallsEqualTo(0, settlement);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageTwiceOnlyCallsReadFromModelOnce, KratosGeoMechanicsFastSuite) {
    auto input_utility = std::make_unique<InputUtilityStub>(parameter_json_settings);
    KratosGeoSettlement settlement(std::move(input_utility));

    RunStage(settlement);
    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
}

}