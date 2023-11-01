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
#include "stub_input_utility.h"
#include "stub_process_info_parser.h"
#include "stub_time_loop_executor.h"

using namespace Kratos;

namespace Kratos::Testing {

const std::string parameter_json_settings = R"(
                                            {
                                                "solver_settings":
                                                {
                                                    "model_part_name":"test",
                                                    "model_import_settings":
                                                    {
                                                        "input_type": "mdpa",
                                                        "input_filename": "mesh_stage1"
                                                    },
                                                    "material_import_settings":
                                                    {
                                                        "materials_filename": "MaterialParameters1.json"
                                                    },
                                                    "problem_domain_sub_model_part_list": [],
                                                    "processes_sub_model_part_list": []
                                                }
                                            }
                                            )"; // these have material_import settings


int RunStage(KratosGeoSettlement& rSettlement) {
    return rSettlement.RunStage("",
                                "",
                                [](const char *) {/* kept empty as a stub method */},
                                [](const double) {/* kept empty as a stub method */},
                                [](const char *) {/* kept empty as a stub method */},
                                []() { return false; });
}

void ExpectNumberOfReadCallsIsEqualToOne(const KratosGeoSettlement &rSettlement)
{
    const auto input_utility_from_settlement = dynamic_cast<const StubInputUtility*>(rSettlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(input_utility_from_settlement, nullptr);
    KRATOS_EXPECT_EQ(input_utility_from_settlement->NumberOfReadCalls(), 1);
}

void ExpectNumberOfMaterialCallsEqualTo(int expectedNumberOfMaterialCalls, const KratosGeoSettlement& rSettlement)
{
    const auto input_utility_from_settlement = dynamic_cast<const StubInputUtility*>(rSettlement.GetInterfaceInputUtility());
    KRATOS_EXPECT_NE(input_utility_from_settlement, nullptr);
    KRATOS_EXPECT_EQ(input_utility_from_settlement->NumberOfMaterialCalls(), expectedNumberOfMaterialCalls);
}

KRATOS_TEST_CASE_IN_SUITE(CreatingKratosGeoSettlementDoesNotThrow, KratosGeoMechanicsFastSuite) {
    bool has_thrown = false;
    try
    {
        KratosGeoSettlement settlement(std::make_unique<StubInputUtility>(parameter_json_settings),
                                       std::make_unique<StubProcessInfoParser>(),
                                       std::make_unique<StubTimeLoopExecutor>());
    }
    catch (...)
    {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(RunStageMakesRelevantCallsOnce, KratosGeoMechanicsFastSuite) {
    KratosGeoSettlement settlement(std::make_unique<StubInputUtility>(parameter_json_settings),
                                   std::make_unique<StubProcessInfoParser>(),
                                   std::make_unique<StubTimeLoopExecutor>());

    const auto exit_status = RunStage(settlement);

    // At the moment the parameter_json_settings do not contain enough info to instantiate the solving strategy.
    // That is why the exit_status is now 1, as with other tests in this suite.
    KRATOS_EXPECT_EQ(exit_status, 1);
    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    ExpectNumberOfMaterialCallsEqualTo(1, settlement);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageDoesNotPerformMaterialCallWhenNotSpecified, KratosGeoMechanicsFastSuite) {
    const std::string parameter_json_string_without_material_import_settings = R"(
                                                                                {
                                                                                    "solver_settings":
                                                                                    {
                                                                                        "model_part_name":"test",
                                                                                        "model_import_settings":
                                                                                        {
                                                                                            "input_type": "mdpa",
                                                                                            "input_filename": "mesh_stage1"
                                                                                        }
                                                                                    }
                                                                                }
                                                                                )";

    KratosGeoSettlement settlement(std::make_unique<StubInputUtility>(parameter_json_string_without_material_import_settings),
                                   std::make_unique<StubProcessInfoParser>(),
                                   std::make_unique<StubTimeLoopExecutor>());

    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
    ExpectNumberOfMaterialCallsEqualTo(0, settlement);
}

KRATOS_TEST_CASE_IN_SUITE(RunStageTwiceOnlyCallsReadFromModelOnce, KratosGeoMechanicsFastSuite) {
    KratosGeoSettlement settlement(std::make_unique<StubInputUtility>(parameter_json_settings),
                                   std::make_unique<StubProcessInfoParser>(),
                                   std::make_unique<StubTimeLoopExecutor>());

    RunStage(settlement);
    RunStage(settlement);

    ExpectNumberOfReadCallsIsEqualToOne(settlement);
}

KRATOS_TEST_CASE_IN_SUITE(ConstructKratosGeoSettlementWithEmptyInputUtilityThrows, KratosGeoMechanicsFastSuite) {
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(KratosGeoSettlement settlement(nullptr,
                                                                     std::make_unique<StubProcessInfoParser>(),
                                                                     std::make_unique<StubTimeLoopExecutor>()),
                                      "Invalid Input Utility")
}

KRATOS_TEST_CASE_IN_SUITE(RunStage_PassesTheCorrectProcessReferences, KratosGeoMechanicsFastSuite) {
    const std::string process_parameters_in_json_settings =
    R"(
    {
        "solver_settings":
        {
            "model_part_name":"test",
            "model_import_settings":
            {
            "input_type": "mdpa",
            "input_filename": "mesh_stage1"
            }
        },

        "processes" :
        [
            { "input" : "true" }
        ]
    }
    )";

    KratosGeoSettlement settlement(std::make_unique<StubInputUtility>(process_parameters_in_json_settings),
                                   std::make_unique<StubProcessInfoParser>(),
                                   std::make_unique<StubTimeLoopExecutor>(StubProcessInfoParser::NumberOfProcesses()));

    RunStage(settlement);
}

}