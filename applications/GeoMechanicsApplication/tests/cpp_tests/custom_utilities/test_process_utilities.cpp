// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_operations/activate_model_part_operation.h"
#include "custom_operations/deactivate_model_part_operation.h"
#include "custom_processes/apply_c_phi_reduction_process.h"
#include "custom_processes/apply_excavation_process.h"
#include "custom_processes/apply_final_stresses_of_previous_stage_to_initial_state.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "custom_utilities/process_utilities.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <iostream>
#include <sstream>

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_SingleModelPart, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    Model model;
    model.CreateModelPart("Main");

    Parameters settings(R"(
        {
            "model_part_name": "Main"
        })");

    // Act
    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    // Assert
    KRATOS_EXPECT_EQ(model_parts.size(), 1);
    KRATOS_EXPECT_EQ(model_parts[0].get().Name(), "Main");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_ListOfModelParts, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    Model model;
    model.CreateModelPart("Part1");
    model.CreateModelPart("Part2");

    Parameters settings(R"(
        {
            "model_part_name_list": ["Part1", "Part2"]
        })");

    // Act
    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    // Assert
    KRATOS_EXPECT_EQ(model_parts.size(), 2);
    KRATOS_EXPECT_EQ(model_parts[0].get().Name(), "Part1");
    KRATOS_EXPECT_EQ(model_parts[1].get().Name(), "Part2");
}

struct NamedFactory {
    std::string                                                    name;
    std::function<void(Kratos::Model&, const Kratos::Parameters&)> factory;
};

class ModelPartsTest : public ::testing::TestWithParam<NamedFactory>
{
};

TEST_P(ModelPartsTest, GetModelPartsFromSettings_BothParametersPresent_Throws)
{
    Model model;

    Parameters settings(R"(
        {
            "model_part_name": "Part1",
            "model_part_name_list": ["Part2"],
            "deactivate_soil_part": false
        })");

    const auto& param = GetParam();

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        param.factory(model, settings),
        "The parameters 'model_part_name' and 'model_part_name_list' are mutually exclusive for " +
            param.name);
}

TEST_P(ModelPartsTest, GetModelPartsFromSettings_MissingParameters_Throws)
{
    Model      model;
    Parameters settings(R"({"deactivate_soil_part": false})");

    const auto& param = GetParam();

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        param.factory(model, settings),
        "Please specify 'model_part_name' or 'model_part_name_list' for " + param.name);
}

TEST_P(ModelPartsTest, GetModelPartsFromSettings_EmptyList_Throws)
{
    Model      model;
    Parameters settings(R"(
        {
            "model_part_name_list": [],
            "deactivate_soil_part": false
        })");

    const auto& param = GetParam();

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        param.factory(model, settings),
        "The parameters 'model_part_name_list' needs to contain at least one model part name for " +
            param.name);
}

static const std::vector<NamedFactory> ProcessesAndOperations = {
    {"ApplyExcavationProcess",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyExcavationProcess>(rModel, rSettings);
}},
    {"ApplyK0ProcedureProcess",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyK0ProcedureProcess>(rModel, rSettings);
}},
    {"ApplyFinalStressesOfPreviousStageToInitialState",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyFinalStressesOfPreviousStageToInitialState>(rModel, rSettings);
}},
    {"ApplyCPhiReductionProcess",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyCPhiReductionProcess>(rModel, rSettings);
}},
    {"ActivateModelPartOperation",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ActivateModelPartOperation>(rModel, rSettings);
}},
    {"DeactivateModelPartOperation", [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::DeactivateModelPartOperation>(rModel, rSettings);
}}};

INSTANTIATE_TEST_SUITE_P(ProcessAndOperationUtilitiesTests, ModelPartsTest, ::testing::ValuesIn(ProcessesAndOperations));

KRATOS_TEST_CASE_IN_SUITE(AddProcessesSubModelPartList_AddsEmptyList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto project_parameters = Parameters{};
    Parameters solver_settings(R"(
        {
            "model_part_name": "Main"
        })");

    // Act
    ProcessUtilities::AddProcessesSubModelPartList(project_parameters, solver_settings);

    // Assert
    KRATOS_EXPECT_TRUE(solver_settings.Has("processes_sub_model_part_list"))
    KRATOS_EXPECT_TRUE(solver_settings["processes_sub_model_part_list"].IsArray())
    KRATOS_EXPECT_TRUE(solver_settings["processes_sub_model_part_list"].size() == 0)
}

KRATOS_TEST_CASE_IN_SUITE(AddProcessesSubModelPartList_AddsFilledList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const Parameters project_parameters(R"(
      {
         "processes":
         {
            "constraints_process_list": [
            {
                "Parameters":    {
                    "model_part_name": "PorousDomain.Sides"
                }
            },{
                "Parameters":    {
                    "model_part_name": "PorousDomain.BottomFixed"
                }
            },{
                "Parameters":    {
                   "model_part_name": "PorousDomain.TopLoad"
                }
            }],
            "loads_process_list": [
            {
                "Parameters":    {
                    "model_part_name": "PorousDomain.Soil"
                }
            },{
                "Parameters":    {
                   "model_part_name": "PorousDomain"
                }
            }],
            "auxiliary_process_list": [{
                "Parameters": {
                    "model_part_name": "PorousDomain.Soil"
                }
            }],
            "arbitrary_process_list": [{
                "Parameters": {
                    "model_part_name": "PorousDomain.Domain"
                }
            }]
         },
         "solver_settings": {
             "model_part_name": "PorousDomain",
             "processes_sub_model_part_list": ["name1", "name2", "name3"]
         }
      })");

    auto solver_settings = project_parameters["solver_settings"];

    std::stringstream buffer;
    std::streambuf*   old = std::cout.rdbuf(buffer.rdbuf()); // NOSONAR

    // Act
    ProcessUtilities::AddProcessesSubModelPartList(project_parameters, solver_settings);

    // Assert
    std::cout.rdbuf(old); // NOSONAR
    std::string output = buffer.str();
    EXPECT_NE(output.find("processes_sub_model_part_list is deprecated"), std::string::npos);

    const auto& r_list = solver_settings["processes_sub_model_part_list"];
    KRATOS_EXPECT_EQ(r_list.size(), 4);

    std::set<std::string, std::less<>> actual_modelpart_names;
    for (const auto& r_name : r_list) {
        actual_modelpart_names.insert(r_name.GetString());
    }
    const std::set<std::string, std::less<>> expected_modelpart_names = {"Sides", "BottomFixed",
                                                                         "TopLoad", "Soil"};
    EXPECT_EQ(actual_modelpart_names, expected_modelpart_names);
}

} // namespace Kratos::Testing