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
#include "custom_processes/apply_vector_constraint_table_process.h"
#include "custom_processes/geo_extrapolate_integration_point_values_to_nodes_process.hpp"
#include "custom_utilities/process_utilities.h"
#include "includes/checks.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetModelPartsFromSettings_SingleModelPart)
{
    // Arrange
    Model model;
    model.CreateModelPart("Main");

    Parameters settings(R"(
        {
            "model_part_name": "Main"
        })");

    // Act
    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    // Assert
    EXPECT_EQ(model_parts.size(), 1);
    EXPECT_EQ(model_parts[0].get().Name(), "Main");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetModelPartsFromSettings_ListOfModelParts)
{
    // Arrange
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
    EXPECT_EQ(model_parts.size(), 2);
    EXPECT_EQ(model_parts[0].get().Name(), "Part1");
    EXPECT_EQ(model_parts[1].get().Name(), "Part2");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetModelPartsFromSettings_CheckForDuplicatedNames)
{
    // Arrange
    Model model;
    model.CreateModelPart("Part1");
    model.CreateModelPart("Part2");

    Parameters settings(R"(
        {
            "model_part_name_list": ["Part1", "Part1"]
        })");

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess"),
                                      "model_part_name_list has duplicated names for TestProcess.")
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
    {"ApplyVectorConstraintTableProcess",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyVectorConstraintTableProcess>(rModel, rSettings);
}},
    {"GeoExtrapolateIntegrationPointValuesToNodesProcess",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::GeoExtrapolateIntegrationPointValuesToNodesProcess>(rModel, rSettings);
}},
    {"ActivateModelPartOperation",
     [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ActivateModelPartOperation>(rModel, rSettings);
}},
    {"DeactivateModelPartOperation", [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::DeactivateModelPartOperation>(rModel, rSettings);
}}};

INSTANTIATE_TEST_SUITE_P(ProcessAndOperationUtilitiesTests, ModelPartsTest, ::testing::ValuesIn(ProcessesAndOperations));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AddProcessesSubModelPartList_AddsEmptyListWhenThereAreNoProcesses)
{
    // Arrange
    const auto project_parameters = Parameters{};
    Parameters solver_settings(R"(
        {
            "model_part_name": "Main"
        })");

    // Act
    ProcessUtilities::AddProcessesSubModelPartListToSolverSettings(project_parameters, solver_settings);

    // Assert
    KRATOS_EXPECT_TRUE(solver_settings.Has("processes_sub_model_part_list"))
    KRATOS_EXPECT_TRUE(solver_settings["processes_sub_model_part_list"].IsArray())
    KRATOS_EXPECT_TRUE(solver_settings["processes_sub_model_part_list"].size() == 0)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       AddProcessesSubModelPartList_AddsFilledListBasedOnConstraintsLoadsAndAuxiliaryProcessesOnly)
{
    // Arrange
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
                    "model_part_name_list": ["PorousDomain.BottomFixed", "PorousDomain.TopLoad"]
                }
            }],
            "loads_process_list": [
            {
                "Parameters":    {
                    "model_part_name": "PorousDomain.Soil"
                }
            }],
            "auxiliary_process_list": [{
                "Parameters": {
                    "model_part_name": "PorousDomain.Soil"
                }
            }],
            "arbitrary_process_list": [{
                "Parameters": {
                    "model_part_name": "PorousDomain"
                }
            }]
         },
         "solver_settings": {
             "model_part_name": "PorousDomain",
             "processes_sub_model_part_list": ["name1", "name2", "name3"]
         }
      })");

    auto solver_settings = project_parameters["solver_settings"];

    // Act
    ProcessUtilities::AddProcessesSubModelPartListToSolverSettings(project_parameters, solver_settings);

    // Assert
    const auto& r_list = solver_settings["processes_sub_model_part_list"];
    KRATOS_EXPECT_EQ(r_list.size(), 4);

    const auto actual_modelpart_names = solver_settings["processes_sub_model_part_list"].GetStringArray();
    const std::vector<std::string> expected_modelpart_names = {"BottomFixed", "Sides", "Soil", "TopLoad"};
    EXPECT_EQ(actual_modelpart_names, expected_modelpart_names);
}

} // namespace Kratos::Testing