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

#include "custom_processes/apply_c_phi_reduction_process.h"
#include "custom_processes/apply_excavation_process.h"
#include "custom_processes/apply_final_stresses_of_previous_stage_to_initial_state.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "custom_utilities/process_utilities.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "custom_operations/activate_model_part_operation.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_SingleModelPart, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    model.CreateModelPart("Main");

    Parameters settings(R"(
        {
            "model_part_name": "Main"
        })");

    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    KRATOS_EXPECT_EQ(model_parts.size(), 1);
    KRATOS_EXPECT_EQ(model_parts[0].get().Name(), "Main");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_ListOfModelParts, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    model.CreateModelPart("Part1");
    model.CreateModelPart("Part2");

    Parameters settings(R"(
        {
            "model_part_name_list": ["Part1", "Part2"]
        })");

    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    KRATOS_EXPECT_EQ(model_parts.size(), 2);
    KRATOS_EXPECT_EQ(model_parts[0].get().Name(), "Part1");
    KRATOS_EXPECT_EQ(model_parts[1].get().Name(), "Part2");
}

struct NamedFactory {
    std::string name;
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

static const std::vector<NamedFactory> kProcessFactories = {
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
    {"ApplyCPhiReductionProcess", [](Model& rModel, const Parameters& rSettings) {
    return std::make_unique<Kratos::ApplyCPhiReductionProcess>(rModel, rSettings);
}}};

static const std::vector<NamedFactory> kOperationFactories = {
    {"ActivateModelPartOperation",
     [](Model& rModel, const Parameters& rSettings) {
         return std::make_unique<Kratos::ActivateModelPartOperation>(rModel, rSettings);
    }}
};

static std::vector<NamedFactory> GetAllFactories()
{
    std::vector<NamedFactory> all = kProcessFactories;
    all.insert(all.end(), kOperationFactories.begin(), kOperationFactories.end());
    return all;
}

INSTANTIATE_TEST_SUITE_P(ProcessAndOperationUtilitiesTests, ModelPartsTest, ::testing::ValuesIn(GetAllFactories()));

} // namespace Kratos::Testing