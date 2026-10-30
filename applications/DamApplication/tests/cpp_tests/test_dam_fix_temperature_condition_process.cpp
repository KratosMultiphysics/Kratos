//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Focused tests of DamFixTemperatureConditionProcess: constant and table-based
//  values, fixity and release, TIME_UNIT_CONVERTER and interval handling.
//
#include <utility>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/table.h"
#include "includes/variables.h"

// Application includes
#include "custom_processes/dam_fix_temperature_condition_process.hpp"

namespace Kratos
{
namespace Testing
{
namespace
{

ModelPart& BuildModel(
    Model& rModel,
    const Variable<double>& rVariable,
    const bool AddDofs)
{
    ModelPart& r_model_part = rModel.CreateModelPart("M");
    r_model_part.AddNodalSolutionStepVariable(rVariable);
    for (int i = 1; i <= 3; ++i) {
        auto p_node = r_model_part.CreateNewNode(
            i, static_cast<double>(i), 0.0, 0.0);
        if (AddDofs) {
            p_node->AddDof(rVariable);
        }
    }
    return r_model_part;
}

void AddTable(
    ModelPart& rModelPart,
    const int TableId,
    const std::vector<std::pair<double, double>>& rRows)
{
    auto p_table = Kratos::make_shared<Table<double, double>>();
    for (const auto& r_row : rRows) {
        p_table->PushBack(r_row.first, r_row.second);
    }
    rModelPart.AddTable(TableId, p_table);
}

void CheckNodalValue(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const double ExpectedValue)
{
    for (auto& r_node : rModelPart.Nodes()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(rVariable), ExpectedValue);
    }
}

void CheckFixity(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const bool ExpectedFixed)
{
    for (auto& r_node : rModelPart.Nodes()) {
        KRATOS_EXPECT_EQ(r_node.IsFixed(rVariable), ExpectedFixed);
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureAssignsConstantValueWithoutFixity, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : false,
        "value" : 12.5,
        "table" : 0
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, TEMPERATURE, 12.5);
    CheckFixity(r_model_part, TEMPERATURE, false);

    // The constant value is re-applied in every solution step.
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = -1.0;
    }
    r_model_part.GetProcessInfo()[TIME] = 2.0;
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 12.5);
    CheckFixity(r_model_part, TEMPERATURE, false);

    process.ExecuteFinalizeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 12.5);
    CheckFixity(r_model_part, TEMPERATURE, false);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureFixesAndReleases, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : true,
        "value" : 10.0,
        "table" : 0
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, TEMPERATURE, 10.0);
    CheckFixity(r_model_part, TEMPERATURE, true);

    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 10.0);
    CheckFixity(r_model_part, TEMPERATURE, true);

    // ExecuteFinalizeSolutionStep releases the prescribed DOFs.
    process.ExecuteFinalizeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 10.0);
    CheckFixity(r_model_part, TEMPERATURE, false);

    // The next step fixes them again.
    process.ExecuteInitializeSolutionStep();
    CheckFixity(r_model_part, TEMPERATURE, true);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureTableUsesConvertedTime, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 2.0;
    AddTable(r_model_part, 1, {{0.0, 1.0}, {2.0, 5.0}});

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : false,
        "value" : 100.0,
        "table" : 1
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    // The initial value is assigned before the solution loop: the table is not
    // read in ExecuteBeforeSolutionLoop.
    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, TEMPERATURE, 100.0);
    CheckFixity(r_model_part, TEMPERATURE, false);

    // TIME / TIME_UNIT_CONVERTER = 2.0 / 2.0 = 1.0 -> table value 3.0
    r_model_part.GetProcessInfo()[TIME] = 2.0;
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 3.0);

    // 4.0 / 2.0 = 2.0 -> table value 5.0
    r_model_part.GetProcessInfo()[TIME] = 4.0;
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 5.0);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureTableFixesAndReleases, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 1.0;
    AddTable(r_model_part, 1, {{0.0, 2.0}, {1.0, 4.0}});

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : true,
        "value" : 7.0,
        "table" : 1
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, TEMPERATURE, 7.0);
    CheckFixity(r_model_part, TEMPERATURE, true);

    r_model_part.GetProcessInfo()[TIME] = 1.0;
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 4.0);
    CheckFixity(r_model_part, TEMPERATURE, true);

    process.ExecuteFinalizeSolutionStep();
    CheckFixity(r_model_part, TEMPERATURE, false);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureIgnoresIntervalParameter, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : false,
        "value" : 3.25,
        "table" : 0,
        "interval" : [5.0, 6.0]
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    // The process does not implement standard IntervalUtility semantics: the
    // value is assigned regardless of the configured interval.
    r_model_part.GetProcessInfo()[TIME] = 1.0;
    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, TEMPERATURE, 3.25);

    r_model_part.GetProcessInfo()[TIME] = 8.0;
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, TEMPERATURE, 3.25);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureWorksForFaceHeatFlux, KratosDamFastSuite)
{
    // Mirrors the usage from ImposeFaceHeatFluxProcess (table branch), which
    // forwards the two assignment callbacks but not ExecuteFinalizeSolutionStep.
    Model model;
    ModelPart& r_model_part = BuildModel(model, FACE_HEAT_FLUX, false);
    r_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER] = 2.0;
    AddTable(r_model_part, 1, {{0.0, 0.0}, {1.0, 50.0}});

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "FACE_HEAT_FLUX",
        "constrained" : false,
        "value" : 0.0,
        "table" : 1
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    process.ExecuteBeforeSolutionLoop();
    CheckNodalValue(r_model_part, FACE_HEAT_FLUX, 0.0);

    r_model_part.GetProcessInfo()[TIME] = 2.0; // 2.0 / 2.0 -> table value 50.0
    process.ExecuteInitializeSolutionStep();
    CheckNodalValue(r_model_part, FACE_HEAT_FLUX, 50.0);
}

KRATOS_TEST_CASE_IN_SUITE(DamFixTemperatureFinalizeAlwaysReleasesDofs, KratosDamFastSuite)
{
    // Historical behavior: ExecuteFinalizeSolutionStep frees the DOF even when
    // constrained is false.
    Model model;
    ModelPart& r_model_part = BuildModel(model, TEMPERATURE, true);

    Parameters parameters(R"({
        "model_part_name" : "M",
        "variable_name" : "TEMPERATURE",
        "constrained" : false,
        "value" : 1.0,
        "table" : 0
    })");
    DamFixTemperatureConditionProcess process(r_model_part, parameters);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(TEMPERATURE);
    }
    CheckFixity(r_model_part, TEMPERATURE, true);

    process.ExecuteFinalizeSolutionStep();
    CheckFixity(r_model_part, TEMPERATURE, false);
}

} // namespace Testing
} // namespace Kratos
