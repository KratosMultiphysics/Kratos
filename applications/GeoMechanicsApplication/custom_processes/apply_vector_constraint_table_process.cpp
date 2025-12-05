// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra,
//                   Anne van de Graaf
//
#include "apply_vector_constraint_table_process.h"
#include "apply_component_table_process.h"
#include "custom_utilities/process_utilities.h"
#include "geo_apply_constant_scalar_value_process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{

ApplyVectorConstraintTableProcess::ApplyVectorConstraintTableProcess(Model& rModel, const Parameters& rSettings)
    : Process(Flags())
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(
        rModel, rSettings, ApplyVectorConstraintTableProcess::Info());
    for (const auto& r_model_part : mrModelParts) {
        const auto parameters_list = CreateParametersForActiveComponents(r_model_part, rSettings);
        for (const auto& parameters : parameters_list) {
            mProcesses.emplace_back(MakeProcessFor(r_model_part, parameters));
        }
    }
}

ApplyVectorConstraintTableProcess::~ApplyVectorConstraintTableProcess() = default;

std::vector<Parameters> ApplyVectorConstraintTableProcess::CreateParametersForActiveComponents(
    const ModelPart& rModelPart, const Parameters& rSettings)
{
    std::vector<Parameters> result;
    for (auto component : ActiveComponents(rSettings)) {
        result.emplace_back(CreateParametersForComponent(rModelPart, rSettings, component));
    }
    return result;
}

std::vector<char> ApplyVectorConstraintTableProcess::ActiveComponents(const Parameters& rSettings)
{
    std::vector<char> result;
    for (auto component : {'X', 'Y', 'Z'}) {
        const auto index = ComponentToIndex(component);
        if (rSettings["active"][index].GetBool()) {
            result.emplace_back(component);
        }
    }
    return result;
}

Parameters ApplyVectorConstraintTableProcess::CreateParametersForComponent(const ModelPart& rModelPart,
                                                                           const Parameters& rSettings,
                                                                           char component)
{
    Parameters result;
    const auto index = ComponentToIndex(component);
    result.AddValue("model_part_name", "\"" + rModelPart.Name() + "\"");
    if (rSettings.Has("is_fixed")) {
        result.AddValue("is_fixed", rSettings["is_fixed"][index]);
    }
    result.AddValue("value", rSettings["value"][index]);
    const auto variable_name = rSettings["variable_name"].GetString() + '_' + component;
    result.AddEmptyValue("variable_name").SetString(variable_name);
    if (rSettings["table"][index].GetInt() != 0) {
        result.AddValue("table", rSettings["table"][index]);
    }
    return result;
}

std::size_t ApplyVectorConstraintTableProcess::ComponentToIndex(char component)
{
    switch (component) {
    case 'X':
        return 0;
    case 'Y':
        return 1;
    case 'Z':
        return 2;

    default:
        KRATOS_ERROR << "ApplyVectorConstraintTableProcess: Unknown component '" << component << "'"
                     << std::endl;
    }
}

ApplyVectorConstraintTableProcess::ProcessUniquePointer ApplyVectorConstraintTableProcess::MakeProcessFor(
    ModelPart& rModelPart, const Parameters& rParameters) const
{
    if (rParameters.Has("table")) {
        return std::make_unique<ApplyComponentTableProcess>(rModelPart, rParameters);
    }

    return std::make_unique<GeoApplyConstantScalarValueProcess>(rModelPart, rParameters);
}

void ApplyVectorConstraintTableProcess::ExecuteInitialize()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitialize();
    }
}

void ApplyVectorConstraintTableProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitializeSolutionStep();
    }
}

void ApplyVectorConstraintTableProcess::ExecuteFinalize()
{
    for (const auto& process : mProcesses) {
        process->ExecuteFinalize();
    }
}

std::string ApplyVectorConstraintTableProcess::Info() const
{
    return "ApplyVectorConstraintTableProcess";
}

} // namespace Kratos