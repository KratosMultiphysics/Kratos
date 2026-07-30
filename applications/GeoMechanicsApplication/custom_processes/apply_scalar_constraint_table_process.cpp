// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf,
//                   Marjan Fathian
//
#include "apply_scalar_constraint_table_process.h"
#include "apply_component_table_process.h"
#include "apply_constant_interpolate_line_pressure_process.h"
#include "apply_constant_phreatic_line_pressure_process.h"
#include "apply_constant_phreatic_surface_pressure_process.h"
#include "apply_hydrostatic_pressure_table_process.h"
#include "apply_phreatic_line_pressure_table_process.h"
#include "apply_phreatic_multi_line_pressure_table_process.h"
#include "apply_phreatic_surface_pressure_table_process.h"
#include "custom_utilities/parameters_utilities.h"
#include "custom_utilities/process_utilities.h"
#include "geo_apply_constant_scalar_value_process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{
using namespace std::string_literals;

ApplyScalarConstraintTableProcess::ApplyScalarConstraintTableProcess(Model& rModel, const Parameters& rProcessSettings)
    : Process(Flags())
{
    const auto model_parts = ProcessUtilities::GetModelPartsFromSettings(
        rModel, rProcessSettings, "ApplyScalarConstraintTableProcess");
    mProcesses.clear();
    mProcesses.reserve(model_parts.size());
    for (const auto& r_model_part : model_parts)
        MakeInternalProcess(r_model_part, rProcessSettings);
}

void ApplyScalarConstraintTableProcess::MakeInternalProcess(ModelPart& rModelPart, const Parameters& rProcessSettings)
{
    auto names_of_settings_to_copy = std::vector<std::string>{"variable_name"};
    ParametersUtilities::AppendParameterNameIfExists("is_fixed", rProcessSettings, names_of_settings_to_copy);

    if (rProcessSettings.Has("fluid_pressure_type")) {
        MakeProcessForFluidPressureType(rModelPart, rProcessSettings, std::move(names_of_settings_to_copy));
    } else {
        MakeScalarConstraintProcess(rModelPart, rProcessSettings, std::move(names_of_settings_to_copy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForFluidPressureType(ModelPart& rModelPart,
                                                                        const Parameters& rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();
    if (fluid_pressure_type == "Uniform") {
        MakeScalarConstraintProcess(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Hydrostatic") {
        MakeProcessForHydrostaticFluidPressure(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Line") {
        MakeProcessForPhreaticLine(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Multi_Line") {
        MakeProcessForPhreaticMultiLine(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Interpolate_Line") {
        MakeProcessForInterpolatedLine(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Surface") {
        MakeProcessForPhreaticSurface(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else {
        KRATOS_ERROR << "Unknown fluid_pressure_type: " << fluid_pressure_type << std::endl;
    }
}

void ApplyScalarConstraintTableProcess::MakeScalarConstraintProcess(ModelPart& rModelPart,
                                                                    const Parameters& rProcessSettings,
                                                                    std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.emplace_back("value");

    InstantiateProcessByTablePresence<ApplyComponentTableProcess, GeoApplyConstantScalarValueProcess>(
        rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForHydrostaticFluidPressure(
    ModelPart& rModelPart, const Parameters& rProcessSettings, std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "reference_coordinate", "specific_weight"});
    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyHydrostaticPressureTableProcess, ApplyConstantHydrostaticPressureProcess>(
        rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticLine(ModelPart& rModelPart,
                                                                   const Parameters& rProcessSettings,
                                                                   std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction", "first_reference_coordinate",
                                  "second_reference_coordinate", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticLinePressureTableProcess, ApplyConstantPhreaticLinePressureProcess>(
        rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticMultiLine(ModelPart& rModelPart,
                                                                        const Parameters& rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction", "x_coordinates",
                                  "y_coordinates", "z_coordinates", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticMultiLinePressureTableProcess, ApplyConstantPhreaticMultiLinePressureProcess>(
        rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticSurface(ModelPart& rModelPart,
                                                                      const Parameters& rProcessSettings,
                                                                      std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "first_reference_coordinate", "second_reference_coordinate",
                                  "third_reference_coordinate", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticSurfacePressureTableProcess, ApplyConstantPhreaticSurfacePressureProcess>(
        rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForInterpolatedLine(ModelPart& rModelPart,
                                                                       const Parameters& rProcessSettings,
                                                                       std::vector<std::string> NamesOfSettingsToCopy)
{
    KRATOS_ERROR_IF(ParametersUtilities::HasTableAttached(rProcessSettings))
        << "No time dependent interpolate line pressure process available" << std::endl;

    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    mProcesses.emplace_back(std::make_unique<ApplyConstantInterpolateLinePressureProcess>(
        rModelPart, PrepareProcessParameters(rModelPart, rProcessSettings, std::move(NamesOfSettingsToCopy))));
}

void ApplyScalarConstraintTableProcess::ExecuteInitialize()
{
    for (auto& r_process : mProcesses)
        r_process->ExecuteInitialize();
}

void ApplyScalarConstraintTableProcess::ExecuteInitializeSolutionStep()
{
    for (auto& r_process : mProcesses)
        r_process->ExecuteInitializeSolutionStep();
}

void ApplyScalarConstraintTableProcess::ExecuteFinalize()
{
    for (auto& r_process : mProcesses)
        r_process->ExecuteFinalize();
}

std::string ApplyScalarConstraintTableProcess::Info() const
{
    return "ApplyScalarConstraintTableProcess"s;
}

void ApplyScalarConstraintTableProcess::AppendOptionalFluidParameters(const Parameters& rProcessSettings,
                                                                      std::vector<std::string>& rNamesOfSettingsToCopy) const
{
    ParametersUtilities::AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings,
                                                     rNamesOfSettingsToCopy);
    ParametersUtilities::AppendParameterNameIfExists("is_seepage", rProcessSettings, rNamesOfSettingsToCopy);
}

template <typename TableProcessType, typename ConstantProcessType>
void ApplyScalarConstraintTableProcess::InstantiateProcessByTablePresence(ModelPart& rModelPart,
                                                                          const Parameters& rProcessSettings,
                                                                          std::vector<std::string>&& rNamesOfSettingsToCopy)
{
    if (ParametersUtilities::HasTableAttached(rProcessSettings)) {
        rNamesOfSettingsToCopy.emplace_back("table");
        mProcesses.emplace_back(std::make_unique<TableProcessType>(
            rModelPart, PrepareProcessParameters(rModelPart, rProcessSettings, std::move(rNamesOfSettingsToCopy))));
    } else {
        mProcesses.emplace_back(std::make_unique<ConstantProcessType>(
            rModelPart, PrepareProcessParameters(rModelPart, rProcessSettings, std::move(rNamesOfSettingsToCopy))));
    }
}

Parameters ApplyScalarConstraintTableProcess::PrepareProcessParameters(const ModelPart& rModelPart,
                                                                       const Parameters& rProcessSettings,
                                                                       std::vector<std::string>&& rNamesOfSettingsToCopy)
{
    auto result = ParametersUtilities::CopyRequiredParameters(rProcessSettings, rNamesOfSettingsToCopy);

    result.AddString("model_part_name", rModelPart.Name());

    return result;
}
} // namespace Kratos
