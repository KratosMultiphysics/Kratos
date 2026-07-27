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
#include "geo_apply_constant_scalar_value_process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{
using namespace std::string_literals;

ApplyScalarConstraintTableProcess::ApplyScalarConstraintTableProcess(ModelPart& rModelPart,
                                                                     const Parameters& rProcessSettings)
    : Process(Flags()), mrModelPart{rModelPart}
{
    MakeInternalProcess(rProcessSettings);
}

void ApplyScalarConstraintTableProcess::MakeInternalProcess(const Parameters& rProcessSettings)
{
    auto names_of_settings_to_copy = std::vector<std::string>{"model_part_name", "variable_name"};
    ParametersUtilities::AppendParameterNameIfExists("is_fixed", rProcessSettings, names_of_settings_to_copy);

    if (rProcessSettings.Has("fluid_pressure_type")) {
        MakeProcessForFluidPressureType(rProcessSettings, std::move(names_of_settings_to_copy));
    } else {
        MakeScalarConstraintProcess(rProcessSettings, std::move(names_of_settings_to_copy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForFluidPressureType(const Parameters& rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();
    if (fluid_pressure_type == "Uniform") {
        MakeScalarConstraintProcess(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Hydrostatic") {
        MakeProcessForHydrostaticFluidPressure(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Line") {
        MakeProcessForPhreaticLine(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Multi_Line") {
        MakeProcessForPhreaticMultiLine(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Interpolate_Line") {
        MakeProcessForInterpolatedLine(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Surface") {
        MakeProcessForPhreaticSurface(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else {
        KRATOS_ERROR << "Unknown fluid_pressure_type: " << fluid_pressure_type << std::endl;
    }
}

void ApplyScalarConstraintTableProcess::MakeScalarConstraintProcess(const Parameters& rProcessSettings,
                                                                    std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.emplace_back("value");

    InstantiateProcessByTablePresence<ApplyComponentTableProcess, GeoApplyConstantScalarValueProcess>(
        rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForHydrostaticFluidPressure(const Parameters& rProcessSettings,
                                                                               std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "reference_coordinate", "specific_weight"});
    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyHydrostaticPressureTableProcess, ApplyConstantHydrostaticPressureProcess>(
        rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticLine(const Parameters& rProcessSettings,
                                                                   std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction", "first_reference_coordinate",
                                  "second_reference_coordinate", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticLinePressureTableProcess, ApplyConstantPhreaticLinePressureProcess>(
        rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticMultiLine(const Parameters& rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction", "x_coordinates",
                                  "y_coordinates", "z_coordinates", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticMultiLinePressureTableProcess, ApplyConstantPhreaticMultiLinePressureProcess>(
        rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticSurface(const Parameters& rProcessSettings,
                                                                      std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "first_reference_coordinate", "second_reference_coordinate",
                                  "third_reference_coordinate", "specific_weight"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    InstantiateProcessByTablePresence<ApplyPhreaticSurfacePressureTableProcess, ApplyConstantPhreaticSurfacePressureProcess>(
        rProcessSettings, std::move(NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::MakeProcessForInterpolatedLine(const Parameters& rProcessSettings,
                                                                       std::vector<std::string> NamesOfSettingsToCopy)
{
    KRATOS_ERROR_IF(ParametersUtilities::HasTableAttached(rProcessSettings))
        << "No time dependent interpolate line pressure process available" << std::endl;

    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(),
                                 {"gravity_direction", "out_of_plane_direction"});

    AppendOptionalFluidParameters(rProcessSettings, NamesOfSettingsToCopy);

    mProcess = std::make_unique<ApplyConstantInterpolateLinePressureProcess>(
        mrModelPart, ParametersUtilities::CopyRequiredParameters(rProcessSettings, NamesOfSettingsToCopy));
}

void ApplyScalarConstraintTableProcess::ExecuteInitialize() { mProcess->ExecuteInitialize(); }

void ApplyScalarConstraintTableProcess::ExecuteInitializeSolutionStep()
{
    mProcess->ExecuteInitializeSolutionStep();
}

void ApplyScalarConstraintTableProcess::ExecuteFinalize() { mProcess->ExecuteFinalize(); }

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
void ApplyScalarConstraintTableProcess::InstantiateProcessByTablePresence(const Parameters& rProcessSettings,
                                                                          std::vector<std::string>&& rNamesOfSettingsToCopy)
{
    if (ParametersUtilities::HasTableAttached(rProcessSettings)) {
        rNamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<TableProcessType>(
            mrModelPart, ParametersUtilities::CopyRequiredParameters(
                             rProcessSettings, std::move(rNamesOfSettingsToCopy)));
    } else {
        mProcess = std::make_unique<ConstantProcessType>(
            mrModelPart, ParametersUtilities::CopyRequiredParameters(
                             rProcessSettings, std::move(rNamesOfSettingsToCopy)));
    }
}
} // namespace Kratos
