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
#include "apply_scalar_constraints_table_process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "apply_component_table_process.hpp"
#include "processes/apply_constant_scalarvalue_process.h"
#include "apply_hydrostatic_pressure_table_process.hpp"
#include "apply_constant_phreatic_line_pressure_process.hpp"
#include "apply_phreatic_line_pressure_table_process.hpp"
#include "apply_constant_interpolate_line_pressure_process.hpp"
#include "apply_constant_phreatic_surface_pressure_process.hpp"
#include "apply_phreatic_surface_pressure_table_process.hpp"

namespace
{

using namespace Kratos;

Parameters ExtractParameters(const Parameters&               rSourceParameters,
                             const std::vector<std::string>& rNamesOfParametersToCopy)
{
    auto result = Parameters{};
    result.CopyValuesFromExistingParameters(rSourceParameters, rNamesOfParametersToCopy);
    return result;
}

void AppendParameterNameIfExists(const std::string&        rParameterName,
                                 const Parameters&         rSourceParameters,
                                 std::vector<std::string>& rResult)
{
    if (rSourceParameters.Has(rParameterName)) {
        rResult.emplace_back(rParameterName);
    }
}

bool HasTableAttached(const Parameters& rSettings)
{
    if (rSettings["table"].IsNumber()) {
        return rSettings["table"].GetInt() != 0;
    }

    KRATOS_ERROR_IF_NOT(rSettings["table"].IsArray()) << "'table' is neither a single number nor an array of numbers";

    const auto& table = rSettings["table"];
    return std::any_of(table.begin(), table.end(),
                       [](const auto& value){return value.GetInt() != 0;});
}

}

namespace Kratos
{

ApplyScalarConstraintsTableProcess::ApplyScalarConstraintsTableProcess(ModelPart&        rModelPart,
                                                                       const Parameters& rProcessSettings)
        : Process(Flags()),
          mrModelPart{rModelPart}
{
    MakeInternalProcess(rProcessSettings);
}

ApplyScalarConstraintsTableProcess::~ApplyScalarConstraintsTableProcess() = default;

void ApplyScalarConstraintsTableProcess::MakeInternalProcess(const Parameters& rProcessSettings)
{
    auto names_of_settings_to_copy = std::vector<std::string>{"model_part_name",
                                                              "variable_name"};
    AppendParameterNameIfExists("is_fixed", rProcessSettings, names_of_settings_to_copy);

    if (rProcessSettings.Has("fluid_pressure_type")) {
        MakeProcessForFluidPressureType(rProcessSettings, std::move(names_of_settings_to_copy));
    }
    else {
        MakeScalarConstraintsProcess(rProcessSettings, std::move(names_of_settings_to_copy));
    }
}

void ApplyScalarConstraintsTableProcess::MakeProcessForFluidPressureType(const Parameters&        rProcessSettings,
                                                                         std::vector<std::string> NamesOfSettingsToCopy)
{
    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();
    if (fluid_pressure_type == "Uniform") {
        MakeScalarConstraintsProcess(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Hydrostatic") {
        MakeProcessForHydrostaticFluidPressure(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Line") {
        MakeProcessForPhreaticLine(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Interpolate_Line") {
        MakeProcessForInterpolatedLine(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else if (fluid_pressure_type == "Phreatic_Surface") {
        MakeProcessForPhreaticSurface(rProcessSettings, std::move(NamesOfSettingsToCopy));
    } else {
        KRATOS_ERROR << "Unknown fluid_pressure_type: " << fluid_pressure_type << std::endl;
    }
}

void ApplyScalarConstraintsTableProcess::MakeScalarConstraintsProcess(const Parameters&        rProcessSettings,
                                                                      std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.emplace_back("value");

    if (HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyComponentTableProcess>(mrModelPart,
                                                                ExtractParameters(rProcessSettings,
                                                                                  NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantScalarValueProcess>(mrModelPart,
                                                                     ExtractParameters(rProcessSettings,
                                                                                       NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintsTableProcess::MakeProcessForHydrostaticFluidPressure(const Parameters&        rProcessSettings,
                                                                                std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "reference_coordinate",
                                                               "specific_weight"});
    AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyHydrostaticPressureTableProcess>(mrModelPart,
                                                                          ExtractParameters(rProcessSettings,
                                                                                            NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantHydrostaticPressureProcess>(mrModelPart,
                                                                             ExtractParameters(rProcessSettings,
                                                                                               NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintsTableProcess::MakeProcessForPhreaticLine(const Parameters&        rProcessSettings,
                                                                    std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "out_of_plane_direction",
                                                               "first_reference_coordinate",
                                                               "second_reference_coordinate",
                                                               "specific_weight"});
    AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyPhreaticLinePressureTableProcess>(mrModelPart,
                                                                           ExtractParameters(rProcessSettings,
                                                                                             NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantPhreaticLinePressureProcess>(mrModelPart,
                                                                              ExtractParameters(rProcessSettings,
                                                                                                NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintsTableProcess::MakeProcessForPhreaticSurface(const Parameters&        rProcessSettings,
                                                                       std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "first_reference_coordinate",
                                                               "second_reference_coordinate",
                                                               "third_reference_coordinate",
                                                               "specific_weight"});
    AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyPhreaticSurfacePressureTableProcess>(mrModelPart,
                                                                              ExtractParameters(rProcessSettings,
                                                                                                NamesOfSettingsToCopy));
    }
    else {
        mProcess = std::make_unique<ApplyConstantPhreaticSurfacePressureProcess>(mrModelPart,
                                                                                 ExtractParameters(rProcessSettings,
                                                                                                   NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintsTableProcess::MakeProcessForInterpolatedLine(const Parameters&        rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    KRATOS_ERROR_IF(HasTableAttached(rProcessSettings)) << "No time dependent interpolate line pressure process available" << std::endl;

    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "out_of_plane_direction"});
    AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    mProcess = std::make_unique<ApplyConstantInterpolateLinePressureProcess>(mrModelPart,
                                                                             ExtractParameters(rProcessSettings,
                                                                                               NamesOfSettingsToCopy));

}

void ApplyScalarConstraintsTableProcess::ExecuteInitialize()
{
    mProcess->ExecuteInitialize();
}

void ApplyScalarConstraintsTableProcess::ExecuteInitializeSolutionStep()
{
    mProcess->ExecuteInitializeSolutionStep();
}

}
