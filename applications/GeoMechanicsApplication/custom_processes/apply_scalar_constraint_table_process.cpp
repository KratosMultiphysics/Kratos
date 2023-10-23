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
#include "custom_utilities/functions_for_parameters.h" 

namespace Kratos
{

ApplyScalarConstraintTableProcess::ApplyScalarConstraintTableProcess(ModelPart&        rModelPart,
                                                                     const Parameters& rProcessSettings)
        : Process(Flags()),
          mrModelPart{rModelPart}
{
    MakeInternalProcess(rProcessSettings);
}

ApplyScalarConstraintTableProcess::~ApplyScalarConstraintTableProcess() = default;

void ApplyScalarConstraintTableProcess::MakeInternalProcess(const Parameters& rProcessSettings)
{
    auto names_of_settings_to_copy = std::vector<std::string>{"model_part_name",
                                                              "variable_name"};
    FunctonsForParameters::AppendParameterNameIfExists("is_fixed", rProcessSettings, names_of_settings_to_copy);

    if (rProcessSettings.Has("fluid_pressure_type")) {
        MakeProcessForFluidPressureType(rProcessSettings, std::move(names_of_settings_to_copy));
    }
    else {
        MakeScalarConstraintProcess(rProcessSettings, std::move(names_of_settings_to_copy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForFluidPressureType(const Parameters&        rProcessSettings,
                                                                        std::vector<std::string> NamesOfSettingsToCopy)
{
    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();
    if (fluid_pressure_type == "Uniform") {
        MakeScalarConstraintProcess(rProcessSettings, std::move(NamesOfSettingsToCopy));
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

void ApplyScalarConstraintTableProcess::MakeScalarConstraintProcess(const Parameters&        rProcessSettings,
                                                                    std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.emplace_back("value");

    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyComponentTableProcess>(mrModelPart,
                                                                FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantScalarValueProcess>(mrModelPart,
                                                                    FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                    NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForHydrostaticFluidPressure(const Parameters&        rProcessSettings,
                                                                               std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "reference_coordinate",
                                                               "specific_weight"});
    FunctonsForParameters::AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    FunctonsForParameters::AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyHydrostaticPressureTableProcess>(mrModelPart,
                                                                          FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                          NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantHydrostaticPressureProcess>(mrModelPart,
                                                                            FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                            NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticLine(const Parameters&        rProcessSettings,
                                                                   std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "out_of_plane_direction",
                                                               "first_reference_coordinate",
                                                               "second_reference_coordinate",
                                                               "specific_weight"});
    FunctonsForParameters::AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    FunctonsForParameters::AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyPhreaticLinePressureTableProcess>(mrModelPart,
                                                                           FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                           NamesOfSettingsToCopy));
    } else {
        mProcess = std::make_unique<ApplyConstantPhreaticLinePressureProcess>(mrModelPart,
                                                                              FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                              NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForPhreaticSurface(const Parameters&        rProcessSettings,
                                                                      std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "first_reference_coordinate",
                                                               "second_reference_coordinate",
                                                               "third_reference_coordinate",
                                                               "specific_weight"});
    FunctonsForParameters::AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    FunctonsForParameters::AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.emplace_back("table");
        mProcess = std::make_unique<ApplyPhreaticSurfacePressureTableProcess>(mrModelPart,
                                                                              FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                              NamesOfSettingsToCopy));
    }
    else {
        mProcess = std::make_unique<ApplyConstantPhreaticSurfacePressureProcess>(mrModelPart,
                                                                                 FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                                 NamesOfSettingsToCopy));
    }
}

void ApplyScalarConstraintTableProcess::MakeProcessForInterpolatedLine(const Parameters&        rProcessSettings,
                                                                       std::vector<std::string> NamesOfSettingsToCopy)
{
    KRATOS_ERROR_IF(FunctonsForParameters::HasTableAttached(rProcessSettings)) << "No time dependent interpolate line pressure process available" << std::endl;

    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), {"gravity_direction",
                                                               "out_of_plane_direction"});
    FunctonsForParameters::AppendParameterNameIfExists("pressure_tension_cut_off", rProcessSettings, NamesOfSettingsToCopy);
    FunctonsForParameters::AppendParameterNameIfExists("is_seepage", rProcessSettings, NamesOfSettingsToCopy);

    mProcess = std::make_unique<ApplyConstantInterpolateLinePressureProcess>(mrModelPart,
                                                                             FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                             NamesOfSettingsToCopy));

}

void ApplyScalarConstraintTableProcess::ExecuteInitialize()
{
    mProcess->ExecuteInitialize();
}

void ApplyScalarConstraintTableProcess::ExecuteInitializeSolutionStep()
{
    mProcess->ExecuteInitializeSolutionStep();
}

}
