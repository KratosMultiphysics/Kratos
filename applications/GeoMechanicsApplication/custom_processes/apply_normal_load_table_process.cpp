// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Anne van de Graaf,
//                   Gennady Markelov
//
#include "apply_normal_load_table_process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/apply_constant_scalarvalue_process.h"
#include "apply_component_table_process.hpp"
#include "apply_boundary_hydrostatic_pressure_table_process.hpp"
#include "apply_boundary_phreatic_line_pressure_table_process.hpp"
#include "apply_boundary_phreatic_surface_pressure_table_process.hpp"
#include "custom_utilities/functions_for_parameters.h" 


namespace Kratos
{

ApplyNormalLoadTableProcess::ApplyNormalLoadTableProcess(ModelPart&        rModelPart,
                                                         const Parameters& rProcessSettings)
    : Process(Flags()),
    mrModelPart{ rModelPart }
{
    MakeInternalProcesses(rProcessSettings);
}

void ApplyNormalLoadTableProcess::MakeInternalProcesses(const Parameters& rProcessSettings)
{
    if (IsNormalComponentActive(rProcessSettings)) {
        MakeProcessForNormalComponent(rProcessSettings);
    }

    if (IsTangentialComponentActive(rProcessSettings)) {
        MakeProcessForTangentialComponent(rProcessSettings);
    }
}

bool ApplyNormalLoadTableProcess::IsNormalComponentActive(const Parameters& rProcessSettings)
{
    constexpr int normalComponentNumber = 0;
    return IsComponentActive(rProcessSettings, normalComponentNumber);
}

bool ApplyNormalLoadTableProcess::IsTangentialComponentActive(const Parameters& rProcessSettings)
{
    constexpr int tangentialComponentNumber = 1;
    return IsComponentActive(rProcessSettings, tangentialComponentNumber);
}

bool ApplyNormalLoadTableProcess::IsComponentActive(const Parameters& rProcessSettings, 
                                                        int componentNumber)
{
    return rProcessSettings["active"][componentNumber].GetBool();
}

void ApplyNormalLoadTableProcess::MakeProcessForNormalComponent(const Parameters& rProcessSettings)
{
    auto NamesOfSettingsToCopy = std::vector<std::string>{ "model_part_name",
                                                            "variable_name" };
    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();

    if (fluid_pressure_type == "Uniform") {
        MakeProcessForUniformFluidPressureType(rProcessSettings, std::move(NamesOfSettingsToCopy));
    }
    else if (fluid_pressure_type == "Hydrostatic") {
        MakeProcessForHydrostaticFluidPressureType(rProcessSettings, std::move(NamesOfSettingsToCopy));
    }
    else if (fluid_pressure_type == "Phreatic_Line") {
        MakeProcessForPhreaticLineFluidPressureType(rProcessSettings, std::move(NamesOfSettingsToCopy));
    }
    else if (fluid_pressure_type == "Phreatic_Surface") {
        MakeProcessForPhreaticSurfaceFluidPressureType(rProcessSettings, std::move(NamesOfSettingsToCopy));
    }
    else {
        KRATOS_ERROR << "Unknown fluid_pressure_type: " << fluid_pressure_type << std::endl;
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForUniformFluidPressureType(const Parameters& rProcessSettings,
                                                                         std::vector<std::string> NamesOfSettingsToCopy)
{  
    auto normal_parameters = FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                      NamesOfSettingsToCopy);
    normal_parameters.AddValue("value", rProcessSettings["value"][0]);

    if (FunctonsForParameters::HasTableAttached(rProcessSettings,0)) {
        normal_parameters.AddValue("table", rProcessSettings["table"][0]);
        mProcesses.push_back(std::make_unique <ApplyComponentTableProcess>(mrModelPart, 
                                                                           normal_parameters));
    } else {
        mProcesses.push_back(std::make_unique <ApplyConstantScalarValueProcess>(mrModelPart, 
                                                                                normal_parameters));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForHydrostaticFluidPressureType(const Parameters& rProcessSettings,
                                                                             std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), { "gravity_direction",
                                                               "reference_coordinate",
                                                               "specific_weight" });
    if (FunctonsForParameters::HasTableAttached(rProcessSettings, 0)) {
        auto normal_parameters = FunctonsForParameters::ExtractParameters(rProcessSettings,
                                                                          NamesOfSettingsToCopy);
        normal_parameters.AddValue("table", rProcessSettings["table"][0]);
        mProcesses.push_back(std::make_unique <ApplyBoundaryHydrostaticPressureTableProcess>(mrModelPart, 
                                                                                             normal_parameters));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryHydrostaticPressureProcess>(mrModelPart, 
                                                                                                FunctonsForParameters::ExtractParameters(rProcessSettings,NamesOfSettingsToCopy)));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForPhreaticLineFluidPressureType(const Parameters& rProcessSettings,
                                                                              std::vector<std::string> NamesOfSettingsToCopy)
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), { "gravity_direction",
                                                               "out_of_plane_direction",
                                                               "first_reference_coordinate",
                                                               "second_reference_coordinate",
                                                               "specific_weight"});
    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), { "table" });
        mProcesses.push_back(std::make_unique <ApplyBoundaryPhreaticLinePressureTableProcess>(mrModelPart, 
                                                                                              FunctonsForParameters::ExtractParameters(rProcessSettings, NamesOfSettingsToCopy)));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryPhreaticLinePressureProcess>(mrModelPart, 
                                                                                                 FunctonsForParameters::ExtractParameters(rProcessSettings, NamesOfSettingsToCopy)));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForPhreaticSurfaceFluidPressureType(const Parameters& rProcessSettings,
                                                                                 std::vector<std::string> NamesOfSettingsToCopy) 
{
    NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), { "gravity_direction",
                                                                "first_reference_coordinate",
                                                                "second_reference_coordinate",
                                                                "third_reference_coordinate",
                                                                "specific_weight" });
    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        NamesOfSettingsToCopy.insert(NamesOfSettingsToCopy.end(), { "table" });
        mProcesses.push_back(std::make_unique <ApplyBoundaryPhreaticSurfacePressureTableProcess>(mrModelPart, 
                                                                                                 FunctonsForParameters::ExtractParameters(rProcessSettings, NamesOfSettingsToCopy)));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryPhreaticSurfacePressureProcess>(mrModelPart, 
                                                                                                    FunctonsForParameters::ExtractParameters(rProcessSettings, NamesOfSettingsToCopy)));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForTangentialComponent(const Parameters& rProcessSettings)
{
    auto tangential_params = Parameters{ "{}" };
    tangential_params.AddValue("model_part_name", rProcessSettings["model_part_name"]);
    tangential_params.AddEmptyValue("variable_name").SetString("TANGENTIAL_CONTACT_STRESS"); // Note : this is not general
    tangential_params.AddValue("value", rProcessSettings["value"][1]);
    if (FunctonsForParameters::HasTableAttached(rProcessSettings, 1)) {
        tangential_params.AddValue("table", rProcessSettings["table"][1]);
        mProcesses.push_back(std::make_unique <ApplyComponentTableProcess>(mrModelPart, 
                                                                           tangential_params));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantScalarValueProcess>(mrModelPart, 
                                                                                tangential_params));
    }
}

ApplyNormalLoadTableProcess::~ApplyNormalLoadTableProcess() = default;

void ApplyNormalLoadTableProcess::ExecuteInitialize()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitialize();
    }
}

void ApplyNormalLoadTableProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitializeSolutionStep();
    }
}

}
