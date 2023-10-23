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
    auto normal_parameters = Parameters{ "{}" };
    normal_parameters.AddValue("model_part_name", rProcessSettings["model_part_name"]);
    normal_parameters.AddValue("variable_name", rProcessSettings["variable_name"]);

    const auto fluid_pressure_type = rProcessSettings["fluid_pressure_type"].GetString();

    if (fluid_pressure_type == "Uniform") {
        MakeProcessForUniformFluidPressureType(rProcessSettings, std::move(normal_parameters));
    }
    else if (fluid_pressure_type == "Hydrostatic") {
        MakeProcessForHydrostaticFluidPressureType(rProcessSettings, std::move(normal_parameters));
    }
    else if (fluid_pressure_type == "Phreatic_Line") {
        MakeProcessForPhreaticLineFluidPressureType(rProcessSettings, std::move(normal_parameters));
    }
    else if (fluid_pressure_type == "Phreatic_Surface") {
        MakeProcessForPhreaticSurfaceFluidPressureType(rProcessSettings, std::move(normal_parameters));
    }
    else {
        KRATOS_ERROR << "Unknown fluid_pressure_type: " << fluid_pressure_type << std::endl;
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForUniformFluidPressureType(const Parameters& rProcessSettings,
                                                                               Parameters normal_parameters)
{
    normal_parameters.AddValue("value", rProcessSettings["value"][0]);

    if (FunctonsForParameters::HasTableAttached(rProcessSettings,0)) {
        normal_parameters.AddValue("table", rProcessSettings["table"][0]);
        mProcesses.push_back(std::make_unique <ApplyComponentTableProcess>(mrModelPart, normal_parameters));
    } else {
        mProcesses.push_back(std::make_unique <ApplyConstantScalarValueProcess>(mrModelPart, normal_parameters));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForHydrostaticFluidPressureType(const Parameters& rProcessSettings,
                                                                                    Parameters normal_parameters)
{
    normal_parameters.AddValue("gravity_direction", rProcessSettings["gravity_direction"]);
    normal_parameters.AddValue("reference_coordinate", rProcessSettings["reference_coordinate"]);
    normal_parameters.AddValue("specific_weight", rProcessSettings["specific_weight"]);
    if (FunctonsForParameters::HasTableAttached(rProcessSettings, 0)) {
        normal_parameters.AddValue("table", rProcessSettings["table"][0]);
        mProcesses.push_back(std::make_unique <ApplyBoundaryHydrostaticPressureTableProcess>(mrModelPart, normal_parameters));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryHydrostaticPressureProcess>(mrModelPart, normal_parameters));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForPhreaticLineFluidPressureType(const Parameters& rProcessSettings,
                                                                                    Parameters normal_parameters)
{
    normal_parameters.AddValue("gravity_direction", rProcessSettings["gravity_direction"]);
    normal_parameters.AddValue("out_of_plane_direction", rProcessSettings["out_of_plane_direction"]);
    normal_parameters.AddValue("first_reference_coordinate", rProcessSettings["first_reference_coordinate"]);
    normal_parameters.AddValue("second_reference_coordinate", rProcessSettings["second_reference_coordinate"]);
    normal_parameters.AddValue("specific_weight", rProcessSettings["specific_weight"]);
    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        normal_parameters.AddValue("table", rProcessSettings["table"]);
        mProcesses.push_back(std::make_unique <ApplyBoundaryPhreaticLinePressureTableProcess>(mrModelPart, normal_parameters));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryPhreaticLinePressureProcess>(mrModelPart, normal_parameters));
    }
}

void ApplyNormalLoadTableProcess::MakeProcessForPhreaticSurfaceFluidPressureType(const Parameters& rProcessSettings,
                                                                                       Parameters normal_parameters) {
    normal_parameters.AddValue("gravity_direction", rProcessSettings["gravity_direction"]);
    normal_parameters.AddValue("first_reference_coordinate", rProcessSettings["first_reference_coordinate"]);
    normal_parameters.AddValue("second_reference_coordinate", rProcessSettings["second_reference_coordinate"]);
    normal_parameters.AddValue("third_reference_coordinate", rProcessSettings["third_reference_coordinate"]);
    normal_parameters.AddValue("specific_weight", rProcessSettings["specific_weight"]);
    if (FunctonsForParameters::HasTableAttached(rProcessSettings)) {
        normal_parameters.AddValue("table", rProcessSettings["table"]);
        mProcesses.push_back(std::make_unique <ApplyBoundaryPhreaticSurfacePressureTableProcess>(mrModelPart, normal_parameters));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantBoundaryPhreaticSurfacePressureProcess>(mrModelPart, normal_parameters));
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
        mProcesses.push_back(std::make_unique <ApplyComponentTableProcess>(mrModelPart, tangential_params));
    }
    else {
        mProcesses.push_back(std::make_unique <ApplyConstantScalarValueProcess>(mrModelPart, tangential_params));
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
