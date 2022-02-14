//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "hrom_visualization_mesh_projection_process.h"
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"

namespace Kratos
{

HRomVisualizationMeshProjectionProcess::HRomVisualizationMeshProjectionProcess(
    Model& rModel,
    Parameters rParameters)
    : Process()
    , mrHRomModelPart(rModel.GetModelPart(rParameters["hrom_model_part_name"].GetString()))
    , mrVisualizationModelPart(rModel.GetModelPart(rParameters["hrom_visualization_model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

    // Parse the ROM input settings
    std::ifstream input_stream(rParameters["rom_settings_filename"].GetString() + ".json", std::ifstream::in);
    Parameters rom_parameters(input_stream);

    // Create an array with pointers to the ROM variables from the provided names
    // Note that these are assumed to be provided in the same order that was used to create the basis
    IndexType i_var = 0;
    const auto& r_rom_var_names = rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray();
    const SizeType n_nodal_dofs = r_rom_var_names.size();
    mRomVariablesList.resize(n_nodal_dofs);
    for (const auto& r_var_name : r_rom_var_names) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_var_name))
            << "Provided variable '" << r_var_name << "' is not in KratosComponents. Note that array-like variables need to be provided componentwise (e.g. DISPLACEMENT_X, DISPLACEMENT_Y)." << std::endl;
        mRomVariablesList[i_var++] = &(KratosComponents<Variable<double>>::Get(r_var_name));
    }

    // We sort the ROM variables alphabetically
    // This is required to match the order of the DOF set
    std::sort(
        mRomVariablesList.begin(),
        mRomVariablesList.end(),
        [](const Variable<double>* pVarA, const Variable<double>* pVarB){return pVarA->Name() < pVarB->Name();});
}

void HRomVisualizationMeshProjectionProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    Parameters default_parameters( R"(
    {
        "hrom_model_part_name" : "",
        "hrom_visualization_model_part_name" : "",
        "rom_settings_filename" : "RomParameters"
    })" );
    rParameters.ValidateAndAssignDefaults(default_parameters);
}

void HRomVisualizationMeshProjectionProcess::ExecuteFinalizeSolutionStep()
{
    // Retrieve the HROM solution from the HROM origin model part
    // Note that we assume this is always saved in the corresponding root model part
    const auto& r_rom_sol_inc = mrHRomModelPart.GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);
    mrVisualizationModelPart.GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT) = r_rom_sol_inc;

    // Perform the HROM solution projection onto the visualization mesh
    // Note that this method assumes that BCs have been already applied
    RomAuxiliaryUtilities::ProjectRomSolutionIncrementToNodes(mRomVariablesList, mrVisualizationModelPart);
}

/* Protected functions ****************************************************/


/* Private functions ****************************************************/

};  // namespace Kratos.
