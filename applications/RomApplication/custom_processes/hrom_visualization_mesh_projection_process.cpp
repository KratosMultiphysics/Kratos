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

    // Set the file name of the input ROM settings
    mRomSettingsFilename = rParameters["rom_settings_filename"].GetString();
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

void HRomVisualizationMeshProjectionProcess::ExecuteInitialize()
{
    // Parse the ROM input settings
    std::ifstream input_stream(mRomSettingsFilename + ".json", std::ifstream::in);
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

    // Set the visualization model part nodal data
    // First, set the origin model part variables container in the visualization model part
    mrVisualizationModelPart.SetNodalSolutionStepVariablesList(mrHRomModelPart.pGetNodalSolutionStepVariablesList());
    // Secondly, update the already existent nodes nodal data in the visualization model part
    mrVisualizationModelPart.SetNodalSolutionStepVariablesList();

    // Set the origin model part buffer size in the visualization one
    mrVisualizationModelPart.SetBufferSize(mrHRomModelPart.GetBufferSize());

    // Set the origin model part ProcessInfo in the visualization model part
    mrVisualizationModelPart.SetProcessInfo(mrHRomModelPart.pGetProcessInfo());

    // Add DOFs to enable fixity
    VariableUtils::AddDofsList(r_rom_var_names, mrVisualizationModelPart);

    // Get and set the ROM basis from ROM settings
    const auto& r_nodal_modes = rom_parameters["nodal_modes"];
    const SizeType n_modes = rom_parameters["rom_settings"]["number_of_rom_dofs"].GetInt();
    Matrix aux_matrix(n_nodal_dofs, n_modes);
    block_for_each(mrVisualizationModelPart.Nodes(), aux_matrix, [&](NodeType& rNode, Matrix& rTLSMatrix){
        // Set the ROM basis matrix from the ROM settings data
        const auto str_id = std::to_string(rNode.Id());
        const auto& r_modes = r_nodal_modes[str_id];
        for (IndexType i = 0; i < n_nodal_dofs; ++i) {
            const auto& r_dof_modes = r_modes[i].GetVector();
            for (IndexType j = 0; j < n_modes; ++j) {
                rTLSMatrix(i,j) = r_dof_modes[j];
            }
        }

        // Save the basis in the ROM_BASIS nodal variable
        rNode.SetValue(ROM_BASIS, rTLSMatrix);
    });

    // Allocate the
}

void HRomVisualizationMeshProjectionProcess::ExecuteBeforeOutputStep()
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
