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
#include "hrom_visualization_mesh_modeler.h"
#include "rom_application_variables.h"

namespace Kratos
{

HRomVisualizationMeshModeler::HRomVisualizationMeshModeler(
    Model& rModel,
    Parameters rParameters)
    : Modeler(rModel, rParameters)
    , mpHRomModelPart(&(rModel.GetModelPart(rParameters["hrom_model_part_name"].GetString())))
    , mpVisualizationModelPart(&(rModel.GetModelPart(rParameters["hrom_visualization_model_part_name"].GetString())))
{
    // Check default settings
    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Set the file name of the input ROM settings
    mRomSettingsFilename = rParameters["rom_settings_filename"].GetString();
}

const Parameters HRomVisualizationMeshModeler::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "echo_level" : 0,
        "hrom_model_part_name" : "",
        "hrom_visualization_model_part_name" : "",
        "rom_settings_filename" : "RomParameters"
    })" );

    return default_parameters;
}

void HRomVisualizationMeshModeler::SetupModelPart()
{
    // Parse the ROM input settings
    std::ifstream input_stream(mRomSettingsFilename + ".json", std::ifstream::in);
    Parameters rom_parameters(input_stream);

    // Create an array with pointers to the ROM variables from the provided names
    // Note that these are assumed to be provided in the same order that was used to create the basis
    const auto& r_rom_var_names = rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray();
    const SizeType n_nodal_dofs = r_rom_var_names.size();
    mRomVariablesList.reserve(n_nodal_dofs);
    for (const auto& r_var_name : r_rom_var_names) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_var_name))
            << "Provided variable '" << r_var_name << "' is not in KratosComponents. Note that array-like variables need to be provided componentwise (e.g. DISPLACEMENT_X, DISPLACEMENT_Y)." << std::endl;
        mRomVariablesList.push_back(&(KratosComponents<Variable<double>>::Get(r_var_name)));
    }

    // Set the visualization model part nodal data
    // First, set the origin model part variables container in the visualization model part
    mpVisualizationModelPart->SetNodalSolutionStepVariablesList(mpHRomModelPart->pGetNodalSolutionStepVariablesList());
    // Secondly, update the already existent nodes nodal data in the visualization model part
    mpVisualizationModelPart->SetNodalSolutionStepVariablesList();

    // Set the origin model part buffer size in the visualization one
    mpVisualizationModelPart->SetBufferSize(mpHRomModelPart->GetBufferSize());

    // Set the origin model part ProcessInfo in the visualization model part
    // Note that we do so also in the visualization submodel parts as this might be required for the visualization mesh BCs
    mpVisualizationModelPart->SetProcessInfo(mpHRomModelPart->pGetProcessInfo());
    for (auto& r_vis_sub_mp : mpVisualizationModelPart->SubModelParts()) {
        r_vis_sub_mp.SetProcessInfo(mpHRomModelPart->pGetProcessInfo());
    }

    // Add DOFs to enable fixity
    VariableUtils::AddDofsList(r_rom_var_names, *mpVisualizationModelPart);

    // Get and set the ROM basis from ROM settings
    const auto& r_nodal_modes = rom_parameters["nodal_modes"];
    const SizeType n_modes = rom_parameters["rom_settings"]["number_of_rom_dofs"].GetInt();
    Matrix aux_matrix(n_nodal_dofs, n_modes);
    block_for_each(mpVisualizationModelPart->Nodes(), aux_matrix, [&](NodeType& rNode, Matrix& rTLSMatrix){
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
}

/* Protected functions ****************************************************/


/* Private functions ****************************************************/

};  // namespace Kratos.
