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
