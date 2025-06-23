// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "convection_diffusion_application_variables.h"
#include "convection_diffusion_testing_utilities.h"


namespace Kratos::Testing
{

void ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(ModelPart &rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Set convection diffusion settings
    auto p_conv_dff_set = Kratos::make_shared<ConvectionDiffusionSettings>();
    p_conv_dff_set->SetDensityVariable(DENSITY);
    p_conv_dff_set->SetDiffusionVariable(CONDUCTIVITY);
    p_conv_dff_set->SetUnknownVariable(TEMPERATURE);
    p_conv_dff_set->SetVolumeSourceVariable(HEAT_FLUX);
    p_conv_dff_set->SetSurfaceSourceVariable(FACE_HEAT_FLUX);
    p_conv_dff_set->SetProjectionVariable(PROJECTED_SCALAR1);
    p_conv_dff_set->SetConvectionVariable(CONVECTION_VELOCITY);
    p_conv_dff_set->SetMeshVelocityVariable(MESH_VELOCITY);
    p_conv_dff_set->SetVelocityVariable(VELOCITY);
    p_conv_dff_set->SetSpecificHeatVariable(SPECIFIC_HEAT);
    p_conv_dff_set->SetReactionVariable(REACTION_FLUX);
    rModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
    rModelPart.AddNodalSolutionStepVariable(HEAT_FLUX);
    rModelPart.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
    rModelPart.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
    rModelPart.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);

    // Create a fake properties container
    auto p_elem_prop = rModelPart.CreateNewProperties(0);
}

} // namespace Kratos::Testing.
