//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "rough_porous_layer_wetting_model.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

RoughPorousLayerWettingModel::RoughPorousLayerWettingModel(ModelPart& rModelPart, Parameters ThisParameters)
 : mrModelPart(rModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_name"        : "rough_porous_layer",
        "layer_thickness"   : 1.0,
        "roughness_factor"  : 0.1
    })");
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    mLayerThickness = ThisParameters["layer_thickness"].GetDouble();
    mRoughnessFactor = ThisParameters["roughness_factor"].GetDouble();
}

RoughPorousLayerWettingModel::RoughPorousLayerWettingModel(ModelPart& rModelPart, double LayerThickness, double RoughnessFactor)
 : mrModelPart(rModelPart)
{
    mLayerThickness = LayerThickness;
    mRoughnessFactor = RoughnessFactor;
}

void RoughPorousLayerWettingModel::ExecuteInitializeSolutionStep()
{
    // Definition of the first node iterator
    auto it_node_begin = mrModelPart.NodesBegin();

    // Execution of the rough porous layer method
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = it_node_begin + i;
        const double topography = it_node->FastGetSolutionStepValue(TOPOGRAPHY);
        const double manning = it_node->FastGetSolutionStepValue(MANNING);
        double& height = it_node->FastGetSolutionStepValue(HEIGHT);
        double& bathymetry = it_node->FastGetSolutionStepValue(BATHYMETRY);
        double& equivalent_manning = it_node->FastGetSolutionStepValue(EQUIVALENT_MANNING);

        // 1. Getting the free surface from the last iteration
        const double free_surface = height - bathymetry;

        // 2. Moving the bathymetry
        double equivalent_topography = topography + mLayerThickness;
        if (equivalent_topography < free_surface) // Wet domain
        {
            bathymetry = -topography;
            equivalent_manning = manning;
        }
        else if (topography < free_surface) // Transition domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning + mRoughnessFactor * (mLayerThickness + topography - free_surface) / mLayerThickness;
        }
        else // Dry domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning + mRoughnessFactor;
        }

        // 3. Updating the water height after moving the bathymetry
        height = free_surface + bathymetry;
    }
}

void RoughPorousLayerWettingModel::ExecuteFinalizeSolutionStep()
{
    ShallowWaterUtilities().IdentifyWetDomain(mrModelPart, FLUID, mLayerThickness);
}

}  // namespace Kratos
