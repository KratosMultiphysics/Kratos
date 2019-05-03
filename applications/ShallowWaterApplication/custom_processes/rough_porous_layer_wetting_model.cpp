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


namespace Kratos
{

RoughPorousLayerWettingModel::RoughPorousLayerWettingModel(ModelPart& rModelPart, Parameters ThisParameters)
 : mrModelPart(rModelPart)
{
    mLayerThickness = ThisParameters["layer_thickness"].GetDouble();
    mRoughnessFactor = ThisParameters["roughness_factor"].GetDouble();
}

RoughPorousLayerWettingModel::RoughPorousLayerWettingModel(ModelPart& rModelPart, double LayerThickness, double RoughnessFactor)
 : mrModelPart(rModelPart)
{
    mLayerThickness = LayerThickness;
    mRoughnessFactor = RoughnessFactor;
}

void RoughPorousLayerWettingModel::Execute()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = mrModelPart.NodesBegin();
        const double free_surface = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        const double topography = it_node->FastGetSolutionStepValue(TOPOGRAPHY);
        const double manning = it_node->FastGetSolutionStepValue(MANNING);
        double& bathymetry = it_node->FastGetSolutionStepValue(BATHYMETRY);
        double& equivalent_manning = it_node->FastGetSolutionStepValue(EQUIVALENT_MANNING);

        double equivalent_topography = topography + mLayerThickness;
        if (equivalent_topography < free_surface) // Wet domain
        {
            bathymetry = -topography;
            equivalent_manning = manning;
        }
        else if ((topography < free_surface) & (free_surface <= equivalent_topography)) // Transition domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning;
        }
        else if (free_surface <= topography) // Dry domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning;
        }
    }
}

}  // namespace Kratos.
