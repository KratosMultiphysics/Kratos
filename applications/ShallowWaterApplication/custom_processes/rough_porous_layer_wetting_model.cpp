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

void RoughPorousLayerWettingModel::Execute()
{
    // Definition of the first node iterator
    auto it_node_begin = mrModelPart.NodesBegin();

    // Execution of the rough porous layer method
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = it_node_begin + i;
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
        else if (topography < free_surface) // Transition domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning;
        }
        else // Dry domain
        {
            bathymetry = -free_surface + mLayerThickness;
            equivalent_manning = manning;
        }
    }
}

void RoughPorousLayerWettingModel::ExecuteFinalizeSolutionStep()
{
    // Definition of the first node iterator
    auto it_elements_begin = mrModelPart.ElementsBegin();

    // Identification of the dry and wet elements
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    {
        auto it_elem = it_elements_begin + i;

        bool wet_element = false;
        for(auto& node : it_elem->GetGeometry())
        {
            const double free_surface = node.FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
            const double topography = node.FastGetSolutionStepValue(TOPOGRAPHY);
            if (free_surface > topography) {
                wet_element = true;  // It means there is almost a wet node
            }
        }

        it_elem->Set(FLUID, wet_element);
    }
}

}  // namespace Kratos
