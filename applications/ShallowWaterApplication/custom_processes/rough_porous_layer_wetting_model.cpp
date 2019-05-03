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
    
}

}  // namespace Kratos.


