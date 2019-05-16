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
#include "shallow_water_application_variables.h"
#include "shallow_water_utilities.h"


namespace Kratos
{

void ShallowWaterUtilities::ComputeFreeSurfaceElevation(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = it_node->FastGetSolutionStepValue(HEIGHT) - it_node->FastGetSolutionStepValue(BATHYMETRY);
    }
}

void ShallowWaterUtilities::ComputeHeightFromFreeSurface(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(HEIGHT) = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) + it_node->FastGetSolutionStepValue(BATHYMETRY);
    }
}

void ShallowWaterUtilities::ComputeVelocity(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(VELOCITY) = it_node->FastGetSolutionStepValue(MOMENTUM) / it_node->FastGetSolutionStepValue(HEIGHT);
    }
}

void ShallowWaterUtilities::ComputeMomentum(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(MOMENTUM) = it_node->FastGetSolutionStepValue(VELOCITY) * it_node->FastGetSolutionStepValue(HEIGHT);
    }
}

void ShallowWaterUtilities::FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(rDestinationVariable) = -it_node->FastGetSolutionStepValue(rOriginVariable);
    }
}

void ShallowWaterUtilities::IdentifySolidBoundary(ModelPart& rSkinModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rSkinModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rSkinModelPart.NodesBegin() + i;
        if (it_node->FastGetSolutionStepValue(TOPOGRAPHY) < SeaWaterLevel)
        {
            it_node->Set(SolidBoundaryFlag, true);
        }
        else
        {
            auto topography_gradient = it_node->FastGetSolutionStepValue(TOPOGRAPHY_GRADIENT);
            auto normal = it_node->FastGetSolutionStepValue(NORMAL);
            double sign = inner_prod(normal, topography_gradient);
            // NOTE: Normal is positive outwards
            // An inwards flow (negative topography_gradient) will produce a positive sign: a SOLID boundary
            if (sign < 0.0)
            {
                it_node->Set(SolidBoundaryFlag, false);
            }
            else
            {
                it_node->Set(SolidBoundaryFlag, true);
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rSkinModelPart.NumberOfConditions()); ++i)
    {
        auto it_cond = rSkinModelPart.ConditionsBegin() + i;
        bool is_solid = true;
        for (auto& node : it_cond->GetGeometry())
        {
            if (node.IsNot(SolidBoundaryFlag)) {
                is_solid = false;
            }
        }
        it_cond->Set(SolidBoundaryFlag, is_solid);
    }
}

}  // namespace Kratos.
