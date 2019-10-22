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
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = it_node->FastGetSolutionStepValue(HEIGHT) - it_node->FastGetSolutionStepValue(BATHYMETRY);
    }
}

void ShallowWaterUtilities::ComputeHeightFromFreeSurface(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(HEIGHT) = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) + it_node->FastGetSolutionStepValue(BATHYMETRY);
    }
}

void ShallowWaterUtilities::ComputeVelocity(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(VELOCITY) = it_node->FastGetSolutionStepValue(MOMENTUM) / it_node->FastGetSolutionStepValue(HEIGHT);
    }
}

void ShallowWaterUtilities::ComputeMomentum(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(MOMENTUM) = it_node->FastGetSolutionStepValue(VELOCITY) * it_node->FastGetSolutionStepValue(HEIGHT);
    }
}

void ShallowWaterUtilities::UpdatePrimitiveVariables(ModelPart& rModelPart)
{
    double dry_height = rModelPart.GetProcessInfo()[DRY_HEIGHT];
    UpdatePrimitiveVariables(rModelPart, dry_height);
}

void ShallowWaterUtilities::UpdatePrimitiveVariables(ModelPart& rModelPart, double Epsilon)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        const double height = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - it_node->FastGetSolutionStepValue(TOPOGRAPHY);
        it_node->FastGetSolutionStepValue(HEIGHT) = height;
        it_node->FastGetSolutionStepValue(VELOCITY) = it_node->FastGetSolutionStepValue(MOMENTUM) / std::max(std::abs(height), Epsilon);
    }
}

void ShallowWaterUtilities::ComputeAccelerations(ModelPart& rModelPart)
{
    double dt_inv = rModelPart.GetProcessInfo()[DELTA_TIME];

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;

        // Free suface derivative or vertical velocity
        auto delta_surface = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION,1);
        it_node->FastGetSolutionStepValue(VELOCITY_Z) = dt_inv * delta_surface;

        // Acceleration
        auto delta_vel = it_node->FastGetSolutionStepValue(VELOCITY) - it_node->FastGetSolutionStepValue(VELOCITY,1);
        it_node->SetValue(ACCELERATION, dt_inv * delta_vel);
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
            // NOTE: The flowstream is opposite to the topography gradient
            // An inwards flow will produce a positive sign: a SOLID boundary
            it_node->Set(SolidBoundaryFlag, (sign >= 0.0));
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

void ShallowWaterUtilities::IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double Thickness)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        const double height = it_node->FastGetSolutionStepValue(HEIGHT);
        it_node->Set(WetFlag, (height > Thickness));
    }

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfElements()); ++i)
    {
        auto it_elem = rModelPart.ElementsBegin() + i;

        bool wet_element = false;
        for(auto& node : it_elem->GetGeometry())
        {
            if (node.Is(WetFlag)) {
                wet_element = true;  // It means there is almost a wet node
                break;
            }
        }

        it_elem->Set(WetFlag, wet_element);
    }
}

void ShallowWaterUtilities::ComputeVisualizationWaterHeight(ModelPart& rModelPart, Flags WetFlag, double SeaWaterLevel)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        if (it_node->Is(WetFlag)) {
            if (it_node->FastGetSolutionStepValue(TOPOGRAPHY) > SeaWaterLevel) {
                it_node->SetValue(WATER_HEIGHT, it_node->FastGetSolutionStepValue(HEIGHT));
            }
            else {
                it_node->SetValue(WATER_HEIGHT, it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - SeaWaterLevel);
            }
        }
        else {
            // This is the undefined value for GiD
            it_node->SetValue(WATER_HEIGHT, std::numeric_limits<float>::lowest());
        }
    }
}

void ShallowWaterUtilities::ComputeVisualizationWaterSurface(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->SetValue(WATER_SURFACE_Z, it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION));
    }
}

}  // namespace Kratos.
