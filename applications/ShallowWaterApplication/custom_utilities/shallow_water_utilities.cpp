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
    const double epsilon = rModelPart.GetProcessInfo()[DRY_HEIGHT];
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        const double height = it_node->FastGetSolutionStepValue(HEIGHT);
        it_node->FastGetSolutionStepValue(VELOCITY) = it_node->FastGetSolutionStepValue(MOMENTUM) / (std::abs(height) + epsilon);
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
        int method = 1;

        auto it_elem = rModelPart.ElementsBegin() + i;
        auto& geom = it_elem->GetGeometry();

        bool is_wet = geom[0].Is(WetFlag);
        bool is_shoreline = false;
        for (size_t j = 1; j < geom.size(); ++j)
        {
            if (geom[j].Is(WetFlag) != is_wet)
                is_shoreline = true;
        }

        if (!is_shoreline)
        {
            it_elem->Set(WetFlag, is_wet);
        }
        else
        {
            if (method == 0) {
                it_elem->Set(WetFlag, false);
            }
            else if (method == 1) {
                it_elem->Set(WetFlag, true);
            }
            else if (method == 2) {
                double height_acc = 0.0;
                for (auto& node : geom)
                {
                    height_acc += node.FastGetSolutionStepValue(VELOCITY_Z);
                }
                it_elem->Set(WetFlag, (height_acc > 0.0));
            }
            else if (method == 3) {
                Geometry<Node<3>>::ShapeFunctionsGradientsType DN_DX(1);
                geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
                array_1d<double,3> height_grad = ZeroVector(3);
                array_1d<double,3> velocity = ZeroVector(3);
                for (size_t j = 0; j < geom.size(); ++j)
                {
                    height_grad[0] += DN_DX[0](j,0) * geom[j].FastGetSolutionStepValue(HEIGHT);
                    height_grad[1] += DN_DX[0](j,1) * geom[j].FastGetSolutionStepValue(HEIGHT);
                    velocity += geom[j].FastGetSolutionStepValue(VELOCITY);
                }
                velocity /= geom.size();

                double run_up = -inner_prod(height_grad, velocity);

                it_elem->Set(WetFlag, (run_up > 0.0));
            }
        }
    }
}

void ShallowWaterUtilities::ResetDryDomain(ModelPart& rModelPart, double Thickness)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        double& height = it_node->FastGetSolutionStepValue(HEIGHT);
        if (height < Thickness)
        {
            height = 0.1 * Thickness;
            it_node->FastGetSolutionStepValue(MOMENTUM) = ZeroVector(3);
        }
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

void ShallowWaterUtilities::NormalizeVector(ModelPart& rModelPart, Variable<array_1d<double,3>>& rVariable)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        auto& vector = it_node->FastGetSolutionStepValue(rVariable);
        const auto modulus = norm_2(vector);
        if (modulus > std::numeric_limits<double>::epsilon())
            vector /= modulus;
    }
}

void ShallowWaterUtilities::SetMinimumValue(ModelPart& rModelPart, const Variable<double>& rVariable, double MinValue)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto& value = (rModelPart.NodesBegin() + i)->FastGetSolutionStepValue(rVariable);
        value = std::max(value, MinValue);
    }
}

}  // namespace Kratos.
