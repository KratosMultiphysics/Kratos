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
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
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
        it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = it_node->FastGetSolutionStepValue(HEIGHT) + it_node->FastGetSolutionStepValue(TOPOGRAPHY);
    }
}

void ShallowWaterUtilities::ComputeHeightFromFreeSurface(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->FastGetSolutionStepValue(HEIGHT) = it_node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - it_node->FastGetSolutionStepValue(TOPOGRAPHY);
    }
}

void ShallowWaterUtilities::ComputeVelocity(ModelPart& rModelPart, bool PerformProjection)
{
    if (PerformProjection) {
        ComputeSmoothVelocity(rModelPart);
    } else {
        const double rel_dry_h = rModelPart.GetProcessInfo()[RELATIVE_DRY_HEIGHT];
        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            const double h = r_node.FastGetSolutionStepValue(HEIGHT);
            const double inv_h = InverseHeight(h, rel_dry_h * r_node.GetValue(NODAL_H));
            r_node.FastGetSolutionStepValue(VELOCITY) = inv_h * r_node.FastGetSolutionStepValue(MOMENTUM);
        });
    }
}

void ShallowWaterUtilities::ComputeSmoothVelocity(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
        r_node.FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
        r_node.SetValue(INTEGRATION_WEIGHT, 0.0);
    });
    Matrix mass_matrix;
    const double rel_dry_h = rModelPart.GetProcessInfo()[RELATIVE_DRY_HEIGHT];
    block_for_each(rModelPart.Elements(), mass_matrix, [&](Element& r_element, Matrix& r_local_mass_matrix){
        auto& r_geom = r_element.GetGeometry();
        const size_t num_nodes = r_geom.size();
        double height = 0.0;
        Vector nodal_discharge_x(num_nodes);
        Vector nodal_discharge_y(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i) {
            height += r_geom[i].FastGetSolutionStepValue(HEIGHT);
            nodal_discharge_x[i] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_X);
            nodal_discharge_y[i] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_Y);
        }
        height /= num_nodes;
        CalculateMassMatrix(r_local_mass_matrix, r_geom);
        r_local_mass_matrix *= InverseHeight(height, rel_dry_h * r_geom.Length());
        Vector nodal_velocity_x(num_nodes);
        Vector nodal_velocity_y(num_nodes);
        nodal_velocity_x = num_nodes * prod(r_local_mass_matrix, nodal_discharge_x);
        nodal_velocity_y = num_nodes * prod(r_local_mass_matrix, nodal_discharge_y);
        for (unsigned int i = 0; i < num_nodes; ++i)
        {
            r_geom[i].SetLock();
            r_geom[i].FastGetSolutionStepValue(VELOCITY_X) += nodal_velocity_x[i];
            r_geom[i].FastGetSolutionStepValue(VELOCITY_Y) += nodal_velocity_y[i];
            r_geom[i].GetValue(INTEGRATION_WEIGHT) += 1.0;
            r_geom[i].UnSetLock();
        }
    });
    block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
        r_node.FastGetSolutionStepValue(VELOCITY) /= r_node.GetValue(INTEGRATION_WEIGHT);
    });
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

void ShallowWaterUtilities::ComputeEnergy(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        const double height = it_node->FastGetSolutionStepValue(HEIGHT);
        const double velocity = norm_2(it_node->FastGetSolutionStepValue(VELOCITY));
        it_node->FastGetSolutionStepValue(INTERNAL_ENERGY) = height + 0.5 * std::pow(velocity, 2);
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
            auto topography_gradient = it_node->GetValue(TOPOGRAPHY_GRADIENT);
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
            height = 0.5 * Thickness;
            it_node->FastGetSolutionStepValue(MOMENTUM) = ZeroVector(3);
        }
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

void ShallowWaterUtilities::SetMeshZCoordinateToZero(ModelPart& rModelPart)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->Z() = 0.0;
    }
}

void ShallowWaterUtilities::SetMeshZCoordinate(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = rModelPart.NodesBegin() + i;
        it_node->Z() = it_node->FastGetSolutionStepValue(rVariable);
    }
}

double ShallowWaterUtilities::InverseHeight(const double Height, const double Epsilon)
{
    const double h4 = std::pow(Height, 4);
    const double epsilon4 = std::pow(Epsilon, 4);
    return std::sqrt(2) * std::max(Height, .0) / std::sqrt(h4 + std::max(h4, epsilon4));
}

void ShallowWaterUtilities::CalculateMassMatrix(Matrix& rMassMatrix, const GeometryType& rGeometry)
{
    const size_t num_nodes = rGeometry.size();
    if (rMassMatrix.size1() != num_nodes) {
        rMassMatrix.resize(num_nodes, num_nodes, false);
    }
    if (num_nodes == 2)
    {
        double one_sixth = 1.0 / 6.0;
        rMassMatrix(0,0) = 2.0 * one_sixth;
        rMassMatrix(0,1) = 1.0 * one_sixth;
        rMassMatrix(1,0) = 1.0 * one_sixth;
        rMassMatrix(1,1) = 2.0 * one_sixth;
    }
    else if (num_nodes == 3)
    {
        double one_twelve = 1.0 / 12.0;
        rMassMatrix(0,0) = 2.0 * one_twelve;
        rMassMatrix(0,1) = 1.0 * one_twelve;
        rMassMatrix(0,2) = 1.0 * one_twelve;
        rMassMatrix(1,0) = 1.0 * one_twelve;
        rMassMatrix(1,1) = 2.0 * one_twelve;
        rMassMatrix(1,2) = 1.0 * one_twelve;
        rMassMatrix(2,0) = 1.0 * one_twelve;
        rMassMatrix(2,1) = 1.0 * one_twelve;
        rMassMatrix(2,2) = 2.0 * one_twelve;
    }
    else if (num_nodes == 4)
    {
        double one_thirty_sixth = 1.0 / 36.0;
        rMassMatrix(0,0) = 4 * one_thirty_sixth;
        rMassMatrix(0,1) = 2 * one_thirty_sixth;
        rMassMatrix(0,2) = 1 * one_thirty_sixth;
        rMassMatrix(0,3) = 2 * one_thirty_sixth;
        rMassMatrix(1,0) = 2 * one_thirty_sixth;
        rMassMatrix(1,1) = 4 * one_thirty_sixth;
        rMassMatrix(1,2) = 2 * one_thirty_sixth;
        rMassMatrix(1,3) = 1 * one_thirty_sixth;
        rMassMatrix(2,0) = 1 * one_thirty_sixth;
        rMassMatrix(2,1) = 2 * one_thirty_sixth;
        rMassMatrix(2,2) = 4 * one_thirty_sixth;
        rMassMatrix(2,3) = 2 * one_thirty_sixth;
        rMassMatrix(3,0) = 2 * one_thirty_sixth;
        rMassMatrix(3,1) = 1 * one_thirty_sixth;
        rMassMatrix(3,2) = 2 * one_thirty_sixth;
        rMassMatrix(3,3) = 4 * one_thirty_sixth;
    }
    else
    {
        KRATOS_ERROR << "ShallowWaterUtilities::MassMatrix. Method implemented for lines, triangles and quadrilaterals" << std::endl;
    }
}

}  // namespace Kratos.
