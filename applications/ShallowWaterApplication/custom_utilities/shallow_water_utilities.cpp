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
#include "shallow_water_utilities.h"


namespace Kratos
{

void ShallowWaterUtilities::ComputeFreeSurfaceElevation(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = rNode.FastGetSolutionStepValue(HEIGHT) + rNode.FastGetSolutionStepValue(TOPOGRAPHY);
    });
}

void ShallowWaterUtilities::ComputeHeightFromFreeSurface(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(HEIGHT) = rNode.FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - rNode.FastGetSolutionStepValue(TOPOGRAPHY);
    });
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
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        noalias(rNode.FastGetSolutionStepValue(MOMENTUM)) = rNode.FastGetSolutionStepValue(VELOCITY) * rNode.FastGetSolutionStepValue(HEIGHT);
    });
}

void ShallowWaterUtilities::ComputeEnergy(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double height = rNode.FastGetSolutionStepValue(HEIGHT);
        const double velocity = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));
        rNode.FastGetSolutionStepValue(INTERNAL_ENERGY) = height + 0.5 * std::pow(velocity, 2);
    });
}

void ShallowWaterUtilities::FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(rDestinationVariable) = -rNode.FastGetSolutionStepValue(rOriginVariable);
    });
}

void ShallowWaterUtilities::IdentifySolidBoundary(ModelPart& rSkinModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag)
{
    block_for_each(rSkinModelPart.Nodes(), [&](NodeType& rNode){
        if (rNode.FastGetSolutionStepValue(TOPOGRAPHY) < SeaWaterLevel)
        {
            rNode.Set(SolidBoundaryFlag, true);
        }
        else
        {
            auto topography_gradient = rNode.GetValue(TOPOGRAPHY_GRADIENT);
            array_1d<double,3> normal = rNode.FastGetSolutionStepValue(NORMAL);
            double sign = inner_prod(normal, topography_gradient);
            // NOTE: Normal is positive outwards
            // NOTE: The flowstream is opposite to the topography gradient
            // An inwards flow will produce a positive sign: a SOLID boundary
            rNode.Set(SolidBoundaryFlag, (sign >= 0.0));
        }
    });

    block_for_each(rSkinModelPart.Conditions(), [&](Condition& rCondition){
        bool is_solid = true;
        for (auto& node : rCondition.GetGeometry())
        {
            if (node.IsNot(SolidBoundaryFlag)) {
                is_solid = false;
            }
        }
        rCondition.Set(SolidBoundaryFlag, is_solid);
    });
}

void ShallowWaterUtilities::IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double RelativeDryHeight)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Set(WetFlag, false);
    });

    IdentifyWetEntities(rModelPart.Elements(), WetFlag, RelativeDryHeight);
    IdentifyWetEntities(rModelPart.Conditions(), WetFlag, RelativeDryHeight);
}

void ShallowWaterUtilities::NormalizeVector(ModelPart& rModelPart, Variable<array_1d<double,3>>& rVariable)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& vector = rNode.FastGetSolutionStepValue(rVariable);
        const auto modulus = norm_2(vector);
        if (modulus > std::numeric_limits<double>::epsilon())
            vector /= modulus;
    });
}

void ShallowWaterUtilities::SetMinimumValue(ModelPart& rModelPart, const Variable<double>& rVariable, double MinValue)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        double& value = rNode.FastGetSolutionStepValue(rVariable);
        value = std::max(value, MinValue);
    });
}

void ShallowWaterUtilities::SetMeshZCoordinateToZero(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Z() = 0.0;
    });
}

void ShallowWaterUtilities::SetMeshZ0CoordinateToZero(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Z0() = 0.0;
    });
}

void ShallowWaterUtilities::SetMeshZCoordinate(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Z() = rNode.FastGetSolutionStepValue(rVariable);
    });
}

double ShallowWaterUtilities::InverseHeight(const double Height, const double Epsilon)
{
    const double h4 = std::pow(Height, 4);
    const double epsilon4 = std::pow(Epsilon, 4);
    return std::sqrt(2) * std::max(Height, .0) / std::sqrt(h4 + std::max(h4, epsilon4));
}

double ShallowWaterUtilities::WetFraction(const double Height, const double Epsilon)
{
    return Height * InverseHeight(Height, Epsilon);
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

template<>
double ShallowWaterUtilities::GetValue<true>(NodeType& rNode, const Variable<double>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
double ShallowWaterUtilities::GetValue<false>(NodeType& rNode, const Variable<double>& rVariable)
{
    return rNode.GetValue(rVariable);
}

template<>
array_1d<double,3> ShallowWaterUtilities::EvaluateHydrostaticForce<ModelPart::ConditionsContainerType>(
    const double Density,
    const double Gravity,
    const double Height,
    const double Area,
    const array_1d<double,3>& rNormal)
{
    return 0.5 * Density * Gravity * Height * Height * Area * rNormal;
}

template<>
array_1d<double,3> ShallowWaterUtilities::EvaluateHydrostaticForce<ModelPart::ElementsContainerType>(
    const double Density,
    const double Gravity,
    const double Height,
    const double Area,
    const array_1d<double,3>& rNormal)
{
    return -Density * Gravity * Height * Area * rNormal;
}


}  // namespace Kratos.
