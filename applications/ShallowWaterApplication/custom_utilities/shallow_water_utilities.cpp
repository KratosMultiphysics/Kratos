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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "shallow_water_utilities.h"
#include "phase_function.h"


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
            const double inv_h = PhaseFunction::InverseHeight(h, rel_dry_h * r_node.GetValue(NODAL_H));
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
        r_local_mass_matrix *= PhaseFunction::InverseHeight(height, rel_dry_h * r_geom.Length());
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

void ShallowWaterUtilities::ComputeLinearizedMomentum(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        noalias(rNode.FastGetSolutionStepValue(MOMENTUM)) = -rNode.FastGetSolutionStepValue(VELOCITY) * rNode.FastGetSolutionStepValue(TOPOGRAPHY);
    });
}

template<bool THistorical>
void ShallowWaterUtilities::ComputeFroude(ModelPart& rModelPart, const double Epsilon)
{
    const double g = rModelPart.GetProcessInfo()[GRAVITY_Z];
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double height = rNode.FastGetSolutionStepValue(HEIGHT);
        const double velocity = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));
        const double inverse_c = std::sqrt(PhaseFunction::InverseHeight(height, Epsilon) / g);
        GetValue<THistorical>(rNode, FROUDE) = velocity * inverse_c;
    });
}

template<bool THistorical>
void ShallowWaterUtilities::ComputeEnergy(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double height = rNode.FastGetSolutionStepValue(HEIGHT);
        const double velocity = norm_2(rNode.FastGetSolutionStepValue(VELOCITY));
        GetValue<THistorical>(rNode, INTERNAL_ENERGY) = height + 0.5 * std::pow(velocity, 2);
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

void ShallowWaterUtilities::FlagWetElements(ModelPart& rModelPart, Flags WetFlag, double RelativeDryHeight)
{
    if (RelativeDryHeight < 0.0) {
        RelativeDryHeight = rModelPart.GetProcessInfo()[RELATIVE_DRY_HEIGHT];
    }
    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const auto& r_geom = rElement.GetGeometry();
        const bool is_wet = IsWet(r_geom, RelativeDryHeight);
        rElement.Set(WetFlag, is_wet);
    });
}

void ShallowWaterUtilities::ExtrapolateElementalFlagToNodes(ModelPart& rModelPart, Flags Flag)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Set(Flag, false);
    });
    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const auto& r_geom = rElement.GetGeometry();
        for (auto& r_node : r_geom)
        {
            if (rElement.Is(Flag))
            {
                if (r_node.IsNot(Flag))
                {
                    r_node.SetLock();
                    r_node.Set(Flag);
                    r_node.UnSetLock();
                }
            }
        }
    });
}

void ShallowWaterUtilities::NormalizeVector(ModelPart& rModelPart, const Variable<array_1d<double,3>>& rVariable)
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

void ShallowWaterUtilities::OffsetMeshZCoordinate(ModelPart& rModelPart, const double Increment)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.Z() += Increment;
    });
}

void ShallowWaterUtilities::SwapYZCoordinates(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
        std::swap(rNode.Y(), rNode.Z());
    });
}

void ShallowWaterUtilities::SwapY0Z0Coordinates(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
        std::swap(rNode.Y0(), rNode.Z0());
    });
}

void ShallowWaterUtilities::StoreNonHistoricalGiDNoDataIfDry(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    const double relative_dry_height = rModelPart.GetProcessInfo()[RELATIVE_DRY_HEIGHT];
    const double length = rModelPart.ElementsBegin()->GetGeometry().Length();
    const double dry_height = relative_dry_height * length;
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double height = rNode.FastGetSolutionStepValue(HEIGHT);
        const bool is_wet = IsWet(height, dry_height);
        const double value = (is_wet) ? rNode.FastGetSolutionStepValue(rVariable) : std::numeric_limits<float>::lowest();
        rNode.SetValue(rVariable, value);
    });
}

template<bool THistorical>
double ShallowWaterUtilities::ComputeL2Norm(ModelPart& rModelPart, const Variable<double>& rVariable)
{
    double l2_norm = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](Element& rElem){
        double partial_l2_norm = 0.0;
        for (auto& r_node : rElem.GetGeometry()) {
            partial_l2_norm += std::pow(GetValue<THistorical>(r_node, rVariable), 2);
        }
        partial_l2_norm *= rElem.GetGeometry().Area();
        partial_l2_norm /= rElem.GetGeometry().size();
        return partial_l2_norm;
    });
    return std::sqrt(l2_norm);
}

template<bool THistorical>
double ShallowWaterUtilities::ComputeL2NormAABB(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    Point& rLow,
    Point& rHigh)
{
    double l2_norm = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](Element& rElem){
        double partial_l2_norm = 0.0;
        if (rElem.GetGeometry().HasIntersection(rLow, rHigh)) {
            for (auto& r_node : rElem.GetGeometry()) {
                partial_l2_norm += std::pow(GetValue<THistorical>(r_node, rVariable), 2);
            }
            partial_l2_norm *= rElem.GetGeometry().Area();
            partial_l2_norm /= rElem.GetGeometry().size();
        }
        return partial_l2_norm;
    });
    return std::sqrt(l2_norm);
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
double& ShallowWaterUtilities::GetValue<true>(NodeType& rNode, const Variable<double>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
double& ShallowWaterUtilities::GetValue<false>(NodeType& rNode, const Variable<double>& rVariable)
{
    return rNode.GetValue(rVariable);
}


bool ShallowWaterUtilities::IsWet(const GeometryType& rGeometry, const double RelativeDryHeight)
{
    double height = 0.0;
    for (const auto& r_node : rGeometry)
    {
        height += r_node.FastGetSolutionStepValue(HEIGHT);
    }
    height /= rGeometry.size();
    return IsWet(rGeometry, height, RelativeDryHeight);
}

bool ShallowWaterUtilities::IsWet(const GeometryType& rGeometry, const double Height, const double RelativeDryHeight)
{
    const double epsilon = RelativeDryHeight * rGeometry.Length();
    return IsWet(Height, epsilon);
}

bool ShallowWaterUtilities::IsWet(const double Height, const double DryHeight)
{
    const double wet_fraction = PhaseFunction::WetFraction(Height, DryHeight);
    const double threshold = 1.0 - 1e-6;
    return (wet_fraction >= threshold);
}

template KRATOS_API(SHALLOW_WATER_APPLICATION) void ShallowWaterUtilities::ComputeFroude<true>(ModelPart&, const double);
template KRATOS_API(SHALLOW_WATER_APPLICATION) void ShallowWaterUtilities::ComputeFroude<false>(ModelPart&, const double);

template KRATOS_API(SHALLOW_WATER_APPLICATION) void ShallowWaterUtilities::ComputeEnergy<true>(ModelPart&);
template KRATOS_API(SHALLOW_WATER_APPLICATION) void ShallowWaterUtilities::ComputeEnergy<false>(ModelPart&);

template KRATOS_API(SHALLOW_WATER_APPLICATION) double ShallowWaterUtilities::ComputeL2Norm<true>(ModelPart&, const Variable<double>&);
template KRATOS_API(SHALLOW_WATER_APPLICATION) double ShallowWaterUtilities::ComputeL2Norm<false>(ModelPart&, const Variable<double>&);

template KRATOS_API(SHALLOW_WATER_APPLICATION) double ShallowWaterUtilities::ComputeL2NormAABB<true>(ModelPart&, const Variable<double>&, Point&, Point&);
template KRATOS_API(SHALLOW_WATER_APPLICATION) double ShallowWaterUtilities::ComputeL2NormAABB<false>(ModelPart&, const Variable<double>&, Point&, Point&);

}  // namespace Kratos.
