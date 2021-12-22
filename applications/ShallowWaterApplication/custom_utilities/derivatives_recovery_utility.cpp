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
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"
#include "derivatives_recovery_utility.h"
#include "includes/global_pointer_variables.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

void DerivativesRecoveryUtility::CalculateDivergence(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep) = 0.0;
    });

    block_for_each(rModelPart.Elements(), [&](Element& rElement){

        auto& r_geom = rElement.GetGeometry();

        array_1d<double, 3> N;
        BoundedMatrix<double, 3, 2> DN_DX;
        double area;

        GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, area);

        double elemental_div = 0;
        array_1d<double,3> gradients_i = ZeroVector(3);
        for (std::size_t i = 0; i < r_geom.size(); ++i) {
            gradients_i[0] = DN_DX(i,0);
            gradients_i[1] = DN_DX(i,1);
            array_1d<double,3> nodal_value = r_geom[i].FastGetSolutionStepValue(rOriginVariable, BufferStep);
            elemental_div += inner_prod(gradients_i, nodal_value);
        }

        double nodal_area = area / 3.0; // The division by the number of nodes comes from the shape functions
        elemental_div *= nodal_area;

        for (auto& r_node : r_geom) {
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(rDestinationVariable, BufferStep) += elemental_div;
            r_node.UnSetLock();
        }
    });

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep) /= rNode.FastGetSolutionStepValue(NODAL_AREA);
    });
}

void DerivativesRecoveryUtility::CalculateGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep) = ZeroVector(3);
    });

    block_for_each(rModelPart.Elements(), [&](Element& rElement){

        auto& r_geom = rElement.GetGeometry();

        array_1d<double, 3> N;
        BoundedMatrix<double, 3, 2> DN_DX;
        double area;

        GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, area);

        array_1d<double,3> elemental_grad = ZeroVector(3);
        for (std::size_t i = 0; i < r_geom.size(); ++i) {
            double nodal_value = r_geom[i].FastGetSolutionStepValue(rOriginVariable, BufferStep);
            elemental_grad[0] += DN_DX(i,0) * nodal_value;
            elemental_grad[1] += DN_DX(i,1) * nodal_value;
        }

        double nodal_area = area / 3.0; // The division by the number of nodes comes from the shape functions
        elemental_grad *= nodal_area;

        for (auto& r_node : r_geom) {
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(rDestinationVariable, BufferStep) += elemental_grad;
            r_node.UnSetLock();
        }
    });

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep) /= rNode.FastGetSolutionStepValue(NODAL_AREA);
    });
}

void DerivativesRecoveryUtility::CalculateLaplacian(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const Variable<double>& rIntermediateVariable,
    const std::size_t BufferStep)
{
    CalculateDivergence(rModelPart, rOriginVariable, rIntermediateVariable, BufferStep);
    CalculateGradient(rModelPart, rIntermediateVariable, rDestinationVariable, BufferStep);
}

void DerivativesRecoveryUtility::CalculateSuperconvergentDivergence(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        double& divergence = rNode.FastGetSolutionStepValue(rDestinationVariable);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(NODAL_WEIGHTS);
        divergence = 0.0;

        for (unsigned int i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
            const array_1d<double,3>& value = neigh_nodes[i_neigh].FastGetSolutionStepValue(rOriginVariable);

            for (unsigned int d = 0; d < TDim; ++d){
                divergence += nodal_weights[3 * i_neigh + d] * value[d];
            }
        }
    });
}

void DerivativesRecoveryUtility::CalculateSuperconvergentGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        array_1d<double,3>& gradient = rNode.FastGetSolutionStepValue(rDestinationVariable);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(NODAL_WEIGHTS);
        gradient = ZeroVector(3);

        for (std::size_t i_neigh = 0; i_neigh < n_neigh; ++i_neigh){
            const double& value = neigh_nodes[i_neigh].FastGetSolutionStepValue(rOriginVariable);

            for (unsigned int d = 0; d < 3; ++d){
                gradient[d] += nodal_weights[3 * i_neigh + d] * value;
            }
        }
    });
}

void DerivativesRecoveryUtility::CalculateSuperconvergentLaplacian(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const Variable<double>& rIntermediateVariable,
    const std::size_t BufferStep)
{
    CalculateSuperconvergentDivergence(rModelPart, rOriginVariable, rIntermediateVariable, BufferStep);
    CalculateSuperconvergentGradient(rModelPart, rIntermediateVariable, rDestinationVariable, BufferStep);
}

}  // namespace Kratos.
