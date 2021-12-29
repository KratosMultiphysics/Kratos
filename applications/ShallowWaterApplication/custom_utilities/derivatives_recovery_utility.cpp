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
#include <unordered_set>


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"
#include "derivatives_recovery_utility.h"
#include "includes/global_pointer_variables.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

#define KRATOS_RECOVERY_FACTORIAL(n) ( \
    (n) == 2 ? 2 : \
    (n) == 3 ? 6 : \
    (n) == 4 ? 24 : \
    (n) == 5 ? 120 : \
    KRATOS_ERROR)

void DerivativesRecoveryUtility::RecoverDivergence(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        double& divergence = rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(FIRST_DERIVATIVE_WEIGHTS);
        divergence = 0.0;

        const array_1d<double,3>& first_value = rNode.FastGetSolutionStepValue(rOriginVariable, BufferStep);
        for (unsigned int d = 0; d < TDim; ++d){
            divergence += nodal_weights[d] * first_value[d];
        }

        for (unsigned int n = 0; n < n_neigh; ++n){
            const array_1d<double,3>& value = neigh_nodes[n].FastGetSolutionStepValue(rOriginVariable, BufferStep);

            for (unsigned int d = 0; d < TDim; ++d){
                divergence += nodal_weights[TDim * (n+1) + d] * value[d];
            }
        }
    });
}

void DerivativesRecoveryUtility::RecoverGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const std::size_t BufferStep)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        array_1d<double,3>& gradient = rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(FIRST_DERIVATIVE_WEIGHTS);
        gradient = ZeroVector(3);

        const double& first_value = rNode.FastGetSolutionStepValue(rOriginVariable, BufferStep);
        for (unsigned int d = 0; d < TDim; ++d) {
            gradient[d] += nodal_weights[d] * first_value;
        }

        for (std::size_t n = 0; n < n_neigh; ++n){
            const double& value = neigh_nodes[n].FastGetSolutionStepValue(rOriginVariable, BufferStep);

            for (unsigned int d = 0; d < TDim; ++d) {
                gradient[d] += nodal_weights[TDim * (n+1) + d] * value;
            }
        }
    });
}

void DerivativesRecoveryUtility::RecoverLaplacian(
    ModelPart& rModelPart,
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    const std::size_t BufferStep)
{
    const std::size_t n_terms = TDim * (TDim + 1) / 2;
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        double& laplacian = rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(SECOND_DERIVATIVE_WEIGHTS);
        laplacian = 0.0;

        const double& first_value = rNode.FastGetSolutionStepValue(rOriginVariable, BufferStep);
        for (unsigned int d = 0; d < TDim; ++d) {
            laplacian += nodal_weights[d] * first_value;
        }

        for (std::size_t n = 0; n < n_neigh; ++n){
            const double& value = neigh_nodes[n].FastGetSolutionStepValue(rOriginVariable, BufferStep);

            for (unsigned int d = 0; d < TDim; ++d) {
                laplacian += nodal_weights[n_terms * (n+1) + d] * value;
            }
        }
    });
}

void DerivativesRecoveryUtility::RecoverLaplacian(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    const std::size_t BufferStep)
{
    const std::size_t n_terms = TDim * (TDim + 1) / 2;
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();

        array_1d<double,3>& laplacian = rNode.FastGetSolutionStepValue(rDestinationVariable, BufferStep);
        const Vector& nodal_weights = rNode.FastGetSolutionStepValue(SECOND_DERIVATIVE_WEIGHTS);
        laplacian = ZeroVector(3);

        const array_1d<double,3>& first_value = rNode.FastGetSolutionStepValue(rOriginVariable, BufferStep);
        for (unsigned int d = 0; d < TDim; ++d) {
            laplacian[d] += nodal_weights[d] * first_value[d];
        }
        if (TDim == 2) {
            laplacian[0] += nodal_weights[TDim] * first_value[1];
            laplacian[1] += nodal_weights[TDim] * first_value[0];
        }
        else {
            laplacian[0] += nodal_weights[TDim] * first_value[1];
            laplacian[1] += nodal_weights[TDim] * first_value[0];

            laplacian[0] += nodal_weights[TDim+1] * first_value[2];
            laplacian[2] += nodal_weights[TDim+1] * first_value[0];

            laplacian[1] += nodal_weights[TDim+2] * first_value[2];
            laplacian[2] += nodal_weights[TDim+2] * first_value[1];
        }

        for (std::size_t n = 0; n < n_neigh; ++n){
            const array_1d<double,3>& value = neigh_nodes[n].FastGetSolutionStepValue(rOriginVariable, BufferStep);

            for (unsigned int d = 0; d < TDim; ++d) {
                laplacian[d] += nodal_weights[n_terms * (n+1) + d] * value[d];
            }
            if (TDim == 2) {
                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim] * first_value[1];
                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim] * first_value[0];
            }
            else {
                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim] * first_value[1];
                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim] * first_value[0];

                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim+1] * first_value[2];
                laplacian[2] += nodal_weights[n_terms * (n+1) + TDim+1] * first_value[0];

                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim+2] * first_value[2];
                laplacian[2] += nodal_weights[n_terms * (n+1) + TDim+2] * first_value[1];
            }
        }
    });
}

void DerivativesRecoveryUtility::ExtendNeighborsPatch(ModelPart& rModelPart)
{
    const std::size_t space_degree = 1;
    const std::size_t required_neighbors = (space_degree + 2) * (space_degree + 3) / 2;
    std::vector<std::unordered_set<int>> second_neighbors_set(rModelPart.NumberOfNodes());
    IndexPartition<int>(rModelPart.NumberOfNodes()).for_each([&](int i){
        auto it_node = rModelPart.NodesBegin() + i;
        auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
        if (neigh_nodes.size() < required_neighbors)
        {
            auto second_neighbors = second_neighbors_set.begin() + i;
            for (auto& r_neigh : neigh_nodes)
            {
                auto& extended_neighbors = r_neigh.GetValue(NEIGHBOUR_NODES);
                for (auto& second_neigh : extended_neighbors) // check if it is a new neighbour
                {
                    if (second_neigh.Id() != it_node->Id())
                    {
                        bool second_neigh_included = false;
                        for (auto& first_neigh : neigh_nodes)
                        {
                            if (second_neigh.Id() == first_neigh.Id())
                            {
                                second_neigh_included = true;
                                break;
                            }
                        }
                        if (!second_neigh_included){
                            second_neighbors->insert(second_neigh.Id());
                        }
                    }
                }
            }
        }
    });
    IndexPartition<int>(rModelPart.NumberOfNodes()).for_each([&](int i){
        auto it_node = rModelPart.NodesBegin() + i;
        auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
        if (neigh_nodes.size() < required_neighbors)
        {
            auto second_neighbors = second_neighbors_set.begin() + i;
            for (auto id : *second_neighbors)
            {
                auto second_neigh = rModelPart.pGetNode(id);
                neigh_nodes.push_back(second_neigh);
            }
        }
    });
}

void DerivativesRecoveryUtility::CalculatePolynomialWeights(ModelPart& rModelPart)
{
    constexpr std::size_t n_poly_terms = KRATOS_RECOVERY_FACTORIAL(TDim+2) / KRATOS_RECOVERY_FACTORIAL(TDim) / 2;
    constexpr std::size_t first_order_terms = KRATOS_RECOVERY_FACTORIAL(TDim+1) / KRATOS_RECOVERY_FACTORIAL(TDim);

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
        auto n_neigh = neigh_nodes.size();
        const double h_inv = 1.0 / DerivativesRecoveryUtility::CalculateMaximumDistance(rNode, neigh_nodes); // we use it as a scaling parameter to improve stability

        Matrix A(n_neigh + 1, n_poly_terms);

        A(0,0) = 1.0; // the current node coordinates
        for (std::size_t i = 1; i < n_poly_terms; ++i)
        {
            A(0,i) = 0.0;
        }

        for (std::size_t n = 0; n < n_neigh; ++n) // the neighbors coordinates
        {
            const auto& r_neigh = neigh_nodes[n];
            const array_1d<double,3> rel_coordinates = h_inv * (r_neigh - rNode);

            A(n+1,0) = 1.0;
            for (std::size_t i = 0; i < TDim; ++i) {
                const auto i_first_order = 1 + i;
                const auto i_second_order = 1 + TDim + i;
                A(n+1,i_first_order) = rel_coordinates[i];
                A(n+1,i_second_order) = rel_coordinates[i] * rel_coordinates[i];
            }
            if (TDim == 2) {
                A(n+1,n_poly_terms-1) = rel_coordinates[0] * rel_coordinates[1];
            }
            else {
                A(n+1,n_poly_terms-3) = rel_coordinates[0] * rel_coordinates[1];
                A(n+1,n_poly_terms-2) = rel_coordinates[0] * rel_coordinates[2];
                A(n+1,n_poly_terms-1) = rel_coordinates[1] * rel_coordinates[2];
            }
        }

        // The least squares projection
        double det;
        Matrix A_pseudo_inv;
        MathUtils<double>::GeneralizedInvertMatrix(A, A_pseudo_inv, det);

        // Each row of the pseudo inverse contributes to a partial derivative evaluated at the current node
        constexpr std::size_t n_first_order_terms = TDim; // x, y...
        constexpr std::size_t n_second_order_terms = n_poly_terms - first_order_terms; // x^2, y^2, xy...
        std::array<int,n_first_order_terms> first_derivative_terms;
        std::array<int,n_second_order_terms> second_derivative_terms;
        if (TDim == 2) {
            first_derivative_terms[0] = 1;
            first_derivative_terms[1] = 2;
            second_derivative_terms[0] = 3;
            second_derivative_terms[1] = 4;
            second_derivative_terms[2] = 5;
        }
        else {
            first_derivative_terms[0] = 1;
            first_derivative_terms[1] = 2;
            first_derivative_terms[2] = 3;
            second_derivative_terms[0] = 4;
            second_derivative_terms[1] = 5;
            second_derivative_terms[2] = 6;
            second_derivative_terms[3] = 7;
            second_derivative_terms[4] = 8;
            second_derivative_terms[5] = 9;
        }

        Vector& first_derivative_weights = rNode.FastGetSolutionStepValue(FIRST_DERIVATIVE_WEIGHTS);
        first_derivative_weights.resize(n_first_order_terms*(n_neigh + 1));
        Vector& second_derivative_weights = rNode.FastGetSolutionStepValue(SECOND_DERIVATIVE_WEIGHTS);
        second_derivative_weights.resize(n_second_order_terms*(n_neigh + 1));

        // Assembly of the partial derivatives contribution for every node (included the current node)
        for (std::size_t n = 0; n < n_neigh + 1; ++n)
        {
            for (std::size_t d = 0; d < n_first_order_terms; ++d)
            {
                auto row = first_derivative_terms[d];
                first_derivative_weights[n_first_order_terms*n + d] = h_inv * A_pseudo_inv(row, n);
            }
            for (std::size_t d = 0; d < n_second_order_terms; ++d)
            {
                auto row = second_derivative_terms[d];
                second_derivative_weights[n_second_order_terms*n + d] = h_inv * h_inv * A_pseudo_inv(row, n);
            }

            // Note, d^2(x^2)/dx^2 = 2  and  d^2(xy)/dxy = 1
            for (std::size_t d = 0; d < TDim; ++d)
            {
                second_derivative_weights[n_second_order_terms*n + d] *= 2.0;
            }
        }
    });
}

double DerivativesRecoveryUtility::CalculateMaximumDistance(
    const NodeType& rNode,
    GlobalPointersVector<NodeType>& rNeighbors)
{
    double max_distance = 0.0;
    for (const auto& r_neigh : rNeighbors) {
        double distance = norm_2(rNode - r_neigh);
        max_distance = std::max(max_distance, distance);
    }
    return max_distance;
}

}  // namespace Kratos.
