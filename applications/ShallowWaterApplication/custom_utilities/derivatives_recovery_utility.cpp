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
#include "includes/checks.h"
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

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::Check(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(FIRST_DERIVATIVE_WEIGHTS, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(SECOND_DERIVATIVE_WEIGHTS, rNode)
    });
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::CheckRequiredNeighborsPatch(ModelPart& rModelPart)
{
    std::size_t space_degree = 1;
    std::size_t required_neighbors = (space_degree + TDim) * (space_degree + TDim + 1) / 2;
    if constexpr (TDim == 3) {
        required_neighbors += 1; //NOTE: the extra node is required for accuracy in the faces for the 3D case
    }
    DerivativesRecoveryUtility<TDim>::ExtendNeighborsPatch(rModelPart, required_neighbors);
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::ExtendNeighborsPatch(
    ModelPart& rModelPart,
    const std::size_t RequiredNeighbors)
{
    std::vector<std::unordered_set<int>> second_neighbors_set(rModelPart.NumberOfNodes());
    IndexPartition<int>(rModelPart.NumberOfNodes()).for_each([&](int i){
        auto it_node = rModelPart.NodesBegin() + i;
        auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
        if (neigh_nodes.size() < RequiredNeighbors)
        {
            auto second_neighbors = second_neighbors_set.begin() + i;
            DerivativesRecoveryUtility<TDim>::FindExtendedNeighbors(*it_node, neigh_nodes, *second_neighbors);
        }
    });
    IndexPartition<int>(rModelPart.NumberOfNodes()).for_each([&](int i){
        auto it_node = rModelPart.NodesBegin() + i;
        auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
        if (neigh_nodes.size() < RequiredNeighbors)
        {
            auto second_neighbors = second_neighbors_set.begin() + i;
            DerivativesRecoveryUtility<TDim>::AppendExtendedNeighbors(rModelPart, neigh_nodes, *second_neighbors);
        }
    });
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::CalculatePolynomialWeights(ModelPart& rModelPart)
{
    CheckRequiredNeighborsPatch(rModelPart);
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        bool is_converged = false;
        int max_iter = 3;
        int iter = 0;
        while (!is_converged && iter < max_iter)
        {
            is_converged = CalculateNodalPolynomialWeights(rNode);
            if (!is_converged)
            {
                auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
                std::unordered_set<int> third_neighbors;
                DerivativesRecoveryUtility<TDim>::FindExtendedNeighbors(rNode, neigh_nodes, third_neighbors);
                DerivativesRecoveryUtility<TDim>::AppendExtendedNeighbors(rModelPart, neigh_nodes, third_neighbors);
            }
            iter++;
        }
    });
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::RecoverDivergence(
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

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::RecoverGradient(
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

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::RecoverLaplacian(
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

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::RecoverLaplacian(
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
        if constexpr (TDim == 2) {
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
            if constexpr (TDim == 2) {
                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim] * value[1];
                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim] * value[0];
            }
            else {
                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim] * value[1];
                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim] * value[0];

                laplacian[0] += nodal_weights[n_terms * (n+1) + TDim+1] * value[2];
                laplacian[2] += nodal_weights[n_terms * (n+1) + TDim+1] * value[0];

                laplacian[1] += nodal_weights[n_terms * (n+1) + TDim+2] * value[2];
                laplacian[2] += nodal_weights[n_terms * (n+1) + TDim+2] * value[1];
            }
        }
    });
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::FindExtendedNeighbors(
    NodeType& rNode,
    GlobalPointersVector<NodeType>& rNeighbors,
    std::unordered_set<int>& rExtendedNeighborsId)
{
    for (auto& r_neigh : rNeighbors)
    {
        auto& extended_neighbors = r_neigh.GetValue(NEIGHBOUR_NODES);
        for (auto& candidate_neigh : extended_neighbors) // check if it is a new neighbour
        {
            if (candidate_neigh.Id() != rNode.Id())
            {
                bool candidate_neigh_included = false;
                for (auto& first_neigh : rNeighbors)
                {
                    if (candidate_neigh.Id() == first_neigh.Id())
                    {
                        candidate_neigh_included = true;
                        break;
                    }
                }
                if (!candidate_neigh_included){
                    rExtendedNeighborsId.insert(candidate_neigh.Id());
                }
            }
        }
    }
}

template<std::size_t TDim>
void DerivativesRecoveryUtility<TDim>::AppendExtendedNeighbors(
    ModelPart& rModelPart,
    GlobalPointersVector<NodeType>& rNeighbors,
    std::unordered_set<int>& rExtendedNeighborsId)
{
    for (auto id : rExtendedNeighborsId)
    {
        auto new_neigh = rModelPart.pGetNode(id);
        rNeighbors.push_back(new_neigh);
    }
}

template<std::size_t TDim>
bool DerivativesRecoveryUtility<TDim>::CalculateNodalPolynomialWeights(NodeType& rNode)
{
    constexpr std::size_t n_poly_terms = KRATOS_RECOVERY_FACTORIAL(TDim+2) / KRATOS_RECOVERY_FACTORIAL(TDim) / 2;
    constexpr std::size_t first_order_terms = KRATOS_RECOVERY_FACTORIAL(TDim+1) / KRATOS_RECOVERY_FACTORIAL(TDim);

    auto& neigh_nodes = rNode.GetValue(NEIGHBOUR_NODES);
    auto n_neigh = neigh_nodes.size();
    const double h_inv = 1.0 / DerivativesRecoveryUtility<TDim>::CalculateMaximumDistance(rNode, neigh_nodes); // we use it as a scaling parameter to improve stability

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
        if constexpr (TDim == 2) {
            A(n+1,n_poly_terms-1) = rel_coordinates[0] * rel_coordinates[1];
        }
        else {
            A(n+1,n_poly_terms-3) = rel_coordinates[0] * rel_coordinates[1];
            A(n+1,n_poly_terms-2) = rel_coordinates[0] * rel_coordinates[2];
            A(n+1,n_poly_terms-1) = rel_coordinates[1] * rel_coordinates[2];
        }
    }

    // The least squares projection
    bool is_invertible;
    Matrix A_pseudo_inv;
    is_invertible = DerivativesRecoveryUtility<TDim>::GeneralizedInvertMatrix(A, A_pseudo_inv);
    if (!is_invertible) {
        return false;
    }

    // Each row of the pseudo inverse contributes to a partial derivative evaluated at the current node
    constexpr std::size_t n_first_order_terms = TDim; // x, y...
    constexpr std::size_t n_second_order_terms = n_poly_terms - first_order_terms; // x^2, y^2, xy...
    std::array<int,n_first_order_terms> first_derivative_terms;
    std::array<int,n_second_order_terms> second_derivative_terms;
    if constexpr (TDim == 2) {
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
    return true;
}

template<std::size_t TDim>
double DerivativesRecoveryUtility<TDim>::CalculateMaximumDistance(
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

template<std::size_t TDim>
bool DerivativesRecoveryUtility<TDim>::GeneralizedInvertMatrix(
    Matrix& rInputMatrix,
    Matrix& rResult)
{
    double det;
    MathUtils<double>::GeneralizedInvertMatrix(rInputMatrix, rResult, det, -1.0);
    const double tolerance = std::numeric_limits<double>::epsilon();
    const bool throw_errors = false;
    return MathUtils<double>::CheckConditionNumber(rInputMatrix, rResult, tolerance, throw_errors);
}

template class DerivativesRecoveryUtility<2>;
template class DerivativesRecoveryUtility<3>;

}  // namespace Kratos.
