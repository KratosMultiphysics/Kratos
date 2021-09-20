//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors: Suneth Warnakulasuriya
//

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/pointer_map_communicator.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "sensitivity_utilities.h"

namespace Kratos
{
void SensitivityUtilities::AddMatrixSubBlock(
    Matrix& rOutput,
    const Matrix& rInput,
    const int RowOffset,
    const int ColOffset)
{
    KRATOS_TRY

    const int rows = rInput.size1();
    const int cols = rInput.size2();

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size1()) < RowOffset + rows)
        << "Output matrix number of rows is smaller than input matrix "
           "number of rows with offset. [ Output Matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2()
        << " ), Input matrix size = ( " << rows << ", " << cols
        << " ), Offsets = ( " << RowOffset << ", " << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size2()) < ColOffset + cols)
        << "Output matrix number of cols is smaller than input matrix "
           "number of cols with offset. [ Output Matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2()
        << " ), Input matrix size = ( " << rows << ", " << cols
        << " ), Offsets = ( " << RowOffset << ", " << ColOffset << " ) ].";

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            rOutput(i + RowOffset, j + ColOffset) += rInput(i, j);
        }
    }

    KRATOS_CATCH("");
}

void SensitivityUtilities::GetMatrixSubBlock(
    Matrix& rOutput,
    const Matrix& rInput,
    const int RowOffset,
    const int Rows,
    const int ColOffset,
    const int Cols)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rInput.size1()) < RowOffset + Rows)
        << "rInput matrix number of rows is smaller than number of rows "
           "with offset. [ Input Matrix size = ( "
        << rInput.size1() << ", " << rInput.size2() << " ), SubBlock size = ( "
        << Rows << ", " << Cols << " ), Offsets = ( " << RowOffset << ", "
        << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rInput.size2()) < ColOffset + Cols)
        << "rInput matrix number of cols is smaller than number of cols "
           "with offset. [ Input Matrix size = ( "
        << rInput.size1() << ", " << rInput.size2() << " ), SubBlock size = ( "
        << Rows << ", " << Cols << " ), Offsets = ( " << RowOffset << ", "
        << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size1()) != Rows ||
                          static_cast<int>(rOutput.size2()) != Cols)
        << "Output matrix type mismatch. [ Output matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2() << " ), SubBlockSize = ( "
        << Rows << ", " << Cols << " ) ].";

    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            rOutput(i, j) = rInput(RowOffset + i, ColOffset + j);
        }
    }

    KRATOS_CATCH("");
}

void SensitivityUtilities::ComputeEntityGeometryNeighbourNodeMap(
    std::unordered_map<int, std::unordered_map<int, int>>& rDerivativeNodesMap,
    const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
    const Geometry<ModelPart::NodeType>& rEntityGeometry,
    const Flags& rFlag,
    const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rEntityGeometry.PointsNumber();

    for (int i_base_node = 0; i_base_node < number_of_nodes; ++i_base_node) {
        const auto& r_base_node = rEntityGeometry[i_base_node];

        if (r_base_node.Is(rFlag) == CheckValue) {
            const int base_node_id = r_base_node.Id();

            std::unordered_map<int, int> derivative_node_map;
            for (int i_deriv_node = 0; i_deriv_node < number_of_nodes; ++i_deriv_node) {
                const int deriv_node_id = rEntityGeometry[i_deriv_node].Id();

                if (base_node_id == deriv_node_id) {
                    derivative_node_map[i_deriv_node] = 0;
                } else {
                    const auto p_itr = rNeighbourNodeIdsMap.find(base_node_id);
                    const auto& r_neighbour_node_indices = p_itr->second;

                    derivative_node_map[i_deriv_node] =
                        r_neighbour_node_indices.size() + 2;
                    for (int j = 0;
                         j < static_cast<int>(r_neighbour_node_indices.size()); ++j) {
                        if (r_neighbour_node_indices[j] == deriv_node_id) {
                            derivative_node_map[i_deriv_node] = j + 1;
                            break;
                        }
                    }

                    KRATOS_ERROR_IF(derivative_node_map[i_deriv_node] ==
                                    static_cast<int>(r_neighbour_node_indices.size() + 2))
                        << "Derivative node id " << deriv_node_id
                        << " not found in neighbour nodes list in node with id "
                        << base_node_id << ".";
                }
            }

            rDerivativeNodesMap[i_base_node] = derivative_node_map;
        }
    }

    KRATOS_CATCH("");
}

template <>
KRATOS_API(KRATOS_CORE) ModelPart::ElementsContainerType& SensitivityUtilities::GetContainer<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
KRATOS_API(KRATOS_CORE) ModelPart::ConditionsContainerType& SensitivityUtilities::GetContainer<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <class TContainerType>
void SensitivityUtilities::AssignEntityDerivativesToNodes(
    ModelPart& rModelPart,
    const int DerivativeDimension,
    const Variable<Matrix>& rDerivativeVariable,
    const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
    const double Weight,
    const Flags& rFlag,
    const bool CheckValue)
{
    KRATOS_TRY

    auto& entity_container = GetContainer<TContainerType>(rModelPart);

    auto& r_nodes = rModelPart.Nodes();

    int value_dimension = 0;
    if (entity_container.size() != 0) {
        value_dimension = entity_container.begin()->GetValue(rDerivativeVariable).size2();
    }

    value_dimension = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(value_dimension);

    VariableUtils().SetFlag(VISITED, false, r_nodes);

    // identify nodes where neighbours are required
    int number_of_visited_nodes = block_for_each<SumReduction<int>>(entity_container, [&](typename TContainerType::value_type& rEntity) {
        if (rEntity.Is(rFlag) == CheckValue) {
            KRATOS_ERROR_IF(!rEntity.Has(rDerivativeVariable))
                << rDerivativeVariable.Name() << " not found in data value container of "
                << rEntity.Info() << ".";

            const Matrix& r_value = rEntity.GetValue(rDerivativeVariable);

            KRATOS_ERROR_IF(value_dimension != static_cast<int>(r_value.size2()))
                << rDerivativeVariable.Name()
                << " matrix value second dimension is not consistent at "
                << rEntity.Info() << " [ required dimension size: " << value_dimension
                << ", obtained matrix: " << r_value << " ].\n";

            auto& r_geometry = rEntity.GetGeometry();
            int number_of_visited_nodes = 0;
            for (auto& r_node : r_geometry) {
                if (r_node.Is(rFlag) == CheckValue) {
                    r_node.SetLock();
                    r_node.Set(VISITED, true);
                    r_node.UnSetLock();
                    ++number_of_visited_nodes;
                }
            }

            return number_of_visited_nodes;
        }

        return 0;
    });

    rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);

    // this is done to identify whether it is required to continue with assigning entity derivatives to nodes
    // if no nodes are flagged to get the entity derivatives (such as in no-slip case) then this can be avoided by checking this.
    number_of_visited_nodes = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(number_of_visited_nodes);

    if (number_of_visited_nodes > 0) {
        KRATOS_ERROR_IF(value_dimension == 0)
            << "Column dimension (representing dimensionality of the value "
                "where derivatives are calculated) of the matrix values are "
                "zero in "
            << rModelPart.Name() << ". Please assign proper matrix values for "
            << rDerivativeVariable.Name() << ".";

        // resizing matrices
        block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
            if (rNode.Is(VISITED)) {
                const int node_id = rNode.Id();

                const auto p_itr = rNeighbourNodeIdsMap.find(node_id);
                KRATOS_ERROR_IF(p_itr == rNeighbourNodeIdsMap.end())
                    << node_id << " is not found in nodal neighbours map.";

                const int number_of_neighbour_nodes = p_itr->second.size();

                // neighbour nodes only contain neighbours, not self. But derivatives are there for neighbour
                // nodes and self node. So we reserve space for number_of_neighbour_nodes + 1 nodes
                // first block always represent self node, rest in the order of rNeighbourNodeIdsMap vector
                rNode.SetValue(rDerivativeVariable,
                                Matrix((number_of_neighbour_nodes + 1) * DerivativeDimension,
                                        value_dimension, 0.0));
            } else {
                // initializing these unused node matrices are required to do non-historical assembly
                // otherwise, this method will fail with TContainerType = ModelPart::ConditionsContainerType
                // at TransferDistributedValues in MPICommunicator::AssembleDynamicMatrixValues
                rNode.SetValue(rDerivativeVariable, Matrix(1, 1, 0.0));
            }
        });

        block_for_each(entity_container, [&](typename TContainerType::value_type& rEntity) {
            if (rEntity.Is(rFlag) == CheckValue) {
                auto& r_geometry = rEntity.GetGeometry();
                const int number_of_nodes = r_geometry.PointsNumber();

                std::unordered_map<int, std::unordered_map<int, int>> derivative_nodes_map;

                // calculate the node mapping
                ComputeEntityGeometryNeighbourNodeMap(
                    derivative_nodes_map, rNeighbourNodeIdsMap, r_geometry, rFlag, CheckValue);

                const Matrix& r_entity_derivatives =
                    rEntity.GetValue(rDerivativeVariable) * Weight;

                // move this variable also to TLS storage
                Matrix nodal_derivative(DerivativeDimension, value_dimension);

                // placing derivatives correctly
                for (int i_base_node = 0; i_base_node < number_of_nodes; ++i_base_node) {
                    auto& r_base_node = r_geometry[i_base_node];

                    if (r_base_node.Is(rFlag) == CheckValue) {
                        const auto& r_derivative_nodes_map =
                            derivative_nodes_map.find(i_base_node)->second;

                        for (int i_deriv_node = 0;
                                i_deriv_node < number_of_nodes; ++i_deriv_node) {
                            GetMatrixSubBlock(nodal_derivative, r_entity_derivatives,
                                                i_deriv_node * DerivativeDimension,
                                                DerivativeDimension, 0, value_dimension);

                            r_base_node.SetLock();
                            AddMatrixSubBlock(
                                r_base_node.GetValue(rDerivativeVariable), nodal_derivative,
                                r_derivative_nodes_map.find(i_deriv_node)->second * DerivativeDimension,
                                0);
                            r_base_node.UnSetLock();
                        }
                    }
                }
            }
        });

        rModelPart.GetCommunicator().AssembleNonHistoricalData(rDerivativeVariable);
    }

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(KRATOS_CORE) void SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ElementsContainerType>
    (
        ModelPart&,
        const int,
        const Variable<Matrix>&,
        const std::unordered_map<int, std::vector<int>>&,
        const double,
        const Flags&,
        const bool
    );

template KRATOS_API(KRATOS_CORE) void SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>
    (
        ModelPart&,
        const int,
        const Variable<Matrix>&,
        const std::unordered_map<int, std::vector<int>>&,
        const double,
        const Flags&,
        const bool
    );

} // namespace Kratos
