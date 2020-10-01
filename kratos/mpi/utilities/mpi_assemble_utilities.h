//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_MPI_ASSEMBLE_UTILITIES)
#define  KRATOS_MPI_ASSEMBLE_UTILITIES

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/assemble_utilities.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class MPIAssembleUtilities : public AssembleUtilities
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MPIAssembleUtilities);

    using BaseType = AssembleUtilities;

    template<class TDataType>
    using TMap = BaseType::TMap<TDataType>;

    using NodeType = BaseType::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    virtual ~MPIAssembleUtilities() = default;

    ///@}
    ///@name Operations
    ///@{

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const override
    {
        BaseType::CheckHistoricalVariable(rModelPart, rVariable);
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateHistoricalNodalValue<int>);
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override
    {
        BaseType::CheckHistoricalVariable(rModelPart, rVariable);
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateHistoricalNodalValue<double>);
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override
    {
        BaseType::CheckHistoricalVariable(rModelPart, rVariable);
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateHistoricalNodalValue<array_1d<double, 3>>);
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const override
    {
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateNonHistoricalNodalValue<int>);
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    }

    void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override
    {
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateNonHistoricalNodalValue<double>);
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    }

    void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override
    {
        MPIAssembleUtilities::MPIAssembleCurrentDataWithNodalValuesMap(
            rModelPart, rVariable, rNodalValuesMap,
            BaseType::UpdateNonHistoricalNodalValue<array_1d<double, 3>>);
        rModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Assembles nodal values given in the map
     *
     * This method can assemble nodal values according to given nodal map.
     * No clearing of nodal values are done, therefore, assemble will add
     * values to existing nodal values.
     *
     * This is the MPI version. Each process can give their own rNodalValuesMap based on their computations.
     * This method finds where given node id (from each rank) belong (the owner) to and updates accordingly.
     * The rNodalValuesMap can have node ids corresponding to LocalMesh/GhostMesh as well as nodes outside
     * LocalMesh/GhostMesh.
     *
     * This only updates the local mesh in each rank, therefore synchronization should be done afterwards.
     *
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rModelPart            Model part where nodes need to update
     * @param rVariable             Variable to store assembled values
     * @param rNodalValuesMap       Nodal values map with node_id and value
     * @param rUpdateFunction       Update function
     */
    template<class TDataType, class TUpdateFunction>
    void static MPIAssembleCurrentDataWithNodalValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rNodalValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        if (!rNodalValuesMap.empty()) {
            const auto& get_keys = [](const TMap<TDataType>& rMap) {
                std::vector<int> keys;
                keys.reserve(rMap.size());
                for (const auto& r_item : rMap) {
                    keys.push_back(r_item.first);
                }

                return keys;
            };

            const auto& update_local_nodes = [&](const std::vector<int>& rNodeIds,
                                                 const TMap<TDataType>& rValuesMap) {
                IndexPartition<int>(rNodeIds.size()).for_each([&](const int Index) {
                    const int node_id = rNodeIds[Index];
                    auto& r_node = rModelPart.GetNode(node_id);
                    rUpdateFunction(r_node, rVariable, rValuesMap.find(node_id)->second);
                });
            };

            // gather all node ids required
            const std::vector<int>& node_ids = get_keys(rNodalValuesMap);

            auto& r_communicator = rModelPart.GetCommunicator();
            const auto& r_data_communicator = r_communicator.GetDataCommunicator();

            const int my_rank = r_data_communicator.Rank();

            const auto& all_rank_node_ids = r_data_communicator.Gatherv(node_ids, 0);
            std::vector<int> all_node_ids;
            if (my_rank == 0) {
                for (const auto& rank : all_rank_node_ids) {
                    for (const auto node_id : rank) {
                        all_node_ids.push_back(node_id);
                    }
                }
                std::sort(all_node_ids.begin(), all_node_ids.end());
                auto last = std::unique(all_node_ids.begin(), all_node_ids.end());
                all_node_ids.erase(last, all_node_ids.end());
            }

            int number_of_communication_nodes = all_node_ids.size();
            // now we get all nodes which requires updating in all_node_ids in all ranks
            r_data_communicator.Broadcast(number_of_communication_nodes, 0);
            if (my_rank != 0) {
                all_node_ids.resize(number_of_communication_nodes);
            }
            r_data_communicator.Broadcast(all_node_ids, 0);

            // identify ranks which these nodes belongs to
            auto& local_mesh = r_communicator.LocalMesh();
            std::vector<int> local_node_id_ranks(number_of_communication_nodes);

            IndexPartition<int>(number_of_communication_nodes).for_each([&](const int Index) {
                local_node_id_ranks[Index] =
                    local_mesh.HasNode(all_node_ids[Index]) ? my_rank : -1;
            });

            const auto& all_rank_nodal_id_ranks =
                r_data_communicator.Gatherv(local_node_id_ranks, 0);
            std::vector<int> all_node_id_ranks(number_of_communication_nodes, -1);
            if (my_rank == 0) {
                IndexPartition<int>(number_of_communication_nodes).for_each([&](const int Index) {
                    for (unsigned int i_rank = 0;
                         i_rank < all_rank_nodal_id_ranks.size(); ++i_rank) {
                        const int node_rank = all_rank_nodal_id_ranks[i_rank][Index];
                        if (node_rank > -1) {
                            if (all_node_id_ranks[Index] == -1) {
                                all_node_id_ranks[Index] = node_rank;
                            } else if (all_node_id_ranks[Index] != node_rank) {
                                KRATOS_ERROR
                                    << "Node id " << all_node_ids[Index]
                                    << " does belong to more than one rank, "
                                    << "which is not allowed.\n";
                            }
                            break;
                        }
                    }

                    KRATOS_ERROR_IF(all_node_id_ranks[Index] == -1)
                        << "Node id " << all_node_ids[Index]
                        << " not found in any of the ranks.\n";
                });
            }

            // now broadcast all the ranks properly to all ranks
            r_data_communicator.Broadcast(all_node_id_ranks, 0);

            // compute send and receive ranks ids
            std::vector<int> send_receive_list;
            TMap<TMap<TDataType>> send_nodal_values_map;
            TMap<TDataType> local_nodal_values_map;

            for (int i_node = 0; i_node < number_of_communication_nodes; ++i_node) {
                const int node_id = all_node_ids[i_node];
                const auto p_itr = rNodalValuesMap.find(node_id);
                if (p_itr != rNodalValuesMap.cend()) {
                    const int node_rank = all_node_id_ranks[i_node];
                    if (node_rank != my_rank) {
                        send_receive_list.push_back(node_rank);
                        auto& nodal_values_map = send_nodal_values_map[node_rank];
                        nodal_values_map[node_id] = p_itr->second;
                    } else {
                        local_nodal_values_map[node_id] = p_itr->second;
                    }
                }
            }

            // compute the communication plan
            std::sort(send_receive_list.begin(), send_receive_list.end());
            const auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(
                send_receive_list, r_data_communicator);

            // update local values
            const auto& local_nodal_ids = get_keys(local_nodal_values_map);
            update_local_nodes(local_nodal_ids, local_nodal_values_map);

            // update values received from remote processes
            for (const int color : colors) {
                if (color >= 0) {
                    const auto& tmp = r_data_communicator.SendRecv(
                        send_nodal_values_map[color], color, color);

                    const auto& local_nodal_ids = get_keys(tmp);
                    update_local_nodes(local_nodal_ids, tmp);
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_MPI_ASSEMBLE_UTILITIES