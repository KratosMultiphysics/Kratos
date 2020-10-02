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
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"

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

            const auto& update_local_nodes = [&](const TMap<TDataType>& rValuesMap) {
                const auto& keys = get_keys(rValuesMap);
                IndexPartition<int>(keys.size()).for_each([&](const int Index) {
                    const int node_id = keys[Index];
                    auto& r_node = rModelPart.GetNode(node_id);
                    rUpdateFunction(r_node, rVariable, rValuesMap.find(node_id)->second);
                });
            };

            auto& r_communicator = rModelPart.GetCommunicator();
            const auto& r_data_communicator = r_communicator.GetDataCommunicator();
            const int my_rank = r_data_communicator.Rank();

            // gather all node ids required
            const std::vector<int>& node_ids = get_keys(rNodalValuesMap);
            auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(
                rModelPart.Nodes(), node_ids, r_data_communicator);

            GlobalPointerCommunicator<NodeType> pointer_comm(
                r_data_communicator, gp_list.ptr_begin(), gp_list.ptr_end());

            KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX)) << "PARTITION_INDEX variable is not found in nodal solution step variables list of "
                                                                                       << rModelPart
                                                                                              .Name()
                                                                                       << ".\n";

            auto partition_index_proxy =
                pointer_comm.Apply([](GlobalPointer<NodeType>& gp) -> int {
                    return gp->FastGetSolutionStepValue(PARTITION_INDEX);
                });

            // compute send and receive ranks ids
            std::vector<int> send_receive_list;
            TMap<TMap<TDataType>> send_nodal_values_map;
            TMap<TDataType> local_nodal_values_map;

            for (unsigned int i = 0; i < node_ids.size(); ++i) {
                auto& gp = gp_list(i);
                const int partition_index = partition_index_proxy.Get(gp);
                const int node_id = node_ids[i];
                if (partition_index != my_rank) {
                    send_receive_list.push_back(partition_index);
                    auto& nodal_values_map = send_nodal_values_map[partition_index];
                    nodal_values_map[node_id] = rNodalValuesMap.find(node_id)->second;
                } else {
                    local_nodal_values_map[node_id] = rNodalValuesMap.find(node_id)->second;
                }
            }

            // compute the communication plan
            std::sort(send_receive_list.begin(), send_receive_list.end());
            const auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(
                send_receive_list, r_data_communicator);

            // update local values
            update_local_nodes(local_nodal_values_map);

            // update values received from remote processes
            for (const int color : colors) {
                if (color >= 0) {
                    const auto& tmp = r_data_communicator.SendRecv(
                        send_nodal_values_map[color], color, color);

                    update_local_nodes(tmp);
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