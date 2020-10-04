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
#define KRATOS_MPI_ASSEMBLE_UTILITIES

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/assemble_utilities.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_MPI_CORE) MPIAssembleUtilities : public AssembleUtilities
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MPIAssembleUtilities);

    using BaseType = AssembleUtilities;

    template <class TDataType>
    using TMap = BaseType::TMap<TDataType>;

    using NodeType = BaseType::NodeType;

    using ElementType = BaseType::ElementType;

    using ConditionType = BaseType::ConditionType;

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
        const TMap<int>& rNodalValuesMap) const override;

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override;

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TMap<int>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TMap<double>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TMap<array_1d<double, 3>>& rNodalValuesMap) const override;

    ///@}

private:
    ///@name Private Operations
    ///@{

    /**
     * @brief Assembles entity values given in the map
     *
     * This method can assemble entity (nodal/elemental/condition) values according to given id map.
     * No clearing of entity values are done, therefore, assemble will add
     * values to existing values.
     *
     * This is the MPI version. Each process can give their own rEntityValuesMap based on their computations.
     * This method finds where given entity id (from each rank) belong (the owner) to and updates accordingly.
     * The rEntityValuesMap can have entity ids corresponding to LocalMesh/GhostMesh as well as entities outside
     * LocalMesh/GhostMesh.
     *
     * This only updates the local mesh in each rank, therefore synchronization should be done afterwards.
     *
     * @tparam TContainerType
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rLocalContainer       Local entities container
     * @param rDataCommunicator     Data communicator
     * @param rVariable             Variable to store assembled values
     * @param rEntityValuesMap      Entity values map with entity_id and value
     * @param rUpdateFunction       Update function
     */
    template <class TContainerType, class TDataType, class TUpdateFunction>
    static void MPIAssembleDataWithEntityValuesMap(
        TContainerType& rLocalContainer,
        const DataCommunicator& rDataCommunicator,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rEntityValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        if (!rEntityValuesMap.empty()) {
            const auto& update_local_entities = [&](const TMap<TDataType>& rValuesMap) {
                const auto& keys = GetKeys(rValuesMap);
                IndexPartition<int>(keys.size()).for_each([&](const int Index) {
                    const int entity_id = keys[Index];
                    auto p_entity = rLocalContainer.find(entity_id);

                    KRATOS_ERROR_IF(p_entity == rLocalContainer.end())
                        << "Entity id " << entity_id << " not found in local entities.\n";

                    rUpdateFunction(*p_entity, rVariable,
                                    rValuesMap.find(entity_id)->second);
                });
            };

            // gather all entity ids, ranks required for communication
            std::vector<int> entity_ids, entity_ranks;
            GetEntityIdsAndRanks(entity_ids, entity_ranks, rLocalContainer,
                                 rDataCommunicator, rEntityValuesMap);

            // compute send and receive ranks ids
            std::vector<int> send_receive_list;
            TMap<TMap<TDataType>> send_entity_values_map;
            TMap<TDataType> local_entity_values_map;

            const int my_rank = rDataCommunicator.Rank();

            for (unsigned int i = 0; i < entity_ids.size(); ++i) {
                const int entity_rank = entity_ranks[i];
                const int entity_id = entity_ids[i];
                if (entity_rank != my_rank) {
                    send_receive_list.push_back(entity_rank);
                    auto& entity_values_map = send_entity_values_map[entity_rank];
                    entity_values_map[entity_id] =
                        rEntityValuesMap.find(entity_id)->second;
                } else {
                    local_entity_values_map[entity_id] =
                        rEntityValuesMap.find(entity_id)->second;
                }
            }

            // compute the communication plan
            std::sort(send_receive_list.begin(), send_receive_list.end());
            const auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(
                send_receive_list, rDataCommunicator);

            // update local values
            update_local_entities(local_entity_values_map);

            // update values received from remote processes
            for (const int color : colors) {
                if (color >= 0) {
                    const auto& tmp = rDataCommunicator.SendRecv(
                        send_entity_values_map[color], color, color);

                    update_local_entities(tmp);
                }
            }
        }

        KRATOS_CATCH("");
    }

    /** Retrieves entity ids and their ranks from all processes in mpi
     *
     * Retrieves all entity ids and their ranks from all processes which needs
     * communicattion for assembly (entity ids are given in the rEntityValuesMap)
     *
     * @tparam TContainerType
     * @tparam TDataType
     * @param rAllEntityIds         Output all entity ids vector
     * @param rAllEntityRanks       Output all entity ids' ranks vector
     * @param rLocalContainer       Local entity container
     * @param rDataCommunicator     Data communicator
     * @param rEntityValuesMap      Entity values map with entity_id and values
     */
    template <class TContainerType, class TDataType>
    static void GetEntityIdsAndRanks(
        std::vector<int>& rAllEntityIds,
        std::vector<int>& rAllEntityRanks,
        const TContainerType& rLocalContainer,
        const DataCommunicator& rDataCommunicator,
        const TMap<TDataType>& rEntityValuesMap)
    {
        KRATOS_TRY

        const int my_rank = rDataCommunicator.Rank();
        const auto& local_ids = GetKeys(rEntityValuesMap);
        const auto& all_rank_ids = rDataCommunicator.Gatherv(local_ids, 0);

        rAllEntityIds.clear();
        if (my_rank == 0) {
            for (const auto& rank : all_rank_ids) {
                for (const auto id : rank) {
                    rAllEntityIds.push_back(id);
                }
            }
            std::sort(rAllEntityIds.begin(), rAllEntityIds.end());
            auto last = std::unique(rAllEntityIds.begin(), rAllEntityIds.end());
            rAllEntityIds.erase(last, rAllEntityIds.end());
        }

        int number_of_communication_ids = rAllEntityIds.size();
        // now we get all nodes which requires updating in rAllEntityIds in all ranks
        rDataCommunicator.Broadcast(number_of_communication_ids, 0);
        if (my_rank != 0) {
            rAllEntityIds.resize(number_of_communication_ids);
        }
        rDataCommunicator.Broadcast(rAllEntityIds, 0);

        // identify ranks which these nodes belongs to
        std::vector<int> local_entity_ranks(number_of_communication_ids);
        IndexPartition<int>(number_of_communication_ids).for_each([&](const int Index) {
            local_entity_ranks[Index] =
                (rLocalContainer.find(rAllEntityIds[Index]) != rLocalContainer.end())
                    ? my_rank
                    : -1;
        });

        const auto& all_rank_entity_ids =
            rDataCommunicator.Gatherv(local_entity_ranks, 0);
        rAllEntityRanks.resize(number_of_communication_ids);
        block_for_each(rAllEntityRanks, [](int& rRank) { rRank = -1; });

        if (my_rank == 0) {
            IndexPartition<int>(number_of_communication_ids).for_each([&](const int Index) {
                for (unsigned int i_rank = 0; i_rank < all_rank_entity_ids.size(); ++i_rank) {
                    const int entity_rank = all_rank_entity_ids[i_rank][Index];
                    if (entity_rank > -1) {
                        if (rAllEntityRanks[Index] == -1) {
                            rAllEntityRanks[Index] = entity_rank;
                        } else if (rAllEntityRanks[Index] != entity_rank) {
                            KRATOS_ERROR
                                << "Entity id " << rAllEntityIds[Index]
                                << " does belong to more than one rank, "
                                << "which is not allowed.\n";
                        }
                        break;
                    }
                }

                KRATOS_ERROR_IF(rAllEntityRanks[Index] == -1)
                    << "Entity id " << rAllEntityIds[Index]
                    << " not found in any of the ranks.\n";
            });
        }

        // now broadcast all the entity ranks properly to all ranks
        rDataCommunicator.Broadcast(rAllEntityRanks, 0);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    std::vector<int> static GetKeys(
        const TMap<TDataType>& rMap)
    {
        std::vector<int> keys;
        keys.reserve(rMap.size());
        for (const auto& r_item : rMap) {
            keys.push_back(r_item.first);
        }

        return keys;
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_MPI_ASSEMBLE_UTILITIES