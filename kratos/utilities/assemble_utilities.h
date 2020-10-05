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

#if !defined(KRATOS_ASSEMBLE_UTILITIES)
#define  KRATOS_ASSEMBLE_UTILITIES

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/global_pointers_unordered_map.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) AssembleUtilities
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AssembleUtilities);

    template<class TDataType>
    using TMap = std::unordered_map<int, TDataType>;

    template<class TEntityType, class TDataType>
    using TGPMap = GlobalPointersUnorderedMap<TEntityType, TDataType>;

    template<class TContainerType, class TDataType>
    using TGPContainerMap = TGPMap<typename TContainerType::value_type, TDataType>;

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    ///@}
    ///@name Operations
    ///@{

    template<class TDataType>
    static void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rNodalValuesMap)
    {
        AssembleCurrentDataWithValuesMap<TDataType>(
            rModelPart, rVariable,
            GetGlobalPointersMap<NodesContainerType, TDataType>(rModelPart, rNodalValuesMap));
    }

    template <class TDataType>
    static void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TGPMap<NodeType, TDataType>& rNodalValuesMap)
    {
        CheckHistoricalVariable(rModelPart, rVariable);
        AssembleDataWithEntityValuesMap(
            rModelPart.Nodes(), rModelPart.GetCommunicator().GetDataCommunicator(),
            rVariable, rNodalValuesMap, UpdateHistoricalNodalValue<TDataType>);
        rModelPart.GetCommunicator().SynchronizeVariable(rVariable);
    }

    template<class TContainerType, class TDataType>
    static void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TMap<TDataType>& rEntityValuesMap)
    {
        AssembleNonHistoricalDataWithValuesMap<TContainerType, TDataType>(
            rModelPart, rVariable,
            GetGlobalPointersMap<TContainerType, TDataType>(rModelPart, rEntityValuesMap));
    }

    template <class TContainerType, class TDataType>
    static void AssembleNonHistoricalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const TGPContainerMap<TContainerType, TDataType>& rEntityValuesMap)
    {
        auto& r_container = GetContainer<TContainerType>(rModelPart);
        auto& r_communicator = rModelPart.GetCommunicator();
        AssembleDataWithEntityValuesMap(
            r_container, r_communicator.GetDataCommunicator(), rVariable, rEntityValuesMap,
            UpdateNonHistoricalValue<typename TContainerType::value_type, TDataType>);
        SynchronizeVariables(r_container, r_communicator, rVariable);
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    template<class TContainerType>
    static TContainerType& GetContainer(
        ModelPart& rModelPart);

    template<class TContainerType, class TDataType>
    static void SynchronizeVariables(
        TContainerType& rContainer,
        Communicator& rCommunicator,
        const Variable<TDataType>& rVariable)
    {
        // do nothing here for ElementsContainerType and ConditionsContainerType
    }

    template <class TDataType>
    static void SynchronizeVariables(
        NodesContainerType& rContainer,         // required for partial template specialization
        Communicator& rCommunicator,
        const Variable<TDataType>& rVariable)
    {
        rCommunicator.SynchronizeNonHistoricalVariable(rVariable);
    }

    template <class TContainerType, class TDataType>
    static TGPContainerMap<TContainerType, TDataType> GetGlobalPointersMap(
        ModelPart& rModelPart,
        const TMap<TDataType>& rValuesMap)
    {
        KRATOS_TRY

        const auto& entity_id_list = GetKeys(rValuesMap);
        const auto& id_gp_map = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(
            GetContainer<TContainerType>(rModelPart), entity_id_list,
            rModelPart.GetCommunicator().GetDataCommunicator());

        TGPContainerMap<TContainerType, TDataType> gp_map;
        for (const auto& r_item : id_gp_map) {
            gp_map[r_item.second] = rValuesMap.find(r_item.first)->second;
        }

        return gp_map;

        KRATOS_CATCH("");
    }

    template<class TDataType>
    static void CheckHistoricalVariable(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!rModelPart.HasNodalSolutionStepVariable(rVariable))
            << rVariable.Name() << " not found in solution step variables list of "
            << rModelPart.Name() << ".\n";

        KRATOS_CATCH("");
    }

    template<class TDataType>
    static void UpdateHistoricalNodalValue(
        NodeType& rNode,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        rNode.FastGetSolutionStepValue(rVariable) += rValue;
    }

    template<class TEntityType, class TDataType>
    static void UpdateNonHistoricalValue(
        TEntityType& rEntity,
        const Variable<TDataType>& rVariable,
        const TDataType& rValue)
    {
        KRATOS_TRY

        rEntity.GetValue(rVariable) += rValue;

        KRATOS_CATCH("");
    }

    template<class TKey, class... TArgs>
    static std::vector<TKey> GetKeys(const std::unordered_map<TKey, TArgs...>& rMap) {
        std::vector<TKey> key_list;
        key_list.reserve(rMap.size());
        for (const auto& r_item : rMap) {
            key_list.push_back(r_item.first);
        }

        return key_list;
    }

    template <class TContainerType, class TDataType, class TUpdateFunction>
    static void AssembleDataWithEntityValuesMap(
        TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator,
        const Variable<TDataType>& rVariable,
        const TGPContainerMap<TContainerType, TDataType>& rValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        if (rDataCommunicator.IsDistributed()) {
            MPIAssembleDataWithEntityValuesMap(rContainer, rDataCommunicator, rVariable,
                                               rValuesMap, rUpdateFunction);
        } else {
            LocalAssembleDataWithEntityValuesMap(rContainer, rVariable,
                                                 rValuesMap, rUpdateFunction);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Assembles entity values given in the map
     *
     * This method can assemble entity(nodal/elemental/condition) values according to given values map.
     * No clearing of entity values are done, therefore, assemble will add
     * values to existing values. (OpenMP version)
     *
     * @tparam TContainerType
     * @tparam TDataType
     * @tparam TUpdateFunction
     * @param rContainer            Entity container
     * @param rVariable             Variable to store assembled values
     * @param rValuesMap            Values map with entity global pointer and value
     * @param rUpdateFunction       Update function
     */
    template<class TContainerType, class TDataType, class TUpdateFunction>
    static void LocalAssembleDataWithEntityValuesMap(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const TGPContainerMap<TContainerType, TDataType>& rValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        const auto& entity_gps = GetKeys(rValuesMap);

        IndexPartition<int>(entity_gps.size()).for_each([&](const int Index) {
            auto entity_gp = entity_gps[Index];
            rUpdateFunction(*entity_gp, rVariable, rValuesMap.find(entity_gp)->second);
        });

        KRATOS_CATCH("");
    }

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
     * @param rEntityValuesMap      Entity values map with entity global pointer and value
     * @param rUpdateFunction       Update function
     */
    template <class TContainerType, class TDataType, class TUpdateFunction>
    static void MPIAssembleDataWithEntityValuesMap(
        TContainerType& rLocalContainer,
        const DataCommunicator& rDataCommunicator,
        const Variable<TDataType>& rVariable,
        const TGPContainerMap<TContainerType, TDataType>& rEntityValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        using gp_map_type = TGPContainerMap<TContainerType, TDataType>;

        const auto& entity_gps = GetKeys(rEntityValuesMap);

        // compute send and receive ranks ids
        std::vector<int> send_receive_list;
        TMap<gp_map_type> send_entity_values_map;
        gp_map_type local_entity_values_map;

        const int my_rank = rDataCommunicator.Rank();

        for (const auto& entity_gp : entity_gps) {
            const int entity_rank = entity_gp.GetRank();
            if (entity_rank != my_rank) {
                send_receive_list.push_back(entity_rank);
                auto& entity_values_map = send_entity_values_map[entity_rank];
                entity_values_map[entity_gp] = rEntityValuesMap.find(entity_gp)->second;
            } else {
                local_entity_values_map[entity_gp] =
                    rEntityValuesMap.find(entity_gp)->second;
            }
        }

        // compute the communication plan
        std::sort(send_receive_list.begin(), send_receive_list.end());
        const auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(
            send_receive_list, rDataCommunicator);

        // update local values
        LocalAssembleDataWithEntityValuesMap(
            rLocalContainer, rVariable, local_entity_values_map, rUpdateFunction);

        // update values received from remote processes
        for (const int color : colors) {
            if (color >= 0) {
                const auto& tmp = rDataCommunicator.SendRecv(
                    send_entity_values_map[color], color, color);

                LocalAssembleDataWithEntityValuesMap(rLocalContainer, rVariable,
                                                     tmp, rUpdateFunction);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_ASSEMBLE_UTILITIES