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

    template <class TEntityType, class TDataType>
    using TGPMap = BaseType::TGPMap<TEntityType, TDataType>;

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
        const TGPMap<NodeType, int>& rNodalValuesMap) const override;

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<NodeType, double>& rNodalValuesMap) const override;

    void AssembleCurrentDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<NodeType, int>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<NodeType, double>& rNodalValuesMap) const override;

    void AssembleNonHistoricalNodalDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<NodeType, array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<ElementType, int>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<ElementType, double>& rNodalValuesMap) const override;

    void AssembleElementDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<ElementType, array_1d<double, 3>>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const TGPMap<ConditionType, int>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const TGPMap<ConditionType, double>& rNodalValuesMap) const override;

    void AssembleConditionDataWithValuesMap(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rVariable,
        const TGPMap<ConditionType, array_1d<double, 3>>& rNodalValuesMap) const override;

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
        const TGPMap<typename TContainerType::value_type, TDataType>& rEntityValuesMap,
        const TUpdateFunction&& rUpdateFunction)
    {
        KRATOS_TRY

        using gp_map_type = TGPMap<typename TContainerType::value_type, TDataType>;

        if (!rEntityValuesMap.empty()) {
            const auto& update_local_entities = [&](const gp_map_type& rValuesMap) {
                const auto& entity_gps = BaseType::GetKeys(rValuesMap);
                IndexPartition<int>(entity_gps.size()).for_each([&](const int Index) {
                    auto entity_gp = entity_gps[Index];
                    rUpdateFunction(*entity_gp, rVariable,
                                    rValuesMap.find(entity_gp)->second);
                });
            };

            const auto& entity_gps = BaseType::GetKeys(rEntityValuesMap);

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
                    entity_values_map[entity_gp] =
                        rEntityValuesMap.find(entity_gp)->second;
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

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_MPI_ASSEMBLE_UTILITIES