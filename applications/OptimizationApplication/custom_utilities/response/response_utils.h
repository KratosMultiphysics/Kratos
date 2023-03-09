//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <map>
#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using NodeIdsType = std::vector<IndexType>;

    ///@}
    ///@name Static operations
    ///@{

    static ModelPart& GetSensitivityModelPartForAdjointSensitivities(
        const std::vector<ModelPart*>& rSensitivityModelParts,
        ModelPart& rAnalysisModelPart,
        const bool AreSensitivityEntityParentsConsidered,
        const bool AreSensitivityEntitesConsidered,
        const bool ForceFindSensitivityEntitiesInAnalysisModelPart = false);

    static ModelPart& GetSensitivityModelPartForDirectSensitivities(
        const std::vector<ModelPart*>& rSensitivityModelParts,
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered);

    ///@}
private:
    ///@name Private classes
    ///@{

    template<class EntityType>
    class ContainerEntityMapReduction
    {
    public:
        using return_type = std::map<IndexType, EntityType*>;
        using value_type = std::vector<std::pair<IndexType, EntityType*>>;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const;

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type& rValue);

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(ContainerEntityMapReduction<EntityType>& rOther);
    };

    ///@}
    ///@name Private static operations
    ///@{

    static std::string GetCombinedModelPartsName(
        const std::string& rPrefix,
        const std::vector<ModelPart*>& rModelParts);

    template<class TContainerType>
    static void AddNeighbourEntitiesToFlaggedNodes(
        TContainerType& rContainer,
        const Variable<GlobalPointersVector<typename TContainerType::value_type>>& rNeighbourEntitiesOutputVariable,
        const Flags& rFlag,
        const bool FlagValue = true);

    template<class TEntityType>
    static void UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(
        std::map<IndexType, TEntityType*>& rOutput,
        const ModelPart::NodesContainerType& rNodes,
        const Variable<GlobalPointersVector<TEntityType>>& rNeighbourEntitiesVariable);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromEntityContainer(
        std::map<IndexType, typename TContainerType::value_type*>& rOutput,
        TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateNodeIdsEntityPtrMapFromEntityContainer(
        std::map<NodeIdsType, typename TContainerType::value_type*>& rOutput,
        TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(
        std::map<IndexType, typename TContainerType::value_type*>& rOutput,
        const std::map<NodeIdsType, typename TContainerType::value_type*>& rNodeIdsEntityPtrMap,
        const TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(
        std::map<IndexType, typename TContainerType::value_type*>& rOutput,
        TContainerType& rContainer,
        const Flags& rFlag,
        const bool FlagValue = true);

    template<class TEntityType>
    static void UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(
        std::map<IndexType, ModelPart::NodeType*>& rOutput,
        const std::map<IndexType, TEntityType*>& rInput);

    ///@}
};

///@}
}