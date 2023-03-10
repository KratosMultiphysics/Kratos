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

    template<class TEntityType>
    using EntityPointerType = typename TEntityType::Pointer;

    template<class TContainerType>
    using ContainerEntityValueType = typename TContainerType::value_type;

    template<class TContainerType>
    using ContainerEntityPointerType = typename ContainerEntityValueType<TContainerType>::Pointer;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief Get the Sensitivity Model Part For Adjoint Sensitivitiy Computation.
     *
     * This does not create or destroy nodes, conditions or elements. This is compatible with OpenMP. In the case of MPI,
     * this will not add neighbour entities from other ranks. Hence, it is required to use Assembling methods for nodal
     * sensitivity computations at the end.
     *
     * This method can be used to obtain combined sensitivity model part as explained below.
     *      1. If AreSensitivityEntityParentsConsidered is made to true, first all nodes of each sensitivity model part
     *         in rSensitivityModelParts will be populated with neighbour elements and conditions from the rAnalysisModelPart.
     *         Then these neighbours are added to the resulting model part. This is useful when adjoint sensitivity
     *         analysis is done on the nodal quantities such as SHAPE_SENSITIVITY, because the nodes requires residual contributions
     *         from all the neighbouring entities. When adding these entities, corresponding nodal container is also filled with
     *         relevant nodes. This requires to have common nodes between rAnalysisModelPart and rSensitivityModelParts.
     *
     *      2. If AreSensitivityEntitesConsidered is made to true, then all the conditions and elements in each sensitivity model part
     *         is added to the resulting model part if they have the same root model part. Otherwise check ForceFindSensitivityEntitiesInAnalysisModelPart = true.
     *         When adding these entities, corresponding nodal container is also filled with relevant nodes. This is useful in the case
     *         when sensitivities are computed using adjoint approach for condition and/or element data values such as DENSITY_SENSITIVITY.
     *         This requires having common conditions and/or elements between rAnalysisModelPart and rSensitivityModelPart.
     *
     *      3. If ForceFindSensitivityEntitiesInAnalysisModelPart is made to true, firstly, condition and element counterparts in rAnalysisModelPart is searched by
     *         matching each entities nodal configurations from each and every entity in each model part in rSensitivityModelParts. Then these condition and
     *         element counterparts in rAnalysisModelPart is added to the resulting model part.  When adding these entities, corresponding nodal container
     *         is also filled with relevant nodes. This is useful in the case when the root model parts differ in the rAnalysisModelPart and rSensitivityModelParts, but
     *         they have common interfaces. This requires having common nodes between rAnalysisModelPart and rSensitivityModelParts and having conditions and/or
     *         elements with the same nodal configuration between rAnalysisModelPart and rSensitivityModelParts.
     *
     * This method creates a root model part with a unique name based on the input arguments. These model part names always starts with "<OPTIMIZATION_APP_AUTO>..."
     * So, they can be removed if required. (such as when remeshing is done after few iterations of optimization.)
     *
     * If the model part is already created, then it will be returning the model part without doing the heavy model part creation steps. So calling this method
     * many times with the same input arguments will only create the resulting model part once. Thereafter, the already created model part is returned.
     *
     * @param rSensitivityModelParts                            List of sensitivity model part pointers.
     * @param rAnalysisModelPart                                Analysis model part.
     * @param AreSensitivityEntityParentsConsidered             true if it is required to have parent conditions and elements from rAnalysisModelPart for nodes in rSensitivityModelParts.
     * @param AreSensitivityEntitesConsidered                   true if it is required to have conditions and elements from rAnalysisModelPart for conditions and elements in rSensitivityModelParts.
     * @param ForceFindSensitivityEntitiesInAnalysisModelPart   true if it is required to have conditions and elements from rAnalysisModelPart for conditions and elements in rSensitivityModelParts when the root model parts does not match.
     * @param EchoLevel                                         Echo level of the operations outputs.
     * @return ModelPart&                                       Resulting model part to carry out sensitivity analysis using the adjoint approach.
     */
    static ModelPart& GetSensitivityModelPartForAdjointSensitivities(
        const std::vector<ModelPart*>& rSensitivityModelParts,
        ModelPart& rAnalysisModelPart,
        const bool AreSensitivityEntityParentsConsidered,
        const bool AreSensitivityEntitesConsidered,
        const bool ForceFindSensitivityEntitiesInAnalysisModelPart = false,
        const IndexType EchoLevel = 0);

    /**
     * @brief Get the Sensitivity Model Part For Direct Sensitivity Computation.
     *
     * This does not create or destroy nodes, conditions or elements. This is compatible with OpenMP.In the case of MPI,
     * this will not add neighbour entities from other ranks. Hence, it is required to use Assembling methods for nodal
     * sensitivity computations at the end.
     *
     * This method can be used to obtain combined sensitivity model part as explained below.
     *      1. If AreNodesConsidered is made to true, then common nodes between each model part in rSensitivityModelParts and each model part
     *         in rEvaluatedModelParts are found, and then they are added to the resulting model part.
     *      2. If AreConditionsConsidered is made to true, then common conditions between each model part in rSensitivityModelParts and each model part
     *         in rEvaluatedModelParts are found, and then they are added to the resulting model part. Thereafter, the relevant nodes are also added
     *         to the resulting model part.
     *      3. If AreConditionsConsidered is made to true, then common elements between each model part in rSensitivityModelParts and each model part
     *         in rEvaluatedModelParts are found, and then they are added to the resulting model part. Thereafter, the relevant nodes are also added
     *         to the resulting model part.
     *
     * This method creates a root model part with a unique name based on the input arguments. These model part names always starts with "<OPTIMIZATION_APP_AUTO>..."
     * So, they can be removed if required. (such as when remeshing is done after few iterations of optimization.)
     *
     * If the model part is already created, then it will be returning the model part without doing the heavy model part creation steps. So calling this method
     * many times with the same input arguments will only create the resulting model part once. Thereafter, the already created model part is returned.
     *
     * @param rSensitivityModelParts        List of sensitivity model part pointers.
     * @param rEvaluatedModelParts          List of evaluated model part pointers.
     * @param AreNodesConsidered            true if it is required to add common nodes to the resulting model part.
     * @param AreConditionsConsidered       true if it is required to add common conditions to the resulting model part.
     * @param AreElementsConsidered         true if it is required to add common elements to the resulting model part.
     * @param EchoLevel                     Echo level of the operations outputs.
     * @return ModelPart&                   Resulting model part to carry out sensitivity analysis using the direct approach.
     */
    static ModelPart& GetSensitivityModelPartForDirectSensitivities(
        const std::vector<ModelPart*>& rSensitivityModelParts,
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered,
        const IndexType EchoLevel = 0);

    ///@}
private:
    ///@name Private classes
    ///@{

    template<class EntityType>
    class ContainerEntityMapReduction
    {
    public:
        using return_type = std::map<IndexType, EntityPointerType<EntityType>>;
        using value_type = std::vector<std::pair<IndexType, EntityPointerType<EntityType>>>;

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
        std::map<IndexType, std::vector<ContainerEntityPointerType<TContainerType>>>& rOutput,
        TContainerType& rContainer,
        const Flags& rFlag,
        const bool FlagValue = true);

    template<class TEntityPointerType>
    static void UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(
        std::map<IndexType, TEntityPointerType>& rOutput,
        const std::map<IndexType, std::vector<TEntityPointerType>>& rNodeIdNeighbourEntityPtrsMap,
        const ModelPart::NodesContainerType& rNodes);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromEntityContainer(
        std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
        TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateNodeIdsEntityPtrMapFromEntityContainer(
        std::map<NodeIdsType, ContainerEntityPointerType<TContainerType>>& rOutput,
        TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(
        std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
        const std::map<NodeIdsType, ContainerEntityPointerType<TContainerType>>& rNodeIdsEntityPtrMap,
        const TContainerType& rContainer);

    template<class TContainerType>
    static void UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(
        std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
        TContainerType& rContainer,
        const Flags& rFlag,
        const bool FlagValue = true);

    template<class TEntityPointerType>
    static void UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(
        std::map<IndexType, ModelPart::NodeType::Pointer>& rOutput,
        const std::map<IndexType, TEntityPointerType>& rInput);

    ///@}
};

///@}
}