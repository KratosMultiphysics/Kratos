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
     * @brief Get Newly Created or Existing Model Parts List With Common Reference Entities Between Reference List And Examined List
     *
     * This method returns newly created or gets an existing sub-model part in each model part in rReferenceModelParts with entities
     * common between that model part and all the model parts in rExaminedModelPartsList with following
     * entities.
     *      1. If AreNodesConsidered is made to true, then common nodes are found and added to the resulting
     *         sub-model part.
     *      2. If AreConditionsConsidered is made to true, then common conditions are found and added to the resulting
     *         sub-model part. All the nodes belonging to added conditions are also added to the sub-model part.
     *      3. If AreElementsConsidered is made to true, then common elements are found and added to the resulting
     *         sub-model part. All the nodes belonging to added elements are also added to the sub-model part.
     *      4. If AreParentsConsidered is made to true, then first common nodes are found. Then neighbour conditions
     *         and elements for those common nodes are found and added to the sub-model part. All nodes belonging to
     *         added conditions and elements are also added to the sub-model part.
     *
     * The created sub-model part is only populated with the entities from the refrence model part. Hence, to find common
     * conditions and elements, each entity from each model part in rExaminedModelPartsList is searched within each of the
     * model parts in rReferenceModelParts.
     *
     * This method creates the sub-model parts with unique name generated for all the input arguments. If it finds an existing
     * sub-model part with the same name, then it is retrieved. Hence this method only creates these model parts with expensive operations
     * only one. If it is required to force create, then the sub-model parts should be removed. All the sub-model parts created
     * by this method have names starting with "<OPTIMIZATION_APP_AUTO>" so it is easier to search and remove them if required.
     *
     * This method keeps the communicator types, process info, properties and tables as same as the rReferenceModelParts.
     * The local mesh, interface mesh and ghost meshes are properly updated. Hence this method is compatible with OpenMP
     * and MPI.
     *
     * This method does not create or destroy nodes, conditions and elements.
     *
     * @param rExaminedModelPartsList       List of input model parts where the entities are searched between.
     * @param rReferenceModelParts          List of reference model parts where common entities are added to output model parts.
     * @param AreNodesConsidered            If true, common nodes from reference model parts are added to output model parts.
     * @param AreConditionsConsidered       If true, common conditions from reference model parts are added to output model parts.
     * @param AreElementsConsidered         If true, common elements from reference model parts are added to output model parts.
     * @param AreParentsConsidered          If true, neighbour conditions and elements from reference model parts for common nodes are added to output model parts.
     * @param EchoLevel                     Echo level for printing info.
     * @return std::vector<ModelPart*>      List of output model parts.
     */
    static std::vector<ModelPart*> GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
        const std::vector<ModelPart*>& rExaminedModelPartsList,
        const std::vector<ModelPart*>& rReferenceModelParts,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered,
        const bool AreParentsConsidered,
        const IndexType EchoLevel = 0);

    ///@}
private:
    ///@name Private classes
    ///@{

    template<class TEntityType, class TMapValueType>
    class ContainerEntityMapReduction
    {
    public:
        using return_type = std::map<IndexType, TMapValueType>;
        using value_type = std::vector<std::pair<IndexType, EntityPointerType<TEntityType>>>;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const;

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type& rValue);

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(ContainerEntityMapReduction<TEntityType, TMapValueType>& rOther);
    };

    ///@}
    ///@name Private static operations
    ///@{

    static std::string GetCombinedModelPartsName(
        const std::string& rPrefix,
        const std::vector<ModelPart*>& rModelParts);

    static std::string GetSensitivityComputationModelPartsInfo(
        const std::vector<ModelPart*>& rSensitivityModelParts,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered,
        const bool AreParentsConsidered);

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

    static void CreateModelPartWithCommonReferenceEntitiesBetweenReferenceAndExamined(
        const std::string& rOutputModelPartName,
        const std::vector<ModelPart*>& rExaminedModelPartsList,
        ModelPart& rReferenceModelPart,
        const bool AreNodesConsidered,
        const bool AreConditionsConsidered,
        const bool AreElementsConsidered,
        const bool AreParentsConsidered,
        const IndexType EchoLevel = 0);

    ///@}
};

///@}
}