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

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "utilities/model_part_operator_utilities.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) ModelPartOperationUtilities {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using CNodePointersType = std::vector<ModelPart::NodeType const*>;

    template<class KeyType, class ValueType>
    using RangedKeyMapType = std::unordered_map<KeyType, ValueType, KeyHasherRange<KeyType>, KeyComparorRange<KeyType>>;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief Checks the validity of main model part and check model parts to be used within operation utilities.
     *
     * This method checks whether all the nodes, and geometries in @ref rCheckModelParts are present in the
     * @ref rMainModelPart.
     *
     * If @ref ThrowError is shown, a meaning full error is shown.
     *
     * @param rMainModelPart        Main model part.
     * @param rCheckModelParts      List of check model parts.
     * @param ThrowError            To throw an error or not.
     * @return true                 If the given model parts are valid.
     * @return false                If the given model parts are invalid.
     */
    static bool CheckValidityOfModelPartsForOperations(
        const ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rCheckModelParts,
        const bool ThrowError = false);

    /**
     * @brief Fill a sub model part with the specified operation.
     *
     * This method finds all the nodes, geometries in @ref rOperands in the @ref rMainModelPart
     * which satisfies the specified @ref TModelPartOperation and adds them to a given @ref rOutputSubModelPart.
     * The @ref TModelPartOperation check is done based on the memory locations of the entities.
     * The corresponding entities from @ref rMainModelPart are used to populate the sub model part.
     * @ref rOutputSubModelPart must be a sub model part of @ref rMainModelPart, but does not
     * necessarily have to be an immediate sub model part.
     *
     * If @ref AddNeighbourEntities is true, then all the neighbours of the newly created model part
     * is found from @ref rMainModelPart and added as well.
     *
     * Please make sure to check validity of model parts using @see CheckValidityOfModelPartsForOperations.
     *
     * @throws If @ref rOutputSubModelPart is not a sub model part of @ref rMainModelPart.
     * @throws If @ref rOutputSubModelPart is not empty.
     *
     * @tparam TModelPartOperation                  Type of the operation.
     * @param rOutputSubModelPart                   Output sub model part.
     * @param rMainModelPart                        Main Model part.
     * @param rOperands                             Model parts to to find satisfying entities w.r.t. @ref TModelPartOperation.
     * @param AddNeighbourEntities                  To add or not the neighbours.
     * @return ModelPart&                           Output sub model part.
     */
    template<class TModelPartOperation>
    static void FillSubModelPart(
        ModelPart& rOutputSubModelPart,
        const ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rOperands,
        const bool AddNeighbourEntities)
    {
        std::vector<ModelPart::NodeType*> output_nodes;
        std::vector<ModelPart::ConditionType*> output_conditions;
        std::vector<ModelPart::ElementType*> output_elements;

        // check whether the given rOutputSubModelPart is a sub model part of rMainModelPart
        // This is done based on pointers to avoid having same names in model parts
        // which are within two different models.
        ModelPart* p_parent_model_part = &rOutputSubModelPart.GetParentModelPart();
        while (p_parent_model_part != &rMainModelPart && p_parent_model_part != &p_parent_model_part->GetParentModelPart()) {
            p_parent_model_part = &p_parent_model_part->GetParentModelPart();
        }

        KRATOS_ERROR_IF_NOT(p_parent_model_part == &rMainModelPart)
            << "The given sub model part " << rOutputSubModelPart.FullName()
            << " is not a sub model part of " << rMainModelPart.FullName() << ".\n";

        // fill vectors
        ModelPartOperation<TModelPartOperation>(output_nodes, output_conditions, output_elements, *p_parent_model_part, rOperands, AddNeighbourEntities);

        // now fill the sub model part
        FillOutputSubModelPart(rOutputSubModelPart, *p_parent_model_part, output_nodes, output_conditions, output_elements);
    }

    /**
     * @brief Checks whether given model parts list has an intersection.
     *
     * This method checks whether @ref rIntersectionModelParts has intersections.
     *
     * The intersection is carried out by finding all the nodes, geometries in @ref rIntersectionModelParts.
     *
     * @param rIntersectionModelParts       Model parts to find intersection.
     * @return true                         If there are entities intersecting.
     * @return false                        If there aren't any entities intersecting.
     */
    static bool HasIntersection(const std::vector<ModelPart*>& rIntersectionModelParts);

    ///@}
private:
    ///@name Private static operations
    ///@{

    template<class TUnaryFunc>
    static void AddNodes(
        std::vector<ModelPart::NodeType*>& rNodesOutput,
        ModelPart::NodesContainerType& rNodesContainer,
        TUnaryFunc&& rIsValidEntity)
    {
        using p_entity_type = ModelPart::NodeType*;

        const auto& result = block_for_each<AccumReduction<p_entity_type>>(rNodesContainer, [&rIsValidEntity](auto& rMainEntity) -> p_entity_type {
            ModelPart::NodeType const* p_entity = &rMainEntity;
            if (rIsValidEntity(p_entity)) {
                return &rMainEntity;
            } else {
                return nullptr;
            }
        });

        std::for_each(result.begin(), result.end(), [&rNodesOutput](auto pEntity) {
            if (pEntity != nullptr) {
                rNodesOutput.push_back(pEntity);
            }
        });
    }

    template<class TContainerType, class TUnaryFunc>
    static void AddEntities(
        std::vector<typename TContainerType::value_type*>& rEntitiesOutput,
        TContainerType& rMainEntityContainer,
        TUnaryFunc&& rIsValidEntity)
    {
        using p_entity_type = typename TContainerType::value_type*;
        const auto& result = block_for_each<AccumReduction<p_entity_type>>(rMainEntityContainer, CNodePointersType(), [&rIsValidEntity](auto& rMainEntity, CNodePointersType& rTLS) -> p_entity_type {
            const auto& r_geometry = rMainEntity.GetGeometry();
            const IndexType number_of_nodes = r_geometry.size();

            if (rTLS.size() != number_of_nodes) {
                rTLS.resize(number_of_nodes);
            }

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                rTLS[i] = &r_geometry[i];
            }

            std::sort(rTLS.begin(), rTLS.end());

            if (rIsValidEntity(rTLS)) {
                return &rMainEntity;
            } else {
                return nullptr;
            }
        });

        std::for_each(result.begin(), result.end(), [&rEntitiesOutput](auto pEntity) {
            if (pEntity != nullptr) {
                rEntitiesOutput.push_back(pEntity);
            }
        });
    }

    template<class TModelPartOperation>
    static void ModelPartOperation(
        std::vector<ModelPart::NodeType*>& rOutputNodes,
        std::vector<ModelPart::ConditionType*>& rOutputConditions,
        std::vector<ModelPart::ElementType*>& rOutputElements,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const*>& rOperands,
        const bool AddNeighbourEntities)
    {
        const IndexType number_of_operation_model_parts = rOperands.size();

        std::vector<std::set<ModelPart::NodeType const*>> set_operation_node_sets(number_of_operation_model_parts);
        std::vector<std::set<CNodePointersType>> set_operation_condition_sets(number_of_operation_model_parts);
        std::vector<std::set<CNodePointersType>> set_operation_element_sets(number_of_operation_model_parts);

        // now iterate through model parts' conditions/elements containers and get available main model part entities.
        for (IndexType i = 0; i < number_of_operation_model_parts; ++i) {
            auto p_model_part = rOperands[i];

            // fill the set_operation nodes
            FillNodesPointerSet(set_operation_node_sets[i], p_model_part->Nodes());

            // fill the set_operation conditions
            FillNodePointersForEntities(set_operation_condition_sets[i], p_model_part->Conditions());

            // fill the set_operation elements
            FillNodePointersForEntities(set_operation_element_sets[i], p_model_part->Elements());
        }

        AddNodes(rOutputNodes, rMainModelPart.Nodes(), [&set_operation_node_sets](auto& rEntity) {
            return TModelPartOperation::IsValid(rEntity, set_operation_node_sets);
        });

        AddEntities(rOutputConditions, rMainModelPart.Conditions(), [&set_operation_condition_sets](auto& rEntity) {
            return TModelPartOperation::IsValid(rEntity, set_operation_condition_sets);
        });

        AddEntities(rOutputElements, rMainModelPart.Elements(), [&set_operation_element_sets](auto& rEntity) {
            return TModelPartOperation::IsValid(rEntity, set_operation_element_sets);
        });

        // now we have all the nodes to find and add neighbour entities.
        if (AddNeighbourEntities) {
            // we need to fill the boundary nodes for elements and conditions
            FillNodesFromEntities<ModelPart::ConditionType>(rOutputNodes, rOutputConditions.begin(), rOutputConditions.end());
            FillNodesFromEntities<ModelPart::ElementType>(rOutputNodes, rOutputElements.begin(), rOutputElements.end());

            // now add all the neighbours
            AddNeighbours(rOutputNodes, rOutputConditions, rOutputElements, rMainModelPart);
        }
    }

    static void FillNodesPointerSet(
        std::set<ModelPart::NodeType const*>& rOutput,
        const ModelPart::NodesContainerType& rNodes);

    template <class TContainerType>
    static void FillNodePointersForEntities(
        std::set<CNodePointersType>& rOutput,
        const TContainerType& rEntityContainer);

    static void FillCommonNodes(
        ModelPart::NodesContainerType& rOutput,
        ModelPart::NodesContainerType& rMainNodes,
        const std::set<ModelPart::NodeType const*>& rNodesSet);

    template<class TEntityType>
    static void FillNodesFromEntities(
        std::vector<ModelPart::NodeType*>& rOutput,
        typename std::vector<TEntityType*>::iterator pEntityBegin,
        typename std::vector<TEntityType*>::iterator pEntityEnd);

    template<class TContainerType>
    static void FindNeighbourEntities(
        std::vector<typename TContainerType::value_type*>& rOutputEntities,
        const Flags& rNodeSelectionFlag,
        TContainerType& rMainEntities);

    static void AddNeighbours(
        std::vector<ModelPart::NodeType*>& rOutputNodes,
        std::vector<ModelPart::ConditionType*>& rOutputConditions,
        std::vector<ModelPart::ElementType*>& rOutputElements,
        ModelPart& rMainModelPart);

    static void FillWithMainNodesFromSearchNodes(
        ModelPart::NodesContainerType& rOutput,
        ModelPart::NodesContainerType& rMainNodes,
        ModelPart::NodesContainerType& rSearchedNodes);

    static void SetCommunicator(
        ModelPart& rOutputModelPart,
        ModelPart& rMainModelPart);

    static void FillOutputSubModelPart(
        ModelPart& rOutputSubModelPart,
        ModelPart& rMainModelPart,
        std::vector<ModelPart::NodeType*>& rOutputNodes,
        std::vector<ModelPart::ConditionType*>& rOutputConditions,
        std::vector<ModelPart::ElementType*>& rOutputElements);

    static void CheckNodes(
        std::vector<IndexType>& rNodeIdsWithIssues,
        const ModelPart::NodesContainerType& rCheckNodes,
        const std::set<ModelPart::NodeType const*>& rMainNodes);

    template<class TContainerType>
    static void CheckEntities(
        std::vector<IndexType>& rEntityIdsWithIssues,
        const TContainerType& rCheckEntities,
        const std::set<CNodePointersType>& rMainEntities);

    ///@}
};

///@}

} // namespace Kratos
