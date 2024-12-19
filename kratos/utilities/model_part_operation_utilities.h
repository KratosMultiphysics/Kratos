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
#include <limits>
#include <algorithm>
#include <iterator>

// External includes

// Project includes
#include "includes/indexed_object.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "utilities/reduction_utilities.h"
#include "utilities/model_part_operator_utilities.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) ModelPartOperationUtilities {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    template<class TValueType>
    using IndexedSetType = std::set<typename TValueType::Pointer, std::less<IndexedObject>>;

    using CNodePointersType = std::vector<IndexType>;

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
     * The @ref TModelPartOperation check is done based on the ids of the entities.
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
        ModelPart::NodesContainerType output_nodes;
        ModelPart::ConditionsContainerType output_conditions;
        ModelPart::ElementsContainerType output_elements;

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

    template<class TModelPartOperation, class TContainer>
    static void AddEntities(
        TContainer& rOutputEntities,
        const TContainer& rMainEntityContainer,
        const std::vector<TContainer const *>& rSetOperands)
    {
        for (auto itr_main = rMainEntityContainer.begin(); itr_main != rMainEntityContainer.end(); ++itr_main) {
            TContainer temp;
            temp.reserve(rMainEntityContainer.size());
            if (TModelPartOperation::IsValid((*itr_main).Id(), rSetOperands)) {
                // since the rMainEntityContainer is a PVS, the following will have a valid
                // hint, hence following insert will use push_back in the back end.
                temp.insert(temp.end(), *(itr_main.base()));
            }

            // since the temp is also of type PVS, there won't be any sorting or making it unique.
            // it will insert as it is to correct places.
            rOutputEntities.insert(temp.begin(), temp.end());
        }
    }

    template<class TContainer>
    static void AddNodesFromEntities(
        ModelPart::NodesContainerType& rOutputNodes,
        const TContainer& rInputEntities);

    static void AddNeighbours(
        ModelPart::NodesContainerType& rOutputNodes,
        ModelPart::ConditionsContainerType& rOutputConditions,
        ModelPart::ElementsContainerType& rOutputElements,
        ModelPart& rMainModelPart);

    template<class TModelPartOperation>
    static void ModelPartOperation(
        ModelPart::NodesContainerType& rOutputNodes,
        ModelPart::ConditionsContainerType& rOutputConditions,
        ModelPart::ElementsContainerType& rOutputElements,
        ModelPart& rMainModelPart,
        const std::vector<ModelPart const *>& rOperands,
        const bool AddNeighbourEntities)
    {
        const IndexType number_of_operation_model_parts = rOperands.size();

        std::vector<ModelPart::NodesContainerType const *> set_operation_nodes(number_of_operation_model_parts);
        std::vector<ModelPart::ConditionsContainerType const *> set_operation_conditions(number_of_operation_model_parts);
        std::vector<ModelPart::ElementsContainerType const *> set_operation_elements(number_of_operation_model_parts);

        // get the nodes sets for set operations
        std::transform(rOperands.begin(), rOperands.end(), set_operation_nodes.begin(), [](const auto& pModelPart) { return &(pModelPart->Nodes()); });

        // get the condition sets for set operations
        std::transform(rOperands.begin(), rOperands.end(), set_operation_conditions.begin(), [](const auto& pModelPart) { return &(pModelPart->Conditions()); });

        // get the element sets for set operations
        std::transform(rOperands.begin(), rOperands.end(), set_operation_elements.begin(), [](const auto& pModelPart) { return &(pModelPart->Elements()); });

        AddEntities<TModelPartOperation, ModelPart::NodesContainerType>(rOutputNodes, rMainModelPart.Nodes(), set_operation_nodes);

        AddEntities<TModelPartOperation, ModelPart::ConditionsContainerType>(rOutputConditions, rMainModelPart.Conditions(), set_operation_conditions);

        AddEntities<TModelPartOperation, ModelPart::ElementsContainerType>(rOutputElements, rMainModelPart.Elements(), set_operation_elements);

        // now we have all the nodes to find and add neighbour entities.
        if (AddNeighbourEntities) {
            // we need to fill the boundary nodes for elements and conditions
            AddNodesFromEntities(rOutputNodes, rOutputConditions);
            AddNodesFromEntities(rOutputNodes, rOutputElements);

            // now add all the neighbours
            AddNeighbours(rOutputNodes, rOutputConditions, rOutputElements, rMainModelPart);
        }
    }

    static void FillOutputSubModelPart(
        ModelPart& rOutputSubModelPart,
        ModelPart& rMainModelPart,
        ModelPart::NodesContainerType& rOutputNodes,
        ModelPart::ConditionsContainerType& rOutputConditions,
        ModelPart::ElementsContainerType& rOutputElements);

    ///@}
};

///@}

} // namespace Kratos
