//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "includes/key_hash.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/parallel_utilities.h"
#include "variable_utils.h"

namespace Kratos
{
void AuxiliarModelPartUtilities::AddElementWithNodes(
    Element::Pointer pNewElement,
    IndexType ThisIndex
    )
{
    const auto& r_geom = pNewElement->GetGeometry();
    std::vector<IndexType> list_of_nodes(r_geom.size());
    for (IndexType i = 0; i < r_geom.size(); ++i) {
        list_of_nodes[i] = r_geom[i].Id();
    }
    mrModelPart.AddNodes(list_of_nodes);
    mrModelPart.AddElement(pNewElement);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::AddElementsWithNodes(
    const std::vector<IndexType>& rElementIds,
    IndexType ThisIndex
    )
{
    KRATOS_TRY
    mrModelPart.AddElements(rElementIds, ThisIndex);
    if(mrModelPart.IsSubModelPart()) { // Does nothing if we are on the root model part, the root model part already contains all the nodes
        // Obtain from the root model part the corresponding list of nodes
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::unordered_set<IndexType> set_of_node_ids;
        for(IndexType i=0; i<rElementIds.size(); ++i) {
            auto it_elem = p_root_model_part->Elements().find(rElementIds[i]);
            if(it_elem!=p_root_model_part->ElementsEnd()) {
                const auto& r_geom = it_elem->GetGeometry();
                for (IndexType j = 0; j < r_geom.size(); ++j) {
                    set_of_node_ids.insert(r_geom[j].Id());
                }
            } else {
                KRATOS_ERROR << "The element with Id " << rElementIds[i] << " does not exist in the root model part";
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_node_ids.begin(), set_of_node_ids.end());
        mrModelPart.AddNodes(list_of_nodes);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            p_current_part->AddNodes(list_of_nodes);
            p_current_part = &(p_current_part->GetParentModelPart());
        }
    } else {
        KRATOS_WARNING("AuxiliarModelPartUtilities") << "Does nothing appart of adding the elements as we are on the root model part, the root model part already contains all the nodes" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::AddConditionWithNodes(
    Condition::Pointer pNewCondition,
    IndexType ThisIndex
    )
{
    const auto& r_geom = pNewCondition->GetGeometry();
    std::vector<IndexType> list_of_nodes(r_geom.size());
    for (IndexType i = 0; i < r_geom.size(); ++i) {
        list_of_nodes[i] = r_geom[i].Id();
    }
    mrModelPart.AddNodes(list_of_nodes);
    mrModelPart.AddCondition(pNewCondition);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::AddConditionsWithNodes(
    const std::vector<IndexType>& rConditionIds,
    IndexType ThisIndex
    )
{
    KRATOS_TRY
    mrModelPart.AddConditions(rConditionIds, ThisIndex);
    if(mrModelPart.IsSubModelPart()) { // Does nothing if we are on the top model part
        // Obtain from the root model part the corresponding list of nodes
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        std::unordered_set<IndexType> set_of_node_ids;
        for(IndexType i=0; i<rConditionIds.size(); ++i) {
            auto it_cond = p_root_model_part->Conditions().find(rConditionIds[i]);
            if(it_cond!=p_root_model_part->ConditionsEnd()) {
                const auto& r_geom = it_cond->GetGeometry();
                for (IndexType j = 0; j < r_geom.size(); ++j) {
                    set_of_node_ids.insert(r_geom[j].Id());
                }
            } else {
                KRATOS_ERROR << "The condition with Id " << rConditionIds[i] << " does not exist in the root model part";
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_node_ids.begin(), set_of_node_ids.end());
        mrModelPart.AddNodes(list_of_nodes);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            p_current_part->AddNodes(list_of_nodes);
            p_current_part = &(p_current_part->GetParentModelPart());
        }
    } else {
        KRATOS_WARNING("AuxiliarModelPartUtilities") << "Does nothing appart of adding the conditions as we are on the root model part, the root model part already contains all the nodes" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::CopySubModelPartStructure(const ModelPart& rModelPartToCopyFromIt, ModelPart& rModelPartToCopyIntoIt)
{
    for (auto& r_sub_model_part : rModelPartToCopyFromIt.SubModelParts()) {
        auto& r_new_sub_model_part = rModelPartToCopyIntoIt.CreateSubModelPart(r_sub_model_part.Name());
        if (r_sub_model_part.NumberOfSubModelParts() > 0) {
            CopySubModelPartStructure(r_sub_model_part, r_new_sub_model_part);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RecursiveEnsureModelPartOwnsProperties(const bool RemovePreviousProperties)
{
    // First we do in this model part
    EnsureModelPartOwnsProperties(RemovePreviousProperties);

    // Now we do in submodelparts
    for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
        AuxiliarModelPartUtilities(r_sub_model_part).RecursiveEnsureModelPartOwnsProperties(RemovePreviousProperties);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::EnsureModelPartOwnsProperties(const bool RemovePreviousProperties)
{
    // First we clear the properties if we want so
    if (RemovePreviousProperties) {
        mrModelPart.GetMesh(0).pProperties()->clear();
    }

    // The list of properties
    std::unordered_set<Properties::Pointer, IndexedObjectPointerHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> list_of_properties;

    // Iterating over the elements
    auto& r_elements_array = mrModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    // Iterating over the conditions
    auto& r_conditions_array = mrModelPart.Conditions();
    const auto it_cond_begin= r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());

    #pragma omp parallel
    {
        // The list of properties
        std::unordered_set<Properties::Pointer, IndexedObjectPointerHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> buffer_list_of_properties;

        #pragma omp for schedule(dynamic, 512) nowait
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            Properties::Pointer p_prop = it_elem->pGetProperties();

            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        #pragma omp for schedule(dynamic, 512) nowait
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            Properties::Pointer p_prop = it_cond->pGetProperties();
            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            list_of_properties.insert(buffer_list_of_properties.begin(),buffer_list_of_properties.end());
        }
    }

    // Add properties to respective model parts
    for (auto p_prop : list_of_properties) {
        if (!mrModelPart.HasProperties(p_prop->Id())) {
            mrModelPart.AddProperties(p_prop);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(
    IndexType ElementId, Flags IdentifierFlag, IndexType ThisIndex)
{
    auto& r_array_nodes = mrModelPart.Nodes(ThisIndex);

    VariableUtils().SetFlag(IdentifierFlag, true, r_array_nodes);

    block_for_each(
        mrModelPart.Elements(ThisIndex), [&IdentifierFlag,ElementId]( ModelPart::ElementType& rElement )
        {
            if (rElement.Id() != ElementId)
                for (auto& r_node : rElement.GetGeometry())
                    r_node.Set(IdentifierFlag, false);
        }
    );

    bool condition_to_remove;
    for (auto& cond : mrModelPart.Conditions(ThisIndex)) {
        condition_to_remove = true;
        for (auto& node : cond.GetGeometry()) {
            if (node.IsNot(IdentifierFlag)) {
                condition_to_remove = false;
                break;
            }
            if (condition_to_remove) cond.Set(IdentifierFlag);
        }
    }

    mrModelPart.RemoveElement(ElementId, ThisIndex);
    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(Element& ThisElement, Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveElementAndBelongings(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(Element::Pointer pThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongings(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(IndexType ElementId, Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(mrModelPart.GetParentModelPart());
        aux_utility.RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);
    } else {
        RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(Element& ThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingsFromAllLevels(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(Element::Pointer pThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingsFromAllLevels(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongings(Flags IdentifierFlag)
{
    //loop over all the meshes
    VariableUtils variable_utils;
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        variable_utils.SetFlag(IdentifierFlag, true, r_array_nodes);

        block_for_each(
            i_mesh->Elements(),
            [&IdentifierFlag](Element& rElement)
            {
                if (rElement.IsNot(IdentifierFlag))
                    for (auto& r_node : rElement.GetGeometry())
                        r_node.Set(IdentifierFlag, false);
            }
        );

        bool condition_to_remove;
        for (auto& cond : i_mesh->Conditions()) {
            condition_to_remove = true;
            for (auto& node : cond.GetGeometry()) {
                if (node.IsNot(IdentifierFlag)) {
                    condition_to_remove = false;
                    break;
                }
                if (condition_to_remove) cond.Set(IdentifierFlag);
            }
        }
    }

    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongingsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveElementsAndBelongings(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)
{
    auto& r_array_nodes = mrModelPart.Nodes(ThisIndex);
    VariableUtils().SetFlag(IdentifierFlag, true, r_array_nodes);

    block_for_each(
        mrModelPart.Conditions(ThisIndex),
        [&IdentifierFlag,ConditionId](ModelPart::ConditionType& rCondition)
        {
            if (rCondition.Id() != ConditionId)
                for (auto& r_node : rCondition.GetGeometry())
                    r_node.Set(IdentifierFlag, false);
        }
    );
    bool element_to_remove;
    for (auto& elem : mrModelPart.Elements(ThisIndex)) {
        element_to_remove = true;
        for (auto& node : elem.GetGeometry()) {
            if (node.IsNot(IdentifierFlag)) {
                element_to_remove = false;
                break;
            }
            if (element_to_remove) elem.Set(IdentifierFlag);
        }
    }

    mrModelPart.RemoveCondition(ConditionId, ThisIndex);
    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(Condition& ThisCondition, Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveConditionAndBelongings(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(Condition::Pointer pThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongings(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(mrModelPart.GetParentModelPart());
        aux_utility.RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);
    } else {
        RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(Condition& ThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingsFromAllLevels(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(Condition::Pointer pThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingsFromAllLevels(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongings(Flags IdentifierFlag)
{
    //loop over all the meshes
    VariableUtils variable_utils;
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        variable_utils.SetFlag(IdentifierFlag, true, r_array_nodes);

        block_for_each(
            i_mesh->Conditions(),
            [&IdentifierFlag](ModelPart::ConditionType& rCondition)
            {
                if (rCondition.IsNot(IdentifierFlag))
                    for (auto& r_node : rCondition.GetGeometry())
                        r_node.Set(IdentifierFlag, false);
            }
        );

        bool element_to_remove;
        for (auto& elem : i_mesh->Elements()) {
            element_to_remove = true;
            for (auto& node : elem.GetGeometry()) {
                if (node.IsNot(IdentifierFlag)) {
                    element_to_remove = false;
                    break;
                }
                if (element_to_remove) elem.Set(IdentifierFlag);
            }
        }
    }

    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongingsFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveConditionsAndBelongings(IdentifierFlag);
}

}  // namespace Kratos.
