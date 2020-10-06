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
    if (mrModelPart.IsSubModelPart()) {
        mrModelPart.GetParentModelPart().AddElement(pNewElement, ThisIndex);
        mrModelPart.GetMesh(ThisIndex).AddElement(pNewElement);
    } else {
        auto existing_element_it = mrModelPart.GetMesh(ThisIndex).Elements().find(pNewElement->Id());
        if( existing_element_it == mrModelPart.GetMesh(ThisIndex).ElementsEnd()) { // Element did not exist
            mrModelPart.GetMesh(ThisIndex).AddElement(pNewElement);
        } else { // Element did exist already
            KRATOS_ERROR_IF(&(*existing_element_it) != (pNewElement.get()))//check if the pointee coincides
                << "Attempting to add pNewElement with Id :" << pNewElement->Id() << ", unfortunately a (different) element with the same Id already exists" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::AddElementsWithNodes(
    std::vector<IndexType> const& rElementIds,
    IndexType ThisIndex
    )
{
    KRATOS_TRY
    if(mrModelPart.IsSubModelPart()) { // Does nothing if we are on the top model part
        // Obtain from the root model part the corresponding list of nodes
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        ModelPart::ElementsContainerType aux;
        aux.reserve(rElementIds.size());
        std::unordered_set<IndexType> set_of_nodes;
        for(IndexType i=0; i<rElementIds.size(); ++i) {
            auto it_elem = p_root_model_part->Elements().find(rElementIds[i]);
            if(it_elem!=p_root_model_part->ElementsEnd()) {
                aux.push_back(*(it_elem.base()));
                const auto& r_geom = it_elem->GetGeometry();
                for (IndexType i = 0; i < r_geom.size(); ++i) {
                    set_of_nodes.insert(r_geom[i].Id());
                }
            } else {
                KRATOS_ERROR << "The element wit_elemh Id " << rElementIds[i] << " does not exist in the root model part";
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_nodes.begin(), set_of_nodes.end());
        mrModelPart.AddNodes(list_of_nodes);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_elem = aux.begin(); it_elem!=aux.end(); it_elem++) {
                p_current_part->Elements().push_back( *(it_elem.base()) );
            }
            p_current_part->AddNodes(list_of_nodes);

            p_current_part->Elements().Unique();
            p_current_part = &(p_current_part->GetParentModelPart());
        }
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
    if (mrModelPart.IsSubModelPart()) {
        mrModelPart.GetParentModelPart().AddCondition(pNewCondition, ThisIndex);
        mrModelPart.GetMesh(ThisIndex).AddCondition(pNewCondition);
    } else {
        auto existing_condition_it = mrModelPart.GetMesh(ThisIndex).Conditions().find(pNewCondition->Id());
        if( existing_condition_it == mrModelPart.GetMesh(ThisIndex).ConditionsEnd()) { // Condition did not exist
            mrModelPart.GetMesh(ThisIndex).AddCondition(pNewCondition);
        } else { // Condition did exist already
            KRATOS_ERROR_IF(&(*existing_condition_it) != (pNewCondition.get()))//check if the pointee coincides
                << "Attempting to add pNewCondition with Id :" << pNewCondition->Id() << ", unfortunately a (different) condition with the same Id already exists" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::AddConditionsWithNodes(
    std::vector<IndexType> const& rConditionIds,
    IndexType ThisIndex
    )
{
    KRATOS_TRY
    if(mrModelPart.IsSubModelPart()) { // Does nothing if we are on the top model part
        // Obtain from the root model part the corresponding list of nodes
        ModelPart* p_root_model_part = &mrModelPart.GetRootModelPart();
        ModelPart::ConditionsContainerType aux;
        aux.reserve(rConditionIds.size());
        std::unordered_set<IndexType> set_of_nodes;
        for(IndexType i=0; i<rConditionIds.size(); ++i) {
            auto it_cond = p_root_model_part->Conditions().find(rConditionIds[i]);
            if(it_cond!=p_root_model_part->ConditionsEnd()) {
                aux.push_back(*(it_cond.base()));
                const auto& r_geom = it_cond->GetGeometry();
                for (IndexType i = 0; i < r_geom.size(); ++i) {
                    set_of_nodes.insert(r_geom[i].Id());
                }
            } else {
                KRATOS_ERROR << "The condition wit_condh Id " << rConditionIds[i] << " does not exist in the root model part";
            }
        }

        // Adding nodes
        std::vector<IndexType> list_of_nodes;
        list_of_nodes.insert(list_of_nodes.end(), set_of_nodes.begin(), set_of_nodes.end());
        mrModelPart.AddNodes(list_of_nodes);

        // Add to all of the leaves
        ModelPart* p_current_part = &mrModelPart;
        while(p_current_part->IsSubModelPart()) {
            for(auto it_cond = aux.begin(); it_cond!=aux.end(); it_cond++) {
                p_current_part->Conditions().push_back( *(it_cond.base()) );
            }
            p_current_part->AddNodes(list_of_nodes);

            p_current_part->Conditions().Unique();
            p_current_part = &(p_current_part->GetParentModelPart());
        }
    }
    KRATOS_CATCH("");
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
    std::unordered_set<Properties::Pointer, IndexedObjecPointertHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> list_of_properties;

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
        std::unordered_set<Properties::Pointer, IndexedObjecPointertHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> buffer_list_of_properties;

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
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(r_array_nodes.size()); ++i) {
        auto it_node = r_array_nodes.begin() + i;
        it_node->Set(IdentifierFlag, true);
    }

    // TODO: Add OMP
    for (auto& elem : mrModelPart.Elements(ThisIndex)) {
        if (elem.Id() != ElementId) {
            for (auto& node : elem.GetGeometry()) {
                node.Set(IdentifierFlag, false);
            }
        }
    }
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
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(r_array_nodes.size()); ++i) {
            auto it_node = r_array_nodes.begin() + i;
            it_node->Set(IdentifierFlag, true);
        }

        // TODO: Add OMP
        for (auto& elem : i_mesh->Elements()) {
            if (elem.IsNot(IdentifierFlag)) {
                for (auto& node : elem.GetGeometry()) {
                    node.Set(IdentifierFlag, false);
                }
            }
        }

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
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(r_array_nodes.size()); ++i) {
        auto it_node = r_array_nodes.begin() + i;
        it_node->Set(IdentifierFlag, true);
    }

    // TODO: Add OMP
    for (auto& cond : mrModelPart.Conditions(ThisIndex)) {
        if (cond.Id() != ConditionId) {
            for (auto& node : cond.GetGeometry()) {
                node.Set(IdentifierFlag, false);
            }
        }
    }
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
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(r_array_nodes.size()); ++i) {
            auto it_node = r_array_nodes.begin() + i;
            it_node->Set(IdentifierFlag, true);
        }

        // TODO: Add OMP
        for (auto& cond : i_mesh->Conditions()) {
            if (cond.IsNot(IdentifierFlag)) {
                for (auto& node : cond.GetGeometry()) {
                    node.Set(IdentifierFlag, false);
                }
            }
        }

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
