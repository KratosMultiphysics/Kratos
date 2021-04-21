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
#include <unordered_set>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/parallel_utilities.h"
#include "variable_utils.h"

namespace Kratos
{
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
