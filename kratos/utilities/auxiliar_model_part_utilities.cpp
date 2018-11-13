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
#include "utilities/auxiliar_model_part_utilities.h"

namespace Kratos
{
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
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(*(mrModelPart.GetParentModelPart()));
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
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(*(mrModelPart.GetParentModelPart()));
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
