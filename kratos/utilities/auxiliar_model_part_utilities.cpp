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
void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodes(
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

    mrModelPart.RemoveElement(ElementId, ThisIndex);

    // TODO: Add OMP
    for(std::size_t i=0; i<r_array_nodes.size(); ++i) {
        auto it_node = r_array_nodes.begin() + i;
        if (it_node->Is(IdentifierFlag)) {
            mrModelPart.RemoveNode(it_node->Id(), ThisIndex);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodes(Element& ThisElement, Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveElementAndBelongingNodes(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodes(Element::Pointer pThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingNodes(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodesFromAllLevels(IndexType ElementId, Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(*(mrModelPart.GetParentModelPart()));
        aux_utility.RemoveElementAndBelongingNodes(ElementId, IdentifierFlag, ThisIndex);
    } else {
        RemoveElementAndBelongingNodes(ElementId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodesFromAllLevels(Element& ThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingNodesFromAllLevels(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingNodesFromAllLevels(Element::Pointer pThisElement, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingNodesFromAllLevels(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongingNodes(Flags IdentifierFlag)
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
    }

    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongingNodesFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveElementsAndBelongingNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodes(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)
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

    mrModelPart.RemoveCondition(ConditionId, ThisIndex);

    // TODO: Add OMP
    for(std::size_t i=0; i<r_array_nodes.size(); ++i) {
        auto it_node = r_array_nodes.begin() + i;
        if (it_node->Is(IdentifierFlag)) {
            mrModelPart.RemoveNode(it_node->Id(), ThisIndex);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodes(Condition& ThisCondition, Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveConditionAndBelongingNodes(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodes(Condition::Pointer pThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingNodes(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodesFromAllLevels(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(*(mrModelPart.GetParentModelPart()));
        aux_utility.RemoveConditionAndBelongingNodes(ConditionId, IdentifierFlag, ThisIndex);
    } else {
        RemoveConditionAndBelongingNodes(ConditionId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodesFromAllLevels(Condition& ThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingNodesFromAllLevels(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingNodesFromAllLevels(Condition::Pointer pThisCondition, Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingNodesFromAllLevels(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongingNodes(Flags IdentifierFlag)
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
    }

    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongingNodesFromAllLevels(Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveConditionsAndBelongingNodes(IdentifierFlag);
}

}  // namespace Kratos.
