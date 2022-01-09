//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "containers/pointer_vector_set.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "rom_auxiliary_utilities.h"

namespace Kratos
{

void RomAuxiliaryUtilities::SetHRomComputingModelPart(
    const Parameters HRomWeights,
    const ModelPart& rOriginModelPart,
    ModelPart& rHRomComputingModelPart)
{
    // Ensure that the provided destination model part is empty
    rHRomComputingModelPart.Clear();

    // Auxiliary containers to save the entities involved in the HROM mesh
    // Note that we use a set for the nodes to make sure that the same node is not added by more than one element/condition
    NodesPointerSet hrom_nodes_set;
    std::vector<Element::Pointer> hrom_elems_vect;
    std::vector<Condition::Pointer> hrom_conds_vect;

    const auto& r_elem_weights = HRomWeights["Elements"];
    hrom_elems_vect.reserve(rOriginModelPart.NumberOfElements());
    for (auto it = r_elem_weights.begin(); it != r_elem_weights.end(); ++it) {
        // Get element from origin mesh
        const IndexType elem_id = stoi(it.name());
        const auto p_elem = rOriginModelPart.pGetElement(elem_id);

        // Add the element to the auxiliary container and to the main HROM model part
        hrom_elems_vect.push_back(p_elem);
        rHRomComputingModelPart.AddElement(p_elem);

        // Add the element nodes to the auxiliary set and to the main HROM model part
        const auto& r_geom = p_elem->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            NodeType::Pointer p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
            rHRomComputingModelPart.AddNode(p_node);
        }
    }
    hrom_elems_vect.shrink_to_fit();

    const auto& r_cond_weights = HRomWeights["Conditions"];
    hrom_conds_vect.reserve(rOriginModelPart.NumberOfConditions());
    for (auto it = r_cond_weights.begin(); it != r_cond_weights.end(); ++it) {
        // Get the condition from origin mesh
        const IndexType cond_id = stoi(it.name());
        auto p_cond = rOriginModelPart.pGetCondition(cond_id);

        // Add the condition to the auxiliary container and to the main HROM model part
        hrom_conds_vect.push_back(p_cond);
        rHRomComputingModelPart.AddCondition(p_cond);

        // Add the condition nodes to the auxiliary set and to the main HROM model part
        const auto& r_geom = p_cond->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
            rHRomComputingModelPart.AddNode(p_node);
        }
    }
    hrom_conds_vect.shrink_to_fit();

    //TODO: ADD MPC'S

    // Create and fill the HROM calculation sub model parts
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(hrom_nodes_set, hrom_elems_vect, hrom_conds_vect, r_orig_sub_mp, rHRomComputingModelPart);
    }
}

void RomAuxiliaryUtilities::RecursiveHRomModelPartCreation(
    const NodesPointerSet& rNodesSet,
    const std::vector<Element::Pointer>& rElementsVector,
    const std::vector<Condition::Pointer>& rConditionsVector,
    const ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    // Emulate the origin submodelpart hierarchy
    auto& r_hrom_sub_mp = rDestinationModelPart.CreateSubModelPart(rOriginModelPart.Name());

    // Add nodes
    std::vector<IndexType> aux_node_ids;
    aux_node_ids.reserve(rOriginModelPart.NumberOfNodes());
    for (const auto& r_node : rOriginModelPart.Nodes()) {
        if (rNodesSet.find(r_node.Id()) != rNodesSet.end()) {
            aux_node_ids.push_back(r_node.Id());
        }
    }
    r_hrom_sub_mp.AddNodes(aux_node_ids);

    // Add elements
    std::vector<IndexType> aux_elem_ids;
    aux_elem_ids.reserve(rOriginModelPart.NumberOfElements());
    for (const auto& r_elem : rOriginModelPart.Elements()) {
        auto is_found = [&r_elem](Element::Pointer p_elem){return r_elem.Id() == p_elem->Id();};
        if (std::find_if(rElementsVector.begin(), rElementsVector.end(), is_found) != rElementsVector.end()) {
            aux_elem_ids.push_back(r_elem.Id());
        }
    }
    r_hrom_sub_mp.AddElements(aux_elem_ids);

    // Add conditions
    std::vector<IndexType> aux_cond_ids;
    aux_cond_ids.reserve(rOriginModelPart.NumberOfConditions());
    for (const auto& r_cond : rOriginModelPart.Conditions()) {
        auto is_found = [&r_cond](Condition::Pointer p_cond){return r_cond.Id() == p_cond->Id();};
        if (std::find_if(rConditionsVector.begin(), rConditionsVector.end(), is_found) != rConditionsVector.end()) {
            aux_cond_ids.push_back(r_cond.Id());
        }
    }
    r_hrom_sub_mp.AddConditions(aux_cond_ids);

    //TODO: ADD MPCs

    // Recursive addition
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(rNodesSet, rElementsVector, rConditionsVector, r_orig_sub_mp, r_hrom_sub_mp);
    }
}

void RomAuxiliaryUtilities::SetHRomVisualizationModelPart(
    const ModelPart& rOriginModelPart,
    ModelPart& rHRomVisualizationModelPart)
{
    KRATOS_ERROR << "TO BE IMPLEMENTED" << std::endl;
}

} // namespace Kratos
