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

    // Auxiliary set to save the nodes involved in the HROM mesh
    NodesPointerSet hrom_nodes_set;
    ElementsPointerSet hrom_elems_set;
    ConditionsPointerSet hrom_conds_set;

    const auto& r_elem_weights = HRomWeights["Elements"];
    for (auto it = r_elem_weights.begin(); it != r_elem_weights.end(); ++it) {
        const IndexType elem_id = stoi(it.name());
        const auto p_elem = rOriginModelPart.pGetElement(elem_id);
        hrom_elems_set.insert(hrom_elems_set.end(), p_elem);
        const auto& r_geom = p_elem->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for(IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
        }
    }

    const auto& r_cond_weights = HRomWeights["Conditions"];
    for (auto it = r_cond_weights.begin(); it != r_cond_weights.end(); ++it) {
        const IndexType cond_id = stoi(it.name());
        const auto p_cond = rOriginModelPart.pGetCondition(cond_id);
        hrom_conds_set.insert(hrom_conds_set.end(), p_cond);
        const auto& r_geom = p_cond->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for(IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
        }
    }

    // Create and fill the HROM calculation model part
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(hrom_nodes_set, hrom_elems_set, hrom_conds_set, r_orig_sub_mp, rHRomComputingModelPart);
    }
}

void RomAuxiliaryUtilities::RecursiveHRomModelPartCreation(
    const NodesPointerSet& rNodesSet,
    const ElementsPointerSet& rElementsSet,
    const ConditionsPointerSet& rConditionsSet,
    const ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    // Emulate the origin submodelpart hierarchy
    auto& r_hrom_sub_mp = rDestinationModelPart.CreateSubModelPart(rOriginModelPart.Name());

    // Add nodes
    SizeType n_orig_nodes = rOriginModelPart.NumberOfNodes();
    #pragma omp parallel for
    for (IndexType i_node = 0; i_node < n_orig_nodes; ++i_node) {
        auto it_node = rOriginModelPart.NodesBegin() + i_node;
        auto p_orig_node = rOriginModelPart.pGetNode(it_node->Id());
        if (rNodesSet.find(p_orig_node) != rNodesSet.end()) {
            #pragma omp critical
            {
                rDestinationModelPart.AddNode(p_orig_node);
            }
        }
    }

    // Add elements
    SizeType n_orig_elems = rOriginModelPart.NumberOfElements();
    #pragma omp parallel for
    for (IndexType i_elem = 0; i_elem < n_orig_elems; ++i_elem) {
        auto it_elem = rOriginModelPart.ElementsBegin() + i_elem;
        auto p_orig_elem = rOriginModelPart.pGetElement(it_elem->Id());
        if (rElementsSet.find(p_orig_elem) != rElementsSet.end()) {
            #pragma omp critical
            {
                rDestinationModelPart.AddElement(p_orig_elem);
            }
        }
    }

    // Add conditions
    SizeType n_orig_conds = rOriginModelPart.NumberOfConditions();
    #pragma omp parallel for
    for (IndexType i_cond = 0; i_cond < n_orig_conds; ++i_cond) {
        auto it_cond = rOriginModelPart.ConditionsBegin() + i_cond;
        auto p_orig_cond = rOriginModelPart.pGetCondition(it_cond->Id());
        if (rConditionsSet.find(p_orig_cond) != rConditionsSet.end()) {
            #pragma omp critical
            {
                rDestinationModelPart.AddCondition(p_orig_cond);
            }
        }
    }

    // Recursive addition
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(rNodesSet, rElementsSet, rConditionsSet, r_orig_sub_mp, r_hrom_sub_mp);
    }
}

void RomAuxiliaryUtilities::SetHRomVisualizationModelPart(
    const ModelPart& rOriginModelPart,
    ModelPart& rHRomVisualizationModelPart)
{
    KRATOS_ERROR << "TO BE IMPLEMENTED" << std::endl;
}

} // namespace Kratos
