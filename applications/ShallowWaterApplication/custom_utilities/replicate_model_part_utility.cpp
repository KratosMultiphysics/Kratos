//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "replicate_model_part_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"


namespace Kratos
{

ReplicateModelPartUtility::ReplicateModelPartUtility(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart
) : mrOriginModelPart(rOriginModelPart)
  , mrDestinationModelPart(rDestinationModelPart)
{
    KRATOS_ERROR_IF(mrOriginModelPart.IsSubModelPart()) << "The origin model part is not a root model part" << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.IsSubModelPart()) << "The destination model part is not a root model part" << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfNodes() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfNodes() << " nodes." << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfElements() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfElements() << " elements." << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfConditions() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfConditions() << " conditions." << std::endl;
}

void ReplicateModelPartUtility::Replicate()
{
    const int num_nodes = static_cast<int>(mrOriginModelPart.NumberOfNodes());
    const int num_elements = static_cast<int>(mrOriginModelPart.NumberOfElements());
    const int num_conditions = static_cast<int>(mrOriginModelPart.NumberOfConditions());

    AssignUniqueModelPartCollectionTagUtility::IndexStringMapType collections;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_nodes_tags;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_elements_tags;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_conditions_tags;

    AssignUniqueModelPartCollectionTagUtility(mrOriginModelPart).ComputeTags(
        origin_nodes_tags,
        origin_conditions_tags,
        origin_elements_tags,
        collections);

    std::unordered_map<IndexType, std::vector<IndexType>> dest_nodes_of_collection;
    std::unordered_map<IndexType, std::vector<IndexType>> dest_elems_of_collection;
    std::unordered_map<IndexType, std::vector<IndexType>> dest_conds_of_collection;

    // Note: we need to manually set the variables list to the new nodes
    VariablesList& dest_variables_list = mrDestinationModelPart.GetNodalSolutionStepVariablesList();
    for (int i = 0; i < num_nodes; i++)
    {
        auto it_node = mrOriginModelPart.NodesBegin() + i;
        auto new_node = it_node->Clone();
        new_node->SetSolutionStepVariablesList(&dest_variables_list);
        mrDestinationModelPart.AddNode(new_node);
        auto id = it_node->Id();
        auto tag = origin_nodes_tags[id];
        dest_nodes_of_collection[tag].push_back(id);
    }

    // Note: we nee to manually create the geometry with the cloned nodes
    for (int i = 0; i < num_elements; i++)
    {
        auto it_elem = mrOriginModelPart.ElementsBegin() + i;
        PointerVector<Node<3>> new_elem_nodes;
        for (auto& node : it_elem->GetGeometry())
        {
            new_elem_nodes.push_back(mrDestinationModelPart.pGetNode(node.Id()));
        }
        auto id = it_elem->Id();
        auto new_elem = it_elem->Clone(id, new_elem_nodes);
        mrDestinationModelPart.AddElement(new_elem);
        auto tag = origin_elements_tags[id];
        dest_elems_of_collection[tag].push_back(id);
    }

    // Note: we need to manually create the geometry with the cloned nodes
    for (int i = 0; i < num_conditions; i++)
    {
        auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
        PointerVector<Node<3>> new_cond_nodes;
        for (auto& node : it_cond->GetGeometry())
        {
            new_cond_nodes.push_back(mrDestinationModelPart.pGetNode(node.Id()));
        }
        auto id = it_cond->Id();
        auto new_cond = it_cond->Clone(id, new_cond_nodes);
        mrDestinationModelPart.AddCondition(new_cond);
        auto tag = origin_conditions_tags[id];
        dest_conds_of_collection[tag].push_back(id);
    }

    // We need to replicate the sub model parts and then, populate them
    auto model_part_names = mrOriginModelPart.GetSubModelPartNames();
    for (auto& name : model_part_names)
    {
        mrDestinationModelPart.CreateSubModelPart(name);
    }
    
    // Finally, we add the nodes, elements and conditions to the sub model parts
    for (auto model_part_collection : collections)
    {
        IndexType tag = model_part_collection.first;
        if (tag != 0) // Tag = 0 means the root model part
        {
            for (auto name : model_part_collection.second)
            {
                ModelPart& model_part = mrDestinationModelPart.GetSubModelPart(name);
                model_part.AddNodes(dest_nodes_of_collection[tag]);
                model_part.AddElements(dest_elems_of_collection[tag]);
                model_part.AddConditions(dest_conds_of_collection[tag]);
            }
        }
    }
}

}  // namespace Kratos.
