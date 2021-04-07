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
    ModelPart& rDestinationModelPart,
    bool ReplicateSubModelParts
) : mrOriginModelPart(rOriginModelPart)
  , mrDestinationModelPart(rDestinationModelPart)
  , mReplicateSubModelParts(ReplicateSubModelParts)
{
    KRATOS_ERROR_IF(mrOriginModelPart.IsSubModelPart() & mReplicateSubModelParts) << "The origin model part is not a root model part" << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.IsSubModelPart()) << "The destination model part is not a root model part" << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfNodes() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfNodes() << " nodes." << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfElements() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfElements() << " elements." << std::endl;
    KRATOS_ERROR_IF(mrDestinationModelPart.NumberOfConditions() != 0) << "The destination model part is not empty. There are " << mrDestinationModelPart.NumberOfConditions() << " conditions." << std::endl;
}

void ReplicateModelPartUtility::Replicate()
{
    IndexType unique_node_id;
    IndexType unique_elem_id;
    IndexType unique_cond_id;
    IndexType unique_prop_id;
    GetMaximumIds(unique_node_id, unique_elem_id, unique_cond_id, unique_prop_id);

    // Setting the nodal solution step variables
    VariablesList& variables_list = mrOriginModelPart.GetNodalSolutionStepVariablesList();
    mrDestinationModelPart.GetNodalSolutionStepVariablesList() = variables_list;

    // Both model parts should share the TIME and the STEP
    mrDestinationModelPart.SetProcessInfo(mrOriginModelPart.pGetProcessInfo());

    // In case we need to replicate the sub model parts
    AssignUniqueModelPartCollectionTagUtility::IndexStringMapType collections;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_nodes_tags;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_elements_tags;
    AssignUniqueModelPartCollectionTagUtility::IndexIndexMapType origin_conditions_tags;

    std::unordered_map<IndexType, std::vector<IndexType>> dest_nodes_of_collection;
    std::unordered_map<IndexType, std::vector<IndexType>> dest_elems_of_collection;
    std::unordered_map<IndexType, std::vector<IndexType>> dest_conds_of_collection;

    if (mReplicateSubModelParts)
    {
        AssignUniqueModelPartCollectionTagUtility(mrOriginModelPart).ComputeTags(
            origin_nodes_tags,
            origin_conditions_tags,
            origin_elements_tags,
            collections);
    }

    // Clone the nodes
    mReplicatedNodesMap.clear();
    ModelPart::NodesContainerType aux_array_with_node_pointers;
    for (IndexType i = 0; i < mrOriginModelPart.NumberOfNodes(); ++i)
    {
        auto it_node = mrOriginModelPart.NodesBegin() + i;
        auto new_node = it_node->Clone();
        new_node->SetId(++unique_node_id);
        new_node->SetSolutionStepVariablesList(&variables_list);
        aux_array_with_node_pointers.push_back(new_node);
        mReplicatedNodesMap[it_node->Id()] = new_node;
        if (mReplicateSubModelParts)
        {
            auto tag = origin_nodes_tags[unique_node_id];
            dest_nodes_of_collection[tag].push_back(unique_node_id);
        }
    }
    mrDestinationModelPart.AddNodes(aux_array_with_node_pointers.begin(), aux_array_with_node_pointers.end());

    // We also need to replicate the properties before creating the new elements
    std::unordered_map<IndexType, Properties::Pointer> replicated_properties_map;
    for (IndexType i = 0; i < mrOriginModelPart.NumberOfProperties(); ++i)
    {
        auto it_prop = mrOriginModelPart.PropertiesBegin() + i;
        auto new_prop = Kratos::make_shared<Properties>(*it_prop);
        new_prop->SetId(++unique_prop_id);
        mrDestinationModelPart.AddProperties(new_prop);
        replicated_properties_map[it_prop->Id()] = new_prop;
    }

    // Note: we need to manually create the geometry with the cloned nodes
    ModelPart::ElementsContainerType aux_array_with_element_pointers;
    for (IndexType i = 0; i < mrOriginModelPart.NumberOfElements(); ++i)
    {
        auto it_elem = mrOriginModelPart.ElementsBegin() + i;
        NodesArrayType new_elem_nodes;
        for (auto& node : it_elem->GetGeometry())
        {
            new_elem_nodes.push_back(mReplicatedNodesMap[node.Id()]);
        }
        auto new_elem = it_elem->Clone(++unique_elem_id, new_elem_nodes);
        new_elem->SetProperties(replicated_properties_map[it_elem->GetProperties().Id()]);
        aux_array_with_element_pointers.push_back(new_elem);
        if (mReplicateSubModelParts)
        {
            auto tag = origin_elements_tags[unique_elem_id];
            dest_elems_of_collection[tag].push_back(unique_elem_id);
        }
    }
    mrDestinationModelPart.AddElements(aux_array_with_element_pointers.begin(), aux_array_with_element_pointers.end());

    // Note: we need to manually create the geometry with the cloned nodes
    ModelPart::ConditionsContainerType aux_array_with_condition_pointers;
    for (IndexType i = 0; i < mrOriginModelPart.NumberOfConditions(); i++)
    {
        auto it_cond = mrOriginModelPart.ConditionsBegin() + i;
        NodesArrayType new_cond_nodes;
        for (auto& node : it_cond->GetGeometry())
        {
            new_cond_nodes.push_back(mReplicatedNodesMap[node.Id()]);
        }
        auto new_cond = it_cond->Clone(++unique_cond_id, new_cond_nodes);
        new_cond->SetProperties(replicated_properties_map[it_cond->GetProperties().Id()]);
        aux_array_with_condition_pointers.push_back(new_cond);
        if (mReplicateSubModelParts)
        {
            auto tag = origin_conditions_tags[unique_cond_id];
            dest_conds_of_collection[tag].push_back(unique_cond_id);
        }
    }
    mrDestinationModelPart.AddConditions(aux_array_with_condition_pointers.begin(), aux_array_with_condition_pointers.end());

    // We need to replicate the sub model parts and then, populate them
    if (mReplicateSubModelParts)
    {
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
    } // if (mReplicateSubModelParts)
}


void ReplicateModelPartUtility::GetMaximumIds(IndexType& rUniqueNodeId, IndexType& rUniqueElemId, IndexType& rUniqueCondId, IndexType& rUniquePropId)
{
    rUniqueNodeId = 0;
    rUniqueElemId = 0;
    rUniqueCondId = 0;
    rUniquePropId = 0;

    for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfNodes()); ++i)
    {
        rUniqueNodeId = std::max(rUniqueNodeId, (mrOriginModelPart.NodesBegin() + i)->Id());
    }

    for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfElements()); ++i)
    {
        rUniqueElemId = std::max(rUniqueNodeId, (mrOriginModelPart.ElementsBegin() + i)->Id());
    }

    for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfConditions()); ++i)
    {
        rUniqueCondId = std::max(rUniqueNodeId, (mrOriginModelPart.ConditionsBegin() + i)->Id());
    }

    for (int i = 0; i < static_cast<int>(mrOriginModelPart.NumberOfProperties()); ++i)
    {
        rUniquePropId = std::max(rUniquePropId, (mrOriginModelPart.PropertiesBegin() + i)->Id());
    }
}

}  // namespace Kratos.
