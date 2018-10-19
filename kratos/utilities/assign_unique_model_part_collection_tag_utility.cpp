//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <set>
#include <unordered_set>

// External includes

// Project includes
#include "utilities/assign_unique_model_part_collection_tag_utility.h"

namespace Kratos
{

/// Default constructor
AssignUniqueModelPartCollectionTagUtility::AssignUniqueModelPartCollectionTagUtility(ModelPart& rModelPart)
    : mrModelPart(rModelPart) {};


/// Destructor
AssignUniqueModelPartCollectionTagUtility::~AssignUniqueModelPartCollectionTagUtility() {};


/// Copmutes the collections and assign the tags
void AssignUniqueModelPartCollectionTagUtility::ComputeTags(
    IndexIndexMapType& rNodeTags,
    IndexIndexMapType& rCondTags,
    IndexIndexMapType& rElemTags,
    IndexStringMapType& rCollections
)
{
    // Initialize and create the auxiliary maps
    IndexIndexSetMapType aux_node_tags, aux_cond_tags, aux_elem_tags;

    // We compute the list of submodelparts and subsubmodelparts
    const StringVectorType& model_part_names = GetRecursiveSubModelPartNames(mrModelPart);

    // Initialize the collections
    IndexType tag = 0;
    for (SizeType i_sub_model_part = 0; i_sub_model_part < model_part_names.size(); ++i_sub_model_part) {
        rCollections[i_sub_model_part].push_back(model_part_names[i_sub_model_part]);

        if (tag > 0) {
            ModelPart& r_sub_model_part = GetRecursiveSubModelPart(mrModelPart, model_part_names[i_sub_model_part]);

            /* Nodes */
            NodesArrayType& nodes_array = r_sub_model_part.Nodes();
            for(SizeType i = 0; i < nodes_array.size(); ++i)
                aux_node_tags[(nodes_array.begin() + i)->Id()].insert(tag);

            /* Conditions */
            ConditionsArrayType& conditions_array = r_sub_model_part.Conditions();
            for(SizeType i = 0; i < conditions_array.size(); ++i)
                aux_cond_tags[(conditions_array.begin() + i)->Id()].insert(tag);

            /* Elements */
            ElementsArrayType& elements_array = r_sub_model_part.Elements();
            for(SizeType i = 0; i < elements_array.size(); ++i)
                aux_elem_tags[(elements_array.begin() + i)->Id()].insert(tag);
        }

        tag++;
    }

    // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously
    std::unordered_map<std::set<IndexType>, IndexType, KeyHasherRange<std::set<IndexType>>, KeyComparorRange<std::set<IndexType>> > combinations;

    /* Nodes */
    for(auto & aux_node_tag : aux_node_tags) {
        const std::set<IndexType>& value = aux_node_tag.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Conditions */
    for(auto & aux_cond_tag : aux_cond_tags) {
        const std::set<IndexType>& value = aux_cond_tag.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Elements */
    for(auto & aux_elem_tag : aux_elem_tags) {
        const std::set<IndexType>& value = aux_elem_tag.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Combinations */
    for(auto & combination : combinations) {
        const std::set<IndexType>& key = combination.first;
        for(IndexType it : key)
            rCollections[tag].push_back(rCollections[it][0]);
        combinations[key] = tag;
        tag++;
    }

    // The final maps are created
    /* Nodes */
    for(auto & aux_node_tag : aux_node_tags) {
        const IndexType key = aux_node_tag.first;
        const std::set<IndexType>& value = aux_node_tag.second;

        if (value.size() == 0)
            rNodeTags[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rNodeTags[key] = *value.begin();
        else // There is a combination
            rNodeTags[key] = combinations[value];
    }

    /* Conditions */
    for(auto & aux_cond_tag : aux_cond_tags) {
        const IndexType key = aux_cond_tag.first;
        const std::set<IndexType>& value = aux_cond_tag.second;

        if (value.size() == 0)
            rCondTags[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rCondTags[key] = *value.begin();
        else // There is a combination
            rCondTags[key] = combinations[value];
    }

    /* Elements */
    for(auto & aux_elem_tag : aux_elem_tags) {
        const IndexType key = aux_elem_tag.first;
        const std::set<IndexType>& value = aux_elem_tag.second;

        if (value.size() == 0)
            rElemTags[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rElemTags[key] = *value.begin();
        else // There is a combination
            rElemTags[key] = combinations[value];
    }

    // Clean up the collections
    for (auto& collection : rCollections) {
        std::unordered_set<std::string> aux_set;
        for (auto& name : collection.second) {
            aux_set.insert(name);
        }
        std::vector<std::string> aux_vector;
        for (auto& name : aux_set) {
            aux_vector.push_back(name);
        }
        collection.second = aux_vector;
    }
}


/// Get the full names of all the submodelparts
std::vector<std::string> AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPartNames(
    ModelPart& rThisModelPart,
    std::string Prefix
)
{
    StringVectorType sub_model_parts_names;

    if (Prefix.empty())
        sub_model_parts_names.push_back(rThisModelPart.Name());
    else
        Prefix += ".";
    
    StringVectorType names = rThisModelPart.GetSubModelPartNames();
    for (auto& name : names)
    {
        ModelPart& sub_model_part = rThisModelPart.GetSubModelPart(name);
        name.insert(0, Prefix);
        sub_model_parts_names.push_back(name);
        auto sub_names = GetRecursiveSubModelPartNames(sub_model_part, name);
        for (auto sub_name : sub_names)
            names.push_back(sub_name);
    }

    return sub_model_parts_names;
}


/// Get a submodelpart given its full name
ModelPart& AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(
    ModelPart& rThisModelPart,
    const std::string& rFullName
)
{
    std::istringstream full_name(rFullName);
    return AuxGetSubModelPart(rThisModelPart, full_name);
}


/// Private function for GetRecursiveSubModelPart
ModelPart& AssignUniqueModelPartCollectionTagUtility::AuxGetSubModelPart(
    ModelPart& rThisModelPart,
    std::istringstream& rFullName
)
{
    std::string name;
    if (std::getline(rFullName, name, '.'))
        return AuxGetSubModelPart(rThisModelPart.GetSubModelPart(name), rFullName);
    else
        return rThisModelPart;
}


/// Debugging purpose
void AssignUniqueModelPartCollectionTagUtility::DebugAssignUniqueModelPartCollectionTag()
{
    IndexIndexMapType node_tags, cond_tags, elem_tags;
    IndexStringMapType collections;

    ComputeTags(node_tags, cond_tags, elem_tags, collections);

    // The collections are the following
    for (auto& collection : collections) {
        KRATOS_INFO("") << "TAG: " << collection.first << " has the following model part collection" << std::endl;
        for (auto& name : collection.second) {
            KRATOS_INFO("") << "\t" << name << std::endl;
        }
    }

    // The nodes belonging to collections are
    for (auto& node_tag : node_tags) {
        KRATOS_INFO("") << "NODE TAG: " << node_tag.second << "\t" << node_tag.first << std::endl;
    }

    // The conditions belonging to collections are
    for (auto& cond_tag : cond_tags) {
        KRATOS_INFO("") << "CONDITION TAG: " << cond_tag.second << "\t" << cond_tag.first << std::endl;
    }

    // The elements belonging to collections are
    for (auto& elem_tag : elem_tags) {
        KRATOS_INFO("") << "ELEMENT TAG: " << elem_tag.second << "\t" << elem_tag.first << std::endl;
    }
}

}  // namespace Kratos.
