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

/***********************************************************************************/
/***********************************************************************************/

/// Destructor
AssignUniqueModelPartCollectionTagUtility::~AssignUniqueModelPartCollectionTagUtility() {};

/***********************************************************************************/
/***********************************************************************************/

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
    const StringVectorType& r_model_part_names = GetRecursiveSubModelPartNames(mrModelPart);

    // Initialize the collections
    IndexType tag = 0;
    for (IndexType i_sub_model_part = 0; i_sub_model_part < r_model_part_names.size(); ++i_sub_model_part) {
        rCollections[i_sub_model_part].push_back(r_model_part_names[i_sub_model_part]);

        if (tag > 0) {
            ModelPart& r_sub_model_part = GetRecursiveSubModelPart(mrModelPart, r_model_part_names[i_sub_model_part]);

            /* Nodes */
            NodesArrayType& r_nodes_array = r_sub_model_part.Nodes();
            const auto it_node_begin = r_nodes_array.begin();
            for(IndexType i_node = 0; i_node < r_nodes_array.size(); ++i_node)
                aux_node_tags[(it_node_begin + i_node)->Id()].insert(tag);

            /* Conditions */
            ConditionsArrayType& r_conditions_array = r_sub_model_part.Conditions();
            const auto it_cond_begin = r_conditions_array.begin();
            for(IndexType i_cond = 0; i_cond < r_conditions_array.size(); ++i_cond)
                aux_cond_tags[(it_cond_begin + i_cond)->Id()].insert(tag);

            /* Elements */
            ElementsArrayType& r_elements_array = r_sub_model_part.Elements();
            const auto it_elem_begin = r_elements_array.begin();
            for(IndexType i_elem = 0; i_elem < r_elements_array.size(); ++i_elem)
                aux_elem_tags[(it_elem_begin + i_elem)->Id()].insert(tag);
        }

        ++tag;
    }

    // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously
   IndexSetIndexMapType combinations;

    /* Nodes */
    for(auto& r_aux_node_tag : aux_node_tags) {
        const std::set<IndexType>& r_value = r_aux_node_tag.second;
        if (r_value.size() > 1) combinations[r_value] = 0;
    }

    /* Conditions */
    for(auto& r_aux_cond_tag : aux_cond_tags) {
        const std::set<IndexType>& r_value = r_aux_cond_tag.second;
        if (r_value.size() > 1) combinations[r_value] = 0;
    }

    /* Elements */
    for(auto& r_aux_elem_tag : aux_elem_tags) {
        const std::set<IndexType>& r_value = r_aux_elem_tag.second;
        if (r_value.size() > 1) combinations[r_value] = 0;
    }

    if (mrModelPart.IsDistributed()) {
        SetParallelModelPartAndSubModelPartCollectionsAndCombinations(rCollections, combinations, tag);
    } else {
        /* Combinations */
        for(auto& combination : combinations) {
            const std::set<IndexType>& r_key_set = combination.first;
            for(IndexType it : r_key_set)
                rCollections[tag].push_back(rCollections[it][0]);
            combinations[r_key_set] = tag;
            ++tag;
        }
    }

    // The final maps are created
    /* Nodes */
    for(auto& r_aux_node_tag : aux_node_tags) {
        const IndexType key = r_aux_node_tag.first;
        const std::set<IndexType>& r_value = r_aux_node_tag.second;

        if (r_value.size() == 0)
            rNodeTags[key] = 0; // Main Model Part
        else if (r_value.size() == 1) // A Sub Model Part
            rNodeTags[key] = *r_value.begin();
        else // There is a combination
            rNodeTags[key] = combinations[r_value];
    }

    /* Conditions */
    for(auto& r_aux_cond_tag : aux_cond_tags) {
        const IndexType key = r_aux_cond_tag.first;
        const std::set<IndexType>& r_value = r_aux_cond_tag.second;

        if (r_value.size() == 0)
            rCondTags[key] = 0; // Main Model Part
        else if (r_value.size() == 1) // A Sub Model Part
            rCondTags[key] = *r_value.begin();
        else // There is a combination
            rCondTags[key] = combinations[r_value];
    }

    /* Elements */
    for(auto& r_aux_elem_tag : aux_elem_tags) {
        const IndexType key = r_aux_elem_tag.first;
        const std::set<IndexType>& r_value = r_aux_elem_tag.second;

        if (r_value.size() == 0)
            rElemTags[key] = 0; // Main Model Part
        else if (r_value.size() == 1) // A Sub Model Part
            rElemTags[key] = *r_value.begin();
        else // There is a combination
            rElemTags[key] = combinations[r_value];
    }

    // Clean up the collections
    for (auto& r_collection : rCollections) {
        std::unordered_set<std::string> aux_set;
        for (auto& r_name : r_collection.second) {
            aux_set.insert(r_name);
        }
        std::vector<std::string> aux_vector;
        for (auto& r_name : aux_set) {
            aux_vector.push_back(r_name);
        }
        r_collection.second = aux_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AssignUniqueModelPartCollectionTagUtility::SetParallelModelPartAndSubModelPartCollectionsAndCombinations(
                                            IndexStringMapType& rCollections,
                                            IndexSetIndexMapType& rCombinations,
                                            IndexType& rTag)
{
    auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
    const IndexType rank = r_data_communicator.Rank();
    const IndexType size = r_data_communicator.Size();

    std::vector<std::vector<IndexType>> local_combination_vector_keys;
    for(auto& key_combination_set : rCombinations) {
        std::vector<IndexType> key_vector;
        for (auto it=key_combination_set.first.begin(); it != key_combination_set.first.end(); ++it)
            key_vector.push_back(*it);
        local_combination_vector_keys.push_back(key_vector);
    }

    std::vector<std::vector<IndexType>> global_combination_vector_keys;
    if (rank>0) {
        r_data_communicator.Send(local_combination_vector_keys, 0);
    }
    else {
        for (IndexType i_rank = 1; i_rank<size; i_rank++) {
            std::vector<std::vector<IndexType>> recv_irank_combinations;
            r_data_communicator.Recv(recv_irank_combinations, i_rank);

            for (auto& key_combination_vector : recv_irank_combinations) {
                // Convert vector to set
                const std::set<IndexType> submodel_part_set(key_combination_vector.begin(), key_combination_vector.end());
                // Check if key exists in the combinations maps, if not, add it
                if (rCombinations.find(submodel_part_set) == rCombinations.end()) {
                    rCombinations[submodel_part_set] = 0; //creating new key, value does not matter
                }
            }
        }
        // Converting back to vector to be broadcasted
        for(auto& r_key_combination_set : rCombinations) {
            std::vector<IndexType> key_vector;
            for (auto it=r_key_combination_set.first.begin(); it != r_key_combination_set.first.end(); ++it)
                key_vector.push_back(*it);
            global_combination_vector_keys.push_back(key_vector);
        }
    }

    r_data_communicator.Broadcast(global_combination_vector_keys, 0);

    rCombinations.clear();

    // Assigning final maps considering all ranks information
    for (auto& r_key_vector : global_combination_vector_keys) {
        const std::set<IndexType> r_key_set(r_key_vector.begin(), r_key_vector.end());
        for(IndexType it : r_key_set)
            rCollections[rTag].push_back(rCollections[it][0]);
        rCombinations[r_key_set] = rTag;
        ++rTag;
    }
}
/***********************************************************************************/
/***********************************************************************************/

/// This functions gets the "colors" from an existing json file
Parameters AssignUniqueModelPartCollectionTagUtility::ReadTagsFromJson(
    const std::string& rFilename,
    IndexStringMapType& rCollections
    )
{
    std::ifstream infile(rFilename + ".json");
    KRATOS_ERROR_IF_NOT(infile.good()) << "Tags file: " << rFilename  + ".json" << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters color_json(buffer.str());
    for (auto it_param = color_json.begin(); it_param != color_json.end(); ++it_param) {
        const std::vector<std::string>& r_sub_model_part_names = it_param->GetStringArray();
        rCollections.insert(std::pair<IndexType,std::vector<std::string>>({std::stoi(it_param.name()), r_sub_model_part_names}));
    }

    return color_json;
}

/***********************************************************************************/
/***********************************************************************************/

/// This functions writes the "colors" to a new json file
Parameters AssignUniqueModelPartCollectionTagUtility::WriteTagsToJson(
    const std::string& rFilename,
    const IndexStringMapType& rCollections
    )
{
    Parameters color_json;
    for (auto& r_color : rCollections) {
        Parameters names_array;
        names_array.AddEmptyArray("r_model_part_names");
        for (auto& r_model_part_name : r_color.second) {
            names_array["r_model_part_names"].Append(r_model_part_name);
        }
        color_json.AddValue(std::to_string(r_color.first), names_array["r_model_part_names"]);
    }

    const std::string& r_json_text = color_json.PrettyPrintJsonString();

    std::filebuf buffer;
    buffer.open(rFilename + ".json",std::ios::out);
    std::ostream os(&buffer);
    os << r_json_text;
    buffer.close();

    return color_json;
}

/***********************************************************************************/
/***********************************************************************************/

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
    for (auto& r_name : names) {
        ModelPart& r_sub_model_part = rThisModelPart.GetSubModelPart(r_name);
        r_name.insert(0, Prefix);
        sub_model_parts_names.push_back(r_name);
        const auto sub_names = GetRecursiveSubModelPartNames(r_sub_model_part, r_name);
        for (auto r_sub_name : sub_names)
            sub_model_parts_names.push_back(r_sub_name);
    }

    return sub_model_parts_names;
}

/***********************************************************************************/
/***********************************************************************************/

/// Get a submodelpart given its full r_name
ModelPart& AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(
    ModelPart& rThisModelPart,
    const std::string& rFullName
)
{
    std::istringstream full_name(rFullName);
    return AuxGetSubModelPart(rThisModelPart, full_name);
}

/***********************************************************************************/
/***********************************************************************************/

/// Private function for GetRecursiveSubModelPart
ModelPart& AssignUniqueModelPartCollectionTagUtility::AuxGetSubModelPart(
    ModelPart& rThisModelPart,
    std::istringstream& rFullName
)
{
    std::string r_name;
    if (std::getline(rFullName, r_name, '.')) {
        if (rThisModelPart.Name() == r_name)
            return AuxGetSubModelPart(rThisModelPart, rFullName);
        else
            return AuxGetSubModelPart(rThisModelPart.GetSubModelPart(r_name), rFullName);
    } else
        return rThisModelPart;
}

/***********************************************************************************/
/***********************************************************************************/

/// Debugging purpose
void AssignUniqueModelPartCollectionTagUtility::DebugAssignUniqueModelPartCollectionTag()
{
    IndexIndexMapType node_tags, cond_tags, elem_tags;
    IndexStringMapType collections;

    ComputeTags(node_tags, cond_tags, elem_tags, collections);

    // The collections are the following
    for (auto& r_collection : collections) {
        KRATOS_INFO("") << "TAG: " << r_collection.first << " has the following model part r_collection" << std::endl;
        for (auto& r_name : r_collection.second) {
            KRATOS_INFO("") << "\t" << r_name << std::endl;
        }
    }

    // The nodes belonging to collections are
    for (auto& r_node_tag : node_tags) {
        KRATOS_INFO("") << "NODE TAG: " << r_node_tag.second << "\t" << r_node_tag.first << std::endl;
    }

    // The conditions belonging to collections are
    for (auto& r_cond_tag : cond_tags) {
        KRATOS_INFO("") << "CONDITION TAG: " << r_cond_tag.second << "\t" << r_cond_tag.first << std::endl;
    }

    // The elements belonging to collections are
    for (auto& r_elem_tag : elem_tags) {
        KRATOS_INFO("") << "ELEMENT TAG: " << r_elem_tag.second << "\t" << r_elem_tag.first << std::endl;
    }
}

}  // namespace Kratos.
