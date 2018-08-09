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
#include "utilities/sub_model_parts_list_utility.h"


namespace Kratos
{
/********************************* CONSTRUCTOR *************************************/
/***********************************************************************************/

SubModelPartsListUtility::SubModelPartsListUtility(ModelPart& rModelPart):mrModelPart(rModelPart) {};


/******************************** DESTRUCTOR ***************************************/
/***********************************************************************************/

SubModelPartsListUtility::~SubModelPartsListUtility() {};

/******************************** PUBLIC METHODS ***********************************/
/***********************************************************************************/

void SubModelPartsListUtility::ComputeSubModelPartsList(
    IndexIntMapType& rNodesColors,
    IndexIntMapType& rCondColors,
    IndexIntMapType& rElemColors,
    IntStringMapType& rColors
    )
{
    // Initialize and create the auxiliary maps
    std::unordered_map<IndexType,std::set<IndexType>> aux_nodes_colors, aux_cond_colors, aux_elem_colors;

    // We compute the list of submodelparts and subsubmodelparts
    const std::vector<std::string>& model_part_names = GetRecursiveSubModelPartNames(mrModelPart);

    // Initialize Colors
    IndexType color = 0;
    for (SizeType i_sub_model_part = 0; i_sub_model_part < model_part_names.size(); ++i_sub_model_part) {
        rColors[i_sub_model_part].push_back(model_part_names[i_sub_model_part]);

        if (color > 0) {
            ModelPart& r_sub_model_part = GetRecursiveSubModelPart(mrModelPart, model_part_names[i_sub_model_part]);

            /* Nodes */
            NodesArrayType& nodes_array = r_sub_model_part.Nodes();
            for(SizeType i = 0; i < nodes_array.size(); ++i)
                aux_nodes_colors[(nodes_array.begin() + i)->Id()].insert(color);

            /* Conditions */
            ConditionsArrayType& conditions_array = r_sub_model_part.Conditions();
            for(SizeType i = 0; i < conditions_array.size(); ++i)
                aux_cond_colors[(conditions_array.begin() + i)->Id()].insert(color);

            /* Elements */
            ElementsArrayType& elements_array = r_sub_model_part.Elements();
            for(SizeType i = 0; i < elements_array.size(); ++i)
                aux_elem_colors[(elements_array.begin() + i)->Id()].insert(color);
        }

        color += 1;
    }

    // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously
    std::unordered_map<std::set<IndexType>, IndexType, KeyHasherRange<std::set<IndexType>>, KeyComparorRange<std::set<IndexType>> > combinations;

    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const std::set<IndexType>& value = aux_nodes_color.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const std::set<IndexType>& value = aux_cond_color.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const std::set<IndexType>& value = aux_elem_color.second;
        if (value.size() > 1) combinations[value] = 0;
    }

    /* Combinations */
    for(auto & combination : combinations) {
        const std::set<IndexType>& key = combination.first;
        for(IndexType it : key)
            rColors[color].push_back(rColors[it][0]);
        combinations[key] = color;
        color += 1;
    }

    // The final maps are created
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const IndexType key = aux_nodes_color.first;
        const std::set<IndexType>& value = aux_nodes_color.second;

        if (value.size() == 0)
            rNodesColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rNodesColors[key] = *value.begin();
        else // There is a combination
            rNodesColors[key] = combinations[value];
    }

    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const IndexType key = aux_cond_color.first;
        const std::set<IndexType>& value = aux_cond_color.second;

        if (value.size() == 0)
            rCondColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rCondColors[key] = *value.begin();
        else // There is a combination
            rCondColors[key] = combinations[value];
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const IndexType key = aux_elem_color.first;
        const std::set<IndexType>& value = aux_elem_color.second;

        if (value.size() == 0)
            rElemColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rElemColors[key] = *value.begin();
        else // There is a combination
            rElemColors[key] = combinations[value];
    }

    // Clean up the colors
    for (auto& color : rColors) {
        std::unordered_set<std::string> aux_set;
        for (auto& name : color.second) {
            aux_set.insert(name);
        }
        std::vector<std::string> aux_vector;
        for (auto& name : aux_set) {
            aux_vector.push_back(name);
        }
        color.second = aux_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> SubModelPartsListUtility::GetRecursiveSubModelPartNames(ModelPart& ThisModelPart)
{
    // We compute the list of submodelparts
    const std::vector<std::string> sub_model_part_names = ThisModelPart.GetSubModelPartNames();

    std::vector<std::string> model_part_names;
    model_part_names.push_back(ThisModelPart.Name());
    for (const auto & sub_model_part_name : sub_model_part_names) {
        model_part_names.push_back(sub_model_part_name);
        ModelPart& r_sub_model_part = ThisModelPart.GetSubModelPart(sub_model_part_name); // We check for sub sub model parts (no more sublevels)
        if (r_sub_model_part.NumberOfSubModelParts() > 0) {
            const std::vector<std::string> sub_sub_model_part_names = r_sub_model_part.GetSubModelPartNames();
            for (const auto& sub_sub_model_part_name : sub_sub_model_part_names) {
                model_part_names.push_back(sub_sub_model_part_name);
            }
        }
    }

    // Check for repeated names on the submodelparts (this is not checked by model_part.h if we work with subsubmodelparts)
    std::sort(model_part_names.begin()+1, model_part_names.end());
    auto last = std::unique(model_part_names.begin()+1, model_part_names.end());
    KRATOS_ERROR_IF_NOT(last == model_part_names.end()) << "ERROR:: Repeated names in subsubmodelparts. Check subsubmodelparts names please" << std::endl;

    return model_part_names;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& SubModelPartsListUtility::GetRecursiveSubModelPart(
    ModelPart& ThisModelPart,
    const std::string& SubModelPartName
    )
{
    // We check if main model_part
    if (ThisModelPart.Name() == SubModelPartName)
        return ThisModelPart;

    // We check a submodelpart
    if (ThisModelPart.HasSubModelPart(SubModelPartName)) // In case we are in a submodelpart
        return ThisModelPart.GetSubModelPart(SubModelPartName);
    else { // In case we are in a subsubmodelpart
        const std::vector<std::string> sub_model_part_names = ThisModelPart.GetSubModelPartNames();
        for (const auto & sub_model_part_name : sub_model_part_names) {
            ModelPart& r_sub_model_part = ThisModelPart.GetSubModelPart(sub_model_part_name); // We check for sub sub model parts (no more sublevels)
            if (r_sub_model_part.HasSubModelPart(SubModelPartName)) {
                return r_sub_model_part.GetSubModelPart(SubModelPartName);
            }
        }
    }

    return ThisModelPart;
}

/***********************************************************************************/
/***********************************************************************************/

void SubModelPartsListUtility::DebugComputeSubModelPartsList()
{
    IndexIntMapType node_colors, cond_colors, elem_colors;
    IntStringMapType colors;

    ComputeSubModelPartsList(node_colors, cond_colors, elem_colors, colors);

    // The colors are the following
    for (auto& color : colors) {
        KRATOS_INFO("") << "COLOR: " << color.first << " has the following model parts" << std::endl;
        for (auto& name : color.second) {
            KRATOS_INFO("") << "\t" << name << std::endl;
        }
    }

    // The nodes belonging to colors are
    for (auto& node_color : node_colors) {
        KRATOS_INFO("") << "NODE COLOR: " << node_color.second << "\t" << node_color.first << std::endl;
    }

    // The conditions belonging to colors are
    for (auto& cond_color : cond_colors) {
        KRATOS_INFO("") << "CONDITION COLOR: " << cond_color.second << "\t" << cond_color.first << std::endl;
    }

    // The elements belonging to colors are
    for (auto& elem_color : elem_colors) {
        KRATOS_INFO("") << "ELEMENT COLOR: " << elem_color.second << "\t" << elem_color.first << std::endl;
    }
}

}  // namespace Kratos.
