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
    std::unordered_map<IndexType,std::set<int>> aux_nodes_colors, aux_cond_colors, aux_elem_colors;

    // We compute the list of submodelparts and subsubmodelparts
    const std::vector<std::string>& model_part_names = GetRecursiveSubModelPartNames(mrModelPart);

    // Initialize Colors
    int color = 0;
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
    std::unordered_map<std::set<int>, int, KeyHasherRange<std::set<int>>, KeyComparorRange<std::set<int>> > combinations;

    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const std::set<int>& value = aux_nodes_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }

    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const std::set<int>& value = aux_cond_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const std::set<int>& value = aux_elem_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }

    /* Combinations */
    for(auto & combination : combinations) {
        const std::set<int>& key = combination.first;
        for(int it : key)
            rColors[color].push_back(rColors[it][0]);
        combinations[key] = color;
        color += 1;
    }

    // The final maps are created
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) {
        const int key = aux_nodes_color.first;
        const std::set<int>& value = aux_nodes_color.second;

        if (value.size() == 0)
            rNodesColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rNodesColors[key] = *value.begin();
        else // There is a combination
            rNodesColors[key] = combinations[value];
    }

    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const int key = aux_cond_color.first;
        const std::set<int>& value = aux_cond_color.second;

        if (value.size() == 0)
            rCondColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rCondColors[key] = *value.begin();
        else // There is a combination
            rCondColors[key] = combinations[value];
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const int key = aux_elem_color.first;
        const std::set<int>& value = aux_elem_color.second;

        if (value.size() == 0)
            rElemColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // A Sub Model Part
            rElemColors[key] = *value.begin();
        else // There is a combination
            rElemColors[key] = combinations[value];
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

std::unordered_map<int,std::vector<ModelPart*>> SubModelPartsListUtility::GetModelPartColorsPointers(
    ModelPart& rThisModelPart,
    const IntStringMapType& rColors
    )
{
    // Initialize output
    std::unordered_map<int,std::vector<ModelPart*>> colors_ptr;

    for (auto color : rColors)
    {
        for (auto name : color.second)
        {
            colors_ptr[color.first].push_back(&SubModelPartsListUtility::GetRecursiveSubModelPart(rThisModelPart, name));
        }
    }

    return colors_ptr;
}

/***********************************************************************************/
/***********************************************************************************/

void SubModelPartsListUtility::IntersectColors(
    IntStringMapType& rColors,
    PairIntMapType& rIntersections
    )
{
    // Check the alphabetic order of the names again
    for (auto& color : rColors)
    {
        std::sort(color.second.begin(), color.second.end());
    }

    // Auxiliary variables
    bool new_color = false;
    int last_color = rColors.size();
    std::unordered_map<std::vector<std::string>, int, KeyHasherRange<std::vector<std::string>>, KeyComparorRange<std::vector<std::string>>> aux_colors;

    // Generate the intersections map
    for (SizeType i = 0; i < rColors.size(); i++)
    {
        for (SizeType j = i; j < rColors.size(); j++)
        {
            std::pair<int,int> intersection_key = {i,j};
            if (i==0)
            {
                rIntersections[intersection_key] = 0;
            }
            else if (i == j)
            {
                rIntersections[intersection_key] = i;
            }
            else
            {
                std::vector<std::string> intersection_names;
                std::set_intersection(rColors[i].begin(), rColors[i].end(),
                                      rColors[j].begin(), rColors[j].end(),
                                      std::back_inserter(intersection_names));

                // Null intersection: Main Model Part
                if (intersection_names.size() == 0)
                {
                    rIntersections[intersection_key] = 0;
                }
                else
                {
                    // Find the intersection color
                    bool intersection_found = false;
                    for (auto color : rColors)
                    {
                        if (color.second == intersection_names)
                        {
                            rIntersections[intersection_key] = color.first;
                            intersection_found = true;
                        }
                    }

                    // It should never happen, just in case
                    if (!intersection_found)
                    {
                        if (aux_colors.count(intersection_names) == 0)
                        {
                            aux_colors[intersection_names] = last_color++;
                            new_color = true;
                        }
                        rIntersections[intersection_key] = aux_colors[intersection_names];
                    }
                }
            }
        }
    }

    // Add the new combinations to the colors
    if (new_color)
    {
        for (auto aux_color : aux_colors)
        {
            rColors[aux_color.second] = aux_color.first;
        }
    }
}

}  // namespace Kratos.
