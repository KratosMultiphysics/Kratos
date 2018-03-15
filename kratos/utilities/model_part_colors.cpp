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
//

// System includes
#include <set>
#include <unordered_set>

// External includes 


// Project includes
#include "utilities/model_part_colors.h"


namespace Kratos
{
// Default constructor
ModelPartColors::ModelPartColors(ModelPart& rModelPart):mrModelPart(rModelPart) {};


// Destructor
ModelPartColors::~ModelPartColors() {};


// Compute colors
void ModelPartColors::ComputeColors(
    std::unordered_map<int,int>& NodesColors,
    std::unordered_map<int,int>& CondColors,
    std::unordered_map<int,int>& ElemColors,
    std::unordered_map<int,std::vector<std::string>>& rColors
    )
{
    // Initialize and create the auxiliar maps
    const std::vector<std::string> sub_model_part_names = mrModelPart.GetSubModelPartNames();
    std::unordered_map<int,std::set<int>> aux_nodes_colors, aux_cond_colors, aux_elem_colors;
    
    std::vector<std::string> model_part_names;
    model_part_names.push_back(mrModelPart.Name());
    for (const auto & sub_model_part_name : sub_model_part_names)
        model_part_names.push_back(sub_model_part_name);
    
    // Initialize Colors
    int color = 0;
    for (SizeType i_sub_model_part = 0; i_sub_model_part < model_part_names.size(); ++i_sub_model_part) {
        rColors[i_sub_model_part].push_back(model_part_names[i_sub_model_part]);
        
        if (color > 0) {            
            ModelPart& r_sub_model_part = mrModelPart.GetSubModelPart(model_part_names[i_sub_model_part]);

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
            NodesColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            NodesColors[key] = *value.begin();
        else // There is a combination
            NodesColors[key] = combinations[value];
    }
    
    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) {
        const int key = aux_cond_color.first;
        const std::set<int>& value = aux_cond_color.second;
        
        if (value.size() == 0)
            CondColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            CondColors[key] = *value.begin();
        else // There is a combination
            CondColors[key] = combinations[value];
    }
    
    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) {
        const int key = aux_elem_color.first;
        const std::set<int>& value = aux_elem_color.second;
        
        if (value.size() == 0)
            ElemColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            ElemColors[key] = *value.begin();
        else // There is a combination
            ElemColors[key] = combinations[value];
    }
}

  
}  // namespace Kratos.


