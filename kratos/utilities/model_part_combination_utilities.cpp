//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/kernel.h"
#include "utilities/model_part_combination_utilities.h"
#include "utilities/parallel_utilities.h"
#include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos
{
ModelPart& ModelPartCombinationUtilities::CombineModelParts(Parameters ThisParameters)
{
    // Ensuring parameters
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    // Getting echo level
    const int echo_level = ThisParameters["echo_level"].GetInt();

    // Retrieve the ModelParts to combine
    std::vector<std::string> model_parts_names;
    const auto& r_model_parts_names = ThisParameters["model_parts_list"];

    // Fill the list
    KRATOS_INFO_IF("ModelPartCombinationUtilities", echo_level > 0) << "The following ModelParts are going to be merged" << std::endl;
    for (auto& r_name : r_model_parts_names) {
        const auto& name = r_name.GetString();
        model_parts_names.push_back(name);
    }
    KRATOS_ERROR_IF(model_parts_names.size() == 0) << "Empty list of ModelParts" << std::endl;

    // Giving information from current ModelParts
    if (echo_level > 0) {
        for (auto& r_name : model_parts_names) {
            auto& r_model_part = mrModel.GetModelPart(r_name);
            KRATOS_INFO("ModelPartCombinationUtilities") << r_model_part << std::endl;
        }
    }

    // Retrieve the new ModelPart name 
    const auto& r_new_model_part_name = ThisParameters["combined_model_part_name"].GetString();

    // Create the new ModelPart
    auto& r_combined_model_part = mrModel.HasModelPart(r_new_model_part_name) ? mrModel.GetModelPart(r_new_model_part_name) : mrModel.CreateModelPart(r_new_model_part_name, ThisParameters["buffer_size"].GetInt());

    // Serial check
    KRATOS_ERROR_IF(mrModel.GetModelPart(mrModel.GetModelPartNames()[0]).IsDistributed()) << "ModelPartCombinationUtilities is only compatible in serial simulations" << std::endl;

    // Before combine the ModelParts we need to check that the submodelparts are not repeated
    CheckSubModelParts(model_parts_names);

    // Reorder Ids before combining model parts
    ReorderIds(model_parts_names);

    // Finally we combine the model parts
    for (auto& r_name : model_parts_names) {
        auto& r_model_part = mrModel.GetModelPart(r_name);   

        // Copy variable list from the first ModelPart
        for (auto& r_var : r_model_part.GetNodalSolutionStepVariablesList()) {
            r_combined_model_part.AddNodalSolutionStepVariable(r_var);
        }

        RecursiveAddEntities(r_combined_model_part, r_model_part);
    }
    // Update database
    r_combined_model_part.SetNodalSolutionStepVariablesList();

    // Print info
    KRATOS_INFO_IF("ModelPartCombinationUtilities", echo_level > 0) << "The merge resulting in the following ModelPart:\n" <<  r_combined_model_part << std::endl;

    // Finally we delete the old model parts
    for (auto& r_name : model_parts_names) {
        mrModel.DeleteModelPart(r_name);
    }

    // Clean up properties
    const int clean_up_properties = ThisParameters["clean_up_properties"].GetBool();
    if (clean_up_properties) {
        for (auto it_prop = r_combined_model_part.PropertiesBegin() + 1; it_prop < r_combined_model_part.PropertiesEnd(); it_prop++) {
            r_combined_model_part.RemovePropertiesFromAllLevels(*(it_prop.base()));
        }
        auto pProperties = *(r_combined_model_part.PropertiesBegin().base());
        pProperties->SetId(0);
        
        // Iterate over conditions
        block_for_each(r_combined_model_part.Conditions(), [&pProperties](Condition& rCondition){
            rCondition.SetProperties(pProperties);
        });

        // Iterate over elements
        block_for_each(r_combined_model_part.Elements(), [&pProperties](Element& rElement){
            rElement.SetProperties(pProperties);
        });
    }

    return r_combined_model_part;
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartCombinationUtilities::CheckSubModelParts(const std::vector<std::string>& rModelPartsNames)
{
    // The list of submodelparts names
    std::unordered_map<std::string, std::size_t> list_of_submodelparts;

    // Interate over the ModelParts
    for (auto& r_name : rModelPartsNames) {
        auto& r_model_part = mrModel.GetModelPart(r_name);

        RecursiveAddOfModelPartsToList(r_model_part, list_of_submodelparts);
    }

    // Check that the submodelparts are not repeated
    bool check = false;
    for (auto& r_map : list_of_submodelparts) {
        if (r_map.second > 1) {
            check = true;
            KRATOS_WARNING("ModelPartCombinationUtilities") << "The submodelpart " << r_map.first << " is defined in different ModelParts and cannot be combined" << std::endl;
        }
    }
    KRATOS_ERROR_IF(check) << "Please check the names of the submodelparts to avoid conflicts in the names" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartCombinationUtilities::ReorderIds(const std::vector<std::string>& rModelPartsNames)
{
    // Prepare the number of entities in model parts
    const std::size_t number_of_model_parts = rModelPartsNames.size();
    std::vector<std::size_t> number_geometries(number_of_model_parts, 0);
    std::vector<std::size_t> number_nodes(number_of_model_parts, 0);
    std::vector<std::size_t> number_conditions(number_of_model_parts, 0);
    std::vector<std::size_t> number_elements(number_of_model_parts, 0);
    std::vector<std::size_t> number_constraints(number_of_model_parts, 0);
    std::vector<std::size_t> number_properties(number_of_model_parts, 0);

    // Interate over the ModelParts
    for (std::size_t i_mp = 0; i_mp < number_of_model_parts - 1; i_mp++) {
        const auto& r_model_part = mrModel.GetModelPart(rModelPartsNames[i_mp]);
        number_geometries[i_mp + 1] = number_geometries[i_mp] + r_model_part.NumberOfGeometries();
        number_nodes[i_mp + 1] = number_nodes[i_mp] + r_model_part.NumberOfNodes();
        number_elements[i_mp + 1] = number_elements[i_mp] + r_model_part.NumberOfElements();
        number_conditions[i_mp + 1] = number_conditions[i_mp] + r_model_part.NumberOfConditions();
        number_constraints[i_mp + 1] = number_constraints[i_mp] + r_model_part.NumberOfMasterSlaveConstraints();
        number_properties[i_mp + 1] = number_properties[i_mp] + r_model_part.NumberOfProperties();
    }

    // Reorder Ids
    for (std::size_t i_mp = 0; i_mp < number_of_model_parts; i_mp++) {
        auto& r_model_part = mrModel.GetModelPart(rModelPartsNames[i_mp]);

        // Iterate over geometries
        auto& r_geometries_array = r_model_part.Geometries();
        const auto it_geom_begin = r_geometries_array.begin();
        const std::size_t geom_initial_index = number_geometries[i_mp];
        IndexPartition<std::size_t>(r_geometries_array.size()).for_each([&it_geom_begin, &geom_initial_index](std::size_t i){
            auto it_geom = it_geom_begin;
            for (std::size_t j = 0; j < i; ++j) ++it_geom;
            it_geom->SetId(geom_initial_index + i + 1);
        });

        // Iterate over nodes
        auto& r_nodes_array = r_model_part.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        const std::size_t node_initial_index = number_nodes[i_mp];
        IndexPartition<std::size_t>(r_nodes_array.size()).for_each([&it_node_begin, &node_initial_index](std::size_t i){
            (it_node_begin + i)->SetId(node_initial_index + i + 1);
        });

        // Iterate over conditions
        auto& r_conditions_array = r_model_part.Conditions();
        const auto it_cond_begin = r_conditions_array.begin();
        const std::size_t cond_initial_index = number_conditions[i_mp];
        IndexPartition<std::size_t>(r_conditions_array.size()).for_each([&it_cond_begin, &cond_initial_index](std::size_t i){
            (it_cond_begin + i)->SetId(cond_initial_index + i + 1);
        });

        // Iterate over elements
        auto& r_elements_array = r_model_part.Elements();
        const auto it_elem_begin = r_elements_array.begin();
        const std::size_t elem_initial_index = number_elements[i_mp];
        IndexPartition<std::size_t>(r_elements_array.size()).for_each([&it_elem_begin, &elem_initial_index](std::size_t i){
            (it_elem_begin + i)->SetId(elem_initial_index + i + 1);
        });

        // Iterate over constraints
        auto& r_constraints_array = r_model_part.MasterSlaveConstraints();
        const auto it_const_begin = r_constraints_array.begin();
        const std::size_t const_initial_index = number_constraints[i_mp];
        IndexPartition<std::size_t>(r_constraints_array.size()).for_each([&it_const_begin, &const_initial_index](std::size_t i){
            (it_const_begin + i)->SetId(const_initial_index + i + 1);
        });

        // Iterate over properties
        auto& r_properties_array = r_model_part.rProperties();
        const auto it_prop_begin = r_properties_array.begin();
        const std::size_t prop_initial_index = number_properties[i_mp];
        IndexPartition<std::size_t>(r_properties_array.size()).for_each([&it_prop_begin, &prop_initial_index](std::size_t i){
            (it_prop_begin + i)->SetId(prop_initial_index + i + 1);
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartCombinationUtilities::RecursiveAddOfModelPartsToList(
    ModelPart& rModelPart,
    std::unordered_map<std::string, std::size_t>& rListModelParts
    )
{
    // Lambda to extend the map of model parts
    auto extend_map = [&rListModelParts](ModelPart& rModelPart) {
        // Check it already exists
        if (rListModelParts.find(rModelPart.Name()) != rListModelParts.end()) {
            rListModelParts[rModelPart.Name()] += 1;
        } else {
            rListModelParts.insert({rModelPart.Name(), 1});
        }
    };

    // Recursively add of ModelParts to the list
    if (rModelPart.NumberOfSubModelParts() > 0) {
        for (auto& r_model_part : rModelPart.SubModelParts()) {
            RecursiveAddOfModelPartsToList(r_model_part, rListModelParts);
        }
    }
    extend_map(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartCombinationUtilities::RecursiveAddEntities(
    ModelPart& rDestinationModelPart,
    ModelPart& rOriginModelPart
    )
{
    // Lambda to transfer the model parts
    auto transfer_lambda = [](ModelPart& rDestinationModelPart, ModelPart& rOriginModelPart) {
        FastTransferBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart).Execute();
        
        // Copy properties
        for (auto it_prop = rOriginModelPart.PropertiesBegin(); it_prop < rOriginModelPart.PropertiesEnd(); it_prop++) {
            if (!rDestinationModelPart.HasProperties(it_prop->Id())) {
                rDestinationModelPart.AddProperties(*(it_prop.base()));
            }
        }
    };

    // Recursively add of ModelParts to the list
    if (rOriginModelPart.NumberOfSubModelParts() > 0) {
        for (auto& r_sub_model_part : rOriginModelPart.SubModelParts()) {
            auto& r_sub_destination_model_part = rDestinationModelPart.CreateSubModelPart(r_sub_model_part.Name());
            RecursiveAddEntities(r_sub_destination_model_part, r_sub_model_part);
        }
    }
    transfer_lambda(rDestinationModelPart, rOriginModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ModelPartCombinationUtilities::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_parts_list"         : [],
        "combined_model_part_name" : "CombinedModelParts",
        "buffer_size"              : 2,
        "echo_level"               : 0,
        "clean_up_properties"      : true
    })" );
    return default_parameters;
}

}  // namespace Kratos.
