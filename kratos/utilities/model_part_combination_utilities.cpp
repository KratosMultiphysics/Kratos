//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "utilities/model_part_combination_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
void ModelPartCombinationUtilities::CombineModelParts(Parameters ThisParameters)
{
    // Ensuring parameters
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    // Retrieve the new ModelPart name 
    //const auto& r_new_model_part_name = ThisParameters["combined_model_part_name"].GetString();

    // Create the new ModelPart
    //auto& r_combined_model_part = mrModel.CreateModelPart(r_new_model_part_name, ThisParameters["buffer_size"].GetInt());

    // Retrieve the ModelParts to combine
    std::vector<std::string> model_parts_names;
    const auto& r_model_parts_names = ThisParameters["model_parts_list"];

    // Fill the list
    for (auto& r_name : r_model_parts_names) {
        model_parts_names.push_back(r_name.GetString());
    }

    // Before combine the ModelParts we need to check that the submodelparts are not repeated
    CheckSubModelParts(model_parts_names);

    // Reorder Ids before combining model parts
    ReorderIds(model_parts_names);

    // Finally we delete the old model parts
    for (auto& r_name : model_parts_names) {
        mrModel.DeleteModelPart(r_name);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartCombinationUtilities::CheckSubModelParts(const std::vector<std::string>& rModelPartsNames)
{
    // The list of submodelparts names
    std::unordered_map<std::string, std::size_t> list_of_submodelparts;

    // // Interate over the ModelParts
    // for (auto& r_name : rModelPartsNames) {
    //     const auto& r_model_part = mrModel.GetModelPart(r_name);
    // }

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
    // // Interate over the ModelParts
    // for (auto& r_name : rModelPartsNames) {
    //     const auto& r_model_part = mrModel.GetModelPart(r_name);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ModelPartCombinationUtilities::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_parts_list"         : [],
        "combined_model_part_name" : "CombinedModelParts",
        "buffer_size"              : 2
    })" );
    return default_parameters;
}

}  // namespace Kratos.
