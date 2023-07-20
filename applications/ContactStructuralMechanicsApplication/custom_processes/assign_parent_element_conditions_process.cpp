// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_processes/assign_parent_element_conditions_process.h"

namespace Kratos
{
void AssignParentElementConditionsProcess::ExecuteInitialize()
{
    KRATOS_TRY;

    // Clear
    mMapFaceIds.clear();

    // The array of elements
    auto& r_elements_array = mrElementsModelPart.Elements();
    for (auto& r_element : r_elements_array) {
        const std::size_t id = r_element.Id();
        auto& r_element_geometry = r_element.GetGeometry();
        const auto faces = r_element_geometry.GenerateBoundariesEntities();
        for (auto& r_face : faces) {
            const std::size_t face_size = r_face.size();
            std::vector<std::size_t> face_indexes(face_size);
            for (std::size_t i_node = 0; i_node < face_size; ++i_node) {
                face_indexes[i_node] = r_face[i_node].Id();
            }
            std::sort(face_indexes.begin(), face_indexes.end());
            mMapFaceIds.insert({face_indexes, id});
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void AssignParentElementConditionsProcess::ExecuteInitializeSolutionStep() 
{
    KRATOS_TRY;

    // We iterate over the conditions
    auto& r_conditions_array = mrConditionsModelPart.Conditions();
    block_for_each(r_conditions_array, [&](Condition& rCondition) {
        auto& r_condition_geometry = rCondition.GetGeometry();
        const std::size_t condition_size = r_condition_geometry.size();
        std::vector<std::size_t> indexes(condition_size);
        for (std::size_t i_node = 0; i_node < condition_size; ++i_node) {
            indexes[i_node] = r_condition_geometry[i_node].Id();
        }
        std::sort(indexes.begin(), indexes.end());
        
        // Find and assign if found
        auto it_find = mMapFaceIds.find(indexes);
        if (it_find != mMapFaceIds.end()) {
            rCondition.SetValue(PARENT_ELEMENT, mrElementsModelPart.pGetElement(it_find->second));
        } else {
            KRATOS_WARNING("AssignParentElementConditionsProcess") << "Condition " << rCondition.Id() << " has not been found in the elements model part" << std::endl;
        }
    });

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

void AssignParentElementConditionsProcess::Execute()
{
    KRATOS_TRY;

    // Execute the initialize of the elements faces database
    ExecuteInitialize();

    // Execute the assign of elements to conditions
    ExecuteInitializeSolutionStep();
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters AssignParentElementConditionsProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "conditions_model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "elements_model_part_name"   : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"                 : 0
    })");

    return default_parameters;
}

} // namespace Kratos.