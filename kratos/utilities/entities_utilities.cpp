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

// System includes

// External includes

// Project includes
#include "utilities/entities_utilities.h"

namespace Kratos
{
namespace EntitiesUtilities
{

void InitializeEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of conditions
    auto& r_conditions_array = rModelPart.Conditions();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());
    const auto it_cond_begin = r_conditions_array.begin();

    // The number of elements
    auto& r_elements_array = rModelPart.Elements();
    const int number_of_elements = static_cast<int>(r_elements_array.size());
    const auto it_elem_begin = r_elements_array.begin();

    // The number of constraints
    auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
    const int number_of_constraints = static_cast<int>(r_constraints_array.size());
    const auto it_const_begin = r_constraints_array.begin();

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Setting to zero the slave dofs
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
            auto it_cond = it_cond_begin + i_cond;

            // Detect if the condition is active or not. If the user did not make any choice the condition
            // It is active by default
            bool condition_is_active = true;
            if (it_cond->IsDefined(ACTIVE))
                condition_is_active = it_cond->Is(ACTIVE);

            if (condition_is_active) {
                it_cond->Initialize(r_current_process_info);
            }
        }

        #pragma omp for schedule(guided, 512)
        for (int i_elem = 0; i_elem < number_of_elements; ++i_elem) {
            auto it_elem = it_elem_begin + i_elem;

            // Detect if the element is active or not. If the user did not make any choice the element
            // It is active by default
            bool element_is_active = true;
            if (it_elem->IsDefined(ACTIVE))
                element_is_active = it_elem->Is(ACTIVE);

            if (element_is_active) {
                it_elem->Initialize(r_current_process_info);
            }
        }

        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = it_const_begin + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                it_const->Initialize(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void InitializeConditions(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of conditions
    auto& r_conditions_array = rModelPart.Conditions();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());
    const auto it_cond_begin = r_conditions_array.begin();

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Initialize
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
            auto it_cond = it_cond_begin + i_cond;

            // Detect if the condition is active or not. If the user did not make any choice the condition
            // It is active by default
            bool condition_is_active = true;
            if (it_cond->IsDefined(ACTIVE))
                condition_is_active = it_cond->Is(ACTIVE);

            if (condition_is_active) {
                it_cond->Initialize(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void InitializeElements(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of elements
    auto& r_elements_array = rModelPart.Elements();
    const int number_of_elements = static_cast<int>(r_elements_array.size());
    const auto it_elem_begin = r_elements_array.begin();

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Initialize
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_elem = 0; i_elem < number_of_elements; ++i_elem) {
            auto it_elem = it_elem_begin + i_elem;

            // Detect if the element is active or not. If the user did not make any choice the element
            // It is active by default
            bool element_is_active = true;
            if (it_elem->IsDefined(ACTIVE))
                element_is_active = it_elem->Is(ACTIVE);

            if (element_is_active) {
                it_elem->Initialize(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void InitializeMasterSlaveConstraints(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of constraints
    auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
    const int number_of_constraints = static_cast<int>(r_constraints_array.size());
    const auto it_const_begin = r_constraints_array.begin();

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Initialize
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = it_const_begin + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                it_const->Initialize(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

} // namespace EntitiesUtilities
} // namespace Kratos
