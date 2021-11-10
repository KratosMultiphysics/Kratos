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
#include "utilities/specifications_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
namespace SpecificationsUtilities
{
void AddMissingVariables(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        AddMissingVariablesFromSpecifications(rModelPart, specifications, element_name);

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                AddMissingVariablesFromSpecifications(rModelPart, specifications, element_name);
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        AddMissingVariablesFromSpecifications(rModelPart, specifications, condition_name);

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                AddMissingVariablesFromSpecifications(rModelPart, specifications, condition_name);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AddMissingVariablesFromSpecifications(
    ModelPart& rModelPart,
    const Parameters SpecificationsParameters,
    const std::string EntityName
    )
{
    KRATOS_TRY

    if (SpecificationsParameters.Has("required_variables")) {
        const std::vector<std::string>& r_variables = SpecificationsParameters["required_variables"].GetStringArray();
        auto& r_variable_list = rModelPart.GetNodalSolutionStepVariablesList();
        for (auto& r_variable_name : r_variables) {
            bool has_variable = false;
            if (KratosComponents<Variable<double> >::Has(r_variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);
                }
            } else if(KratosComponents<Variable<bool> >::Has(r_variable_name)) {
                const Variable<bool>& r_variable = KratosComponents<Variable<bool>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else if(KratosComponents<Variable<int> >::Has(r_variable_name)) {
                const Variable<int>& r_variable = KratosComponents<Variable<int>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(r_variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(r_variable_name)) {
                const Variable<array_1d<double, 6>>& r_variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else if(KratosComponents<Variable<Vector > >::Has(r_variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else if(KratosComponents<Variable<Matrix> >::Has(r_variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>().Get(r_variable_name);
                has_variable = rModelPart.HasNodalSolutionStepVariable(r_variable);

                // If variable is missign is added to the model part
                if (!has_variable) {
                    r_variable_list.Add(r_variable);;
                }
            } else {
                KRATOS_ERROR << "Value type for \"" << r_variable_name << "\" not defined" << std::endl;
            }
            KRATOS_WARNING_IF("SpecificationsUtilities", !has_variable) << "Variable:" << r_variable_name << " is not in the model part. Required by entity: " << EntityName << ". Added to the model part" << std::endl;
        }

        // Set the variable list in all the nodes
        auto& r_nodes_array = rModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (r_nodes_array.size()); ++k) {
            auto it_node = it_node_begin + k;
            it_node->SetSolutionStepVariablesList(&r_variable_list);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AddMissingDofs(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        AddMissingDofsFromSpecifications(rModelPart, specifications, element_name);

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                AddMissingDofsFromSpecifications(rModelPart, specifications, element_name);
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        AddMissingDofsFromSpecifications(rModelPart, specifications, condition_name);

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                AddMissingDofsFromSpecifications(rModelPart, specifications, condition_name);
            }
        }
    }


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AddMissingDofsFromSpecifications(
    ModelPart& rModelPart,
    const Parameters SpecificationsParameters,
    const std::string EntityName
    )
{
    KRATOS_TRY

    // First node iterator
    const auto it_node_begin = rModelPart.Nodes().begin();

    if (SpecificationsParameters.Has("required_dofs")) {
        const std::vector<std::string>& r_variables = SpecificationsParameters["required_dofs"].GetStringArray();
        for (auto& r_variable_name : r_variables) {
            bool has_dof = false;
            if (KratosComponents<Variable<double> >::Has(r_variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(r_variable_name);
                has_dof = it_node_begin->HasDofFor(r_variable);

                // If variable is missign is added to the model part
                if (!has_dof) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else {
                KRATOS_ERROR << "Value type for \"" << r_variable_name << "\" not defined" << std::endl;
            }
            KRATOS_WARNING_IF("SpecificationsUtilities", !has_dof) << "Variable:" << r_variable_name << " is not in the model part nodes. Required by entity: " << EntityName << ". Added to the model part nodes" << std::endl;
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DetermineFlagsUsed(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("flags_used")) {
            const std::vector<std::string>& r_flags_used = specifications["flags_used"].GetStringArray();
            for (auto& r_flag_name : r_flags_used) {
                KRATOS_INFO("SpecificationsUtilities") << "Flag:" << r_flag_name << " is used by the element: " << element_name << std::endl;
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                if (specifications.Has("flags_used")) {
                    const std::vector<std::string>& r_flags_used = specifications["flags_used"].GetStringArray();
                    for (auto& r_flag_name : r_flags_used) {
                        KRATOS_INFO("SpecificationsUtilities") << "Flag:" << r_flag_name << " is used by the element: " << element_name << std::endl;
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("flags_used")) {
            const std::vector<std::string>& r_flags_used = specifications["flags_used"].GetStringArray();
            for (auto& r_flag_name : r_flags_used) {
                KRATOS_INFO("SpecificationsUtilities") << "Flag:" << r_flag_name << " is used by the condition: " << condition_name << std::endl;
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                if (specifications.Has("flags_used")) {
                    const std::vector<std::string>& r_flags_used = specifications["flags_used"].GetStringArray();
                    for (auto& r_flag_name : r_flags_used) {
                        KRATOS_INFO("SpecificationsUtilities") << "Flag:" << r_flag_name << " is used by the condition: " << condition_name << std::endl;
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> DetermineTimeIntegration(ModelPart& rModelPart)
{
    KRATOS_TRY

    std::unordered_map<std::string, std::size_t> time_integrations_map;
    std::vector<std::string> time_integrations;

    // Auxiliar classes
    std::size_t counter = 0;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        ++counter;
        if (specifications.Has("time_integration")) {
            const std::vector<std::string> time_integrations = specifications["time_integration"].GetStringArray();
            for (auto& r_time_integration : time_integrations) {
                time_integrations_map.insert(std::pair<std::string, std::size_t>(r_time_integration, 1));
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                ++counter;
                if (specifications.Has("time_integration")) {
                    const std::vector<std::string> time_integrations = specifications["time_integration"].GetStringArray();
                    for (auto& r_time_integration : time_integrations) {
                        if (time_integrations_map.find(r_time_integration) != time_integrations_map.end()) {
                            time_integrations_map[r_time_integration] += 1;
                        } else {
                            time_integrations_map.insert(std::pair<std::string, std::size_t>(r_time_integration, 1));
                        }
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        ++counter;
        if (specifications.Has("time_integration")) {
            const std::vector<std::string> time_integrations = specifications["time_integration"].GetStringArray();
            for (auto& r_time_integration : time_integrations) {
                if (time_integrations_map.find(r_time_integration) != time_integrations_map.end()) {
                    time_integrations_map[r_time_integration] += 1;
                } else {
                    time_integrations_map.insert(std::pair<std::string, std::size_t>(r_time_integration, 1));
                }
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                ++counter;
                if (specifications.Has("time_integration")) {
                    const std::vector<std::string> time_integrations = specifications["time_integration"].GetStringArray();
                    for (auto& r_time_integration : time_integrations) {
                        if (time_integrations_map.find(r_time_integration) != time_integrations_map.end()) {
                            time_integrations_map[r_time_integration] += 1;
                        } else {
                            time_integrations_map.insert(std::pair<std::string, std::size_t>(r_time_integration, 1));
                        }
                    }
                }
            }
        }
    }

    // We create the final method
    for (auto& r_time_integration : time_integrations_map) {
        if (r_time_integration.second == counter) {
            time_integrations.push_back(r_time_integration.first);
        }
    }

    return time_integrations;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string DetermineFramework(ModelPart& rModelPart)
{
    KRATOS_TRY

    std::string framework = "NONE";

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("framework")) {
            const std::string& r_framework = specifications["framework"].GetString();
            if (framework == "NONE") {
                framework = r_framework;
            } else {
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                KRATOS_WARNING("SpecificationsUtilities") << "Framework:" << r_framework << " is used by the element: " << element_name << ". But also the framework: " << framework << " is used" << std::endl;
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("framework")) {
                    const std::string& r_framework = specifications["framework"].GetString();
                    if (framework == "NONE") {
                        framework = r_framework;
                    } else {
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "Framework:" << r_framework << " is used by the element: " << element_name << ". But also the framework: " << framework << " is used" << std::endl;
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("framework")) {
            const std::string& r_framework = specifications["framework"].GetString();
            if (framework == "NONE") {
                framework = r_framework;
            } else {
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                KRATOS_WARNING("SpecificationsUtilities") << "Framework:" << r_framework << " is used by the condition: " << condition_name << ". But also the framework: " << framework << " is used" << std::endl;
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("framework")) {
                    const std::string& r_framework = specifications["framework"].GetString();
                    if (framework == "NONE") {
                        framework = r_framework;
                    } else {
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "Framework:" << r_framework << " is used by the condition: " << condition_name << ". But also the framework: " << framework << " is used" << std::endl;
                    }
                }
            }
        }
    }

    return framework;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool DetermineSymmetricLHS(ModelPart& rModelPart)
{
    KRATOS_TRY

    bool symmetric = true;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("symmetric_lhs")) {
            const bool symmetric_lhs = specifications["symmetric_lhs"].GetBool();
            if (symmetric != symmetric_lhs) {
                if (symmetric) {
                    symmetric = symmetric_lhs;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines non-symmetric LHS" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines symmetric LHS, but at least one element requires non-symmetric LHS. Non-symmetric LHS will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("symmetric_lhs")) {
                    const bool symmetric_lhs = specifications["symmetric_lhs"].GetBool();
                    if (symmetric != symmetric_lhs) {
                        if (symmetric) {
                            symmetric = symmetric_lhs;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines non-symmetric LHS" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines symmetric LHS, but at least one element requires non-symmetric LHS. Non-symmetric LHS will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("symmetric_lhs")) {
            const bool symmetric_lhs = specifications["symmetric_lhs"].GetBool();
            if (symmetric != symmetric_lhs) {
                if (symmetric) {
                    symmetric = symmetric_lhs;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines non-symmetric LHS" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines symmetric LHS, but at least one condition requires non-symmetric LHS. Non-symmetric LHS will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("symmetric_lhs")) {
                    const bool symmetric_lhs = specifications["symmetric_lhs"].GetBool();
                    if (symmetric != symmetric_lhs) {
                        if (symmetric) {
                            symmetric = symmetric_lhs;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines non-symmetric LHS" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines symmetric LHS, but at least one condition requires non-symmetric LHS. Non-symmetric LHS will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    return symmetric;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool DeterminePositiveDefiniteLHS(ModelPart& rModelPart)
{
    KRATOS_TRY

    bool positive_definite = true;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("positive_definite_lhs")) {
            const bool positive_definite_lhs = specifications["positive_definite_lhs"].GetBool();
            if (positive_definite != positive_definite_lhs) {
                if (positive_definite) {
                    positive_definite = positive_definite_lhs;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines non-positive definite" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines positive definite, but at least one element requires non-positive definite. Non-positive definite will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("positive_definite_lhs")) {
                    const bool positive_definite_lhs = specifications["positive_definite_lhs"].GetBool();
                    if (positive_definite != positive_definite_lhs) {
                        if (positive_definite) {
                            positive_definite = positive_definite_lhs;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines non-positive definite" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines positive definite, but at least one element requires non-positive definite. Non-positive definite will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("positive_definite_lhs")) {
            const bool positive_definite_lhs = specifications["positive_definite_lhs"].GetBool();
            if (positive_definite != positive_definite_lhs) {
                if (positive_definite) {
                    positive_definite = positive_definite_lhs;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines non-positive definite" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines positive definite, but at least one condition requires non-positive definite. Non-positive definite will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("positive_definite_lhs")) {
                    const bool positive_definite_lhs = specifications["positive_definite_lhs"].GetBool();
                    if (positive_definite != positive_definite_lhs) {
                        if (positive_definite) {
                            positive_definite = positive_definite_lhs;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines non-positive definite" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines positive definite, but at least one condition requires non-positive definite. Non-positive definite will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    return positive_definite;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool DetermineIfCompatibleGeometries(ModelPart& rModelPart)
{
    KRATOS_TRY

    bool compatible_geometries = true;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        // Getting geometry type
        const auto& r_geometry = it_elem_begin->GetGeometry();
        const auto& r_geometry_type = r_geometry.GetGeometryType();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("compatible_geometries")) {
            const std::vector<std::string> compatible_geometries = specifications["compatible_geometries"].GetStringArray();
            bool compatible_geometry = false;
            for (auto& r_geometry_name : compatible_geometries) {
                if (r_geometry_type == string_geometry_map[r_geometry_name]) {
                    compatible_geometry = true;
                }
            }
            if (!compatible_geometry) {
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "is not compatible with the following geometries:\n";
                for (auto& r_geometry_name : compatible_geometries) {
                    KRATOS_WARNING("")  << r_geometry_name << "\n";
                }
                return false;
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);

                // Getting geometry type
                const auto& r_geometry = it_elem_current->GetGeometry();
                const auto& r_geometry_type = r_geometry.GetGeometryType();

                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                if (specifications.Has("compatible_geometries")) {
                    const std::vector<std::string> compatible_geometries = specifications["compatible_geometries"].GetStringArray();
                    bool compatible_geometry = false;
                    for (auto& r_geometry_name : compatible_geometries) {
                        if (r_geometry_type == string_geometry_map[r_geometry_name]) {
                            compatible_geometry = true;
                        }
                    }
                    if (!compatible_geometry) {
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "is not compatible with the following geometries:\n";
                        for (auto& r_geometry_name : compatible_geometries) {
                            KRATOS_WARNING("")  << r_geometry_name << "\n";
                        }
                        return false;
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condition has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        // Getting geometry type
        const auto& r_geometry = it_cond_begin->GetGeometry();
        const auto& r_geometry_type = r_geometry.GetGeometryType();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("compatible_geometries")) {
            const std::vector<std::string> compatible_geometries = specifications["compatible_geometries"].GetStringArray();
            bool compatible_geometry = false;
            for (auto& r_geometry_name : compatible_geometries) {
                if (r_geometry_type == string_geometry_map[r_geometry_name]) {
                    compatible_geometry = true;
                }
            }
            if (!compatible_geometry) {
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                KRATOS_WARNING("SpecificationsUtilities") << "The element: " << condition_name << "is not compatible with the following geometries:\n";
                for (auto& r_geometry_name : compatible_geometries) {
                    KRATOS_WARNING("")  << r_geometry_name << "\n";
                }
                return false;
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);

                // Getting geometry type
                const auto& r_geometry = it_cond_current->GetGeometry();
                const auto& r_geometry_type = r_geometry.GetGeometryType();

                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                if (specifications.Has("compatible_geometries")) {
                    const std::vector<std::string> compatible_geometries = specifications["compatible_geometries"].GetStringArray();
                    bool compatible_geometry = false;
                    for (auto& r_geometry_name : compatible_geometries) {
                        if (r_geometry_type == string_geometry_map[r_geometry_name]) {
                            compatible_geometry = true;
                        }
                    }
                    if (!compatible_geometry) {
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "The element: " << condition_name << "is not compatible with the following geometries:\n";
                        for (auto& r_geometry_name : compatible_geometries) {
                            KRATOS_WARNING("")  << r_geometry_name << "\n";
                        }
                        return false;
                    }
                }
            }
        }
    }

    return compatible_geometries;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool DetermineIfRequiresTimeIntegration(ModelPart& rModelPart)
{
    KRATOS_TRY

    bool requires_time_integration = true;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("element_integrates_in_time")) {
            const bool element_integrates_in_time = specifications["element_integrates_in_time"].GetBool();
            if (requires_time_integration != element_integrates_in_time) {
                if (requires_time_integration) {
                    requires_time_integration = element_integrates_in_time;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines does not require time integration" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines requires time integration, but at least one element does not require time integration. Explicit formulation will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("element_integrates_in_time")) {
                    const bool element_integrates_in_time = specifications["element_integrates_in_time"].GetBool();
                    if (requires_time_integration != element_integrates_in_time) {
                        if (requires_time_integration) {
                            requires_time_integration = element_integrates_in_time;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines does not require time integration" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "defines requires time integration, but at least one element does not require time integration. Explicit formulation will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("element_integrates_in_time")) {
            const bool element_integrates_in_time = specifications["element_integrates_in_time"].GetBool();
            if (requires_time_integration != element_integrates_in_time) {
                if (requires_time_integration) {
                    requires_time_integration = element_integrates_in_time;
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines does not require time integration" << std::endl;
                } else {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines requires time integration, but at least one condition does not require time integration. Explicit formulation will be considered" << std::endl;
                }
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("element_integrates_in_time")) {
                    const bool element_integrates_in_time = specifications["element_integrates_in_time"].GetBool();
                    if (requires_time_integration != element_integrates_in_time) {
                        if (requires_time_integration) {
                            requires_time_integration = element_integrates_in_time;
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines does not require time integration" << std::endl;
                        } else {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "defines requires time integration, but at least one condition does not require time integration. Explicit formulation will be considered" << std::endl;
                        }
                    }
                }
            }
        }
    }

    return requires_time_integration;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

bool CheckCompatibleConstitutiveLaws(ModelPart& rModelPart)
{
    KRATOS_TRY

    bool compatible_cl = true;

    // Auxiliar things
    std::vector<ConstitutiveLaw::Pointer> cl_vector;
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    ConstitutiveLaw::Features features;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        // Getting CL vectors
        cl_vector.clear();
        it_elem_begin->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, cl_vector, r_process_info);

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("compatible_constitutive_laws")) {
            const std::vector<std::string> compatible_constitutive_laws_type = specifications["compatible_constitutive_laws"]["type"].GetStringArray();
            const std::vector<std::string> compatible_constitutive_laws_dimension = specifications["compatible_constitutive_laws"]["dimension"].GetStringArray();
            const Vector compatible_constitutive_laws_strain_size = specifications["compatible_constitutive_laws"]["strain_size"].GetVector();
            const std::size_t number_of_cl_compatible = compatible_constitutive_laws_strain_size.size();
            // If the vector is full of nullptr we just return true
            if (cl_vector.size() > 0) 
                if (cl_vector[0] == nullptr) 
                    return true; 
            for (auto& p_cl : cl_vector) {
                bool check = false;
                for (std::size_t i = 0; i < number_of_cl_compatible; ++i) {
                    const std::string& r_type = compatible_constitutive_laws_type[i];
                    const std::string& r_dimension = compatible_constitutive_laws_dimension[i];
                    const std::size_t strain_size = compatible_constitutive_laws_strain_size[i];

                    // Getting features
                    p_cl->GetLawFeatures(features);

                    // Getting flag
                    const Flags* p_flag_considered = r_type == "ThreeDimensional" ?  &ConstitutiveLaw::THREE_DIMENSIONAL_LAW : r_type == "PlaneStrain" ? &ConstitutiveLaw::PLANE_STRAIN_LAW : r_type == "PlaneStress" ? &ConstitutiveLaw::PLANE_STRESS_LAW : &ConstitutiveLaw::AXISYMMETRIC_LAW; // NOTE: "Axisymmetric"

                    if (p_cl->GetStrainSize() == strain_size &&
                        p_cl->WorkingSpaceDimension() == string_dimension_map[r_dimension] &&
                        features.mOptions.Is( *p_flag_considered )
                    ) {
                        check = true;
                    }
                }
                if (!check) {
                    CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                    KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << " is considering a not compatible CL" << std::endl;
                    return false;
                }
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                if (specifications.Has("compatible_constitutive_laws")) {
                    const std::vector<std::string> compatible_constitutive_laws_type = specifications["compatible_constitutive_laws"]["type"].GetStringArray();
                    const std::vector<std::string> compatible_constitutive_laws_dimension = specifications["compatible_constitutive_laws"]["dimension"].GetStringArray();
                    const Vector compatible_constitutive_laws_strain_size = specifications["compatible_constitutive_laws"]["strain_size"].GetVector();
                    const std::size_t number_of_cl_compatible = compatible_constitutive_laws_strain_size.size();
                    // If the vector is full of nullptr we just return true
                    if (cl_vector.size() > 0) 
                        if (cl_vector[0] == nullptr) 
                            return true; 
                    for (auto& p_cl : cl_vector) {
                        bool check = false;
                        for (std::size_t i = 0; i < number_of_cl_compatible; ++i) {
                            const std::string& r_type = compatible_constitutive_laws_type[i];
                            const std::string& r_dimension = compatible_constitutive_laws_dimension[i];
                            const std::size_t strain_size = compatible_constitutive_laws_strain_size[i];

                            // Getting features
                            p_cl->GetLawFeatures(features);

                            // Getting flag
                            const Flags* p_flag_considered = r_type == "ThreeDimensional" ?  &ConstitutiveLaw::THREE_DIMENSIONAL_LAW : r_type == "PlaneStrain" ? &ConstitutiveLaw::PLANE_STRAIN_LAW : r_type == "PlaneStress" ? &ConstitutiveLaw::PLANE_STRESS_LAW : &ConstitutiveLaw::AXISYMMETRIC_LAW; // NOTE: "Axisymmetric"

                            if (p_cl->GetStrainSize() == strain_size &&
                                p_cl->WorkingSpaceDimension() == string_dimension_map[r_dimension] &&
                                features.mOptions.Is( *p_flag_considered )
                            ) {
                                check = true;
                            }
                        }
                        if (!check) {
                            CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                            KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << " is considering a not compatible CL" << std::endl;
                            return false;
                        }
                    }
                }
            }
        }
    }

    return compatible_cl;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int CheckGeometricalPolynomialDegree(ModelPart& rModelPart)
{
    KRATOS_TRY

    int polynomial_degree_of_geometry = -1;

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();

        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("required_polynomial_degree_of_geometry")) {
            const int required_polynomial_degree_of_geometry = specifications["required_polynomial_degree_of_geometry"].GetInt();
            if (required_polynomial_degree_of_geometry > polynomial_degree_of_geometry) {
                polynomial_degree_of_geometry = required_polynomial_degree_of_geometry;
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "requires a polynomial degree of geometry of: " << required_polynomial_degree_of_geometry << std::endl;
            }
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;

            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("required_polynomial_degree_of_geometry")) {
                    const int required_polynomial_degree_of_geometry = specifications["required_polynomial_degree_of_geometry"].GetInt();
                    if (required_polynomial_degree_of_geometry > polynomial_degree_of_geometry) {
                        polynomial_degree_of_geometry = required_polynomial_degree_of_geometry;
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "The element: " << element_name << "requires a polynomial degree of geometry of: " << required_polynomial_degree_of_geometry << std::endl;
                    }
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();

        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("required_polynomial_degree_of_geometry")) {
            const int required_polynomial_degree_of_geometry = specifications["required_polynomial_degree_of_geometry"].GetInt();
            if (required_polynomial_degree_of_geometry > polynomial_degree_of_geometry) {
                polynomial_degree_of_geometry = required_polynomial_degree_of_geometry;
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "requires a polynomial degree of geometry of: " << required_polynomial_degree_of_geometry << std::endl;
            }
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;

            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("required_polynomial_degree_of_geometry")) {
                    const int required_polynomial_degree_of_geometry = specifications["required_polynomial_degree_of_geometry"].GetInt();
                    if (required_polynomial_degree_of_geometry > polynomial_degree_of_geometry) {
                        polynomial_degree_of_geometry = required_polynomial_degree_of_geometry;
                        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);
                        KRATOS_WARNING("SpecificationsUtilities") << "The condition: " << condition_name << "requires a polynomial degree of geometry of: " << required_polynomial_degree_of_geometry << std::endl;
                    }
                }
            }
        }
    }

    // Final message 
    KRATOS_WARNING_IF("SpecificationsUtilities", polynomial_degree_of_geometry > 0) << "Finally, the polynomial degree of geometry considered is: " << polynomial_degree_of_geometry << std::endl;

    return polynomial_degree_of_geometry;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Parameters GetDocumention(ModelPart& rModelPart)
{
    KRATOS_TRY

    Parameters documentation = Parameters(R"({})");

    // Define specifications
    Parameters specifications;

    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    const auto& r_elements_array = rModelPart.Elements();
    if (r_elements_array.size() > 0) {
        std::string element_name;
        const auto it_elem_begin = r_elements_array.begin();
        specifications = it_elem_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
        if (specifications.Has("documentation")) {
            const Parameters& r_documentation = specifications["documentation"];
            documentation.AddValue(element_name, r_documentation);
        }

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;
            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_begin, element_name);
                if (specifications.Has("documentation")) {
                    const Parameters& r_documentation = specifications["documentation"];
                    documentation.AddValue(element_name, r_documentation);
                }
            }
        }
    }

    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        std::string condition_name;
        const auto it_cond_begin = r_conditions_array.begin();
        specifications = it_cond_begin->GetSpecifications();
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
        if (specifications.Has("documentation")) {
            const Parameters& r_documentation = specifications["documentation"];
            documentation.AddValue(condition_name, r_documentation);
        }

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;
            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_begin, condition_name);
                if (specifications.Has("documentation")) {
                    const Parameters& r_documentation = specifications["documentation"];
                    documentation.AddValue(condition_name, r_documentation);
                }
            }
        }
    }

    return documentation;

    KRATOS_CATCH("")
}

} // namespace SpecificationsUtilities
} // namespace Kratos
