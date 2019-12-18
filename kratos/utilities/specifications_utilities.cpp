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
        const auto it_elem_begin = r_elements_array.begin();
        const auto elements_components = KratosComponents<Element>::GetComponents();
        
        specifications = it_elem_begin->GetSpecifications();
        AddMissingVariablesFromSpecifications(rModelPart, specifications);
            
        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;
            
            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                AddMissingVariablesFromSpecifications(rModelPart, specifications);
            }
        }
    }
    
    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        const auto it_cond_begin = r_conditions_array.begin();
        const auto conditions_components = KratosComponents<Condition>::GetComponents();
        
        specifications = it_cond_begin->GetSpecifications();
        AddMissingVariablesFromSpecifications(rModelPart, specifications);
            
        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;
            
            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                AddMissingVariablesFromSpecifications(rModelPart, specifications);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AddMissingVariablesFromSpecifications(
    ModelPart& rModelPart,
    const Parameters SpecificationsParameters 
    )
{
    KRATOS_TRY

    if (SpecificationsParameters.Has("required_variables")) {
        const std::vector<std::string>& r_variables = SpecificationsParameters["required_variables"].GetStringArray();
        for (auto& r_variable_name : r_variables) {
            bool variable_is_missing = false;
            if (KratosComponents<Variable<double> >::Has(r_variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<bool> >::Has(r_variable_name)) {
                const Variable<bool>& r_variable = KratosComponents<Variable<bool>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<int> >::Has(r_variable_name)) {
                const Variable<int>& r_variable = KratosComponents<Variable<int>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(r_variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(r_variable_name)) {
                const Variable<array_1d<double, 6>>& r_variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<Vector > >::Has(r_variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else if(KratosComponents<Variable<Matrix> >::Has(r_variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>().Get(r_variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
                
                // If variable is missign is added to the model part
                if (variable_is_missing) {
                    rModelPart.AddNodalSolutionStepVariable(r_variable);
                }
            } else {
                KRATOS_ERROR << "Value type for \"" << r_variable_name << "\" not defined";
            }
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
        const auto it_elem_begin = r_elements_array.begin();
        const auto elements_components = KratosComponents<Element>::GetComponents();
        
        specifications = it_elem_begin->GetSpecifications();
        AddMissingDofsFromSpecifications(rModelPart, specifications);
            
        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;
            
            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                specifications = it_elem_current->GetSpecifications();
                AddMissingDofsFromSpecifications(rModelPart, specifications);
            }
        }
    }
    
    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condent has changed
    const auto& r_conditions_array = rModelPart.Conditions();
    if (r_conditions_array.size() > 0) {
        const auto it_cond_begin = r_conditions_array.begin();
        const auto conditions_components = KratosComponents<Condition>::GetComponents();
        
        specifications = it_cond_begin->GetSpecifications();
        AddMissingDofsFromSpecifications(rModelPart, specifications);
            
        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < r_conditions_array.size(); i++) {
            const auto it_cond_previous = it_cond_begin + i - 1;
            const auto it_cond_current = it_cond_begin + i;
            
            if(!GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                specifications = it_cond_current->GetSpecifications();
                AddMissingDofsFromSpecifications(rModelPart, specifications);
            }
        }
    }


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AddMissingDofsFromSpecifications(
    ModelPart& rModelPart,
    const Parameters SpecificationsParameters 
    )
{
    KRATOS_TRY

    // First node iterator
    const auto it_node_begin = rModelPart.Nodes().begin();
    
    if (SpecificationsParameters.Has("required_dofs")) {
        const std::vector<std::string>& r_variables = SpecificationsParameters["required_dofs"].GetStringArray();
        for (auto& r_variable_name : r_variables) {
            bool dof_is_missing = false;
            if (KratosComponents<Variable<double> >::Has(r_variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(r_variable_name);
                dof_is_missing = it_node_begin->HasDofFor(r_variable);
                
                // If variable is missign is added to the model part
                if (dof_is_missing) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else if (KratosComponents<Component3VarType>::Has(r_variable_name)) {
                const Component3VarType& r_variable = KratosComponents<Component3VarType>().Get(r_variable_name);
                dof_is_missing = it_node_begin->HasDofFor(r_variable);
                
                // If variable is missign is added to the model part
                if (dof_is_missing) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else if (KratosComponents<Component4VarType>::Has(r_variable_name)) {
                const Component4VarType& r_variable = KratosComponents<Component4VarType>().Get(r_variable_name);
                dof_is_missing = it_node_begin->HasDofFor(r_variable);
                
                // If variable is missign is added to the model part
                if (dof_is_missing) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else if (KratosComponents<Component6VarType>::Has(r_variable_name)) {
                const Component6VarType& r_variable = KratosComponents<Component6VarType>().Get(r_variable_name);
                dof_is_missing = it_node_begin->HasDofFor(r_variable);
                
                // If variable is missign is added to the model part
                if (dof_is_missing) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else if (KratosComponents<Component9VarType>::Has(r_variable_name)) {
                const Component9VarType& r_variable = KratosComponents<Component9VarType>().Get(r_variable_name);
                dof_is_missing = it_node_begin->HasDofFor(r_variable);
                
                // If variable is missign is added to the model part
                if (dof_is_missing) {
                    VariableUtils().AddDof(r_variable, rModelPart);
                }
            } else {
                KRATOS_ERROR << "Value type for \"" << r_variable_name << "\" not defined";
            }
        }
    }

    KRATOS_CATCH("")
}

} // namespace SpecificationsUtilities
} // namespace Kratos
