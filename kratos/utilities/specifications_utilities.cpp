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
    const auto& r_elements_array = rModelPart.Elements()
    if (r_elements_array.size() > 0) {
        
        const auto it_elem_begin = r_elements_array.begin();
        const auto elements_components = KratosComponents<Element>::GetComponents();
        
        specifications = it_element->GetSpecifications();
        
        Parameters variables = MaterialData["Variables"];
        for (auto iter = variables.begin(); iter != variables.end(); ++iter) {
            const Parameters value = variables.GetValue(iter.name());

            const std::string variable_name = iter.name();
            
            bool variable_is_missing = false;
            if (KratosComponents<Variable<double> >::Has(variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
                const Variable<bool>& r_variable = KratosComponents<Variable<bool>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
                const Variable<int>& r_variable = KratosComponents<Variable<int>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
                const Variable<array_1d<double, 6>>& r_variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
                variable_is_missing = rModelPart.HasNodalSolutionStepVariable(r_variable);
            } else {
                KRATOS_ERROR << "Value type for \"" << variable_name << "\" not defined";
            }
        }
        
        // Now we iterate over all the elements
        for(std::size_t i = 1; i < r_elements_array.size(); i++) {
            const auto it_elem_previous = it_elem_begin + i - 1;
            const auto it_elem_current = it_elem_begin + i;
            
            if(!GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
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


    KRATOS_CATCH("")
}

} // namespace SpecificationsUtilities
} // namespace Kratos
