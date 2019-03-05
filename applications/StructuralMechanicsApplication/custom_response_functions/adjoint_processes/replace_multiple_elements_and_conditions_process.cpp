// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes

// Project includes
#include "replace_multiple_elements_and_conditions_process.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{


void ReplaceMultipleElementsAndConditionsProcess::Execute()
{
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();

    bool throw_not_found_error = mSettings["throw_error"].GetBool();

    // replace elements
    std::map<std::string, const Element*> elements_table;
    Parameters element_name_table = mSettings["element_name_table"];
    for(Parameters::iterator it=element_name_table.begin(); it!=element_name_table.end(); ++it){
        elements_table[it.name()] = &KratosComponents<Element>::Get(element_name_table[it.name()].GetString());
    }

    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(r_root_model_part.Elements().size()); i++) {
        auto it_elem = r_root_model_part.ElementsBegin() + i;

        std::string name_current_element;
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem, name_current_element);

        auto it_reference_element = elements_table.find(name_current_element);

        if (it_reference_element == elements_table.end()) {
            // This error is thrown in a parallel region and can not get catched or even printed properly!
            KRATOS_ERROR_IF(throw_not_found_error) << name_current_element
                << " was not defined in the replacement table!" << std::endl;
            continue; // skip if no error should be thrown
        }

        auto p_element = it_reference_element->second->Create(it_elem->Id(),
            it_elem->pGetGeometry(),
            it_elem->pGetProperties()
        );

        // Deep copy elemental data and flags
        p_element->Data() = it_elem->Data();
        p_element->Set(Flags(*it_elem));

        (*it_elem.base()) = p_element;
    }


    // replace conditions
    std::map<std::string, const Condition*> conditions_table;
    Parameters condition_name_table = mSettings["condition_name_table"];
    for(Parameters::iterator it=condition_name_table.begin(); it!=condition_name_table.end(); ++it){
        conditions_table[it.name()] = &KratosComponents<Condition>::Get(condition_name_table[it.name()].GetString());
    }

    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(r_root_model_part.Conditions().size()); i++) {
        auto it_cond = r_root_model_part.ConditionsBegin() + i;

        std::string name_current_condition;
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond, name_current_condition);

        auto it_reference_condition = conditions_table.find(name_current_condition);

        if (it_reference_condition == conditions_table.end()) {
            // This error is thrown in a parallel region and can not get catched or even printed properly!
            KRATOS_ERROR_IF(throw_not_found_error) << name_current_condition
                << " was not defined in the replacement table!" << std::endl;
            continue; // skip if no error should be thrown
        }

        auto p_condition = it_reference_condition->second->Create(it_cond->Id(),
            it_cond->pGetGeometry(),
            it_cond->pGetProperties()
        );

        // Deep copy conditional data and flags
        p_condition->Data() = it_cond->Data();
        p_condition->Set(Flags(*it_cond));

        (*it_cond.base()) = p_condition;
    }


    // Change the submodelparts
    for (auto& i_sub_model_part : r_root_model_part.SubModelParts()) {
        UpdateSubModelPart( i_sub_model_part, r_root_model_part );
    }
}


void ReplaceMultipleElementsAndConditionsProcess::UpdateSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart
    )
{
    // Change the model part itself
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(rModelPart.Elements().size()); i++) {
        auto it_elem = rModelPart.ElementsBegin() + i;

        (*it_elem.base()) = rRootModelPart.Elements()(it_elem->Id());
    }

    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(rModelPart.Conditions().size()); i++) {
        auto it_cond = rModelPart.ConditionsBegin() + i;

        (*it_cond.base()) = rRootModelPart.Conditions()(it_cond->Id());
    }

    // Change the submodelparts
    for (auto& i_sub_model_part : rModelPart.SubModelParts()) {
        UpdateSubModelPart( i_sub_model_part, rRootModelPart );
    }
}

}  // namespace Kratos.



