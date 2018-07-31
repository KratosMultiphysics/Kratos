//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                    
//

// System includes

// External includes

// Project includes
#include "processes/replace_elements_and_condition_process.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::Execute()
{
    ModelPart& r_root_model_part = ObtainRootModelPart( mrModelPart );
    
    const Element& rReferenceElement = KratosComponents<Element>::Get(mSettings["element_name"].GetString());
    const Condition& rReferenceCondition = KratosComponents<Condition>::Get(mSettings["condition_name"].GetString());
    
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(r_root_model_part.Elements().size()); i++) {
        auto it_elem = r_root_model_part.ElementsBegin() + i;
        
        auto p_element = rReferenceElement.Create(it_elem->Id(), it_elem->pGetGeometry(), it_elem->pGetProperties());
        
        // Deep copy elemental data and flags
        p_element->Data() = it_elem->Data();
        p_element->Set(Flags(*it_elem));
        
        (*it_elem.base()) = p_element;
    }
    
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(r_root_model_part.Conditions().size()); i++) {
        auto it_cond = r_root_model_part.ConditionsBegin() + i;
        
        auto p_condition = rReferenceCondition.Create(it_cond->Id(), it_cond->pGetGeometry(), it_cond->pGetProperties());
        
        // Deep copy elemental data and flags
        p_condition->Data() = it_cond->Data();
        p_condition->Set(Flags(*it_cond));
        
        (*it_cond.base()) = p_condition;

    }      
    
    // Change the submodelparts
    for (auto& i_sub_model_part : r_root_model_part.SubModelParts())
        UpdateSubModelPart( i_sub_model_part, r_root_model_part );

}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& ReplaceElementsAndConditionsProcess::ObtainRootModelPart( ModelPart& rModelPart )
{
    if (rModelPart.IsSubModelPart())
        return ObtainRootModelPart(*rModelPart.GetParentModelPart());
    else
        return rModelPart;
}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::UpdateSubModelPart(
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
    for (auto& i_sub_model_part : rModelPart.SubModelParts())
        UpdateSubModelPart( i_sub_model_part, rRootModelPart );

}

}  // namespace Kratos.



