// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//                  


// Project includes
#include "utilities/compare_elements_and_conditions_utility.h"
#include "replace_elements_and_conditions_for_adjoint_problem_process.h"

namespace Kratos
{


    ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part, 
                              Parameters Settings
                                   ) : Process(Flags()) , mr_model_part(model_part), mSettings( Settings)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "add_string": "NAME_OF_ADD_STRING",
                "add_before_in_element_name": "ADD_STRING_BEFORE",
                "add_before_in_condition_name": "ADD_STRING_BEFORE",
                "elements_conditions_to_ignore": "NAME_OF_EXEPTION",
                "from_primal_to_adjoint": true
            }  )" );    
        //now validate agains defaults -- this also ensures no type mismatch*/
        Settings.ValidateAndAssignDefaults(default_parameters); 
        
        KRATOS_CATCH("")
    }
    

    /// Destructor.
    ReplaceElementsAndConditionsForAdjointProblemProcess::~ReplaceElementsAndConditionsForAdjointProblemProcess()
    {

    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::Execute() 
    {
        ModelPart& r_root_model_part = ObtainRootModelPart( mr_model_part );

        std::string sub_name_element = mSettings["add_before_in_element_name"].GetString();
        std::string sub_name_condition = mSettings["add_before_in_condition_name"].GetString();
        std::string adding_string = mSettings["add_string"].GetString();
        std::string ignore_string = mSettings["elements_conditions_to_ignore"].GetString();
        bool  from_primal_to_adjoint = mSettings["from_primal_to_adjoint"].GetBool();
        
    #pragma omp parallel for                              
        for(int i=0; i< (int)r_root_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_root_model_part.ElementsBegin() + i;

            std::string element_name; 
            CompareElementsAndConditionsUtility::GetRegisteredName(*it, element_name);
 
            if(!(element_name == ignore_string))
            {
                if(from_primal_to_adjoint) 
                {
                    std::string::size_type position_ele = 0;
                    std::string::size_type found_ele;
                    found_ele = element_name.find(sub_name_element, position_ele);
                    element_name.insert(found_ele, adding_string);
                }
                else
                {
                   auto pos = element_name.find(adding_string);
                   element_name.erase(pos,adding_string.length());   
                }

                KRATOS_ERROR_IF_NOT( KratosComponents< Element >::Has( element_name ) )
                    << "Element name not found in KratosComponents< Element > -- name is " << element_name << std::endl;
                const Element& rReferenceElement = KratosComponents<Element>::Get(element_name); 

                Element::Pointer p_element = rReferenceElement.Create(it->Id(), it->GetGeometry().Points(), it->pGetProperties());
 
                //deep copy elemental data
                p_element->Data() = it->Data();
            
                (*it.base()) = p_element;
            }
        }
        
    #pragma omp parallel for                              
        for(int i=0; i< (int)r_root_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_root_model_part.ConditionsBegin() + i;

            std::string condition_name; 
            CompareElementsAndConditionsUtility::GetRegisteredName(*it, condition_name);

            if(!(condition_name == ignore_string))
            {
                if(from_primal_to_adjoint) 
                {
                    std::string::size_type position_cond = 0;
                    std::string::size_type found_cond;
                    found_cond = condition_name.find(sub_name_condition, position_cond);
                    condition_name.insert(found_cond, adding_string);
                }
                else
                {
                    auto pos = condition_name.find(adding_string);
                    condition_name.erase(pos,adding_string.length());
                }
                KRATOS_ERROR_IF_NOT( KratosComponents< Condition >::Has( condition_name ) )
                    << "Condition name not found in KratosComponents< Condition > -- name is " << condition_name;
                const Condition& rReferenceCondition = KratosComponents<Condition>::Get(condition_name); 

                Condition::Pointer p_condition = rReferenceCondition.Create(it->Id(),it->GetGeometry().Points(), it->pGetProperties());
                
                //deep copy elemental data
                p_condition->Data() = it->Data();
            
                (*it.base()) = p_condition;
            }
        }      
        
        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_root_model_part.SubModelPartsBegin(); i_sub_model_part != r_root_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, r_root_model_part );


    }

    /// Turn back information as a string.
    std::string ReplaceElementsAndConditionsForAdjointProblemProcess::Info() const
    {
        return "ReplaceElementsAndConditionsForAdjointProblemProcess";
    }

    /// Print information about this object.
    void ReplaceElementsAndConditionsForAdjointProblemProcess::PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << "ReplaceElementsAndConditionsForAdjointProblemProcess";
    }

    /// Print object's data.
    void ReplaceElementsAndConditionsForAdjointProblemProcess::PrintData(std::ostream& rOStream) const 
    {
    }
  
    ModelPart& ReplaceElementsAndConditionsForAdjointProblemProcess::ObtainRootModelPart( ModelPart& r_model_part )
    {
        if (r_model_part.IsSubModelPart())
            return ObtainRootModelPart(*r_model_part.GetParentModelPart());
        else
            return r_model_part;
    }
    
    void ReplaceElementsAndConditionsForAdjointProblemProcess::UpdateSubModelPart(ModelPart& r_model_part, ModelPart& r_root_model_part)
    {
        //change the model part itself
        #pragma omp parallel for
        for(int i=0; i< (int)r_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_model_part.ElementsBegin() + i;
            
            (*it.base()) = r_root_model_part.Elements()(it->Id());
        }
        
        #pragma omp parallel for
        for(int i=0; i< (int)r_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;
            
            (*it.base()) = r_root_model_part.Conditions()(it->Id());
        }
            
        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_model_part.SubModelPartsBegin(); i_sub_model_part != r_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, r_root_model_part );
    }

}  // namespace Kratos.




