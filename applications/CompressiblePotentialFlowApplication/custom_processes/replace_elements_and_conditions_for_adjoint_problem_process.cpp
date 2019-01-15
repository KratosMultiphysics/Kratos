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
#include "custom_elements/incompressible_adjoint_potential_flow_element.h"
#include "custom_conditions/incompressible_adjoint_potential_wall_condition.h"


namespace Kratos
{


    ReplaceElementsAndConditionsAdjointProcess::ReplaceElementsAndConditionsAdjointProcess(ModelPart& model_part) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    void ReplaceElementsAndConditionsAdjointProcess::Execute()
    {
        KRATOS_ERROR_IF(mrModelPart.IsSubModelPart()) << "The replacement process can only be done for the root model part!" << std::endl;

        // if ( (!mrModelPart.GetProcessInfo().Has(IS_ADJOINT)) || (!mrModelPart.GetProcessInfo()[IS_ADJOINT]) )
        // {
        std::cout<<"starting replace to adjoint"<< std::endl;
        this->ReplaceToAdjoint();
        std::cout<<"replace to adjoint completed"<< std::endl;
            // mrModelPart.GetProcessInfo()[IS_ADJOINT] = true;
        // }
        // else
        // {
        //     this->ReplaceToPrimal();
        //     mrModelPart.GetProcessInfo()[IS_ADJOINT] = false;
        // }
    }

    // void ReplaceElementsAndConditionsAdjointProcess::ReplaceToPrimal()
    // {
    //     #pragma omp parallel for
    //     for(int i=0; i<static_cast<int>(mrModelPart.NumberOfElements()); ++i)
    //     {
    //         const auto it = mrModelPart.ElementsBegin() + i;

    //         AdjointFiniteDifferencingBaseElement::Pointer p_adjoint_element = dynamic_pointer_cast<AdjointFiniteDifferencingBaseElement>(*it.base());
    //         if (p_adjoint_element != nullptr)
    //         {
    //             (*it.base()) = p_adjoint_element->pGetPrimalElement();
    //         }
    //     }

    //     #pragma omp parallel for
    //     for(int i=0; i<static_cast<int>(mrModelPart.NumberOfConditions()); ++i)
    //     {
    //         const auto it = mrModelPart.ConditionsBegin() + i;

    //         AdjointSemiAnalyticBaseCondition::Pointer p_adjoint_condition = dynamic_pointer_cast<AdjointSemiAnalyticBaseCondition>(*it.base());
    //         if (p_adjoint_condition != nullptr)
    //         {
    //             (*it.base()) = p_adjoint_condition->pGetPrimalCondition();
    //         }
    //     }

    //     //change the sons
    //     for (ModelPart::SubModelPartIterator i_sub_model_part = mrModelPart.SubModelPartsBegin(); i_sub_model_part != mrModelPart.SubModelPartsEnd(); i_sub_model_part++)
    //         UpdateSubModelPart( *i_sub_model_part, mrModelPart );
    // }

    void ReplaceElementsAndConditionsAdjointProcess::ReplaceToAdjoint()
    {

                    // constexpr int NumNodes = it->GetGeometry().size();
                    // constexpr int Dim = it->GetGeometry().Dimension();
        const Element& rReferenceElement = IncompressibleAdjointPotentialFlowElement<2,3>();
        const Condition& rReferenceCondition = IncompressibleAdjointPotentialWallCondition<2,2>();

        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            const auto it = mrModelPart.ElementsBegin() + i;

            std::string element_name;
            const bool replace_element = GetAdjointElementName(*it, element_name);
            if(replace_element)
            {
                if (element_name == "IncompressibleAdjointPotentialFlowElement2D3N")
                {   
                    // auto p_element = rReferenceElement.Create(*it.base());
                    Element::Pointer p_element= Kratos::make_shared<IncompressibleAdjointPotentialFlowElement<2,3>>(*it.base());
                    // Deep copy elemental data and flags                
                    p_element->Data() = it->Data();
                    p_element->Set(Flags(*it));
                    (*it.base()) = p_element;

                }                
                else
                    KRATOS_ERROR << "Unknown adjoint element: " << element_name << std::endl;
            }
        }
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.NumberOfConditions()); ++i)
        {
            const auto it = mrModelPart.ConditionsBegin() + i;

            std::string condition_name;
            const bool replace_condition = GetAdjointConditionName(*it, condition_name);

            if(replace_condition)
            {
                if (condition_name == "IncompressibleAdjointPotentialWallCondition2D2N")
                {
                    auto p_condition = rReferenceCondition.Create(it->Id(), it->pGetGeometry(), it->pGetProperties());
                    p_condition->Data() = it->Data();
                    p_condition->Set(Flags(*it));
                    (*it.base()) = p_condition;
                }
                else
                    KRATOS_ERROR << "Unknown adjoint condition: " << condition_name << std::endl;
            }
        }
        std::cout << "All elements and conditions replaced"<<std::endl;
        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = mrModelPart.SubModelPartsBegin(); i_sub_model_part != mrModelPart.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, mrModelPart );
    }

    /// Turn back information as a string.
    std::string ReplaceElementsAndConditionsAdjointProcess::Info() const
    {
        return "ReplaceElementsAndConditionsAdjointProcess";
    }

    /// Print information about this object.
    void ReplaceElementsAndConditionsAdjointProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ReplaceElementsAndConditionsAdjointProcess";
    }

    /// Print object's data.
    void ReplaceElementsAndConditionsAdjointProcess::PrintData(std::ostream& rOStream) const
    {
    }

    void ReplaceElementsAndConditionsAdjointProcess::UpdateSubModelPart(ModelPart& r_model_part, ModelPart& r_root_model_part)
    {
        //change the model part itself
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(r_model_part.NumberOfElements()); ++i)
        {
            const auto it = r_model_part.ElementsBegin() + i;

            (*it.base()) = r_root_model_part.Elements()(it->Id());
        }

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(r_model_part.NumberOfConditions()); ++i)
        {
            const auto it = r_model_part.ConditionsBegin() + i;

            (*it.base()) = r_root_model_part.Conditions()(it->Id());
        }

        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_model_part.SubModelPartsBegin(); i_sub_model_part != r_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, r_root_model_part );
    }

    bool ReplaceElementsAndConditionsAdjointProcess::GetAdjointElementName(const Element& rElement, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_element;
        CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name_current_element);
        
        // Add here all new adjoint elements or elements which should be ignored by the replacement process
        if(name_current_element == "IncompressiblePotentialFlowElement2D3N")
            rName = "IncompressibleAdjointPotentialFlowElement2D3N";
        else
        {
            KRATOS_ERROR << "It is not possible to replace the " << name_current_element <<
             " because there is no equivalent adjoint/primal element available." << std::endl;
        }

        return replacement_necessary;

        KRATOS_CATCH("")
    }

    bool ReplaceElementsAndConditionsAdjointProcess::GetAdjointConditionName(const Condition& rCondition, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_condition;

        CompareElementsAndConditionsUtility::GetRegisteredName(rCondition, name_current_condition);


        // Add here all new adjoint conditions or conditions which should be ignored by the replacement process
        if(name_current_condition == "IncompressiblePotentialWallCondition2D2N")
            rName = "IncompressibleAdjointPotentialWallCondition2D2N";
        else
        {
            KRATOS_ERROR << "It is not possible to replace the " << name_current_condition <<
             " because there is no equivalent adjoint/primal condition available." << std::endl;
        }


        return replacement_necessary;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.




