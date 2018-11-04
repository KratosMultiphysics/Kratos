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
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_shell_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "custom_response_functions/adjoint_conditions/adjoint_semi_analytic_base_condition.h"
#include "custom_response_functions/adjoint_conditions/adjoint_semi_analytic_point_load_condition.h"

namespace Kratos
{


    ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::Execute()
    {
        KRATOS_ERROR_IF(mrModelPart.IsSubModelPart()) << "The replacement process can only be done for the root model part!" << std::endl;

        if ( (!mrModelPart.GetProcessInfo().Has(IS_ADJOINT)) || (!mrModelPart.GetProcessInfo()[IS_ADJOINT]) )
        {
            this->ReplaceToAdjoint();
            mrModelPart.GetProcessInfo()[IS_ADJOINT] = true;
        }
        else
        {
            this->ReplaceToPrimal();
            mrModelPart.GetProcessInfo()[IS_ADJOINT] = false;
        }
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceToPrimal()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            const auto it = mrModelPart.ElementsBegin() + i;

            AdjointFiniteDifferencingBaseElement::Pointer p_adjoint_element = dynamic_pointer_cast<AdjointFiniteDifferencingBaseElement>(*it.base());
            if (p_adjoint_element != nullptr)
            {
                (*it.base()) = p_adjoint_element->pGetPrimalElement();
            }
        }

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.NumberOfConditions()); ++i)
        {
            const auto it = mrModelPart.ConditionsBegin() + i;

            AdjointSemiAnalyticBaseCondition::Pointer p_adjoint_condition = dynamic_pointer_cast<AdjointSemiAnalyticBaseCondition>(*it.base());
            if (p_adjoint_condition != nullptr)
            {
                (*it.base()) = p_adjoint_condition->pGetPrimalCondition();
            }
        }

        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = mrModelPart.SubModelPartsBegin(); i_sub_model_part != mrModelPart.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, mrModelPart );
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceToAdjoint()
    {

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            const auto it = mrModelPart.ElementsBegin() + i;

            std::string element_name;
            const bool replace_element = GetAdjointElementName(*it, element_name);

            if(replace_element)
            {
                if (element_name == "AdjointFiniteDifferencingShellElement")
                {
                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferencingShellElement>(*it.base() );

                    (*it.base()) = p_element;
                }
                else if (element_name == "AdjointFiniteDifferenceCrBeamElement")
                {
                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferenceCrBeamElement>(*it.base() );

                    (*it.base()) = p_element;
                }
                else if (element_name == "AdjointFiniteDifferenceTrussElement")
                {
                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferenceTrussElement>(*it.base() );

                    (*it.base()) = p_element;
                }
                else if (element_name == "AdjointFiniteDifferenceTrussLinearElement")
                {
                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferenceTrussElementLinear>(*it.base() );

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
                if (condition_name == "AdjointSemiAnalyticPointLoadCondition")
                {
                    Condition::Pointer p_condition = Kratos::make_shared<AdjointSemiAnalyticPointLoadCondition>(*it.base() );

                    (*it.base()) = p_condition;
                }
                else
                    KRATOS_ERROR << "Unknown adjoint condition: " << condition_name << std::endl;
            }
        }

        //change the sons
        for (ModelPart::SubModelPartIterator i_sub_model_part = mrModelPart.SubModelPartsBegin(); i_sub_model_part != mrModelPart.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, mrModelPart );
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

    void ReplaceElementsAndConditionsForAdjointProblemProcess::UpdateSubModelPart(ModelPart& r_model_part, ModelPart& r_root_model_part)
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

    bool ReplaceElementsAndConditionsForAdjointProblemProcess::GetAdjointElementName(const Element& rElement, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_element;
        CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name_current_element);

        // Add here all new adjoint elements or elements which should be ignored by the replacement process
        if(name_current_element == "CrLinearBeamElement3D2N")
            rName = "AdjointFiniteDifferenceCrBeamElement";
        else if(name_current_element == "ShellThinElement3D3N")
            rName = "AdjointFiniteDifferencingShellElement";
        else if(name_current_element == "TrussElement3D2N")
            rName = "AdjointFiniteDifferenceTrussElement";
        else if(name_current_element == "TrussLinearElement3D2N")
            rName = "AdjointFiniteDifferenceTrussLinearElement";
        else
        {
            KRATOS_ERROR << "It is not possible to replace the " << name_current_element <<
             " because there is no equivalent adjoint/primal element available." << std::endl;
        }

        return replacement_necessary;

        KRATOS_CATCH("")
    }

    bool ReplaceElementsAndConditionsForAdjointProblemProcess::GetAdjointConditionName(const Condition& rCondition, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_condition;
        CompareElementsAndConditionsUtility::GetRegisteredName(rCondition, name_current_condition);

        // Add here all new adjoint conditions or conditions which should be ignored by the replacement process
        if(name_current_condition == "PointLoadCondition2D1N")
            rName = "AdjointSemiAnalyticPointLoadCondition";
        else if(name_current_condition == "PointLoadCondition3D1N")
            rName = "AdjointSemiAnalyticPointLoadCondition";
        else if(name_current_condition == "ShapeOptimizationCondition3D3N")
            replacement_necessary = false;
        else if(name_current_condition == "ShapeOptimizationCondition3D4N")
            replacement_necessary = false;
        else
        {
            KRATOS_ERROR << "It is not possible to replace the " << name_current_condition <<
             " because there is no equivalent adjoint/primal condition available." << std::endl;
        }


        return replacement_necessary;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.




