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
#include "../adjoint_elements/adjoint_finite_difference_base_element.h"
#include "../adjoint_elements/adjoint_finite_difference_shell_element.h"
#include "../adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "../adjoint_conditions/adjoint_semi_analytic_point_load_condition.h"

namespace Kratos
{


    ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::Execute()
    {
        if (mr_model_part.IsSubModelPart())
            KRATOS_ERROR << "The replacement process can only be done for the root model part!" << std::endl;

        if ( (!mr_model_part.GetProcessInfo().Has(IS_ADJOINT)) or (!mr_model_part.GetProcessInfo()[IS_ADJOINT]) )
        {
            KRATOS_WATCH("ReplaceToAdjoint")
            this->ReplaceToAdjoint();
            mr_model_part.GetProcessInfo()[IS_ADJOINT] = true;
        }
        else
        {
            KRATOS_WATCH("ReplaceToPrimal")
            this->ReplaceToPrimal();
            mr_model_part.GetProcessInfo()[IS_ADJOINT] = false;
        }
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceToPrimal()
    {
        // TODO
    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceToAdjoint()
    {

        #pragma omp parallel for
        for(int i=0; i< (int)mr_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = mr_model_part.ElementsBegin() + i;

            std::string element_name;
            bool replace_element = GetAdjointElementName(*it, element_name);

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
                else if (element_name == "AdjointFiniteDifferencingBaseElement")
                {
                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferencingBaseElement>(*it.base() );

                    (*it.base()) = p_element;
                }
                else
                    KRATOS_ERROR << "Unknown adjoint element: " << element_name << std::endl;
            }
        }

        #pragma omp parallel for
        for(int i=0; i< (int)mr_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = mr_model_part.ConditionsBegin() + i;

            std::string condition_name;
            bool replace_condition = GetAdjointConditionName(*it, condition_name);

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
        for (ModelPart::SubModelPartIterator i_sub_model_part = mr_model_part.SubModelPartsBegin(); i_sub_model_part != mr_model_part.SubModelPartsEnd(); i_sub_model_part++)
            UpdateSubModelPart( *i_sub_model_part, mr_model_part );
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




