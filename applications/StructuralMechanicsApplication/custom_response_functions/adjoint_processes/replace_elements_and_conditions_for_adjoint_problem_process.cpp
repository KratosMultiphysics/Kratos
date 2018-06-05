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

namespace Kratos
{


    ReplaceElementsAndConditionsForAdjointProblemProcess::ReplaceElementsAndConditionsForAdjointProblemProcess(ModelPart& model_part) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
        
        KRATOS_CATCH("")
    }


    /// Destructor.
    ReplaceElementsAndConditionsForAdjointProblemProcess::~ReplaceElementsAndConditionsForAdjointProblemProcess()
    {

    }

    void ReplaceElementsAndConditionsForAdjointProblemProcess::Execute()
    {
        ModelPart& r_root_model_part = ObtainRootModelPart( mr_model_part );

    #pragma omp parallel for
        for(int i=0; i< (int)r_root_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator it = r_root_model_part.ElementsBegin() + i;

            std::string element_name;
            bool replace_element = GetNewElementName(*it, element_name);

            if(replace_element)
            {
                if (element_name == "AdjointFiniteDifferencingBaseElement")
                {
                    KRATOS_ERROR_IF_NOT( KratosComponents< Element >::Has( element_name ) )
                        << "Element name not found in KratosComponents< Element > -- name is " << element_name << std::endl;

                    Element::Pointer p_element = Kratos::make_shared<AdjointFiniteDifferencingBaseElement>(
                        it->Id(), it->pGetGeometry(), it->pGetProperties(), *it.base() );

                    //deep copy elemental data
                    p_element->Data() = it->Data();

                    (*it.base()) = p_element;
                }
                else
                {
                    KRATOS_ERROR_IF_NOT( KratosComponents< Element >::Has( element_name ) )
                        << "Element name not found in KratosComponents< Element > -- name is " << element_name << std::endl;
                    const Element& rReferenceElement = KratosComponents<Element>::Get(element_name);

                    Element::Pointer p_element = rReferenceElement.Create(it->Id(), it->GetGeometry().Points(), it->pGetProperties());

                    //deep copy elemental data
                    p_element->Data() = it->Data();

                    (*it.base()) = p_element;
                }
            }
        }

    #pragma omp parallel for
        for(int i=0; i< (int)r_root_model_part.Conditions().size(); i++)
        {
            ModelPart::ConditionsContainerType::iterator it = r_root_model_part.ConditionsBegin() + i;

            std::string condition_name;
            bool replace_condition = GetNewConditionName(*it, condition_name);

            if(replace_condition)
            {
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

    bool ReplaceElementsAndConditionsForAdjointProblemProcess::GetNewElementName(const Element& rElement, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_element;
        CompareElementsAndConditionsUtility::GetRegisteredName(rElement, name_current_element);

        // Add here all new adjoint elements or elements which should be ignored by the replacement process
        if(name_current_element == "CrLinearBeamElement3D2N")
            rName = "CrLinearBeamAdjointElement3D2N";
        else if(name_current_element == "CrLinearBeamAdjointElement3D2N")
            rName = "CrLinearBeamElement3D2N";
        else if(name_current_element == "ShellThinElement3D3N")
            rName = "AdjointFiniteDifferencingBaseElement";
        else if(name_current_element == "AdjointFiniteDifferencingBaseElement")
            rName = "ShellThinElement3D3N";
        else
            KRATOS_ERROR << "It is not possible to replace the " << name_current_element <<
             " because there is no equivalent adjoint/primal element available." << std::endl;

        return replacement_necessary;

        KRATOS_CATCH("")
    }

    bool ReplaceElementsAndConditionsForAdjointProblemProcess::GetNewConditionName(const Condition& rCondition, std::string& rName)
    {
        KRATOS_TRY

        bool replacement_necessary = true;
        std::string name_current_condition;
        CompareElementsAndConditionsUtility::GetRegisteredName(rCondition, name_current_condition);

        // Add here all new adjoint conditions or conditions which should be ignored by the replacement process
        if(name_current_condition == "PointLoadCondition2D1N")
            rName = "PointLoadAdjointCondition2D1N";
        else if(name_current_condition == "PointLoadCondition3D1N")
            rName = "PointLoadAdjointCondition3D1N";
        else if(name_current_condition == "PointLoadAdjointCondition2D1N")
            rName = "PointLoadCondition2D1N";
        else if(name_current_condition == "PointLoadAdjointCondition3D1N")
            rName = "PointLoadCondition3D1N";
        else if(name_current_condition == "ShapeOptimizationCondition")
            replacement_necessary = false;
        else
            KRATOS_ERROR << "It is not possible to replace the " << name_current_condition <<
             " because there is no equivalent adjoint/primal condition available." << std::endl;


        return replacement_necessary;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.




