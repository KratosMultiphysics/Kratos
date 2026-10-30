// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    E. Wehrle
//                   based on max stress nad previous MR to aggregation from Daniel Baumgaertner, https://github.com/dbaumgaertner
//

// System includes

// External includes

// Project includes
#include "adjoint_KS_max_stress_response_function.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    AdjointKSMaxStressResponseFunction::AdjointKSMaxStressResponseFunction(ModelPart& rAdjointModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rAdjointModelPart, ResponseSettings),
      mrAdjointModelPart(rAdjointModelPart),
      mCriticalPartName(ResponseSettings["critical_part_name"].GetString())
    {
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment(ResponseSettings["stress_treatment"].GetString());

        pKS = ResponseSettings["aggregation_penalty"].GetDouble();

        if(ResponseSettings.Has("echo_level"))
            mEchoLevel = ResponseSettings["echo_level"].GetInt();

        KRATOS_ERROR_IF(mStressTreatment != StressTreatment::Mean)
            << "AdjointKSMaxStressResponseFunction::AdjointKSMaxStressResponseFunction: Specified stress treatment not supported: " << ResponseSettings["stress_type"].GetString() << std::endl;

        ModelPart& adjoint_aggregation_part = mrAdjointModelPart.GetSubModelPart(mCriticalPartName);
        for(auto& elem : adjoint_aggregation_part.Elements())
            mAggregatedElementIds.push_back(elem.Id());

        // Initialize
        for(auto& elem : adjoint_aggregation_part.Elements())
        {
            elem.SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType));
        }
    }

    // primal analysis for max stress approximated via modified Kreisselmeier-Steinhauser aggregation
    AdjointKSMaxStressResponseFunction::~AdjointKSMaxStressResponseFunction(){}

    double AdjointKSMaxStressResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;

        ModelPart& primal_agglomeration_part = rPrimalModelPart.GetSubModelPart(mCriticalPartName);

        IndexType elem_id_at_max = 0;
        max_mean_stress = 0.0;

        for(auto& elem : primal_agglomeration_part.Elements())
        {
            Vector element_stress;
            const ProcessInfo &r_process_info = rPrimalModelPart.GetProcessInfo();
            StressCalculation::CalculateStressOnGP(elem,
                                                   mTracedStressType,
                                                   element_stress,
                                                   r_process_info);

            const SizeType stress_vec_size = element_stress.size();
            double mean_stress = 0.0;

            for(IndexType i = 0; i < stress_vec_size; ++i)
                mean_stress += element_stress[i];
            mean_stress /= stress_vec_size;

            if(mean_stress > max_mean_stress)
            {
                max_mean_stress = mean_stress;
                elem_id_at_max = elem.Id();
            }
            mean_stress_vector[elem.Id()] = mean_stress;
        }

        // reloop over elements, faster as vector operation?
        for(auto& elem : primal_agglomeration_part.Elements())
        {
            KS_exp_sum += std::exp(pKS*(mean_stress_vector[elem.Id()]-max_mean_stress));
        }
        double KS_max_mean_stress = max_mean_stress+1/pKS*std::log(KS_exp_sum);

        KRATOS_INFO_IF("AdjointKSMaxStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with max stress value = " << elem_id_at_max << std::endl;
        KRATOS_INFO_IF("AdjointKSMaxStressResponseFunction::CalculateValue", mEchoLevel > 0) << "max mean stress approximated with KS aggregation = " << KS_max_mean_stress << std::endl;

        // Set traced element in adjoint model part corresponding to the found element in the primal part
        mpTracedElementInAdjointPart = mrAdjointModelPart.pGetElement(elem_id_at_max);
        mpTracedElementInAdjointPart->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        return KS_max_mean_stress;

        KRATOS_CATCH("");
    }

    // adjoint equation for max stress approximated via modified Kreisselmeier-Steinhauser aggregation
    void AdjointKSMaxStressResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                               const Matrix& rResidualGradient,
                                                               Vector& rResponseGradient,
                                                               const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        //KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
        //    << "AdjointKSMaxStressResponseFunction::CalculateGradient: Traced element not initialized. First call \"CalculateValue\"!" << std::endl;


        rResponseGradient.clear();

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {

            Matrix stress_displacement_derivative;

            Element::Pointer pElem = mrAdjointModelPart.pGetElement(rAdjointElement.Id());
            pElem->Calculate(STRESS_DISP_DERIV_ON_GP,
                             stress_displacement_derivative,
                             rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative,
                                              rResponseGradient);

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                << "AdjointKSMaxStressResponseFunction::CalculateGradient: Size of stress displacement derivative does not fit!" << std::endl;

            //double mean_stress = mean_stress_vector[rAdjointElement.Id()];
            rResponseGradient *= (-1) *std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum;
        }
        else
        {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                         const Variable<double>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement,
                                                                   rVariable.Name(),
                                                                   rSensitivityMatrix,
                                                                   rSensitivityGradient,
                                                                   rProcessInfo);
            rSensitivityGradient *= std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum;
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("")
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                         const Variable<double>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                         const Variable<array_1d<double,
                                                                         3>>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;


         if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement,
                                                                   rVariable.Name(),
                                                                   rSensitivityMatrix,
                                                                   rSensitivityGradient,
                                                                   rProcessInfo);

            rSensitivityGradient *= std::exp(pKS*(mean_stress_vector[rAdjointElement.Id()]-max_mean_stress))/KS_exp_sum;

        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                         const Variable<array_1d<double,
                                                                         3>>& rVariable,
                                                                         const Matrix& rSensitivityMatrix,
                                                                         Vector& rSensitivityGradient,
                                                                         const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                                                                              const std::string& rVariableName,
                                                                                              const Matrix& rSensitivityMatrix,
                                                                                              Vector& rSensitivityGradient,
                                                                                              const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP,
                                  stress_design_variable_derivative,
                                  rProcessInfo);

        this->ExtractMeanStressDerivative(stress_design_variable_derivative,
                                          rSensitivityGradient);

        KRATOS_ERROR_IF(rSensitivityGradient.size() != rSensitivityMatrix.size1()) << "Size of partial stress design variable derivative does not fit!" << std::endl;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, "");

        KRATOS_CATCH("");
    }

    void AdjointKSMaxStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix,
                                                                         Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();
        double stress_derivative_value = 0.0;

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
        {
            for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
                stress_derivative_value += rStressDerivativesMatrix(deriv_it, stress_it);

            stress_derivative_value /= num_of_stress_positions;

            rResponseGradient[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.
