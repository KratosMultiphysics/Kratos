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

// System includes

// External includes

// Project includes
#include "adjoint_min_stress_response_function.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    // --------------------------------------------------------------------------
    AdjointMinStressResponseFunction::AdjointMinStressResponseFunction(ModelPart& rAdjointModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rAdjointModelPart, ResponseSettings),
      mrAdjointModelPart(rAdjointModelPart),
      mCriticalPartName(ResponseSettings["critical_part_name"].GetString())
    {
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment(ResponseSettings["stress_treatment"].GetString());

        if(ResponseSettings.Has("echo_level"))
            mEchoLevel = ResponseSettings["echo_level"].GetInt();

        KRATOS_ERROR_IF(mStressTreatment != StressTreatment::Mean)
            << "AdjointMinStressResponseFunction::AdjointMinStressResponseFunction: Specified stress treatment not supported: " << ResponseSettings["stress_type"].GetString() << std::endl;
    }

    // --------------------------------------------------------------------------
    AdjointMinStressResponseFunction::~AdjointMinStressResponseFunction(){}

    // --------------------------------------------------------------------------
    double AdjointMinStressResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;

        ModelPart& primal_agglomeration_part = rPrimalModelPart.GetSubModelPart(mCriticalPartName);

        // Find min mean value and corresponding element in primal model part
        double min_mean_stress = 1.0e12;
        IndexType elem_id_at_min = 0;

        for(auto& elem : primal_agglomeration_part.Elements())
        {
            Vector element_stress;
            StressCalculation::CalculateStressOnGP(elem, mTracedStressType, element_stress, rPrimalModelPart.GetProcessInfo());

            const SizeType stress_vec_size = element_stress.size();
            double mean_stress = 0.0;
            for(IndexType i = 0; i < stress_vec_size; ++i)
                mean_stress += element_stress[i];
            mean_stress /= stress_vec_size;

            if(mean_stress < min_mean_stress)
            {
                min_mean_stress = mean_stress;
                elem_id_at_min = elem.Id();
            }
        }

        KRATOS_INFO_IF("AdjointMinStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with min stress value = " << elem_id_at_min << std::endl;
        KRATOS_INFO_IF("AdjointMinStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Min mean stress value = " << min_mean_stress << std::endl;

        // Set traced element in adjoint model part corresponding to the found element in the primal part
        mpTracedElementInAdjointPart = mrAdjointModelPart.pGetElement(elem_id_at_min);
        mpTracedElementInAdjointPart->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        return min_mean_stress;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                             const Matrix& rResidualGradient,
                                                             Vector& rResponseGradient,
                                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointMinStressResponseFunction::CalculateGradient: Traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            Matrix stress_displacement_derivative;

            mpTracedElementInAdjointPart->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                << "AdjointMinStressResponseFunction::CalculateGradient: Size of stress displacement derivative does not fit!" << std::endl;

            rResponseGradient *= (-1);
        }
        else
        {
            if(rResponseGradient.size() != rResidualGradient.size1())
                rResponseGradient.resize(rResidualGradient.size1(), false);

            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<double>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointMinStressResponseFunction::CalculateGradient: traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("")
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                       const Variable<double>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<array_1d<double, 3>>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointMinStressResponseFunction::CalculateGradient: traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                                                       const Variable<array_1d<double, 3>>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                      const std::string& rVariableName,
                                      const Matrix& rSensitivityMatrix,
                                      Vector& rSensitivityGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        rAdjointElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
        this->ExtractMeanStressDerivative(stress_design_variable_derivative, rSensitivityGradient);

        KRATOS_ERROR_IF(rSensitivityGradient.size() != rSensitivityMatrix.size1()) << "Size of partial stress design variable derivative does not fit!" << std::endl;

        rAdjointElement.SetValue(DESIGN_VARIABLE_NAME, "");

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointMinStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
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

