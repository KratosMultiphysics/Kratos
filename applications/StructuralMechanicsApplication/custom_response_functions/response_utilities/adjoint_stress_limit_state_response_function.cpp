// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//

// System includes

// External includes

// Project includes
#include "adjoint_stress_limit_state_response_function.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    // --------------------------------------------------------------------------
    AdjointStressLimitStateResponseFunction::AdjointStressLimitStateResponseFunction(ModelPart& rAdjointModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rAdjointModelPart, ResponseSettings),
      mrAdjointModelPart(rAdjointModelPart),
      mCriticalPartName(ResponseSettings["critical_part_name"].GetString())
    {
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment(ResponseSettings["stress_treatment"].GetString());

        if(ResponseSettings.Has("echo_level"))
            mEchoLevel = ResponseSettings["echo_level"].GetInt();

        KRATOS_ERROR_IF(mStressTreatment != StressTreatment::Mean)
            << "AdjointStressLimitStateResponseFunction::AdjointStressLimitStateResponseFunction: Specified stress treatment not supported: " << ResponseSettings["stress_type"].GetString() << std::endl;
    }

    // --------------------------------------------------------------------------
    AdjointStressLimitStateResponseFunction::~AdjointStressLimitStateResponseFunction(){}

    // --------------------------------------------------------------------------
    double AdjointStressLimitStateResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;

        ModelPart& primal_agglomeration_part = rPrimalModelPart.GetSubModelPart(mCriticalPartName);

        // Find max mean value and corresponding element in primal model part
        double max_stress_utilization = std::numeric_limits<double>::max();
        IndexType elem_id_at_max = 0;

        for(auto& elem : primal_agglomeration_part.Elements())
        {
            KRATOS_ERROR_IF_NOT(elem.Has(YIELD_STRESS)) << 
                                "YIELD_STRESS not available for element with ID " << elem.Id() << ".\n";
            const auto& yield_stress = elem.GetValue(YIELD_STRESS);

            Vector element_stress;
            const ProcessInfo &r_process_info = rPrimalModelPart.GetProcessInfo();
            StressCalculation::CalculateStressOnGP(elem, mTracedStressType, element_stress, r_process_info);

            const SizeType stress_vec_size = element_stress.size();
            double mean_stress = 0.0;
            for(IndexType i = 0; i < stress_vec_size; ++i)
                mean_stress += element_stress[i];
            mean_stress /= stress_vec_size;

            double stress_limit_state_value = 0.0;
            stress_limit_state_value = yield_stress - mean_stress;

            if(max_stress_utilization > stress_limit_state_value)
            {
                max_stress_utilization = stress_limit_state_value;
                elem_id_at_max = elem.Id();
            }
        }

        KRATOS_INFO_IF("AdjointStressLimitStateResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with max stress utilization = " << elem_id_at_max << std::endl;
        KRATOS_INFO_IF("AdjointStressLimitStateResponseFunction::CalculateValue", mEchoLevel > 0) << "Max stress utilization = " << max_stress_utilization << std::endl;

        // Set traced element in adjoint model part corresponding to the found element in the primal part
        mpTracedElementInAdjointPart = mrAdjointModelPart.pGetElement(elem_id_at_max);
        mpTracedElementInAdjointPart->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        return max_stress_utilization;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointStressLimitStateResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                             const Matrix& rResidualGradient,
                                                             Vector& rResponseGradient,
                                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointStressLimitStateResponseFunction::CalculateGradient: Traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            Matrix stress_displacement_derivative;

            mpTracedElementInAdjointPart->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                << "AdjointStressLimitStateResponseFunction::CalculateGradient: Size of stress displacement derivative does not fit!" << std::endl;

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
    void AdjointStressLimitStateResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<double>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointStressLimitStateResponseFunction::CalculateGradient: traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            if(rVariable == YIELD_STRESS) {
                rSensitivityGradient.resize(1, false);
                rSensitivityGradient[0] = 1.0;
            } 
            else {
                this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, rProcessInfo);
            }
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("")
    }

    // --------------------------------------------------------------------------
    void AdjointStressLimitStateResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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
    void AdjointStressLimitStateResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<array_1d<double, 3>>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mpTracedElementInAdjointPart == nullptr)
            << "AdjointStressLimitStateResponseFunction::CalculateGradient: traced element not initialized. First call \"CalculateValue\"!" << std::endl;

        if(rAdjointElement.Id() == mpTracedElementInAdjointPart->Id())
        {
            if(rVariable == YIELD_STRESS) {
                rSensitivityGradient.resize(1, false);
                rSensitivityGradient[0] = 1.0;
            } 
            else {
                this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, rProcessInfo);
            }
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointStressLimitStateResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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
    void AdjointStressLimitStateResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
                                      const std::string& rVariableName,
                                      const Matrix& rSensitivityMatrix,
                                      Vector& rSensitivityGradient,
                                      const ProcessInfo& rProcessInfo)
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
    void AdjointStressLimitStateResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
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

            stress_derivative_value *= -1.0; // Multiply by "-1" because: LS = R - S(p) -> dLS/dp = -dS/dp and LS = R - S(u) -> dLS/du = -dS/du

            rResponseGradient[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }

} // namespace Kratos.

