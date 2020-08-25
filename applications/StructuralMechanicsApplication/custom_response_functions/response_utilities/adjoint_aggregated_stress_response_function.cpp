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
#include "adjoint_aggregated_stress_response_function.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
    // --------------------------------------------------------------------------
    AdjointAggregatedStressResponseFunction::AdjointAggregatedStressResponseFunction(ModelPart& rAdjointModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rAdjointModelPart, ResponseSettings),
      mrAdjointModelPart(rAdjointModelPart),
      mCriticalPartName(ResponseSettings["critical_part_name"].GetString()),
      mRho(ResponseSettings["rho"].GetDouble())
    {
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment(ResponseSettings["stress_treatment"].GetString());

        if(ResponseSettings.Has("echo_level"))
            mEchoLevel = ResponseSettings["echo_level"].GetInt();

        KRATOS_ERROR_IF(mStressTreatment != StressTreatment::Mean)
            << "AdjointAggregatedStressResponseFunction::AdjointAggregatedStressResponseFunction: Specified stress treatment not supported: " << ResponseSettings["stress_type"].GetString() << std::endl;

        // Collect relevant elements
        ModelPart& adjoint_aggregation_part = mrAdjointModelPart.GetSubModelPart(mCriticalPartName);
        for(auto& elem : adjoint_aggregation_part.Elements())
            mAggregatedElementIds.push_back(elem.Id());

        // Initialize
        for(auto& elem : adjoint_aggregation_part.Elements())
        {
            elem.SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType));
            mKSPrefactors[elem.Id()] = 0.0;
        }
    }

    // --------------------------------------------------------------------------
    AdjointAggregatedStressResponseFunction::~AdjointAggregatedStressResponseFunction(){}

    // --------------------------------------------------------------------------
    double AdjointAggregatedStressResponseFunction::CalculateValue(ModelPart& rPrimalModelPart)
    {
        KRATOS_TRY;

        ModelPart& primal_aggregation_part = rPrimalModelPart.GetSubModelPart(mCriticalPartName);

        // Find max mean value and corresponding element in primal model part
        std::vector<double> mean_stresses;
        mean_stresses.reserve(primal_aggregation_part.NumberOfElements());

        double max_mean_stress = 0.0;
        IndexType elem_id_at_max = 0;

        // Determine max value
        for(auto& elem : primal_aggregation_part.Elements())
        {
            Vector element_stress;
            StressCalculation::CalculateStressOnGP(elem, mTracedStressType, element_stress, rPrimalModelPart.GetProcessInfo());

            const SizeType stress_vec_size = element_stress.size();
            double mean_stress = 0.0;
            for(IndexType i = 0; i < stress_vec_size; ++i)
                mean_stress += element_stress[i];

            mean_stress /= stress_vec_size;
            mean_stresses.push_back(mean_stress);

            if(mean_stress > max_mean_stress)
            {
                max_mean_stress = mean_stress;
                elem_id_at_max = elem.Id();
            }
        }

        KRATOS_INFO_IF("  AdjointAggregatedStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Id of element with max stress value = " << elem_id_at_max << std::endl;
        KRATOS_INFO_IF("  AdjointAggregatedStressResponseFunction::CalculateValue", mEchoLevel > 0) << "Max mean stress value = " << max_mean_stress << std::endl;

        if(std::abs(mStressScalingFactor) < 1e-13)
            mStressScalingFactor = max_mean_stress;

        // Aggegration according Kreiselmeier-Steinhauser Function
        int elem_index = 0;
        mSumPrefactors = 0.0;

        for(auto& elem : primal_aggregation_part.Elements())
        {
            double KS_value_contribution = std::exp(mRho * mean_stresses[elem_index] / mStressScalingFactor);
            mKSPrefactors[elem.Id()] = KS_value_contribution;
            mSumPrefactors += KS_value_contribution;
            elem_index++;
        }
        mArePrefactorsInitialized = true;

        double KS_value = 1/mRho * std::log(mSumPrefactors);
        return KS_value;

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointAggregatedStressResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                             const Matrix& rResidualGradient,
                                                             Vector& rResponseGradient,
                                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mArePrefactorsInitialized == false) << "AdjointAggregatedStressResponseFunction::CalculateGradient: Prefactors missing. First calculate value before calculating gradients!" << std::endl;

        rResponseGradient.clear();

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            Element::Pointer pElem = mrAdjointModelPart.pGetElement(rAdjointElement.Id());

            Matrix stress_displacement_derivative;
            pElem->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);

            KRATOS_ERROR_IF(rResponseGradient.size() != rResidualGradient.size1())
                << "Size of stress displacement derivative does not fit!" << std::endl;

            rResponseGradient *= (-1) * mKSPrefactors[rAdjointElement.Id()] / (mStressScalingFactor * mSumPrefactors);
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
    void AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<double>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(mArePrefactorsInitialized == false) << "AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity: Prefactors missing. First calculate value before calculating gradients!" << std::endl;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
            rSensitivityGradient *= mKSPrefactors[rAdjointElement.Id()] / (mStressScalingFactor * mSumPrefactors);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("")
    }

    // --------------------------------------------------------------------------
    void AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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
    void AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                                                       const Variable<array_1d<double, 3>>& rVariable,
                                                                       const Matrix& rSensitivityMatrix,
                                                                       Vector& rSensitivityGradient,
                                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(mArePrefactorsInitialized == false) << "AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity: Prefactors missing. First calculate value before calculating gradients!" << std::endl;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        {
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rAdjointElement, rVariable.Name(), rSensitivityMatrix,
                                                                    rSensitivityGradient, process_info);
            rSensitivityGradient *= mKSPrefactors[rAdjointElement.Id()] / (mStressScalingFactor * mSumPrefactors);
        }
        else
            rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AdjointAggregatedStressResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
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
    void AdjointAggregatedStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rAdjointElement,
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
    void AdjointAggregatedStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
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

    // --------------------------------------------------------------------------
    void AdjointAggregatedStressResponseFunction::FinalizeSolutionStep()
    {
        KRATOS_TRY;

        mArePrefactorsInitialized = false;

        KRATOS_CATCH("");
    }

} // namespace Kratos.

