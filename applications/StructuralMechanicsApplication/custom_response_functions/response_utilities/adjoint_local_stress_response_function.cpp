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
#include "adjoint_local_stress_response_function.h"
#include "local_stress_response_function.h"

namespace Kratos
{
    AdjointLocalStressResponseFunction::AdjointLocalStressResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        // Get traced element
        const int id_of_traced_element = ResponseSettings["traced_element_id"].GetInt();
        mpTracedElement = rModelPart.pGetElement(id_of_traced_element);

        // Tell traced element the stress type
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        mpTracedElement->SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        // Get info how and where to treat the stress
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );

        if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node)
        {
            mIdOfLocation = ResponseSettings["stress_location"].GetInt();
            KRATOS_ERROR_IF(mIdOfLocation < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
        }
    }

    AdjointLocalStressResponseFunction::~AdjointLocalStressResponseFunction(){}


    double AdjointLocalStressResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        LocalStressResponseFunction zero_order_response(rModelPart, mResponseSettings);

        return zero_order_response.CalculateValue();

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateGradient(const Element& rAdjointElem, const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            Matrix stress_displacement_derivative;

            if(mStressTreatment == StressTreatment::Mean)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
                this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);
            }
            else if(mStressTreatment == StressTreatment::GaussPoint)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
                this->ExtractGaussPointStressDerivative(stress_displacement_derivative, rResponseGradient);
            }
            else if(mStressTreatment == StressTreatment::Node)
            {
                mpTracedElement->Calculate(STRESS_DISP_DERIV_ON_NODE, stress_displacement_derivative, rProcessInfo);
                this->ExtractNodeStressDerivative(stress_displacement_derivative, rResponseGradient);
            }

            KRATOS_ERROR_IF(rResponseGradient.size() != rAdjointMatrix.size1())
                 << "Size of stress displacement derivative does not fit!" << std::endl;

            rResponseGradient *= (-1);
        }
        else
        {
            if(rResponseGradient.size() != rAdjointMatrix.size1())
                rResponseGradient.resize(rAdjointMatrix.size1(), false);

            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<double>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            this->CalculateElementContributionToSensitivityGradient(rAdjointElem, rVariable.Name(), rDerivativesMatrix,
                                                                    rResponseGradient, rProcessInfo);
        }
        else
        {
            if (rResponseGradient.size() != 0)
                rResponseGradient.resize(0, false);
        }

        KRATOS_CATCH("")
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rDerivativesMatrix,
                                     Vector& rResponseGradient,
                                     ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != 0)
            rResponseGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if(rAdjointElem.Id() == mpTracedElement->Id())
        {
            this->CalculateElementContributionToSensitivityGradient(rAdjointElem, rVariable.Name(), rDerivativesMatrix,
                                                                    rResponseGradient, rProcessInfo);
        }
        else
        {
            if (rResponseGradient.size() != 0)
                rResponseGradient.resize(0, false);
        }

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != 0)
            rResponseGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::CalculateElementContributionToSensitivityGradient(Element& rAdjointElem,
                                      const std::string& rVariableName,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rResponseGradient,
                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        if(mStressTreatment == StressTreatment::Mean)
        {
            rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_design_variable_derivative, rResponseGradient);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_design_variable_derivative, rResponseGradient);
        }
        else if(mStressTreatment == StressTreatment::Node)
        {
            rAdjointElem.Calculate(STRESS_DESIGN_DERIVATIVE_ON_NODE, stress_design_variable_derivative, rProcessInfo);
            this->ExtractNodeStressDerivative(stress_design_variable_derivative, rResponseGradient);
        }

        KRATOS_ERROR_IF(rResponseGradient.size() != rDerivativesMatrix.size1())
             << "Size of partial stress design variable derivative does not fit!" << std::endl;

        rAdjointElem.SetValue(DESIGN_VARIABLE_NAME, "");

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
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

    void AdjointLocalStressResponseFunction::ExtractNodeStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= mIdOfLocation ) << "Chosen node is not available. The element has only " <<
                            num_of_stress_positions << " nodes."<< std::endl;

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rResponseGradient[deriv_it] = rStressDerivativesMatrix(deriv_it, (mIdOfLocation-1));

        KRATOS_CATCH("");
    }

    void AdjointLocalStressResponseFunction::ExtractGaussPointStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= mIdOfLocation ) <<
                "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                num_of_stress_positions  << "!"<< std::endl;

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rResponseGradient[deriv_it] = rStressDerivativesMatrix(deriv_it, (mIdOfLocation-1));

        KRATOS_CATCH("");
    }

} // namespace Kratos.

