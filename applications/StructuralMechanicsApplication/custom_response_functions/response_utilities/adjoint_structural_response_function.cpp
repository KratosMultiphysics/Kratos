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
#include "utilities/openmp_utils.h"
#include "adjoint_structural_response_function.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

    /// Constructor.
    AdjointStructuralResponseFunction::AdjointStructuralResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
      : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        // Set gradient mode
        const std::string& gradient_mode = ResponseSettings["gradient_mode"].GetString();

        // Mode 1: semi-analytic sensitivities
        if (gradient_mode == "semi_analytic")
        {
            mGradientMode = 1;
            double delta = ResponseSettings["step_size"].GetDouble();
            mDelta = delta;
        }
        else
            KRATOS_ERROR << "Specified gradient_mode not recognized. The only option is: semi_analytic. Specified gradient_mode: " <<  gradient_mode << std::endl;


        KRATOS_CATCH("");
    }

    /// Destructor.
    AdjointStructuralResponseFunction::~AdjointStructuralResponseFunction()
    {
    }

    ModelPart& AdjointStructuralResponseFunction::GetModelPart()
    {
      return mrModelPart;
    }

    ModelPart& AdjointStructuralResponseFunction::GetModelPart() const
    {
      return mrModelPart;
    }

    void AdjointStructuralResponseFunction::Initialize()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        if(mGradientMode == 1)
        {
            VariableUtils().SetNonHistoricalVariable(PERTURBATION_SIZE, mDelta, r_model_part.Elements());
            VariableUtils().SetNonHistoricalVariable(PERTURBATION_SIZE, mDelta, r_model_part.Conditions());
        }

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. primal.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateGradient(const Element& rAdjointElem,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateFirstDerivativesGradient(const Element& rAdjointElem,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointElem      the adjoint element.
     * @param[in]     rAdjointMatrix    the transposed gradient of the
     *                                  element's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient the gradient of the response function.
     * @param[in]     rProcessInfo      the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSecondDerivativesGradient(const Element& rAdjointElem,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                                    const Matrix& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rAdjointMatrix.size1())
            rResponseGradient.resize(rAdjointMatrix.size1(), false);

        rResponseGradient.clear();

        KRATOS_CATCH("");
    }


    /// Calculate the scalar valued response function
    double AdjointStructuralResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        return 0.0;
    }


    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
         KRATOS_TRY;

         KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

         KRATOS_CATCH("");
    }

    /// Calculate the local gradient of response function w.r.t. the sensitivity variable.
    /**
     * @param[in]     rAdjointElem       the adjoint element.
     * @param[in]     rVariable          the sensitivity variable.
     * @param[in]     rDerivativesMatrix the transposed gradient of the element's
     *                                   residual w.r.t. the sensitivity variable.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointStructuralResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "This should be implemented in the derived class." << std::endl;

        KRATOS_CATCH("");
    }
};


