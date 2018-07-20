//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#ifndef ADJOINT_RESPONSE_FUNCTION_H
#define ADJOINT_RESPONSE_FUNCTION_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "containers/variable.h"


// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/// A base class for adjoint response functions.
class AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointResponseFunction()
    {
    }

    /// Destructor.
    virtual ~AdjointResponseFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize()
    {
    }

    virtual void InitializeSolutionStep()
    {
    }

    virtual void FinalizeSolutionStep()
    {
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. primal.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the local gradient w.r.t. primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. primal.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                                   const Matrix& rResidualGradient,
                                                   Vector& rResponseGradient,
                                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the local gradient w.r.t. first derivatives of primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. first derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                                   const Matrix& rResidualGradient,
                                                   Vector& rResponseGradient,
                                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointElement    the adjoint element.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   element's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                                    const Matrix& rResidualGradient,
                                                    Vector& rResponseGradient,
                                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the local gradient w.r.t. second derivatives of primal solution.
    /**
     * @param[in]     rAdjointCondition  the adjoint condition.
     * @param[in]     rResidualGradient  the transposed gradient of the
     *                                   condition's residual w.r.t. second derivatives.
     * @param[out]    rResponseGradient  the gradient of the response function.
     * @param[in]     rProcessInfo       the current process info.
     */
    virtual void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                                    const Matrix& rResidualGradient,
                                                    Vector& rResponseGradient,
                                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointElement       the adjoint element.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      element's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointCondition     the adjoint condition.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      condition's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointElement       the adjoint element.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      element's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the partial sensitivity w.r.t. design variable.
    /**
     * @param[in]     rAdjointCondition     the adjoint condition.
     * @param[in]     rVariable             the design variable.
     * @param[in]     rSensitivityMatrix    the transposed gradient of the
     *                                      condition's residual w.r.t. design variable.
     * @param[out]    rSensitivityGradient  the gradient of the response function.
     * @param[in]     rProcessInfo          the current process info.
     */
    virtual void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Calling base class response function method.\n";
    }

    /// Calculate the scalar valued response function.
    virtual double CalculateValue() = 0;

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* ADJOINT_RESPONSE_FUNCTION_H defined */
