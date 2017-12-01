//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_RESPONSE_FUNCTION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A base class for scalar response functions.
class ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~ResponseFunction() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
    }

    virtual void FinalizeSolutionStep(ModelPart& rModelPart)
    {
    }

    virtual void Check(ModelPart const& rModelPart)
    {
    }

    virtual void Clear()
    {
    }

    /// Calculate the response function value.
    virtual double CalculateValue(ModelPart const& rModelPart)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        return 0.0;
        KRATOS_CATCH("");
    }

    // Update the sensitivities of the response function.
    virtual void UpdateSensitivities(ModelPart& rModelPart) = 0;

    /// Calculate the local gradient w.r.t. the primal variable.
    /**
     * @param[in]     rElement          local adjoint element.
     * @param[in]     rAdjointMatrix    transposed gradient of the residual
     *                                  w.r.t. primal variable.
     * @param[out]    rResponseGradient gradient of the response function w.r.t.
     *                                  primal variable.
     * @param[in]     rProcessInfo      current process info.
     */
    virtual void CalculateGradient(Element const& rElement,
                                   Matrix const& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculateGradient(Condition const& rCondition,
                                   Matrix const& rAdjointMatrix,
                                   Vector& rResponseGradient,
                                   ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. the first derivative of the primal variable.
    /**
     * @param[in]     rElement          local adjoint element.
     * @param[in]     rAdjointMatrix    transposed gradient of the residual
     *                                  w.r.t. first derivative of the primal variable.
     * @param[out]    rResponseGradient gradient of the response function w.r.t.
     *                                  first derivative of the primal variable.
     * @param[in]     rProcessInfo      current process info.
     */
    virtual void CalculateFirstDerivativesGradient(Element const& rElement,
                                                   Matrix const& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculateFirstDerivativesGradient(Condition const& rCondition,
                                                   Matrix const& rAdjointMatrix,
                                                   Vector& rResponseGradient,
                                                   ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    /// Calculate the local gradient w.r.t. the second derivative of the primal variable.
    /**
     * @param[in]     rElement          local adjoint element.
     * @param[in]     rAdjointMatrix    transposed gradient of the residual
     *                                  w.r.t. second derivative of the primal variable.
     * @param[out]    rResponseGradient local gradient of the response function
     *                                  w.r.t. the second derivative of the primal variable.
     * @param[in]     rProcessInfo      current process info.
     */
    virtual void CalculateSecondDerivativesGradient(Element const& rElement,
                                                    Matrix const& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    virtual void CalculateSecondDerivativesGradient(Condition const& rCondition,
                                                    Matrix const& rAdjointMatrix,
                                                    Vector& rResponseGradient,
                                                    ProcessInfo const& rProcessInfo) const
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the base class method." << std::endl;

        KRATOS_CATCH("");
    }

    ///@}
}; // class ResponseFunction

///@} // Kratos Classes
} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION_H_INCLUDED defined */
