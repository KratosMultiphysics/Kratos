//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_RESPONSE_FUNCTION)
#define KRATOS_RESPONSE_FUNCTION

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A base class for response functions.
class ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
      : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~ResponseFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ModelPart& GetModelPart()
    {
      return mrModelPart;
    }

    const ModelPart& GetModelPart() const
    {
      return mrModelPart;
    }

    virtual void Initialize()
    {
    }

    virtual void InitializeSolutionStep()
    {
    }

    virtual void FinalizeSolutionStep()
    {
    }

    virtual void Check()
    {
    }

    /// Calculate the local gradient w.r.t. velocity
    /**
     * @param[in]     rElem            the local adjoint element.
     * @param[in]     rAdjointMatrix   the transposed gradient of the local
     *                                 element's residual w.r.t. velocity.
     * @param[out]    rRHSContribution the gradient of the response function.
     * @param[in,out] rProcessInfo     the current process info.
     */
    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
                                                      const Matrix& rAdjointMatrix,
                                                      Vector& rRHSContribution,
                                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    /// Calculate the local gradient w.r.t. acceleration
    /**
     * @param[in]     rElem            the local adjoint element.
     * @param[in]     rAdjointMatrix   the transposed gradient of the local
     *                                 element's residual w.r.t. acceleration.
     * @param[out]    rRHSContribution the gradient of the response function.
     * @param[in,out] rProcessInfo     the current process info.
     */
    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
                                                          const Matrix& rAdjointMatrix,
                                                          Vector& rRHSContribution,
                                                          ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    /// Calculate the local gradient w.r.t. nodal coordinates
    /**
     * @param[in]     rElem              the local adjoint element.
     * @param[in]     rDerivativesMatrix the transposed gradient of the local
     *                                   element's residual w.r.t. nodal
     *                                   coordinates.
     * @param[out]    rRHSContribution   the gradient of the response function.
     * @param[in,out] rProcessInfo       the current process info.
     */
    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    /// Calculate the scalar valued response function
    virtual double Calculate(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0.0;

        KRATOS_CATCH("")
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;

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

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION defined */
