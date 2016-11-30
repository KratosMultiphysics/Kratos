//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined(KRATOS_AITKEN_CONVERGENCE_ACCELERATOR)
#define  KRATOS_AITKEN_CONVERGENCE_ACCELERATOR

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/math_utils.h"

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "convergence_accelerator.hpp"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

///@name  Enum's
///@{
///@}

///@name  Functions
///@{
///@}

///@name Kratos Classes
///@{

/** @brief Aitken acceleration scheme 
 */
template<class TSpace>
class AitkenConvergenceAccelerator: public ConvergenceAccelerator<TSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( AitkenConvergenceAccelerator );
    
    typedef ConvergenceAccelerator<TSpace>                                 BaseType;
    
    typedef typename BaseType::Pointer                              BaseTypePointer;
    
    typedef typename BaseType::VectorType                                VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * Aitken convergence accelerator
     */
    AitkenConvergenceAccelerator(double rOmegaInitial = 0.825)
    {
        mOmega_0 = rOmegaInitial;
    }

    /** 
     * Copy Constructor.
     */
    AitkenConvergenceAccelerator(const AitkenConvergenceAccelerator& rOther)
    {
        mOmega_0 = rOther.mOmega_0;
    }

    /**
     * Clone
     */
    virtual BaseTypePointer Clone()
    {
        return BaseTypePointer( new AitkenConvergenceAccelerator(*this) );
    }

    /** 
     * Destructor.
     */
    virtual ~AitkenConvergenceAccelerator
    () {}

    ///@}
    
    ///@name Operators
    ///@{
    ///@}
    
    ///@name Operations
    ///@{

    /**
     * Initialize the internal iteration counter
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration = 1;

        KRATOS_CATCH( "" );
    }

    /**
     * Performs the solution update 
     * The correction is done as u_i+1 = u_i + w_i+1*r_i+1 where w_i+1 is de relaxation parameter computed using the Aitken's formula.
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess
     * @param rCorrectedValue: Computed corrected value
     */
    void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        mResidualVector_1 = rResidualVector;

        if (mConvergenceAcceleratorIteration == 1)
        {
            rIterationGuess += mOmega_0*mResidualVector_1;
        }
        else
        {
            double aux1 = inner_prod(mResidualVector_0,(mResidualVector_1-mResidualVector_0));
            double aux2 = inner_prod((mResidualVector_1-mResidualVector_0),(mResidualVector_1-mResidualVector_0));
            
            mOmega_1 = -mOmega_0*(aux1/aux2);
            
            rIterationGuess += mOmega_1*mResidualVector_1;
            mOmega_0 = mOmega_1;
        }
                        
        KRATOS_CATCH( "" );
    }
    
    /**
     * Updates the Aitken iteration values for the next non-linear iteration
     */
    void FinalizeNonLinearIteration() override
    {
        KRATOS_TRY;

        mResidualVector_0 = mResidualVector_1;
        mConvergenceAcceleratorIteration += 1;

        KRATOS_CATCH( "" );
    }
    
    /**
     * Reset the convergence accelerator iterations counter
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration = 1;

        KRATOS_CATCH( "" );
    }

    ///@}
    
    ///@name Access
    ///@{
    ///@}
    
    ///@name Inquiry
    ///@{
    ///@}
    
    ///@name Input and output
    ///@{
    ///@}
    
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    
    ///@name Protected member Variables
    ///@{

    unsigned int mConvergenceAcceleratorIteration;

    double mOmega_0;
    double mOmega_1;
    
    VectorType mResidualVector_0;
    VectorType mResidualVector_1;
    
    ///@}
    
    ///@name Protected Operators
    ///@{
    ///@}
    
    ///@name Protected Operations
    ///@{
    ///@}
    
    ///@name Protected  Access
    ///@{
    ///@}
    
    ///@name Protected Inquiry
    ///@{
    ///@}
    
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    
    ///@name Member Variables
    ///@{
    ///@}
    
    ///@name Private Operators
    ///@{
    ///@}
    
    ///@name Private Operations
    ///@{
    ///@}
    
    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{
    ///@}

    ///@name Private Inquiry
    ///@{
    ///@}
    
    ///@name Un accessible methods
    ///@{
    ///@}
    
}; /* Class AitkenConvergenceAccelerator */

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_AITKEN_CONVERGENCE_ACCELERATOR defined */
