//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined(KRATOS_CONSTANT_RELAXATION_CONVERGENCE_ACCELERATOR)
#define  KRATOS_CONSTANT_RELAXATION_CONVERGENCE_ACCELERATOR

/* System includes */

/* External includes */

/* Project includes */
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
class ConstantRelaxationConvergenceAccelerator: public ConvergenceAccelerator<TSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ConstantRelaxationConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSpace>                                 BaseType;

    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::VectorType                                VectorType;
    typedef typename BaseType::VectorPointerType                  VectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * Aitken convergence accelerator
     */
    ConstantRelaxationConvergenceAccelerator(Parameters& rConvAcceleratorParameters)
    {
        Parameters default_parameters(R"(
        {
            "solver_type": "constant_relaxation",
            "w": 0.5
        }
        )");

        rConvAcceleratorParameters.ValidateAndAssignDefaults(default_parameters);

        mOmega = rConvAcceleratorParameters["w"].GetDouble();
    }


    ConstantRelaxationConvergenceAccelerator(double rOmega = 0.5)
    {
        mOmega = rOmega;
    }

    /**
     * Copy Constructor.
     */
    ConstantRelaxationConvergenceAccelerator(const ConstantRelaxationConvergenceAccelerator& rOther)
    {
        mOmega = rOther.mOmega;
    }

    /**
     * Clone
     */
    virtual BaseTypePointer Clone()
    {
        return BaseTypePointer( new ConstantRelaxationConvergenceAccelerator(*this) );
    }

    /**
     * Destructor.
     */
    virtual ~ConstantRelaxationConvergenceAccelerator
    () {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * Performs the solution update
     * The correction is done as u_i+1 = u_i + w_i+1*r_i+1 where w_i+1 is de relaxation parameter computed using the Aitken's formula.
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess
     */
    void UpdateSolution(const VectorType& rResidualVector,
                        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        TSpace::UnaliasedAdd(rIterationGuess, mOmega, rResidualVector);

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

    double mOmega;

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

}; /* Class ConstantRelaxationConvergenceAccelerator */

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_CONSTANT_RELAXATION_CONVERGENCE_ACCELERATOR defined */
