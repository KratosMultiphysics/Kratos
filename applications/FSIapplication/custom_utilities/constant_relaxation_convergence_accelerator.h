//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
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

/**
 * @brief Constant relaxation convergence accelerator
 * This utility corrects the iteration guess with a constant relaxation factor
 * @tparam TSpace Linear algebra space type
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
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Constructor with json string settings
     * @param rConvAcceleratorParameters json string encapsulating the settings
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

    /**
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Constructor with given relaxation factor
     * @param rOmega relaxation factor
     */
    ConstantRelaxationConvergenceAccelerator(double rOmega = 0.5)
    {
        mOmega = rOmega;
    }

    /**
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Explicitly deleted constant relaxation copy constructor
     * @param rOther constant relaxation convergence accelerator to be copied
     */
    ConstantRelaxationConvergenceAccelerator(const ConstantRelaxationConvergenceAccelerator& rOther) = delete;

    /**
     * @brief Destroy the Constant Relaxation Convergence Accelerator object
     * Default constant relaxation constructor
     */
    virtual ~ConstantRelaxationConvergenceAccelerator() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Performs the solution update
     * The correction is performed as u_{k+1} = u_{k} + w * r_{k+1}, being u_{k+1} the current
     * iteration relaxed solution, u_{k} the previous iteration relaxed solution and r_{k+1}
     * the current iteration residual, which is computed as r_{k+1} = \hat{u}_{k+1} - u_{k},
     * being \hat{u}_{k+1} the current iteration unrelaxed solution.
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess
     */
    void UpdateSolution(const VectorType& rResidualVector,
                        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        TSpace::UnaliasedAdd(rIterationGuess, mOmega, rResidualVector);

        KRATOS_CATCH("");
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
///@name Input and output
///@{

///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_CONSTANT_RELAXATION_CONVERGENCE_ACCELERATOR defined */
