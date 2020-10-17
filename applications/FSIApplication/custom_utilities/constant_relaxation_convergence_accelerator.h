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
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"


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
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class ConstantRelaxationConvergenceAccelerator: public ConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ConstantRelaxationConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSparseSpace, TDenseSpace>              BaseType;

    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::VectorType                                VectorType;
    typedef typename BaseType::VectorPointerType                  VectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Get the Default Settings object
     * This method returns the default parameters for this convergence accelerator.
     * Note that it is required to be static since it is called during
     * the construction of the object so no instantation exists yet.
     * @return Parameters Default parameters json string
     */
    static Parameters GetDefaultParameters()
    {
        Parameters default_settings(R"(
        {
            "solver_type": "constant_relaxation",
            "w": 0.5
        }
        )");

        return default_settings;
    }

    /**
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Constructor with json string settings
     * @param rConvAcceleratorParameters json string encapsulating the settings
     */
    ConstantRelaxationConvergenceAccelerator(Parameters rConvAcceleratorParameters)
    : BaseType(),
      mOmega([] (Parameters x) -> double {x.ValidateAndAssignDefaults(GetDefaultParameters()); return x["w"].GetDouble();} (rConvAcceleratorParameters))
    {
    }

    /**
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Constructor with given relaxation factor
     * @param rOmega relaxation factor
     */
    ConstantRelaxationConvergenceAccelerator(const double rOmega = 0.5)
    : BaseType(),
      mOmega(rOmega)
    {
    }

    /**
     * @brief Construct a new Constant Relaxation Convergence Accelerator object
     * Explicitly deleted constant relaxation copy constructor
     * @param rOther constant relaxation convergence accelerator to be copied
     */
    ConstantRelaxationConvergenceAccelerator(const ConstantRelaxationConvergenceAccelerator &rOther) = delete;

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
    void UpdateSolution(const VectorType &rResidualVector,
                        VectorType &rIterationGuess) override
    {
        KRATOS_TRY;

        TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega, rResidualVector);

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
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const double mOmega;

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
