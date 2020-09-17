//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_CONVERGENCE_ACCELERATOR )
#define  KRATOS_CONVERGENCE_ACCELERATOR

// System includes
#include <set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "input_output/logger.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Base class for convergence accelerators
 * This class is intended to be the base of any convergence accelerator in Kratos
 * @tparam TSpace Linear algebra space
 */
template<class TSpace>
class ConvergenceAccelerator
{

public:

    ///@name Type Definitions
    ///@{

    typedef typename TSpace::VectorType                             VectorType;
    typedef typename TSpace::MatrixType                             MatrixType;

    typedef typename TSpace::VectorPointerType               VectorPointerType;
    typedef typename TSpace::MatrixPointerType               MatrixPointerType;

    // Counted pointer of ConvergenceAccelerator
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceAccelerator);

    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor
    ConvergenceAccelerator() = default;

    // Deleted copy constructor
    ConvergenceAccelerator(const ConvergenceAccelerator& Other) = delete;

    // Default destructor
    virtual ~ConvergenceAccelerator() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the convergence accelerator
     * Operations that are required to be done once before the resolution of the problem
     */
    virtual void Initialize()
    {
    }

    /**
     * @brief Operations before solving the solution step
     * Performs all the required operations that should be done (for each step) before solving the solution step.
     */
    virtual void InitializeSolutionStep()
    {
    }

    /**
     * @brief Operations before each non-linear iteration
     * Performs all the required operations that should be done before each non-linear iteration.
     */
    virtual void InitializeNonLinearIteration()
    {
    }

    /**
     * @brief Corrects the given iteration guess
     * Computes the correction of the iteration guess according to the provided residual vector
     * Note that the correction is performed over the given iteration guess, so the return value is already modified
     * @param rResidualVector Residual vector from which the correction is computed
     * @param rIterationGuess Iteration guess to be corrected
     */
    virtual void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess)
    {
    }

    /**
     * @brief Operations after each non-linear iteration
     * Performs all the required operations that should be done at the end of each non-linear iteration.
     */
    virtual void FinalizeNonLinearIteration()
    {
    }

    /**
     * @brief Operations after solving the solution step
     * Performs all the required operations that should be done (for each step) after solving the solution step
     */
    virtual void FinalizeSolutionStep()
    {
    }

    /**
     * @brief Finalize the convergence accelerator
     * Perform all the operations required after the resolution of the problem
     */
    virtual void Finalize()
    {
    }

    /**
     * @brief Clear the internal storage
     * Clears all the internal data (e.g. previous observation matrices)
     */
    virtual void Clear()
    {
    }

    /**
     * @brief Set the Echo Level object
     * Set the echo level of the convergence accelerator
     * @param Level Echo level value
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief Get the Echo Level object
     * Get the echo level of the convergence accelerator
     * @return int Echo level value
     */
    int GetEchoLevel()
    {
        return mEchoLevel;
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    // Level of echo for the convergence accelerator
    int mEchoLevel;

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
}; // Class ConvergenceAccelerator

///@} // Kratos classes

///@} // Application group

} /// namespace Kratos

#endif // KRATOS_CONVERGENCE_ACCELERATOR  defined
