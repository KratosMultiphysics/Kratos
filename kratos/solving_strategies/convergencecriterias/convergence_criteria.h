//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_BASE_CONVERGENCE_CRITERIA_H )
#define  KRATOS_BASE_CONVERGENCE_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/model_part.h"

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
 * @class ConvergenceCriteria
 * @ingroup KratosCore
 * @brief This is the base class to define the  different convergence criterion considered
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Riccardo Rossi
*/
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ConvergenceCriteria
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteria);

    /// Data type definition
    typedef typename TSparseSpace::DataType TDataType;
    /// Matrix type definition
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    /// Vector type definition
    typedef typename TSparseSpace::VectorType TSystemVectorType;
    /// Local system matrix type definition
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// DoF array type definition
    typedef ModelPart::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    explicit ConvergenceCriteria()
    {
        mActualizeRHSIsNeeded = false;
        mConvergenceCriteriaIsInitialized = false;
        SetEchoLevel(1);
    }

    /** Copy constructor.
     */
    explicit ConvergenceCriteria( ConvergenceCriteria const& rOther)
      :mActualizeRHSIsNeeded(rOther.mActualizeRHSIsNeeded)
      ,mConvergenceCriteriaIsInitialized(rOther.mConvergenceCriteriaIsInitialized)
      ,mEchoLevel(rOther.mEchoLevel)
    {
    }

    /** Destructor.
     */
    virtual ~ConvergenceCriteria()
    {
    }

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get component wise element components
     * @warning Must be defined on the derived classes
     * @return The RHS element components
     */
    virtual std::vector<TSystemVectorType>&  GetRHS_Element_Components()
    {
        KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * @brief Get component wise element variables
     * @warning Must be defined on the derived classes
     * @return The RHS element variables
     */
    virtual std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Element_Variables()
    {
        KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * @brief Get component wise condition components
     * @warning Must be defined on the derived classes
     * @return The RHS condition components
     */
    virtual std::vector<TSystemVectorType>&  GetRHS_Condition_Components()
    {
        KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * @brief Get component wise condition variables
     * @warning Must be defined on the derived classes
     * @return The RHS condition variables
     */
    virtual std::vector< Variable< LocalSystemVectorType > >&  GetRHS_Condition_Variables()
    {
        KRATOS_ERROR <<"Asking for Global Components to the CONVERGENCE CRITERION base class which is not component wise and not contains this member variable" << std::endl;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel()
    {
        return mEchoLevel;
    }

    /**
     * @brief This method sets the flag mActualizeRHSIsNeeded
     * @param ActualizeRHSIsNeeded The flag that tells if actualize RHS is needed
     */
    void SetActualizeRHSFlag(bool ActualizeRHSIsNeeded)
    {
        mActualizeRHSIsNeeded = ActualizeRHSIsNeeded;
    }

    /**
     * @brief This method gets the flag mActualizeRHSIsNeeded
     * @return mActualizeRHSIsNeeded: The flag that tells if actualize RHS is needed
     */
    bool GetActualizeRHSflag()
    {
        return mActualizeRHSIsNeeded;
    }

    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    virtual bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
        return true;
    }

    /**
     * @brief Criterias that need to be called after getting the solution
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     * @return true if convergence is achieved, false otherwise
     */
    virtual bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
        return true;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the problem. (unused)
     */
    virtual void Initialize(ModelPart& rModelPart)
    {
        mConvergenceCriteriaIsInitialized = true;
    }

    /**
     * @brief This function returns if the convergence criteria is initialized
     * @return mConvergenceCriteriaIsInitialized, true if initialized, false otherwise
     */
    virtual bool IsInitialized()
    {
        return mConvergenceCriteriaIsInitialized;
    }

    /**
     * @brief This function initializes the solution step
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief This function initializes the non-linear iteration
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief This function finalizes the solution step
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    virtual void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        )
    {
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @warning Must be defined on the derived classes
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
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

    bool mActualizeRHSIsNeeded = false;             /// This "flag" is set in order to know if it is necessary to actualize the RHS
    bool mConvergenceCriteriaIsInitialized = false; /// This "flag" is set in order to know if it is convergence criteria is initialized

    int mEchoLevel; /// The echo level

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ConvergenceCriteria */
} /* namespace Kratos.*/

#endif /* KRATOS_BASE_CONVERGENCE_CRITERIA_H  defined */

