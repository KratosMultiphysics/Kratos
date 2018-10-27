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

#if !defined(KRATOS_AND_CRITERIA_H)
#define  KRATOS_AND_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "factories/convergence_criteria_factory.h"

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
 * @class And_Criteria
 * @ingroup KratosCore
 * @brief This convergence criteria checks simultaneously two convergence criteria (both must be satisfied)
 * @details It takes two different convergence criteria in order to work
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class And_Criteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of And_Criteria */
    KRATOS_CLASS_POINTER_DEFINITION(And_Criteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;

    typedef ConvergenceCriteriaFactory< TSparseSpace, TDenseSpace > ConvergenceCriteriaFactoryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @details It takes two different convergence criteria in order to work
     */
    explicit And_Criteria(Kratos::Parameters Settings)
        :BaseType()
    {
        // We check if the criterias are defined
        if (Settings.Has("first_criterion_settings") && Settings.Has("second_criterion_settings")) {
            mpFirstCriterion = ConvergenceCriteriaFactoryType().Create(Settings["first_criterion_settings"]);
            mpSecondCriterion = ConvergenceCriteriaFactoryType().Create(Settings["second_criterion_settings"]);
        } else { // Displacement criteria and residual criteria will be combined
            Settings["convergence_criterion"].SetString("displacement_criterion");
            mpFirstCriterion = ConvergenceCriteriaFactoryType().Create(Settings);
            Settings["convergence_criterion"].SetString("residual_criterion");
            mpSecondCriterion = ConvergenceCriteriaFactoryType().Create(Settings);
        }
    }

    /**
     * @brief Default constructor.
     * @details It takes two different convergence criteria in order to work
     * @param pFirstCriterion The first convergence criteria
     * @param pSecondCriterion The second convergence criteria
     */
    explicit And_Criteria(
        ConvergenceCriteriaPointerType pFirstCriterion,
        ConvergenceCriteriaPointerType pSecondCriterion
        ) :BaseType(),
           mpFirstCriterion(pFirstCriterion),
           mpSecondCriterion(pSecondCriterion)
    {
    }

    /**
     * @brief Copy constructor.
     * @param rOther The other And_Criteria to be copied
     */
    explicit And_Criteria(And_Criteria const& rOther)
        :BaseType(rOther),
         mpFirstCriterion(rOther.mpFirstCriterion),
         mpSecondCriterion(rOther.mpSecondCriterion)
     {
     }

    /** Destructor.
    */
    ~And_Criteria () override {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing extra informations
     */
    void SetEchoLevel(int Level) override
    {
        BaseType::SetEchoLevel(Level);
        mpFirstCriterion->SetEchoLevel(Level);
        mpSecondCriterion->SetEchoLevel(Level);
    }

    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        const bool first_criterion_result  = mpFirstCriterion ->PreCriteria(rModelPart,rDofSet,A,Dx,b);
        const bool second_criterion_result = mpSecondCriterion ->PreCriteria(rModelPart,rDofSet,A,Dx,b);

        return (first_criterion_result && second_criterion_result);
    }

    /**
     * @brief Criteria that need to be called after getting the solution
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        const bool first_criterion_result  = mpFirstCriterion ->PostCriteria(rModelPart,rDofSet,A,Dx,b);
        const bool second_criterion_result = mpSecondCriterion ->PostCriteria(rModelPart,rDofSet,A,Dx,b);

        return (first_criterion_result && second_criterion_result);
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        mpFirstCriterion->Initialize(rModelPart);
        mpSecondCriterion->Initialize(rModelPart);
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        mpFirstCriterion->InitializeSolutionStep(rModelPart,rDofSet,A,Dx,b);
        mpSecondCriterion->InitializeSolutionStep(rModelPart,rDofSet,A,Dx,b);
    }

    /**
     * @brief This function initializes the non linear iteration
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void InitializeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        mpFirstCriterion->InitializeNonLinearIteration(rModelPart,rDofSet,A,Dx,b);
        mpSecondCriterion->InitializeNonLinearIteration(rModelPart,rDofSet,A,Dx,b);
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        mpFirstCriterion->FinalizeSolutionStep(rModelPart,rDofSet,A,Dx,b);
        mpSecondCriterion->FinalizeSolutionStep(rModelPart,rDofSet,A,Dx,b);
    }

    /**
     * @brief This function finalizes the non linear iteration
     * @param rModelPart ModelPart containing the problem.
     * @param rDofSet Container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        mpFirstCriterion->FinalizeNonLinearIteration(rModelPart,rDofSet,A,Dx,b);
        mpSecondCriterion->FinalizeNonLinearIteration(rModelPart,rDofSet,A,Dx,b);
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart ModelPart containing the problem.
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const int check1 = mpFirstCriterion->Check(rModelPart);
        const int check2 = mpSecondCriterion->Check(rModelPart);

        return check1 + check2;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "And_Criteria";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

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

    ConvergenceCriteriaPointerType mpFirstCriterion;  /// The pointer to the first convergence criterion
    ConvergenceCriteriaPointerType mpSecondCriterion; /// The pointer to the second convergence criterion

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

}; /* Class And_Criteria */

///@}

///@name Type Definitions */
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_AND_CRITERIA_H  defined */

