// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Mohamed Khalil
//                   Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY

/* System Includes */

/* External Includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/contact_utilities.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// TODO: Extend the descriptions

namespace Kratos {

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
    
/** \brief  Short class definition.
This class 
*/

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
         
class ResidualBasedNewtonRaphsonContactStrategy :
    public ResidualBasedNewtonRaphsonStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonContactStrategy );

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    
    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     * @param SplitFactor: Value by one we split the time step when the problem does not converge
     * @param MaxNumberSplits: Maximim number of splits of the time step
     */
    
    ResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        double SplitFactor = 10.0,
        unsigned int MaxNumberSplits = 3
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        // Splitting values
        mSplitFactor = SplitFactor;
        mMaxNumberSplits = MaxNumberSplits;

        KRATOS_CATCH("");
    }

    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     * @param SplitFactor: Value by one we split the time step when the problem does not converge
     * @param MaxNumberSplits: Maximim number of splits of the time step
     */
    
    ResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        double SplitFactor = 10,
        unsigned int MaxNumberSplits = 3
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        // Splitting values
        mSplitFactor = SplitFactor;
        mMaxNumberSplits = MaxNumberSplits;

        KRATOS_CATCH("");
    }
    
    /******************** OPERATIONS ACCESSIBLE FROM THE INPUT: ************************/
    /***********************************************************************************/
 
    /**
     * This is the initialization of the iteration cycle
     * @return is_converged: True when the problem has converged
     * @return ResidualIsUpdated: True when the residual has been updated
     * @return iteration_number: Number of iterations done 
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     * @param mb: The RHS of the problem
     */
    
    void InitiliazeCycle(
        bool& is_converged,
        bool & ResidualIsUpdated,
        unsigned int& iteration_number,
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb
        )
    {
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
      
        CoutSolvingProblem();
        
//         BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);

        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        else
        {
            TSparseSpace::SetToZero(mDx); // mDx=0.00;
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }

        if (this->GetEchoLevel() == 3) // If it is needed to print the debug info
        {
            std::cout << "SystemMatrix = " << mA << std::endl;
            std::cout << "solution obtained = " << mDx << std::endl;
            std::cout << "RHS  = " << mb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) // Print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), mA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), mb);
        }
        
        // Updating the results stored in the database
        UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb, BaseType::MoveMeshFlag());

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

        if (is_converged == true)
        {
            // Initialisation of the convergence criteria
            rDofSet = pBuilderAndSolver->GetDofSet();
            BaseType::mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);

                pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }
    }
 
    /**
     * This is the inner iteration cycle considered in the SolveSolutionStep
     * @return is_converged: True when the problem has converged
     * @return ResidualIsUpdated: True when the residual has been updated
     * @return iteration_number: Number of iterations done 
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     * @param mb: The RHS of the problem
     */
    
    void IterationCycle(
        bool& is_converged,
        bool & ResidualIsUpdated,
        unsigned int& iteration_number,
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb
        )
    {
        // Iteration Cicle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < BaseType::mMaxIterationNumber)
        {
            // Setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number; 
            
            CoutSolvingProblem();
        
            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            // Call the linear system solver to find the correction mDx for the
            // It is not called if there is no system to solve
            if (SparseSpaceType::Size(mDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
                {
                    if( BaseType::GetKeepSystemConstantDuringIterations() == false)
                    {
                        //mA = 0.00;
                        TSparseSpace::SetToZero(mA);
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(mDx);
                        TSparseSpace::SetToZero(mb);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(mDx);
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                }
            }
            else
            {
                std::cout << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Updating the results stored in the database
            UpdateDatabase(pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb, BaseType::MoveMeshFlag());

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(mb);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
                    ResidualIsUpdated = true;
                    //std::cout << "mb is calculated" << std::endl;
                }
                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
            }
        }
    }
    
        
    /**
     * Here the database is updated
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     */
    void UpdateDatabase(
        typename TSchemeType::Pointer& pScheme,
        typename TBuilderAndSolverType::Pointer& pBuilderAndSolver,
        DofsArrayType& rDofSet,
        TSystemMatrixType& mA,
        TSystemVectorType& mDx,
        TSystemVectorType& mb,
        const bool MoveMesh
    )
    {
        rDofSet = pBuilderAndSolver->GetDofSet();
        pScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        // Move the mesh if needed
        if (MoveMesh == true) 
        {
            BaseType::MoveMesh();
        }
    }
    
    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    
    bool SolveSolutionStep()
    {
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();

        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        TSystemMatrixType& mA  = *BaseType::mpA;
        TSystemVectorType& mDx = *BaseType::mpDx;
        TSystemVectorType& mb  = *BaseType::mpb;

        //  Initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        unsigned int split_number     = 1;
        bool is_converged         = false;
        bool ResidualIsUpdated    = false;
        
        InitiliazeCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);
       
        IterationCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);

        if (iteration_number >= BaseType::mMaxIterationNumber)
        {
//             double original_delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]; // We save the delta time to restore later
            
//             // We iterate until we reach the convergence or we split more than desired
//             while (is_converged == false && split_number <= mMaxNumberSplits)
//             {   
//                 // Expliting time step as a way to try improve the convergence
//                 split_number += 1;
//                 iteration_number = 1;
//                 
//                 double aux_time       = BaseType::GetModelPart().GetProcessInfo()[TIME];
//                 double aux_delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]; // FIXME: The DELTA_TIME is set to 0 for some reason!!!!
//                 double current_time   = aux_time - aux_delta_time;
//                 
//                 BaseType::GetModelPart().GetProcessInfo()[TIME]       =   current_time; // Restore time to the previous one
//                 aux_delta_time /= mSplitFactor;
//                 BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME] = aux_delta_time; // Change delta time
//                 
//                 CoutSplittingTime(aux_delta_time);
//                 
//                 unsigned int aux_cout = 0;
//                 while (is_converged == false && BaseType::GetModelPart().GetProcessInfo()[TIME] <= aux_time && iteration_number < mMaxIterationNumber)
//                 {      
//                     iteration_number = 1;
//                     current_time += aux_delta_time;
//          
//                     aux_cout += 1;
//                     if (aux_cout > 1) // We avoid to restore the database if we create the new step just after the first iteration
//                     {
//                         BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] += 1;
//                     }
// 
//                     BaseType::GetModelPart().GetProcessInfo()[TIME] = current_time; // Increase the time in the new delta time        
//                     BaseType::GetModelPart().CloneTimeStep(current_time);
//                     
//                     // We repeat the predict with the new DELTA_TIME
//                     Predict();
//                     
//                     InitiliazeCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);
//                     IterationCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb); 
//                 
//                     // Plots a warning if the maximum number of iterations is exceeded
//                     if (is_converged == false  && iteration_number >= mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
//                     {
//                         MaxIterationsExceeded();
//                     }
//                 }
//                 
//                 if (is_converged == true)
//                 {
//                     // Restoring original DELTA_TIME
//                     BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME] = original_delta_time;
//                 }
//             }
            
            // Plots a warning if the maximum number of iterations and splits are exceeded
            if (is_converged == false  && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                MaxIterationsAndSplitsExceeded();
            }
            
        }

        // Recalculate residual if needed (note that some convergence criteria need it to be recalculated)
        if (ResidualIsUpdated == false)
        {
            // NOTE: The following part will be commented because it is time consuming
      
            //    TSparseSpace::SetToZero(mb);
            //    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
        }

        // Calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        
        return is_converged;
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

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    /**
    * Number by one the delta time is split
    */
    double mSplitFactor;
    
    /**
    * Maximum number of splits
    */
    unsigned int mMaxNumberSplits;

    ///@}
    ///@name Protected Operators
    ///@{

    void CoutSolvingProblem()
    {
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "STEP: " << BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] << "\t NON LINEAR ITERATION: " << BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] << "\t TIME: " << BaseType::GetModelPart().GetProcessInfo()[TIME] << "\t DELTA TIME: " << BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]  << std::endl;
        }
    }
    
    void CoutSplittingTime(const double aux_delta_time)
    {
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "**** Max. iter. exceeded: SPLITTING TIME STEP *****" << std::endl;
            std::cout << "***\t\t COMING BACK TO TIME: " << BaseType::GetModelPart().GetProcessInfo()[TIME] << "\t\t ***" << std::endl;
            std::cout << "***\t\t NEW TIME STEP: "<< aux_delta_time << "\t\t ***" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }
    
    void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }
    
    void MaxIterationsAndSplitsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "********** Max number of splits exceeded **********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }
    
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
    ///@{

    /** 
     * Copy constructor.
     */
    
    ResidualBasedNewtonRaphsonContactStrategy(const ResidualBasedNewtonRaphsonContactStrategy& Other)
    {
    };

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; /* Class ResidualBasedNewtonRaphsonContactStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY */
