// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
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
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            StrategyBaseType;
    
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef TSparseSpace                                                          SparseSpaceType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;

    typedef typename BaseType::DofsArrayType                                        DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;
    
    typedef ModelPart::NodesContainerType                                          NodesArrayType;
    
    typedef ModelPart::ConditionsContainerType                                ConditionsArrayType;
    
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
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag )
    {
        KRATOS_TRY;

        // Splitting values
        mSplitFactor = SplitFactor;
        mMaxNumberSplits = MaxNumberSplits;

        KRATOS_CATCH("");
    }

    /** 
     * Destructor.
     */
    
    virtual ~ResidualBasedNewtonRaphsonContactStrategy()
    {
    }
    
    /******************** OPERATIONS ACCESSIBLE FROM THE INPUT: ************************/
    /***********************************************************************************/
    
    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    
    bool SolveSolutionStep() override
    {
        bool is_converged = BaseType::SolveSolutionStep();
        
//         // Pointers needed in the solution
//         typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
//         typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
// 
//         TSystemMatrixType& A  = *BaseType::mpA;
//         TSystemVectorType& Dx = *BaseType::mpDx;
//         TSystemVectorType& b  = *BaseType::mpb;
// 
//         //initializing the parameters of the Newton-Raphson cicle
//         unsigned int iteration_number = 1;
//         BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
// //         BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
//         bool is_converged = false;
//         bool ResidualIsUpdated = false;
//         pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
//         is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
// 
//         //function to perform the building and the solving phase.
//         if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
//         {
//             TSparseSpace::SetToZero(A);
//             TSparseSpace::SetToZero(Dx);
//             TSparseSpace::SetToZero(b);
// 
//             pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
//         }
//         else
//         {
//             TSparseSpace::SetToZero(Dx); //Dx=0.00;
//             TSparseSpace::SetToZero(b);
// 
//             pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
//         }
//         
//         // Debugging info
//         BaseType::EchoInfo(iteration_number);
//         
//         // Updating the results stored in the database
//         BaseType::UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());
// 
//         // Calculate reactions
//         CalculateContactReactions(pBuilderAndSolver,pScheme, b);
//         
//         pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
// 
//         if (is_converged == true)
//         {
//             //initialisation of the convergence criteria
//             BaseType::mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
// 
//             if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
//             {
//                 TSparseSpace::SetToZero(b);
// 
//                 pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
//             }
// 
//             is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
//         }
// 
//         //Iteration Cicle... performed only for NonLinearProblems
//         while (is_converged == false && iteration_number++ < BaseType::mMaxIterationNumber)
//         {
//             // Setting the number of iteration
//             BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
// 
//             pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
// 
//             is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
// 
//             //call the linear system solver to find the correction mDx for the
//             //it is not called if there is no system to solve
//             if (SparseSpaceType::Size(Dx) != 0)
//             {
//                 if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
//                 {
//                     if( BaseType::GetKeepSystemConstantDuringIterations() == false)
//                     {
//                         //A = 0.00;
//                         TSparseSpace::SetToZero(A);
//                         TSparseSpace::SetToZero(Dx);
//                         TSparseSpace::SetToZero(b);
// 
//                         pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
//                     }
//                     else
//                     {
//                         TSparseSpace::SetToZero(Dx);
//                         TSparseSpace::SetToZero(b);
// 
//                         pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
//                     }
//                 }
//                 else
//                 {
//                     TSparseSpace::SetToZero(Dx);
//                     TSparseSpace::SetToZero(b);
// 
//                     pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
//                 }
//             }
//             else
//             {
//                 std::cout << "ATTENTION: no free DOFs!! " << std::endl;
//             }
//             
//             // Debugging info
//             BaseType::EchoInfo(iteration_number);
//         
//             // Updating the results stored in the database
//             BaseType::UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());
// 
//             // Calculate reactions
//             CalculateContactReactions(pBuilderAndSolver,pScheme, b);
//             
//             pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
//             
//             ResidualIsUpdated = false;
// 
//             if (is_converged == true)
//             {
//                 if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
//                 {
//                     TSparseSpace::SetToZero(b);
// 
//                     pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
//                     ResidualIsUpdated = true;
//                     //std::cout << "mb is calculated" << std::endl;
//                 }
// 
//                 is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
//             }
//         }
// 
// 
//         // Plots a warning if the maximum number of iterations is exceeded
//         if (iteration_number >= BaseType::mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
//         {
//             MaxIterationsExceeded();
//         }
// 
//         // Recalculate residual if needed (note that some convergence criteria need it to be recalculated)
//         if (ResidualIsUpdated == false)
//         {
// //             TSparseSpace::SetToZero(mb);
// //             pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), mb);
//         }
// 
//         // Calculate reactions if required
//         CalculateContactReactions(pBuilderAndSolver,pScheme, b);
//         if (BaseType::mCalculateReactionsFlag == true)
//         {
//             pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), A, Dx, b);
//         }

        // TODO: Finish this!!!
        // Plots a warning if the maximum number of iterations is exceeded
        if (is_converged == false)
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
//                 while (is_converged == false && BaseType::GetModelPart().GetProcessInfo()[TIME] <= aux_time && iteration_number < BaseType::mMaxIterationNumber)
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
//                     BaseType::InitiliazeCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);
//                     BaseType::IterationCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb); 
//                 
//                     // Plots a warning if the maximum number of iterations is exceeded
//                     if (is_converged == false  && iteration_number >= BaseType::mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
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
    
    double mSplitFactor;           // Number by one the delta time is split
    unsigned int mMaxNumberSplits; // Maximum number of splits

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * This method calculates the reactions concerning the contact (residual of the contact)
     */
    
    void CalculateContactReactions(
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TSchemeType::Pointer pScheme,
        TSystemVectorType& b
        )
    {
        // We recalulate the RHS (NOTE: EXPENSIVE)
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
        
        // Now we iterate over all the nodes
        NodesArrayType& pNode = BaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(SLAVE) == true)
            {                
                for(auto itDoF = itNode->GetDofs().begin() ; itDoF != itNode->GetDofs().end() ; itDoF++)
                {
                    const int j = (itDoF)->EquationId();
                    (itDoF)->GetSolutionStepReactionValue() = -b[j];
                }
            }
        }
    }
    
    /**
     * This method prints information after solving the problem
     */
    
    void CoutSolvingProblem()
    {
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "STEP: " << BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] << "\t NON LINEAR ITERATION: " << BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] << "\t TIME: " << BaseType::GetModelPart().GetProcessInfo()[TIME] << "\t DELTA TIME: " << BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]  << std::endl;
        }
    }
    
    /**
     * This method prints information after split the increment of time
     */
        
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
    
    /**
     * This method prints information after reach the max number of interations
     */
    
    void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }
    
    /**
     * This method prints information after reach the max number of interations and splits
     */
        
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
