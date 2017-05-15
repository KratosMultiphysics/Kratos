// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY

/* System Includes */

/* External Includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/kratos_parameters.h"
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
        Parameters ThisParameters =  Parameters(R"({})")
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        Parameters DefaultParameters = Parameters(R"(
        {
            "split_factor"                     : 10.0,
            "max_number_splits"                : 3,
            "rescale_factor"                   : false,
            "path_following_penalty"           : false
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mSplitFactor = ThisParameters["split_factor"].GetDouble();
        mMaxNumberSplits = ThisParameters["max_number_splits"].GetInt();
        mRecalculateFactor = ThisParameters["rescale_factor"].GetBool();
        mPenaltyPathFollowing = ThisParameters["path_following_penalty"].GetBool();
        mInitialPenaltyParameter = 0.0;

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
        Parameters ThisParameters =  Parameters(R"({})")
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag )
    {
        KRATOS_TRY;

        Parameters DefaultParameters = Parameters(R"(
        {
            "split_factor"                     : 10.0,
            "max_number_splits"                : 3,
            "rescale_factor"                   : false,
            "path_following_penalty"           : false
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mSplitFactor = ThisParameters["split_factor"].GetDouble();
        mMaxNumberSplits = ThisParameters["max_number_splits"].GetInt();
        mRecalculateFactor = ThisParameters["rescale_factor"].GetBool();
        mPenaltyPathFollowing = ThisParameters["path_following_penalty"].GetBool();
        mInitialPenaltyParameter = 0.0;

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
        
//         // We update the contact reactions
//         TSystemVectorType& b = *BaseType::mpb;
//         CalculateContactReactions(b);
        
//         // TODO: Finish this!!!
//         Plots a warning if the maximum number of iterations is exceeded
//         if (is_converged == false)
//         {
//             double original_delta_time = StrategyBaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]; // We save the delta time to restore later
//             
//             // We iterate until we reach the convergence or we split more than desired
//             while (is_converged == false && split_number <= mMaxNumberSplits)
//             {   
//                 // Expliting time step as a way to try improve the convergence
//                 split_number += 1;
//                 iteration_number = 1;
//                 
//                 double aux_time       = StrategyBaseType::GetModelPart().GetProcessInfo()[TIME];
//                 double aux_delta_time = StrategyBaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]; // FIXME: The DELTA_TIME is set to 0 for some reason!!!!
//                 double current_time   = aux_time - aux_delta_time;
//                 
//                 StrategyBaseType::GetModelPart().GetProcessInfo()[TIME]       =   current_time; // Restore time to the previous one
//                 aux_delta_time /= mSplitFactor;
//                 StrategyBaseType::GetModelPart().GetProcessInfo()[DELTA_TIME] = aux_delta_time; // Change delta time
//                 
//                 CoutSplittingTime(aux_delta_time);
//                 
//                 unsigned int aux_cout = 0;
//                 while (is_converged == false && StrategyBaseType::GetModelPart().GetProcessInfo()[TIME] <= aux_time && iteration_number < BaseType::mMaxIterationNumber)
//                 {      
//                     iteration_number = 1;
//                     current_time += aux_delta_time;
//          
//                     aux_cout += 1;
//                     if (aux_cout > 1) // We avoid to restore the database if we create the new step just after the first iteration
//                     {
//                         StrategyBaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] += 1;
//                     }
// 
//                     StrategyBaseType::GetModelPart().GetProcessInfo()[TIME] = current_time; // Increase the time in the new delta time        
//                     StrategyBaseType::GetModelPart().CloneTimeStep(current_time);
//                     
//                     // We repeat the predict with the new DELTA_TIME
//                     Predict();
//                     
//                     BaseType::InitiliazeCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb);
//                     BaseType::IterationCycle(is_converged, ResidualIsUpdated, iteration_number, pScheme, pBuilderAndSolver, rDofSet, mA, mDx, mb); 
//                 
//                     // Plots a warning if the maximum number of iterations is exceeded
//                     if (is_converged == false  && iteration_number >= BaseType::mMaxIterationNumber && StrategyBaseType::GetModelPart().GetCommunicator().MyPID() == 0)
//                     {
//                         MaxIterationsExceeded();
//                     }
//                 }
//                 
//                 if (is_converged == true)
//                 {
//                     // Restoring original DELTA_TIME
//                     StrategyBaseType::GetModelPart().GetProcessInfo()[DELTA_TIME] = original_delta_time;
//                 }
//             }
//             
//             // Plots a warning if the maximum number of iterations and splits are exceeded
//             if (is_converged == false  && StrategyBaseType::GetModelPart().GetCommunicator().MyPID() == 0)
//             {
//                 MaxIterationsAndSplitsExceeded();
//             }
//         }

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
    
    double mSplitFactor;             // Number by one the delta time is split
    unsigned int mMaxNumberSplits;   // Maximum number of splits
    bool mRecalculateFactor;         // To check if we recalculate or not the scale factor
    bool mPenaltyPathFollowing;      // To check if we recalculate or not the penalty parameter
    double mInitialPenaltyParameter; // The initial penalty parameter

    ///@}
    ///@name Protected Operators
    ///@{
    
    /**
     * Performs all the required operations that should be done (for each step) 
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
        
    void InitializeSolutionStep() override
    {
        BaseType::InitializeSolutionStep();
        
        // Now we rescale the scale factor
        if (mRecalculateFactor == true && StrategyBaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] == 1)
        {
            RescaleFactor();
        }
    }
    
    /**
     * We rescale the scale factor in function of the norm of the RHS
     */
    
    void RescaleFactor()
    {
        // We get the scale factor
        double& ScaleFactor = StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR]; 
        if (ScaleFactor == 0.0)
        {
            KRATOS_ERROR << "You don't have any value assigned to SCALE_FACTOR" << std::endl;
        }
        
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        // We recalculate the RHS
        TSystemVectorType& b = *BaseType::mpb;
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
        
        // We initialize the values of the contact and non contact norm
        double AuxContact = 0.0;
        double AuxNonContact = 0.0;
        
        // Now we iterate over all the nodes
        NodesArrayType& pNode = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
    
            for(auto itDoF = itNode->GetDofs().begin() ; itDoF != itNode->GetDofs().end() ; itDoF++)
            {
                const int j = (itDoF)->EquationId();
                std::size_t CurrVar = (itDoF)->GetVariable().Key();
                
                if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                {          
                    #pragma omp atomic
                    AuxNonContact += b[j] * b[j];
                }
                else // Corresponding with contact
                {
                    #pragma omp atomic
                    AuxContact += b[j] * b[j];
                }
            }
        }
        
        ScaleFactor *= std::sqrt(AuxNonContact/AuxContact);
        
        if (StrategyBaseType::mEchoLevel > 0)
        {
            std::cout << "The new SCALE_FACTOR is: " << ScaleFactor << std::endl;
        }
    }
    
//     /**
//      * We recalculate the penalty parameter using the path following method for the frictionless case
//      */
//     
//     void CalculatePenaltyPathFollowingFrictionless()
//     {
//         // We get the penalty parameter
//         double& PenaltyParameter = StrategyBaseType::GetModelPart().GetProcessInfo()[PENALTY_PARAMETER]; 
//         
//         if (PenaltyParameter == 0.0)
//         {
//             KRATOS_ERROR << "You don't have any value assigned to PENALTY_PARAMETER" << std::endl;
//         }
//         
// //         const double Tolerance = 1.0e-4;
//         
//         // We get the scale factor
//         const double ScaleFactor = StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR]; 
// 
//         // We initialize the values for the path following
//         double Vfunction0      = 0.0;
//         double Vfunction       = 0.0;
//         double DeltaVfunction  = 0.0;
//         
//         // Now we iterate over all the nodes
//         NodesArrayType& pNode = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
//         auto numNodes = pNode.end() - pNode.begin();
//         
//         #pragma omp parallel for
//         for(unsigned int i = 0; i < numNodes; i++)  // TODO: ADDtangent contact
//         {
//             auto itNode = pNode.begin() + i;
//     
//             if (itNode->Is(ACTIVE) == true)
//             {
//                 const double NormalPressure = (itNode)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
//                 const double WeightedGap    = (itNode)->FastGetSolutionStepValue(WEIGHTED_GAP);    
//                 
//                 const double AugmentedNormalPressure  = (ScaleFactor * NormalPressure + PenaltyParameter * WeightedGap);
//                 const double AugmentedNormalPressure0 = (ScaleFactor * NormalPressure + mInitialPenaltyParameter * WeightedGap);
//                                 
//                 // Initial penalized values
//                 #pragma omp atomic
//                 Vfunction0 +=  AugmentedNormalPressure0 * AugmentedNormalPressure0/(2.0 * mInitialPenaltyParameter);
//                 
//                 // Current penalized values 
//                 #pragma omp atomic
//                 Vfunction +=  AugmentedNormalPressure * AugmentedNormalPressure/(2.0 * PenaltyParameter);
//                 #pragma omp atomic
//                 DeltaVfunction += AugmentedNormalPressure * AugmentedNormalPressure/(2.0 * PenaltyParameter * PenaltyParameter) + WeightedGap/(2.0 * PenaltyParameter) * AugmentedNormalPressure;
//             }
//         }
//         
//         const double E = PenaltyParameter * PenaltyParameter * DeltaVfunction * 1.0/(Vfunction - Vfunction0 - PenaltyParameter * DeltaVfunction);
//         const double C2 = (E * (E + PenaltyParameter) * (Vfunction - Vfunction0))/PenaltyParameter;
//         const double C1 = Vfunction0 + C2/E;
// 
//         const double m = C1 - C2/(E + PenaltyParameter); // NOTE: We need to calculate the new one!!!!
//         const double Beta = C1 - m;
//         PenaltyParameter = C2/Beta - E;
// 
//         if (StrategyBaseType::mEchoLevel > 0)
//         {
//             std::cout << "The path following parameters: " << std::endl;
//             std::cout << "E: " << E << " C1: " << C1 << " C2: " << C2 << " m: " << m << " Beta: " << Beta << std::endl;
//             std::cout << "The new PENALTY_PARAMETER is: " << PenaltyParameter << std::endl;
//         }
//     }
//     
//     /**
//      * We recalculate the penalty parameter using the path following method for the frictional case
//      */
//     
//     void CalculatePenaltyPathFollowingFrictional()
//     {
//         // TODO: Finish me!!!
//     }
    
    /**
     * Here the database is updated
     */
     
    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh
    ) override
    {
        BaseType::UpdateDatabase(A,Dx,b,MoveMesh);
        
        // First we activate all the nodes in the active conditions
//         ReactivateNodes();
//         CalculateContactReactions(b); // TODO: Remove this
    }
    
    /**
     * We activate all yhe nodes in the active conditions
     */
    
    void ReactivateNodes()
    {
        // Now we iterate over all the conditions and we set all the nodes as ACTIVE in the active conditions 
        ConditionsArrayType& pCond = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Conditions();
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            
            if (itCond->Is(ACTIVE) == true)
            {
                for (unsigned int iNode = 0; iNode < itCond->GetGeometry().size(); iNode++)
                {
                    itCond->GetGeometry()[iNode].SetLock();
                    itCond->GetGeometry()[iNode].Set(ACTIVE, true);
                    itCond->GetGeometry()[iNode].UnSetLock();
                }
            }
        }
    }
    
    /**
     * This method calculates the reactions concerning the contact (residual of the contact)
     * @param b: The residual vector
     */
    
    void CalculateContactReactions(TSystemVectorType& b)
    {       
        // First we activate all the nodes in the active conditions
        ReactivateNodes();
        
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        const double ScaleFactor = (StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR] > 0.0) ? StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR] : 1.0;
        
        // We recalculate the RHS (NOTE: EXPENSIVE)
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
        
        // Now we iterate over all the nodes
        NodesArrayType& pNode = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {                
                for(auto itDoF = itNode->GetDofs().begin() ; itDoF != itNode->GetDofs().end() ; itDoF++)
                {
                    const int j = (itDoF)->EquationId();
                    
                    if (((itDoF)->GetReaction().Name()).find("WEIGHTED") != std::string::npos) // Corresponding with contact
                    {                        
                        (itDoF)->GetSolutionStepReactionValue() = b[j]/ScaleFactor;
                    }
                    else if ((itDoF)->GetReaction().Name() != "NONE") // The others
                    {                        
                        (itDoF)->GetSolutionStepReactionValue() = -b[j];
                    }
                }
            }
        }
        
//         // We recalculate the penalty parameter in case we want
//         if (mPenaltyPathFollowing == true)
//         {
//             if (StrategyBaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] == 1)
//             {
//                 mInitialPenaltyParameter = StrategyBaseType::GetModelPart().GetProcessInfo()[PENALTY_PARAMETER];
//             }
// 
//             if (StrategyBaseType::GetModelPart().Is(SLIP) == true)
//             {
//                 CalculatePenaltyPathFollowingFrictional();
//             }
//             else
//             {
//                 CalculatePenaltyPathFollowingFrictionless();
//             }
//         }
    }
    
    /**
     * This method prints information after solving the problem
     */
    
    void CoutSolvingProblem()
    {
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "STEP: " << StrategyBaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] << "\t NON LINEAR ITERATION: " << StrategyBaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] << "\t TIME: " << StrategyBaseType::GetModelPart().GetProcessInfo()[TIME] << "\t DELTA TIME: " << StrategyBaseType::GetModelPart().GetProcessInfo()[DELTA_TIME]  << std::endl;
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
            std::cout << "***\t\t COMING BACK TO TIME: " << StrategyBaseType::GetModelPart().GetProcessInfo()[TIME] << "\t\t ***" << std::endl;
            std::cout << "***\t\t NEW TIME STEP: "<< aux_delta_time << "\t\t ***" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }
    
    /**
     * This method prints information after reach the max number of interations
     */
    
    void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && StrategyBaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
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
        if (this->GetEchoLevel() != 0 && StrategyBaseType::GetModelPart().GetCommunicator().MyPID() == 0 )
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
