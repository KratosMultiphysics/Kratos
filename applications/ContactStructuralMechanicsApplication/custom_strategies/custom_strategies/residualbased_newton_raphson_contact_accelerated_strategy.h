// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY

/* System Includes */

/* External Includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/contact_utilities.h"
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"

// Convergence accelerators
#include "../FSIapplication/custom_utilities/convergence_accelerator.hpp"

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
         
class ResidualBasedNewtonRaphsonContactAcceleratedStrategy :
    public ResidualBasedNewtonRaphsonContactStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonContactAcceleratedStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>                   StrategyBaseType;
    
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>  NewtonBaseType;
    
    typedef ResidualBasedNewtonRaphsonContactStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                               TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                       TDataType;

    typedef TSparseSpace                                                                 SparseSpaceType;

    typedef typename BaseType::TSchemeType                                                   TSchemeType;
    
    typedef ConvergenceAccelerator<TDenseSpace>                              TConvergenceAcceleratorType;

    typedef typename BaseType::DofsArrayType                                               DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                       TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                       TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                               LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                               LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                         TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                         TSystemVectorPointerType;
    
    typedef ModelPart::NodesContainerType                                                 NodesArrayType;
    
    typedef ModelPart::ConditionsContainerType                                       ConditionsArrayType;
    
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
     * @param ConvergenceAccelerator: The convergence accelerator to use
     */
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = true,
        typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator = nullptr,
        Parameters ThisParameters =  Parameters(R"({})"
        )
        : ResidualBasedNewtonRaphsonContactStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        Parameters DefaultParameters = Parameters(R"(
        {
            "max_number_accel_iterations" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        mMaxNumberAccelerationIterations = ThisParameters["max_number_accel_iterations"].GetInt();
        
        // Saving convergence accelerator
        mpConvergenceAccelerator = pConvergenceAccelerator;

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
     * @param ConvergenceAccelerator: The convergence accelerator to use
     */
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        unsigned int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = true,
        typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator = nullptr
    )
        : ResidualBasedNewtonRaphsonContactStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag )
    {
        KRATOS_TRY;

        // Saving convergence accelerator
        mpConvergenceAccelerator = pConvergenceAccelerator;
        
        KRATOS_CATCH("");
    }

    /** 
     * Destructor.
     */
    
    virtual ~ResidualBasedNewtonRaphsonContactAcceleratedStrategy()
    {
    }
    
    /******************** OPERATIONS ACCESSIBLE FROM THE INPUT: ************************/
    /***********************************************************************************/
    
    /**
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    
    bool SolveSolutionStep() override
    {
        bool IsConverged = BaseType::SolveSolutionStep();
    
        TSystemMatrixType& A = *NewtonBaseType::mpA;
        TSystemVectorType& Dx = *NewtonBaseType::mpDx;
        TSystemVectorType& b = *NewtonBaseType::mpb;
        
        ResidualIsUpdated = false;
            
        ApplyAcceleration(A, Dx, b, StrategyBaseType::MoveMeshFlag(), IsConverged, ResidualIsUpdated);
        
        return IsConverged;
    }
    
    /**
     * Initialization of member variables and prior operations
     */
    
    void Initialize() override
    {
        KRATOS_TRY;
        
        BaseType::Initialize();
        
        mpConvergenceAccelerator->Initialize();

        KRATOS_CATCH("");
    }

    /**
     * Clears the internal storage
     */
    
    void Clear() override
    {
        KRATOS_TRY;
        
        BaseType::Clear();

        KRATOS_CATCH("");
    }

   /**
    * Performs all the required operations that should be done (for each step) before solving the solution step.
    * A member variable should be used as a flag to make sure this function is called only once per step.
    */
   
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep();
        
        mpConvergenceAccelerator->InitializeSolutionStep();

        KRATOS_CATCH("");
    }
    
    /**
     * Performs all the required operations that should be done (for each step) after solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
    
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;
        
        BaseType::FinalizeSolutionStep();
        
        mpConvergenceAccelerator->FinalizeSolutionStep();

        KRATOS_CATCH("");
    }
 
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
        BaseType::UpdateDatabase(A, Dx, b, MoveMesh);
        
        // TODO: Separate residual by DoF type
    }

    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    
    int Check() override
    {
        KRATOS_TRY;

        StrategyBaseType::Check();

//         StrategyBaseType::mpConvergenceCriteria->Check(StrategyBaseType::GetModelPart());

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
    
    typename TConvergenceAcceleratorType::Pointer mpConvergenceAccelerator;
    unsigned int mMaxNumberAccelerationIterations;

    ///@}
    ///@name Protected Operators
    ///@{
 
    /**
     * Here the convergence accelerator is applied
     */
    
    void ApplyAcceleration( 
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh,
        bool& IsConverged,
        bool& ResidualIsUpdated
    )
    {
        typename TSchemeType::Pointer pScheme = StrategyBaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = StrategyBaseType::GetBuilderAndSolver();
        
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "  Applying the convergence accelerator" << std::endl;
        }
        
        while (IsConverged == false && iAccel++ < mMaxNumberAccelerationIterations)
        {
            pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
            
            if (this->GetEchoLevel() != 0)
            {
                std::cout << "  CONV_IT: " << iAccel << " RESIDUAL NORM: " << norm_2(b) << std::endl;
            }
        
            Vector UpdatedX( Dx.size() ); //UpdatedX // TODO: Divide residual
            UpdatedX.clear();
            DofsArrayType DofsArray = pBuilderAndSolver->GetDofSet();
            for(typename DofsArrayType::iterator iDof = DofsArray.begin() ; iDof != DofsArray.end() ; ++iDof)
            {
                if (iDof->IsFree() == true)
                {
                    UpdatedX[ iDof->EquationId() ] = iDof->GetSolutionStepValue();
                }
            }

            // Calculate the new displacement
            Vector tmp = UpdatedX;
            mpConvergenceAccelerator->InitializeNonLinearIteration();  
            pScheme->InitializeNonLinIteration(StrategyBaseType::GetModelPart(), A, Dx, b);
//             IsConverged = StrategyBaseType::mpConvergenceCriteria->PreCriteria(StrategyBaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            mpConvergenceAccelerator->UpdateSolution(Dx, UpdatedX);
            mpConvergenceAccelerator->FinalizeNonLinearIteration();   
            
            // Update residual variables
            Dx = UpdatedX - tmp;
             
            BaseType::UpdateDatabase(A, Dx, b, MoveMesh);
            
            pScheme->FinalizeNonLinIteration(StrategyBaseType::GetModelPart(), A, Dx, b);

            ResidualIsUpdated = false;

            if (IsConverged == true)
            {
                if (StrategyBaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {                    
                    TSparseSpace::SetToZero(b);
                                
                    pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
                    ResidualIsUpdated = true;
                }
                                    
                IsConverged = StrategyBaseType::mpConvergenceCriteria->PostCriteria(StrategyBaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            }
        }
        
        if (this->GetEchoLevel() != 0)
        {
            std::cout << "  Finishing the convergence accelerator" << std::endl;
        }
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
    
    ResidualBasedNewtonRaphsonContactAcceleratedStrategy(const ResidualBasedNewtonRaphsonContactAcceleratedStrategy& Other)
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

}; /* Class ResidualBasedNewtonRaphsonContactAcceleratedStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_ACCELERATED_STRATEGY */
