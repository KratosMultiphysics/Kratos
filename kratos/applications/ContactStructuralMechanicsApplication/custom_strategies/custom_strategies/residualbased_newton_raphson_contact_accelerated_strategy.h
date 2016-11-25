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
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
    public ResidualBasedNewtonRaphsonStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonContactAcceleratedStrategy );

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef TSparseSpace                                                          SparseSpaceType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;
    
    typedef ConvergenceAccelerator<TDenseSpace>                       TConvergenceAcceleratorType;

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
        typename TConvergenceAcceleratorType::Pointer pConvergenceAccelerator = nullptr
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

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
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag )
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
     * Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
 
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();

        TSystemMatrixType& A  = *BaseType::mpA;
        TSystemVectorType& Dx = *BaseType::mpDx;
        TSystemVectorType& b  = *BaseType::mpb;

        //initializing the parameters of the Newton-Raphson cicle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        //                      BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
        bool is_converged = false;
        bool ResidualIsUpdated = false;
        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

        //function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(A);
            TSparseSpace::SetToZero(Dx);
            TSparseSpace::SetToZero(b);

            pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx); //Dx=0.00;
            TSparseSpace::SetToZero(b);

            pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }

        if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            //                          std::cout << "After first system solution" << std::endl;
            std::cout << "SystemMatrix = " << A << std::endl;
            std::cout << "solution obtained = " << Dx << std::endl;
            std::cout << "RHS  = " << b << std::endl;
        }
        if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), A, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), b);
        }
        
        // Updating the results stored in the database
        BaseType::UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());

        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

        if (is_converged == true)
        {
            //initialisation of the convergence criteria
            BaseType::mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(b);

                pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
        }

        // Acceleration is applied
        ResidualIsUpdated = false;
        ApplyAcceleration( A, Dx, b, BaseType::MoveMeshFlag(), is_converged, ResidualIsUpdated);
        
        //Iteration Cicle... performed only for NonLinearProblems
        while (is_converged == false &&
                iteration_number++<BaseType::mMaxIterationNumber)
        {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(Dx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false )
                {
                    if( BaseType::GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(A);
                        TSparseSpace::SetToZero(Dx);
                        TSparseSpace::SetToZero(b);

                        pBuilderAndSolver->BuildAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(Dx);
                        TSparseSpace::SetToZero(b);

                        pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(Dx);
                    TSparseSpace::SetToZero(b);

                    pBuilderAndSolver->BuildRHSAndSolve(pScheme, BaseType::GetModelPart(), A, Dx, b);
                }
            }
            else
            {
                std::cout << "ATTENTION: no free DOFs!! " << std::endl;
            }

            
            // Updating the results stored in the database
            BaseType::UpdateDatabase(A, Dx, b, BaseType::MoveMeshFlag());

            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {

                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(b);

                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
                    ResidualIsUpdated = true;
                    //std::cout << "b is calculated" << std::endl;
                }

                is_converged =BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            }
            
            // Acceleration is applied
            ResidualIsUpdated = false;
            ApplyAcceleration( A, Dx, b, BaseType::MoveMeshFlag(), is_converged, ResidualIsUpdated);
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            MaxIterationsExceeded();

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (ResidualIsUpdated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(b);
            //    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
        }

        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
        {
            pBuilderAndSolver->CalculateReactions(pScheme, BaseType::GetModelPart(), A, Dx, b);
        }

        return is_converged;
    }
 
    /**
     * Here the convergence accelerator is applied
     */
    void ApplyAcceleration( 
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh,
        bool& is_converged,
        bool& ResidualIsUpdated
    )
    {
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
//         unsigned int accel_it = 0;
//         unsigned int MaxNumberAccel = 5;
//         while (is_converged == false && accel_it++ < MaxNumberAccel)
        if (is_converged == false)// && is_iteration == true)
        {
//             if (accel_it == 0)
            if (this->GetEchoLevel() != 0)
            {
                std::cout << "  Applying the convergence accelerator" << std::endl;
            }
//             std::cout << "  CONV_IT: " << accel_it << " RESIDUAL NORM: " << norm_2(b) << std::endl;
           
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
        
            Vector updated_x( Dx.size() ); //updated_x
            updated_x.clear();
            for(typename DofsArrayType::iterator i_dof = this->GetBuilderAndSolver()->GetDofSet().begin() ; i_dof != this->GetBuilderAndSolver()->GetDofSet().end() ; ++i_dof)
            {
                if (i_dof->IsFree() == true)
                {
                    updated_x[ i_dof->EquationId() ] = i_dof->GetSolutionStepValue();
                }
            }
            
            // We calculate the norm of b
            const double normA = SparseSpaceType::TwoNorm(A);

            // Calculate the new displacement
            Vector tmp = updated_x;
            mpConvergenceAccelerator->InitializeNonLinearIteration();  
            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            mpConvergenceAccelerator->UpdateSolution(b/normA, updated_x);
            mpConvergenceAccelerator->FinalizeNonLinearIteration();   
            
            // Update residual variables
            Dx = updated_x - tmp;
            
            BaseType::UpdateDatabase(A, Dx, b, MoveMesh);
            
            pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);

            ResidualIsUpdated = false;

            if (is_converged == true)
            {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {                    
                    TSparseSpace::SetToZero(b);
                                
                    pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b);
                    ResidualIsUpdated = true;
                }
                                    
                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), pBuilderAndSolver->GetDofSet(), A, Dx, b);
            }
            
            if (this->GetEchoLevel() != 0)
//             if (is_converged == true || accel_it == MaxNumberAccel)
            {
                std::cout << "  Finishing the convergence accelerator" << std::endl;
            }
        }
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
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    
//     int Check() override
//     {
//         KRATOS_TRY;
// 
//         BaseType::Check();
// 
//         BaseType::mpConvergenceCriteria->Check(BaseType::GetModelPart());
// 
//         return 0;
// 
//         KRATOS_CATCH("");
//     }

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
