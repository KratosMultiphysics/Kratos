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
     * Here the database is updated
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pBuilderAndSolver: The builder and solver considered
     * @param rDofSet: The set of degrees of freedom
     * @param mA: The LHS of the problem
     * @return mDx: The solution to the problem
     */
    void UpdateDatabase( 
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh,
        bool& is_converged,
        bool& ResidualIsUpdated,
        const bool is_iteration
    ) override
    {
        BaseType::UpdateDatabase(A, Dx, b, MoveMesh, is_converged, ResidualIsUpdated, is_iteration);
        
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
//         unsigned int accel_it = 0;
//         unsigned int MaxNumberAccel = 50;
//         while (is_converged == false && accel_it++ < MaxNumberAccel)
        if (is_converged == false)// && is_iteration == true)
        {
//             if (accel_it == 0)
            if (this->GetEchoLevel() != 0)
            {
                std::cout << "Applying the convergence accelerator" << std::endl;
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
            
            // Calculate the new displacement
            Vector tmp = updated_x;
            mpConvergenceAccelerator->InitializeNonLinearIteration();  
            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), A, Dx, b);
//             mpConvergenceAccelerator->UpdateSolution(b, Dx);
            mpConvergenceAccelerator->UpdateSolution(b, updated_x);
            mpConvergenceAccelerator->FinalizeNonLinearIteration();   
        
            // Update residual variables
            Dx = updated_x - tmp;
            
            BaseType::UpdateDatabase(A, Dx, b, MoveMesh, is_converged, ResidualIsUpdated, true);
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
