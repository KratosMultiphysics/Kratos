//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:    BSD License
//              Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_RESIDUALBASED_DEM_COUPLED_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_RESIDUALBASED_DEM_COUPLED_NEWTON_RAPHSON_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"

#include "custom_processes/update_dem_kinematics_process.h"
#include "custom_processes/transfer_nodal_forces_to_fem.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
 * @class ResidualBasedDEMCoupledNewtonRaphsonStrategy
 * @ingroup KratosFEMDEM App
 * @brief This is the base Newton Raphson strategy coupled with the DEM strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Alejandro Cornejo
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedDEMCoupledNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    typedef ExplicitSolverStrategy ExplicitSolverStrategyType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedDEMCoupledNewtonRaphsonStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedDEMCoupledNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;
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

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedDEMCoupledNewtonRaphsonStrategy() : BaseType()
    {
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit ResidualBasedDEMCoupledNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        ExplicitSolverStrategyType::Pointer pDEMStrategy,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart,
                   pScheme,
                   pNewConvergenceCriteria,
                   pNewBuilderAndSolver,
                   MaxIterations,
                   CalculateReactions,
                   ReformDofSetAtEachStep,
                   MoveMeshFlag),
          mpDEMStrategy(pDEMStrategy)
    {
    }


    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ResidualBasedDEMCoupledNewtonRaphsonStrategy() override
    {
        auto p_builder_and_solver = this->GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->Clear();
        }

        this->mpA.reset();
        this->mpDx.reset();
        this->mpb.reset();

        this->Clear();
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/


    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;
        BaseType::Initialize();
        mpDEMStrategy->Initialize();
        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;
        BaseType::InitializeSolutionStep();
        mpDEMStrategy->InitializeSolutionStep();
        TransferNodalForcesToFem(this->GetModelPart(), false).Execute();
        UpdateDemKinematicsProcess(this->GetModelPart()).Execute();
        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;
        BaseType::FinalizeSolutionStep();
        mpDEMStrategy->FinalizeSolutionStep();
        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // We compute the contact forces with the DEM
        auto update_dem_kinematics_process = UpdateDemKinematicsProcess(this->GetModelPart());
        update_dem_kinematics_process.Execute();
        mpDEMStrategy->SolveSolutionStep();
        auto transfer_process = TransferNodalForcesToFem(this->GetModelPart(), false);
        transfer_process.Execute();
        update_dem_kinematics_process.Execute();

        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *this->mpA;
        TSystemVectorType& rDx = *this->mpDx;
        TSystemVectorType& rb  = *this->mpb;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        this->mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = this->mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (this->mUseOldStiffnessInFirstIteration){
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(p_scheme, r_model_part, rA, rDx, rb, BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx);  // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        this->EchoInfo(iteration_number);

        // Updating the results stored in the database
        this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        this->mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (is_converged) {
            if (this->mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }
            is_converged = this->mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < this->mMaxIterationNumber) {
            // We compute the contact forces with the DEM
            update_dem_kinematics_process.Execute();
            mpDEMStrategy->SolveSolutionStep();
            auto transfer_process_damped = TransferNodalForcesToFem(this->GetModelPart(), true);
            transfer_process_damped.Execute();
            update_dem_kinematics_process.Execute();

            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            this->mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = this->mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0) {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false) {
                    if (this->GetKeepSystemConstantDuringIterations() == false) {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    } else {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                } else {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            this->EchoInfo(iteration_number);

            // Updating the results stored in the database
            this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            this->mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true) {
                if (this->mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                    TSparseSpace::SetToZero(rb);
                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }
                is_converged = this->mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= this->mMaxIterationNumber) {
            this->MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedDEMCoupledNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << this->mMaxIterationNumber << " iterations" << std::endl;
        }

        // Calculate reactions if required
        if (this->mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        return is_converged;
    }

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        const DataCommunicator &r_comm = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator();
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (this->mInitializeWasPerformed == false)
            BaseType::Initialize();

        //initialize solution step
        BaseType::InitializeSolutionStep();
        mpDEMStrategy->InitializeSolutionStep();

        TSystemMatrixType& rA  = *this->mpA;
        TSystemVectorType& rDx = *this->mpDx;
        TSystemVectorType& rb  = *this->mpb;

        DofsArrayType& r_dof_set = this->GetBuilderAndSolver()->GetDofSet();

        this->GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        // Applying constraints if needed
        auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        const int local_number_of_constraints = r_constraints_array.size();
        const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        if (global_number_of_constraints != 0) {
            const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

            const auto it_const_begin = r_constraints_array.begin();

            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(local_number_of_constraints); ++i)
                (it_const_begin + i)->ResetSlaveDofs(r_process_info);

            #pragma omp parallel for
            for (int i=0; i<static_cast<int>(local_number_of_constraints); ++i)
                 (it_const_begin + i)->Apply(r_process_info);

            // The following is needed since we need to eventually compute time derivatives after applying
            // Master slave relations
            TSparseSpace::SetToZero(rDx);
            this->GetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true)
            BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Operators

    ///@{

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
        return "ResidualBasedDEMCoupledNewtonRaphsonStrategy";
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
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

  protected:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename ExplicitSolverStrategy::Pointer mpDEMStrategy = nullptr;

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

    /**
     * Copy constructor.
     */

    ResidualBasedDEMCoupledNewtonRaphsonStrategy(const ResidualBasedDEMCoupledNewtonRaphsonStrategy &Other){};

    ///@}

}; /* Class ResidualBasedDEMCoupledNewtonRaphsonStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_RESIDUALBASED_DEM_COUPLED_NEWTON_RAPHSON_STRATEGY  defined */
