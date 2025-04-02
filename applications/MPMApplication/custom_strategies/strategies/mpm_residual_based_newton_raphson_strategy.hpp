//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY )
#define      KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY

/* System includes */

/* External includes */
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "mpm_application_variables.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/// Short class definition.

/**
 * @class MPMResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosMPM
 * @brief Newton Raphson strategy suited for MPM simulations
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is achieved) using a Newton Raphson algorithm
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class MPMResidualBasedNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(MPMResidualBasedNewtonRaphsonStrategy);

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructors.
     */

     MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        bool MoveMeshFlag = false
    ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rModelPart, MoveMeshFlag)
    {
    }

    MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = true,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria,
            MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
    }

    MPMResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = true,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
            rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver,
            MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
    }

    /** Destructor.
     */
    virtual ~MPMResidualBasedNewtonRaphsonStrategy()
    {
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();

        TSystemMatrixType& rA = *(this->mpA);
        TSystemVectorType& rDx = *(this->mpDx);
        TSystemVectorType& rb = *(this->mpb);
        DofsArrayType& r_dof_set = p_builder_and_solver->GetDofSet();

        // Initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;

        p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);
        is_converged = this->mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "PreCriteria:"
            << "\tIs_converged: " << is_converged << "\tmRebuildLevel: " << BaseType::mRebuildLevel
            << "\tmStiffnessMatrixIsBuilt: " << BaseType::mStiffnessMatrixIsBuilt << std::endl;

        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "SetToZero the matrix and vectors of the system" << std::endl;

            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build and Solve" << std::endl;

            p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx); // rDx=0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
            KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "BuildRHSAndSolve" << std::endl;
        }


        if (this->GetEchoLevel() == 3) // If it is needed to print the debug info
        {
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "solution obtained = " << rDx << std::endl;
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) // Print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char*)(matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << iteration_number << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char*)(matrix_market_vectname.str()).c_str(), rb);
        }

        // Update results
        r_dof_set = p_builder_and_solver->GetDofSet();
        p_scheme->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        if (is_converged == true)
        {
            // Initialisation of the convergence criteria
            this->mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

            if (this->mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                // due to the assembly procedure the friction related values need to be reset
                const bool friction_active = BaseType::GetModelPart().GetProcessInfo()[FRICTION_ACTIVE];
                if (friction_active){
                    for (Node &curr_node : BaseType::GetModelPart().Nodes()){
                    curr_node.SetValue(FRICTION_ASSIGNED, false);
                    curr_node.FastGetSolutionStepValue(STICK_FORCE).clear();
                    }
                }

                TSparseSpace::SetToZero(rb);
                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
            }

            is_converged = this->mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }
        KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3 && !is_converged) << "Starting Nonlinear iteration" << std::endl;

        // Iteration Loop
        while (is_converged == false &&
            iteration_number++ < this->mMaxIterationNumber)
        {
            SolveSolutionStepIteration(p_scheme, p_builder_and_solver, rA, rDx, rb,
                                                                       r_dof_set, iteration_number, is_converged);
        }

        // Plot a warning if the maximum number of iterations is exceeded
        if (iteration_number >= this->mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            if (this->GetEchoLevel() > 1) this->MaxIterationsExceeded();
        }


        // recompute first timestep if friction is active
        const bool friction_active = BaseType::GetModelPart().GetProcessInfo()[FRICTION_ACTIVE];
        const bool is_initial_loop = (BaseType::GetModelPart().GetProcessInfo()[STEP] ==  1);

        bool mCalculateReactionsFlag = true;


        //calculate reactions if required
        if (mCalculateReactionsFlag == true){
            // due to the assembly procedure the friction related values need to be reset
            if (friction_active){
                for (Node &curr_node : BaseType::GetModelPart().Nodes()){
                curr_node.SetValue(FRICTION_ASSIGNED, false);
                curr_node.FastGetSolutionStepValue(STICK_FORCE).clear();
                }
            }

            p_builder_and_solver->CalculateReactions(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }

        if (is_initial_loop && friction_active){

            BaseType::GetModelPart().GetProcessInfo()[INITIAL_LOOP_COMPLETE] = true;
            // Transfer STICK_FORCE into a 'previous' timestep
            // [since friction algorithm requires values of normal forces at the 'previous' timestep]
            for (Node &curr_node : BaseType::GetModelPart().Nodes()){
                curr_node.FastGetSolutionStepValue(STICK_FORCE, 1) = curr_node.FastGetSolutionStepValue(STICK_FORCE, 0);
            }

            is_converged=false;

	    KRATOS_INFO("") << "\n";
            KRATOS_INFO("MPMNewtonRaphsonStrategy") << "Running another Newton-Rhapson loop for friction" << std::endl;


            // Reset iteration_number and run another loop
            iteration_number = 1;
            while (is_converged == false &&
                   iteration_number++ < this->mMaxIterationNumber)
            {
                SolveSolutionStepIteration(p_scheme, p_builder_and_solver, rA, rDx, rb,
                                           r_dof_set, iteration_number, is_converged);
            }

            // Plot a warning if the maximum number of iterations is exceeded
            if (iteration_number >= this->mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                if (this->GetEchoLevel() > 1) this->MaxIterationsExceeded();
            }

            if (mCalculateReactionsFlag == true){
                // due to the assembly procedure the friction related values need to be reset
                if (friction_active){
                    for (Node &curr_node : BaseType::GetModelPart().Nodes()){
                    curr_node.SetValue(FRICTION_ASSIGNED, false);
                    curr_node.FastGetSolutionStepValue(STICK_FORCE).clear();
                    }
                }
                
                p_builder_and_solver->CalculateReactions(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
            }
        }

        return is_converged;
    }

private:
    /**
     * @brief Auxilliary function to perform a single Newton-Rhapson iteration. For use in the SolveSolutionStep
     */
    void inline SolveSolutionStepIteration(typename TSchemeType::Pointer p_scheme,
                                    typename TBuilderAndSolverType::Pointer p_builder_and_solver,
                                    TSystemMatrixType& rA,
                                    TSystemVectorType& rDx,
                                    TSystemVectorType& rb,
                                    DofsArrayType& r_dof_set,
                                    const unsigned int iteration_number,
                                    bool& is_converged)
    { // Setting the number of iteration
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);
        is_converged = this->mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        // Call the linear system solver to find the correction rDx. It is not called if there is no system to solve
        if (SparseSpaceType::Size(rDx) != 0)
        {
            if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
            {
                KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Iteration Number: " << iteration_number << std::endl;

                if (this->GetKeepSystemConstantDuringIterations() == false)        
		{
                    TSparseSpace::SetToZero(rA);
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build and Solve" << std::endl;
                    p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build RHS and Solve" << std::endl;
                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                }
            }
            else
            {
                TSparseSpace::SetToZero(rDx);
                TSparseSpace::SetToZero(rb);

                KRATOS_INFO_IF("MPMNewtonRaphsonStrategy", this->GetEchoLevel() >= 3) << "Build RHS and Solve" << std::endl;
                p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
            }
        }

        std::cout<<"Number of non-linear iterations: "<<iteration_number<<"\n";

        // Plot a warning if the maximum number of iterations is exceeded
        if (iteration_number >= this->mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            if (this->GetEchoLevel() > 1) this->MaxIterationsExceeded();


            std::cout<<"MAX ITERATIONS EXCEEDED!!!!!!!!!!!!!!! "<<"\n";
        }
        else
        {
            KRATOS_WARNING("MPMNewtonRaphsonStrategy") << "ATTENTION: no free DOFs!! " << std::endl;

        }

        // Updating the results stored in the database
        r_dof_set = p_builder_and_solver->GetDofSet();

        p_scheme->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        // If converged
        if (is_converged == true)
        {
            if (this->mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                // due to the assembly procedure the friction related values need to be reset
                const bool friction_active = BaseType::GetModelPart().GetProcessInfo()[FRICTION_ACTIVE];
                if (friction_active){
                    for (Node &curr_node : BaseType::GetModelPart().Nodes()){
                    curr_node.SetValue(FRICTION_ASSIGNED, false);
                    curr_node.FastGetSolutionStepValue(STICK_FORCE).clear();
                    }
                }
                
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);

            }

            is_converged = this->mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }
    }

}; /* Class MPMResidualBasedNewtonRaphsonStrategy */

} /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_NEWTON_RAPHSON_STRATEGY  defined */
