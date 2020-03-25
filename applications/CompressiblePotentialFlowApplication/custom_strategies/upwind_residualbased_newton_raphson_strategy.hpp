//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

#if !defined(KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY)
#define  KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
template<class TSparseSpace,
class TDenseSpace, // = DenseSpace<double>,
class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>
class UpwindResidualBasedNewtonRaphsonStrategy
: public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(UpwindResidualBasedNewtonRaphsonStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// DoF array type definition
    typedef ModelPart::DofsArrayType DofsArrayType;

    /// Data type definition
    typedef typename TSparseSpace::DataType TDataType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
    * Constructor.
    */

    // constructor without Builder and Solver
    UpwindResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
    )
    : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme,
        pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
        MoveMeshFlag)
    {
    }


    /**
     * Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    UpwindResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme,
        pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
        MoveMeshFlag)
    {
        BaseType::SetEchoLevel(0);
    }

    // Destructor
    ~UpwindResidualBasedNewtonRaphsonStrategy() = default;

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        std::ofstream outfile;
        BaseType::SetEchoLevel(0);
        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename BaseType::TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        typename BaseType::TSystemMatrixType& rA  = *BaseType::mpA;
        typename BaseType::TSystemVectorType& rDx = *BaseType::mpDx;
        typename BaseType::TSystemVectorType& rb  = *BaseType::mpb;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // TSparseSpace::SetToZero(rb);
        // p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        // double absolute_residual_norm0 = CalculateResidualNorm(r_model_part, r_dof_set, rb);
        // KRATOS_INFO("Absolute norm = ") << absolute_residual_norm0 << std::endl;

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            //setting up the Vectors involved to the correct size
            BuiltinTimer system_matrix_resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(
                p_scheme, BaseType::mpA, BaseType::mpDx, BaseType::mpb, r_model_part);
            KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0)
                << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }else {
            TSparseSpace::SetToZero(rDx); //Dx=0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        BaseType::EchoInfo(iteration_number);

        // TSparseSpace::SetToZero(rb);
        // p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        // double absolute_residual_norm = CalculateResidualNorm(r_model_part, r_dof_set, rb);
        // double relaxation_factor = ComputeRelaxationFactor(absolute_residual_norm);
        // KRATOS_INFO("Absolute norm = ") << absolute_residual_norm << std::endl;
        // KRATOS_INFO("Relaxation factor = ") << relaxation_factor << std::endl;
        // KRATOS_WATCH(rDx(692))
        // rDx *= relaxation_factor;
        // KRATOS_WATCH(rDx(692))
        double displacement_norm = CalculateFinalCorrectionNorm(r_dof_set, rDx);
        KRATOS_INFO("DISPLACEMENT Absolute norm = ") << displacement_norm << std::endl;

        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        TSparseSpace::SetToZero(rb);
        p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        double absolute_residual_norm = CalculateResidualNorm(r_dof_set, rb);
        KRATOS_INFO("RESIDUAL Absolute norm = ") << absolute_residual_norm << std::endl;

        outfile.open("/media/inigo/10740FB2740F9A1C/2d_results_test/plots/convergence_results.dat", std::ios_base::app);
        outfile << "   " << r_model_part.GetProcessInfo()[STEP]
        << "      " << r_model_part.GetProcessInfo()[FREE_STREAM_MACH]
        << "      " << absolute_residual_norm
        << "      " << displacement_norm;
        outfile.close();

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        // TSparseSpace::SetToZero(rb);
        // p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        // absolute_residual_norm = CalculateResidualNorm(r_model_part, r_dof_set, rb);
        // KRATOS_INFO("Absolute norm = ") << absolute_residual_norm << std::endl;

        if (is_converged) {
            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        // TSparseSpace::SetToZero(rb);
        // p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        // absolute_residual_norm = CalculateResidualNorm(r_model_part, r_dof_set, rb);
        // KRATOS_INFO("Absolute norm = ") << absolute_residual_norm << std::endl;

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < BaseType::mMaxIterationNumber)
        {
            std::cout << std::endl;
            KRATOS_INFO("ITERATION NUMBER ") << iteration_number << std::endl;
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (BaseType::SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (BaseType::GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        //setting up the Vectors involved to the correct size
                        BuiltinTimer system_matrix_resize_time;
                        p_builder_and_solver->ResizeAndInitializeVectors(
                            p_scheme, BaseType::mpA, BaseType::mpDx, BaseType::mpb, r_model_part);
                        KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0)
                            << system_matrix_resize_time.ElapsedSeconds() << std::endl;

                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                    else
                    {
                        KRATOS_WATCH(BaseType::mMaxIterationNumber)
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                }
                else
                {
                    KRATOS_WATCH(iteration_number)
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(iteration_number);

            // absolute_residual_norm = CalculateResidualNorm(r_model_part, r_dof_set, rb);
            // relaxation_factor = ComputeRelaxationFactor(absolute_residual_norm);
            // KRATOS_INFO("Absolute norm = ") << absolute_residual_norm << std::endl;
            // KRATOS_INFO("Relaxation factor = ") << relaxation_factor << std::endl;
            // KRATOS_WATCH(rDx(692))
            // rDx *= relaxation_factor;
            // KRATOS_WATCH(rDx(692))

            // Updating the results stored in the database
            UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            BaseType::MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("NR-Strategy", BaseType::GetEchoLevel() > -1)
                << "Convergence achieved after " << iteration_number << " / "
                << BaseType::mMaxIterationNumber << " iterations" << std::endl;
            outfile.open("/media/inigo/10740FB2740F9A1C/2d_results_test/plots/convergence_results.dat", std::ios_base::app);
            outfile << "            " << iteration_number;
            outfile.close();
        }

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        return is_converged;
    }
protected:

virtual void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh) override
    {
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        p_scheme->Update(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        // Move the mesh if needed
        if (MoveMesh == true)
            BaseType::MoveMesh();
    }



private:

    /**
     * @brief This method computes the norm of the residual
     * @details It checks if the dof is fixed
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rResidualSolutionNorm The norm of the residual
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rb RHS vector (residual + reactions)
     */
    double CalculateResidualNorm(
        DofsArrayType& rDofSet,
        const TSystemVectorType& rb
        )
    {
        // Initialize
        TDataType residual_solution_norm = TDataType();
        SizeType dof_num = 0;

        // Auxiliar values
        TDataType residual_dof_value = 0.0;
        const auto it_dof_begin = rDofSet.begin();
        const int number_of_dof = static_cast<int>(rDofSet.size());

        // Loop over Dofs
        #pragma omp parallel for firstprivate(residual_dof_value) reduction(+:residual_solution_norm, dof_num)
        for (int i = 0; i < number_of_dof; i++) {
            auto it_dof = it_dof_begin + i;

            if (!it_dof->IsFixed()) {
                const IndexType dof_id = it_dof->EquationId();
                residual_dof_value = TSparseSpace::GetValue(rb,dof_id);
                residual_solution_norm += std::pow(residual_dof_value, 2);
                dof_num++;
            }
        }

        return std::sqrt(residual_solution_norm) / (double)dof_num;
    }

    /**
     * @brief This method computes the final norm
     * @details It checks if the dof is fixed
     * @param rDofNum The number of DoFs
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param Dx Vector of results (variations on nodal variables)
     */
    TDataType CalculateFinalCorrectionNorm(
        DofsArrayType& rDofSet,
        const TSystemVectorType& Dx
        )
    {
        // Initialize
        TDataType final_correction_norm = TDataType();
        SizeType dof_num = 0;

        // Loop over Dofs
        #pragma omp parallel for reduction(+:final_correction_norm,dof_num)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
            auto it_dof = rDofSet.begin() + i;

            IndexType dof_id;
            TDataType variation_dof_value;

            if (it_dof->IsFree()) {
                dof_id = it_dof->EquationId();
                variation_dof_value = Dx[dof_id];
                final_correction_norm += std::pow(variation_dof_value, 2);
                dof_num++;
            }
        }

        return std::sqrt(final_correction_norm) / (double)dof_num;
    }

    double ComputeRelaxationFactor(const double& rAbsoluteResidualNorm)
    {
        if(rAbsoluteResidualNorm < 1e-4){
            return 1.0;
        }
        else{
            return 0.5;
        }
    }


}; /* Class UpwindResidualBasedNewtonRaphsonStrategy */
} /* namespace Kratos. */

#endif /* KRATOS_UPWIND_RESIDUALBASED_NEWTON_RAPHSON_STRATEGY defined */
