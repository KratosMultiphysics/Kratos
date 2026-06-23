// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// Project includes
#include "factories/linear_solver_factory.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsQuasiNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsQuasiNewtonRaphsonStrategy);

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpace, TDenseSpace>;
    using BaseType   = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mInitializeWasPerformed;
    using MotherType::mMaxIterationNumber;
    using MotherType::mpA; // Tangent matrix
    using MotherType::mpb; // Residual vector of iteration i
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpConvergenceCriteria;
    using MotherType::mpDx; // Delta x of iteration i
    using MotherType::mpScheme;

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    GeoMechanicsQuasiNewtonRaphsonStrategy(ModelPart&                    model_part,
                                           typename TSchemeType::Pointer pScheme,
                                           typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                           typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                           Parameters& rParameters,
                                           int         MaxIterations          = 30,
                                           bool        CalculateReactions     = false,
                                           bool        ReformDofSetAtEachStep = false,
                                           bool        MoveMeshFlag           = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        // only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
                {
                    "min_iteration":    2,
                    "number_cycles":    5,
                    "increase_factor":  2.0,
                    "reduction_factor": 0.5,
                    "max_piping_iterations": 50,
                    "desired_iterations": 4,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "rebuild_level": 2,
                    "quasi_newton_type": "broyden",
                    "quasi_newton_raphson_restart_interval": 5
                }  )");

        // Validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // Select the quasi-Newton update scheme (low-rank secant approximation of the Jacobian/its inverse).
        const std::string quasi_newton_type = rParameters["quasi_newton_type"].GetString();
        if (quasi_newton_type == "bfgs") {
            mQuasiNewtonType = QuasiNewtonType::BFGS;
        } else if (quasi_newton_type == "broyden") {
            mQuasiNewtonType = QuasiNewtonType::Broyden;
        } else {
            KRATOS_ERROR << "Unknown 'quasi_newton_type': '" << quasi_newton_type
                         << "'. Valid options are 'broyden' and 'bfgs'." << std::endl;
        }


		mRestartInterval = rParameters["quasi_newton_raphson_restart_interval"].GetInt();

        // it is required to use a direct solver without scaling for the quasi-newton strategy, as the strategy depends on reusing the factorization of the initial stiffness matrix.
        auto linear_solver_settings = Parameters(R"(
            {
                "solver_type": "LinearSolversApplication.sparse_lu",
                "scaling": false
            }  )");

        mpLinearSolver = LinearSolverFactoryType().Create(linear_solver_settings);
        
    }

    bool SolveSolutionStep() override
    {
        KRATOS_TRY

        ModelPart&                    r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme     = MotherType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = MotherType::GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        // Are there master-slave constraints active? If so, we must operate in the
        // reduced (constrained) space: A_hat = T^T A T, b_hat = T^T b, and recover dx = T dx_hat.
        const bool HasConstraints = !r_model_part.MasterSlaveConstraints().empty();

        unsigned int iteration_number                      = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // ---- Build K0 ONCE and LU-factor it ONCE ----
        TSparseSpace::SetToZero(rDx);
        TSparseSpace::SetToZero(rb);

        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            BuildAndConstrainK0(rA, rDx, rb); // builds A & b, applies constraints + Dirichlet, factorizes
            BaseType::mStiffnessMatrixIsBuilt = true;
        } else {
            BuildReducedResidual(rb);
        }

        // ---- Broyden low-rank storage (all vectors live in the REDUCED space): K_k = K0 + sum_i U_i V_i^T ----
        std::vector<TSystemVectorType> numerator_list;   // u_i
        std::vector<TSystemVectorType> denominator_list; // v_i
        std::vector<TSystemVectorType> K0inv_Ui_list;    // z_i = K0^{-1} u_i (cached)

        // ---- BFGS low-rank storage (all vectors live in the REDUCED space) ----
        std::vector<TSystemVectorType> bfgs_s_list;   // delta_u_i (steps)
        std::vector<TSystemVectorType> bfgs_y_list;   // delta_g_i (residual/gradient changes)
        std::vector<double>            bfgs_rho_list; // 1 / (delta_g_i . delta_u_i)

        const std::size_t n = rb.size();
        TSystemVectorType rb_new(n);
        TSystemVectorType K0s(n);
        TSystemVectorType u_col(n);
        TSystemVectorType z_col(n);


        for (; iteration_number <= mMaxIterationNumber; ++iteration_number) {
            if (iteration_number > 1 && (iteration_number - 1) % mRestartInterval == 0) {
                // rebuild & refactorize K0 (reduced), clear Broyden memory
                BuildAndConstrainK0(rA, rDx, rb);
                numerator_list.clear();
                denominator_list.clear();
                K0inv_Ui_list.clear();
                bfgs_s_list.clear();
                bfgs_y_list.clear();
                bfgs_rho_list.clear();
            }

            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            if (iteration_number > 1) {
                p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
                mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            if (is_converged) {
                break;
            }

            // Quasi-Newton solve in reduced space: rDx_hat = H_k * rb_hat   (H_k approximates K_k^{-1})
            if (mQuasiNewtonType == QuasiNewtonType::BFGS) {
                this->BfgsSolve(rA, rDx, rb, bfgs_s_list, bfgs_y_list, bfgs_rho_list);
            } else {
                this->ShermanMorrisSolve(rA, rDx, rb, numerator_list, denominator_list, K0inv_Ui_list);
            }

            bool do_line_search = false;
            if (do_line_search) {
                LineSearch(rA, rDx, rb, rb_new, n, HasConstraints);
            } else {
                // Apply the reduced increment to the database (mapped to full space if MPC), then rebuild residual.
                UpdateDatabaseReduced(rA, rDx, rb, HasConstraints);
                BuildReducedResidual(rb_new);
            }

            // ---- Store the new secant pair / low-rank update (reduced space) ----
            if (mQuasiNewtonType == QuasiNewtonType::BFGS) {
                // BFGS curvature pair:  s = delta_u = rDx,  y = delta_g = rb - rb_new
                // The pair is only stored if the curvature condition y . s > 0 holds, which
                // guarantees the (implicit) inverse Hessian approximation stays positive definite.
                TSystemVectorType y_col(n);
                noalias(y_col) = rb - rb_new;
                const double y_dot_s = TSparseSpace::Dot(y_col, rDx);
                if (y_dot_s > std::numeric_limits<double>::epsilon()) {

                    const std::size_t max_memory = 5;   // e.g. 10
                    if (bfgs_s_list.size() >= max_memory) {
                        bfgs_s_list.erase(bfgs_s_list.begin());
                        bfgs_y_list.erase(bfgs_y_list.begin());
                        bfgs_rho_list.erase(bfgs_rho_list.begin());
                    }
                    //numerator_list.push_back(u_col);
                    //denominator_list.push_back(v_col);
                    //K0inv_Ui_list.push_back(z_col);


                    bfgs_s_list.push_back(rDx);
                    bfgs_y_list.push_back(y_col);
                    bfgs_rho_list.push_back(1.0 / y_dot_s);
                }
            } else {
                // ---- Build the new Broyden rank-1 column (reduced space) ----
                //   s     = rDx (reduced step actually taken)
                //   K_k s = K0 s + U (V^T s)
                //   y     = rb - rb_new
                //   u_col = y - K_k s
                //   v_col = s / (s . s)
                //   z_col = K0^{-1} u_col
                const double s_dot_s = TSparseSpace::Dot(rDx, rDx);
                if (s_dot_s > std::numeric_limits<double>::epsilon()) {
                    const std::size_t k = numerator_list.size();
                    TSparseSpace::Mult(rA, rDx, K0s); // K0 * s   (reduced A)
                    if (k > 0) {
                        Vector VtS(k);
                        for (std::size_t i = 0; i < k; ++i) {
                            VtS[i] = TSparseSpace::Dot(denominator_list[i], rDx);
                        }
                        for (std::size_t i = 0; i < k; ++i) {
                            noalias(K0s) += VtS[i] * numerator_list[i];
                        }
                    }

                    TSystemVectorType y_col(n);
                    noalias(y_col) = rb - rb_new;

                    noalias(u_col) = y_col - K0s;

                    TSystemVectorType v_col(n);
                    noalias(v_col) = rDx / s_dot_s;

                    TSparseSpace::SetToZero(z_col);
                    TSystemVectorType u_solve(n);
                    TSparseSpace::Copy(u_col, u_solve);
                    mpLinearSolver->PerformSolutionStep(rA, z_col, u_solve);


                    const std::size_t max_memory = 5;   // e.g. 10
                    if (numerator_list.size() >= max_memory) {
                        numerator_list.erase(numerator_list.begin());
                        denominator_list.erase(denominator_list.begin());
                        K0inv_Ui_list.erase(K0inv_Ui_list.begin());
                    }
                    numerator_list.push_back(u_col);
                    denominator_list.push_back(v_col);
                    K0inv_Ui_list.push_back(z_col);


                    //numerator_list.push_back(u_col);
                    //denominator_list.push_back(v_col);
                    //K0inv_Ui_list.push_back(z_col);
                }
            }

            // Advance: R := R_new
            TSparseSpace::Copy(rb_new, rb);

            // Per-iteration housekeeping
            MotherType::EchoInfo(iteration_number);
            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        }

        if (!is_converged) {
            MotherType::MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("GeoMechanicsQuasiNewtonRaphsonStrategy[Broyden]", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / " << mMaxIterationNumber
                << " Broyden iterations (rank-" << numerator_list.size() << " update)" << std::endl;
        }

        if (mCalculateReactionsFlag) {
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
        }

        return is_converged;

        KRATOS_CATCH("")
    }

private:
    /// The available quasi-Newton update schemes.
    enum class QuasiNewtonType { Broyden, BFGS };

    QuasiNewtonType                mQuasiNewtonType = QuasiNewtonType::Broyden;
	unsigned int 			  mRestartInterval = 10; // rebuild K0 every N iterations to refresh curvature
    typename TLinearSolver::Pointer mpLinearSolver;

    /// Builds A and b, applies master-slave constraints (A <- T^T A T, b <- T^T b) and
    /// Dirichlet conditions, then factorizes the resulting (reduced) K0. Leaves rA holding
    /// the reduced operator used throughout the Broyden iteration.
    void BuildAndConstrainK0(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb)
    {
        auto       p_builder_and_solver = MotherType::GetBuilderAndSolver();
        auto       p_scheme             = MotherType::GetScheme();
        ModelPart& r_model_part         = BaseType::GetModelPart();

        TSparseSpace::SetToZero(rA);
        TSparseSpace::SetToZero(rb);

        // Build raw A and b together (ApplyConstraints needs the raw b to form T^T b consistently).
        p_builder_and_solver->Build(p_scheme, r_model_part, rA, rb);

        if (!r_model_part.MasterSlaveConstraints().empty()) {
            p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA, rb); // A <- T^T A T ; b <- T^T b ; builds mT
        }
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

        mpLinearSolver->InitializeSolutionStep(rA, rDx, rb); // factorize the reduced K0
    }

    /// Builds the residual and reduces it to the constrained space: b_hat = T^T b (slaves zeroed).
    void BuildReducedResidual(TSystemVectorType& rb)
    {
        auto       p_builder_and_solver = MotherType::GetBuilderAndSolver();
        auto       p_scheme             = MotherType::GetScheme();
        ModelPart& r_model_part         = BaseType::GetModelPart();

        TSparseSpace::SetToZero(rb);

        // note that BuildRHS also applies the Dirichlet conditions on the RHS
        p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        if (!r_model_part.MasterSlaveConstraints().empty()) {
            p_builder_and_solver->ApplyRHSConstraints(p_scheme, r_model_part, rb);
        }
    }

    void UpdateDatabaseReduced(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, bool HasConstraints)
    {
        auto p_builder_and_solver = MotherType::GetBuilderAndSolver();

        if (HasConstraints) {
            auto&             rT = p_builder_and_solver->GetConstraintRelationMatrix();
            TSystemVectorType dx_full(rDx.size());
            TSparseSpace::Mult(rT, rDx, dx_full); // dx_full = T * dx_reduced
            MotherType::UpdateDatabase(rA, dx_full, rb, BaseType::MoveMeshFlag());
        } else {
            MotherType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());
        }
    }

    void LineSearch(TSystemMatrixType& rA,
                    TSystemVectorType& rDx,
                    TSystemVectorType& rb,
                    TSystemVectorType& rb_new,
                    unsigned int       n,
                    bool               HasConstraints)
    {
        // rDx (reduced) is the proposed full Broyden step direction. Find alpha that reduces ||R||.
        TSystemVectorType s_dir(n);
        TSparseSpace::Copy(rDx, s_dir); // keep the reduced direction

        const double r0 = TSparseSpace::TwoNorm(rb); // ||R|| at alpha = 0

        // alpha = 1/2 (incremental: apply +0.5 s)
        noalias(rDx) = 0.5 * s_dir;
        UpdateDatabaseReduced(rA, rDx, rb, HasConstraints);
        BuildReducedResidual(rb_new);
        const double rh = TSparseSpace::TwoNorm(rb_new);

        // alpha = 1 (apply the other +0.5 s)
        noalias(rDx) = 0.5 * s_dir;
        UpdateDatabaseReduced(rA, rDx, rb, HasConstraints);
        BuildReducedResidual(rb_new);
        const double rf = TSparseSpace::TwoNorm(rb_new);

        // parabola through (0, r0), (1/2, rh), (1, rf) -> vertex
        double       alpha = 1.0;
        const double denom = (2.0 * r0 - 4.0 * rh + 2.0 * rf);
        if (std::abs(denom) > std::numeric_limits<double>::epsilon()) {
            alpha = (3.0 * r0 - 4.0 * rh + rf) / denom; // = -b/2a
        }
        alpha = std::clamp(alpha, 0.1, 1.0);

        if (alpha != 1.0) {
            // currently at alpha=1; move to alpha
            noalias(rDx) = (alpha - 1.0) * s_dir;
            UpdateDatabaseReduced(rA, rDx, rb, HasConstraints);
            BuildReducedResidual(rb_new);
        }
        // the ACTUAL reduced step taken is alpha*s_dir -> used for the Broyden update
        noalias(rDx) = alpha * s_dir;
    }

    /// Solves rDx = H_k * rb using the L-BFGS two-loop recursion, with the seed
    /// inverse Hessian H_0 = K0^{-1} applied through the cached factorization of the
    /// (reduced) initial stiffness matrix. The stored pairs satisfy the secant condition
    /// H_{i+1} * delta_g_i = delta_u_i.
    void BfgsSolve(TSystemMatrixType&              rA,
                   TSystemVectorType&              rDx,
                   TSystemVectorType&              rb,
                   std::vector<TSystemVectorType>& rSList,
                   std::vector<TSystemVectorType>& rYList,
                   std::vector<double>&            rRhoList)
    {
        const std::size_t k = rSList.size();
        Vector            alpha(k);

        TSystemVectorType q(rb.size());
        TSparseSpace::Copy(rb, q); // q = rb

        // First loop: newest to oldest
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(k) - 1; i >= 0; --i) {
            alpha[i] = rRhoList[i] * TSparseSpace::Dot(rSList[i], q);
            noalias(q) -= alpha[i] * rYList[i];
        }

        // Seed: rDx = H_0 * q = K0^{-1} q (reuses the cached factorization)
        TSparseSpace::SetToZero(rDx);
        mpLinearSolver->PerformSolutionStep(rA, rDx, q);

        // Second loop: oldest to newest
        for (std::size_t i = 0; i < k; ++i) {
            const double beta = rRhoList[i] * TSparseSpace::Dot(rYList[i], rDx);
            noalias(rDx) += (alpha[i] - beta) * rSList[i];
        }
    }

    /// Solves rDx = (K0 + U V^T)^{-1} rb (Sherman-Morrison-Woodbury).
    void ShermanMorrisSolve(TSystemMatrixType&              rA,
                            TSystemVectorType&              rDx,
                            TSystemVectorType&              rb,
                            std::vector<TSystemVectorType>& rUVectorList,
                            std::vector<TSystemVectorType>& rVVectorList,
                            std::vector<TSystemVectorType>& rZVectorList)
    {
        mpLinearSolver->PerformSolutionStep(rA, rDx, rb); // rDx = K0^{-1} rb

        const std::size_t k = rUVectorList.size();
        if (k > 0) {
            Vector VtW(k);
            for (std::size_t i = 0; i < k; ++i) {
                VtW[i] = TSparseSpace::Dot(rVVectorList[i], rDx);
            }
            Matrix M(k, k);
            for (std::size_t i = 0; i < k; ++i) {
                for (std::size_t j = 0; j < k; ++j) {
                    M(i, j) = (i == j ? 1.0 : 0.0) + TSparseSpace::Dot(rVVectorList[i], rZVectorList[j]);
                }
            }
            Vector alpha;
            SolveSmallDense(M, alpha, VtW);
            for (std::size_t j = 0; j < k; ++j) {
                noalias(rDx) -= alpha[j] * rZVectorList[j];
            }
        }
    }

    /**
     * @brief Solve a small dense k x k linear system M * x = b using
     *        Gaussian elimination with partial pivoting. Returns false if M is singular.
     *        Note: rM and rRhs are taken by value (mutated locally).
     */
    bool SolveSmallDense(Matrix rM, Vector& rSolution, Vector rRhs)
    {
        const std::size_t k = rRhs.size();
        rSolution.resize(k, false);
        if (k == 0) return true;

        for (std::size_t i = 0; i < k; ++i) {
            std::size_t pivot     = i;
            double      pivot_val = std::abs(rM(i, i));
            for (std::size_t r = i + 1; r < k; ++r) {
                const double v = std::abs(rM(r, i));
                if (v > pivot_val) {
                    pivot_val = v;
                    pivot     = r;
                }
            }
            if (pivot_val < std::numeric_limits<double>::epsilon()) {
                return false; // singular
            }
            if (pivot != i) {
                for (std::size_t c = i; c < k; ++c)
                    std::swap(rM(i, c), rM(pivot, c));
                std::swap(rRhs[i], rRhs[pivot]);
            }
            const double diag = rM(i, i);
            for (std::size_t r = i + 1; r < k; ++r) {
                const double factor = rM(r, i) / diag;
                if (factor == 0.0) continue;
                for (std::size_t c = i; c < k; ++c) {
                    rM(r, c) -= factor * rM(i, c);
                }
                rRhs[r] -= factor * rRhs[i];
            }
        }
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(k) - 1; i >= 0; --i) {
            double s = rRhs[i];
            for (std::size_t c = i + 1; c < k; ++c) {
                s -= rM(i, c) * rSolution[c];
            }
            rSolution[i] = s / rM(i, i);
        }
        return true;
    }

}; // Class GeoMechanicsQuasiNewtonRaphsonStrategy

} // namespace Kratos
