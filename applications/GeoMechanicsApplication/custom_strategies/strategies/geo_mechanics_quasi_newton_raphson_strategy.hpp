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
                    "quasi_newton_raphson_restart_interval": 50,
                    "quasi_newton_raphson_max_rank" : 10

                }  )");

        // Validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // Select the quasi-Newton update scheme (low-rank secant approximation of the Jacobian/its inverse).
        const std::string quasi_newton_type = rParameters["quasi_newton_type"].GetString();
        if (quasi_newton_type == "lbfgs") {
            mQuasiNewtonType = QuasiNewtonType::LBFGS;
        } else if (quasi_newton_type == "broyden") {
            mQuasiNewtonType = QuasiNewtonType::Broyden;
        } else {
            KRATOS_ERROR << "Unknown 'quasi_newton_type': '" << quasi_newton_type
                         << "'. Valid options are 'broyden' and 'lbfgs'." << std::endl;
        }

        mRestartInterval = rParameters["quasi_newton_raphson_restart_interval"].GetInt();
        mMaxRank         = rParameters["quasi_newton_raphson_max_rank"].GetInt();

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

        TSparseSpace::SetToZero(rDx);
        TSparseSpace::SetToZero(rb);

        LBFGSRankStorage   lbfgs_rank_storage;
        BroydenRankStorage broyden_rank_storage;

        bool clear_storage = false;
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            clear_storage = true;
        } else {
            BuildReducedResidual(rb);
        }

        TSystemVectorType rb_new(rb.size());

        for (; iteration_number <= mMaxIterationNumber; ++iteration_number) {
            // check if we need to clear the low-rank storage and rebuild A_0
            if (iteration_number > 1 && (iteration_number - 1) % mRestartInterval == 0) {
                clear_storage = true;
            }

            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            if (iteration_number > 1) {
                p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
                mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
            }

            // Quasi-Newton solve: rDx = H_k * rb   (H_k approximates rA_k^{-1})
            if (mQuasiNewtonType == QuasiNewtonType::LBFGS) {
                if (clear_storage) {
                    BuildAndConstrainA0(rA, rDx, rb, lbfgs_rank_storage);
                    clear_storage = false;
                }
                this->LBfgsSolve(rA, rDx, rb, lbfgs_rank_storage);
            } else if (mQuasiNewtonType == QuasiNewtonType::Broyden) {
                if (clear_storage) {
                    BuildAndConstrainA0(rA, rDx, rb, broyden_rank_storage);
                    clear_storage = false;
                }
                this->ShermanMorrisSolve(rA, rDx, rb, broyden_rank_storage);
            } else {
                KRATOS_ERROR
                    << "Unknown 'quasi_newton_type'. Valid options are 'broyden' and 'lbfgs'." << std::endl;
            }

            UpdateDatabaseReduced(rA, rDx, rb, HasConstraints);
            BuildReducedResidual(rb_new);

            TSystemVectorType delta_b = rb - rb_new;

            if (mQuasiNewtonType == QuasiNewtonType::LBFGS) {
                this->UpdateLBFGSRank(lbfgs_rank_storage, rDx, delta_b);

            } else if (mQuasiNewtonType == QuasiNewtonType::Broyden) {
                this->UpdateBroydenRank(broyden_rank_storage, rA, rDx, delta_b);
            } else {
                KRATOS_ERROR
                    << "Unknown 'quasi_newton_type'. Valid options are 'broyden' and 'lbfgs'." << std::endl;
            }

            TSparseSpace::Copy(rb_new, rb);
            MotherType::EchoInfo(iteration_number);
            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            if (is_converged) {
                break;
            }
        }

        if (!is_converged) {
            MotherType::MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("GeoMechanicsQuasiNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / " << mMaxIterationNumber
                << std::endl;
        }

        if (mCalculateReactionsFlag) {
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
        }

        return is_converged;

        KRATOS_CATCH("")
    }

private:
    /// The available quasi-Newton update schemes.
    enum class QuasiNewtonType { Broyden, LBFGS };

    QuasiNewtonType mQuasiNewtonType = QuasiNewtonType::Broyden;
    unsigned int    mRestartInterval = 100; // rebuild K0 every N iterations to refresh curvature
    unsigned int    mMaxRank = 10; // maximum number of low-rank updates to store (Broyden or BFGS)
    typename TLinearSolver::Pointer mpLinearSolver;

    struct RankStorage {
        virtual void Clear()   = 0;
        virtual ~RankStorage() = default;
    };

    struct BroydenRankStorage : RankStorage {
        std::vector<TSystemVectorType> u_list; // u_i
        std::vector<TSystemVectorType> v_list; // v_i
        std::vector<TSystemVectorType> z_list; // z_i = rA0^{-1} u_i (cached)

        void Clear() override
        {
            u_list.clear();
            v_list.clear();
            z_list.clear();
        }
    };

    struct LBFGSRankStorage : RankStorage {
        std::vector<TSystemVectorType> dx_list;  // step differences
        std::vector<TSystemVectorType> db_list;  // gradient/residual differences
        std::vector<double>            rho_list; // curvature scalars

        void Clear() override
        {
            dx_list.clear();
            db_list.clear();
            rho_list.clear();
        }
    };

    void UpdateBroydenRank(BroydenRankStorage&      rBroydenRankStorage,
                           TSystemMatrixType&       rA_0,
                           const TSystemVectorType& rDx,
                           const TSystemVectorType& rDb)
    {
        // ---- Build the new Broyden rank-1 column  ----
        //   rA_k rDx = rA_0 rDx + U (V^T rDx)
        //   db_hat = rA_k rDx
        //   u_col = rDb - y_hat
        //   v_col = rDx / (rDx . rDx)
        //   z_col = rA_0^{-1} u_col

        TSystemVectorType db_hat(rDx.size());

        const double d_x_squared = TSparseSpace::Dot(rDx, rDx);
        if (d_x_squared > std::numeric_limits<double>::epsilon()) {
            const std::size_t k = rBroydenRankStorage.u_list.size();
            TSparseSpace::Mult(rA_0, rDx, db_hat);
            if (k > 0) {
                Vector v_list_dot_dx(k);
                for (std::size_t i = 0; i < k; ++i) {
                    v_list_dot_dx[i] = TSparseSpace::Dot(rBroydenRankStorage.v_list[i], rDx);
                }
                for (std::size_t i = 0; i < k; ++i) {
                    noalias(db_hat) += v_list_dot_dx[i] * rBroydenRankStorage.u_list[i];
                }
            }
            TSystemVectorType u_col = rDb - db_hat;
            TSystemVectorType v_col = rDx / d_x_squared;

            TSystemVectorType z_col(rDx.size());
            mpLinearSolver->PerformSolutionStep(rA_0, z_col, u_col);

            if (rBroydenRankStorage.u_list.size() >= mMaxRank) {
                rBroydenRankStorage.u_list.erase(rBroydenRankStorage.u_list.begin());
                rBroydenRankStorage.v_list.erase(rBroydenRankStorage.v_list.begin());
                rBroydenRankStorage.z_list.erase(rBroydenRankStorage.z_list.begin());
            }
            rBroydenRankStorage.u_list.push_back(u_col);
            rBroydenRankStorage.v_list.push_back(v_col);
            rBroydenRankStorage.z_list.push_back(z_col);
        }
    }

    void UpdateLBFGSRank(LBFGSRankStorage& rLBFGSRankStorage, const TSystemVectorType& rDx, const TSystemVectorType& rDb)
    {
        const double db_dot_dx = TSparseSpace::Dot(rDb, rDx);
        if (db_dot_dx > std::numeric_limits<double>::epsilon()) {
            if (rLBFGSRankStorage.dx_list.size() >= mMaxRank) {
                rLBFGSRankStorage.dx_list.erase(rLBFGSRankStorage.dx_list.begin());
                rLBFGSRankStorage.db_list.erase(rLBFGSRankStorage.db_list.begin());
                rLBFGSRankStorage.rho_list.erase(rLBFGSRankStorage.rho_list.begin());
            }
            rLBFGSRankStorage.dx_list.push_back(rDx);
            rLBFGSRankStorage.db_list.push_back(rDb);
            rLBFGSRankStorage.rho_list.push_back(1.0 / db_dot_dx);
        }
    }

    /// Builds A and b, applies master-slave constraints
    /// Dirichlet conditions, then factorizes the resulting (reduced) rA_0.
    void BuildAndConstrainA0(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rb, RankStorage& rRankStorage)
    {
        rRankStorage.Clear();

        auto       p_builder_and_solver = MotherType::GetBuilderAndSolver();
        auto       p_scheme             = MotherType::GetScheme();
        ModelPart& r_model_part         = BaseType::GetModelPart();

        TSparseSpace::SetToZero(rA_0);
        TSparseSpace::SetToZero(rb);

        // Build raw A and b together (ApplyConstraints needs the raw b to form T^T b consistently).
        p_builder_and_solver->Build(p_scheme, r_model_part, rA_0, rb);

        if (!r_model_part.MasterSlaveConstraints().empty()) {
            p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA_0, rb); // A <- T^T A T ; b <- T^T b ; builds mT
        }
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA_0, rDx, rb);

        mpLinearSolver->InitializeSolutionStep(rA_0, rDx, rb); // factorize rA_0

        BaseType::mStiffnessMatrixIsBuilt = true;
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

    void LBfgsSolve(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rb, const LBFGSRankStorage& rLBFGSRankStorage)
    {
        const std::size_t k = rLBFGSRankStorage.dx_list.size();
        Vector            alpha(k);

        TSystemVectorType q(rb.size());
        TSparseSpace::Copy(rb, q); // q = rb

        // First loop: newest to oldest
        for (auto i = k; i-- > 0;) {
            alpha[i] = rLBFGSRankStorage.rho_list[i] * TSparseSpace::Dot(rLBFGSRankStorage.dx_list[i], q);
            noalias(q) -= alpha[i] * rLBFGSRankStorage.db_list[i];
        }

        // Seed: rDx = H_0 * q = rA_0^{-1} q (reuses the cached factorization)
        TSparseSpace::SetToZero(rDx);
        mpLinearSolver->PerformSolutionStep(rA_0, rDx, q);

        // Second loop: oldest to newest
        for (std::size_t i = 0; i < k; ++i) {
            const double beta =
                rLBFGSRankStorage.rho_list[i] * TSparseSpace::Dot(rLBFGSRankStorage.db_list[i], rDx);
            noalias(rDx) += (alpha[i] - beta) * rLBFGSRankStorage.dx_list[i];
        }
    }

    /// Solves rDx = (rA_0 + U V^T)^{-1} rb (Sherman-Morrison-Woodbury).
    void ShermanMorrisSolve(TSystemMatrixType&        rA_0,
                            TSystemVectorType&        rDx,
                            TSystemVectorType&        rb,
                            const BroydenRankStorage& rBroydenRankStorage)
    {
        mpLinearSolver->PerformSolutionStep(rA_0, rDx, rb); // rDx = rA_0^{-1} rb

        // add the low rank correction, rDx = rDx - Z * (I + V^T Z)^{-1} (V^T rDx), where Z = rA_0^{-1} U
        const std::size_t k = rBroydenRankStorage.u_list.size();
        if (k > 0) {
            Vector v_list_dot_dx(k);
            for (std::size_t i = 0; i < k; ++i) {
                v_list_dot_dx[i] = TSparseSpace::Dot(rBroydenRankStorage.v_list[i], rDx);
            }
            // M = I + V^T Z, where Z = rA_0^{-1} U
            Matrix M(k, k);
            for (std::size_t i = 0; i < k; ++i) {
                for (std::size_t j = 0; j < k; ++j) {
                    M(i, j) = (i == j ? 1.0 : 0.0) + TSparseSpace::Dot(rBroydenRankStorage.v_list[i],
                                                                       rBroydenRankStorage.z_list[j]);
                }
            }

            // alpha = M^{-1} (V^T rDx)
            Vector alpha;
            SolveSmallDense(M, alpha, v_list_dot_dx);
            for (std::size_t j = 0; j < k; ++j) {
                noalias(rDx) -= alpha[j] * rBroydenRankStorage.z_list[j];
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
