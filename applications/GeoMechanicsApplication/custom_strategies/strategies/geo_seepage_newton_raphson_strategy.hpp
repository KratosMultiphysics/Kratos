// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include <string>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "custom_utilities/seepage_boundary_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

// A Newton-Raphson strategy that can switch seepage boundary conditions between Dirichlet and
// Neumann while iterating.
//
// Applying a switch is the condition's job: the scheme calls
// Condition::InitializeNonLinearIteration on every non-linear iteration, and GeoSeepageCondition
// fixes or frees its nodes there. This strategy owns the two things the condition cannot do:
// deciding when to switch, and making sure the solver does not declare convergence on an iteration
// whose boundary configuration has just changed underneath it.
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoSeepageNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoSeepageNewtonRaphsonStrategy);

    using BaseType   = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;

    // The base class is a dependent template base, so every inherited data member used by the
    // copied SolveSolutionStep must be pulled into scope explicitly.
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mMaxIterationNumber;
    using MotherType::mNonconvergedSolutionsMatrix;
    using MotherType::mpA;  // Tangent matrix
    using MotherType::mpb;  // Residual vector of iteration i
    using MotherType::mpConvergenceCriteria;
    using MotherType::mpDx; // Delta x of iteration i
    using MotherType::mStoreNonconvergedSolutionsFlag;
    using MotherType::mUseOldStiffnessInFirstIteration;

    GeoSeepageNewtonRaphsonStrategy(ModelPart&                                 rModelPart,
                                    typename TSchemeType::Pointer              pScheme,
                                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                    typename TBuilderAndSolverType::Pointer    pNewBuilderAndSolver,
                                    int  MaxIterations          = 30,
                                    bool CalculateReactions     = false,
                                    bool ReformDofSetAtEachStep = false,
                                    bool MoveMeshFlag           = false)
        : MotherType(rModelPart,
                     pScheme,
                     pNewConvergenceCriteria,
                     pNewBuilderAndSolver,
                     MaxIterations,
                     CalculateReactions,
                     ReformDofSetAtEachStep,
                     MoveMeshFlag)
    {
    }

    [[nodiscard]] std::string Info() const override { return "GeoSeepageNewtonRaphsonStrategy"; }

    void Initialize() override
    {
        KRATOS_TRY

        MotherType::Initialize();
        mSeepageNodes = Geo::SeepageBoundaryUtilities::CollectSeepageNodes(BaseType::GetModelPart());
        for (auto* p_node : mSeepageNodes)
        {
            p_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
            p_node->Fix(WATER_PRESSURE);
        }
        KRATOS_INFO_IF("GeoSeepageNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
            << "Found " << mSeepageNodes.size() << " seepage nodes" << std::endl;

        KRATOS_CATCH("")
    }

    // NOTE: This is a deliberate copy of
    // ResidualBasedNewtonRaphsonStrategy::SolveSolutionStep
    // (kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h, lines
    // 919-1105). The ONLY intended differences are the blocks marked "SEEPAGE SEAM". Please keep it
    // that way: a diff against the core method should show nothing else. It is copied rather than
    // hooked because seam 3 has to run after PostCriteria, and no existing virtual method of the
    // base class sits at that point.
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart&                              r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer           p_scheme     = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        auto&                                   r_dof_set = p_builder_and_solver->GetDofSet();
        std::vector<Vector>                     NonconvergedSolutions;

        if (mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            this->GetCurrentSolution(r_dof_set, initial);
            NonconvergedSolutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        // initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number                      = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated                           = false;

        // ---- SEEPAGE SEAM 0 ------------------------------------------------------------------
        // Declared once here, and only re-assigned at the seams below. Declaring it inside both
        // the first-iteration block and the loop body would shadow one with the other.
        bool any_switched = false;
        // --------------------------------------------------------------------------------------

        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // ---- SEEPAGE SEAM 1: decide the switch, then force a stiffness rebuild ----------------
        // Deciding before the build means this iteration's BuildAndSolve reassembles with the new
        // fixity, so the block builder applies the switched node's Dirichlet/Neumann state at once.
        any_switched = UpdateSeepageBoundaryConditions();
        if (any_switched) this->SetStiffnessMatrixIsBuilt(false);
        // --------------------------------------------------------------------------------------

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (mUseOldStiffnessInFirstIteration) {
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(
                    p_scheme, r_model_part, rA, rDx, rb, BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx); // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        this->EchoInfo(iteration_number);

        // Updating the results stored in the database
        this->UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            Vector first;
            this->GetCurrentSolution(r_dof_set, first);
            NonconvergedSolutions.push_back(first);
        }

        if (is_converged) {
            if (mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        // ---- SEEPAGE SEAM 3 ------------------------------------------------------------------
        // A switch was just applied and solved for. The settled solution might itself warrant a
        // further switch, so convergence may only be declared on an iteration that switches
        // nothing. Force at least one more iteration here.
        if (any_switched) is_converged = false;
        // --------------------------------------------------------------------------------------

        // Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false && iteration_number++ < mMaxIterationNumber) {
            // setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            // ---- SEEPAGE SEAM 1: decide the switch, then force a stiffness rebuild ------------
            any_switched = UpdateSeepageBoundaryConditions();
            if (any_switched) this->SetStiffnessMatrixIsBuilt(false);
            // ----------------------------------------------------------------------------------

            // call the linear system solver to find the correction mDx for the
            // it is not called if there is no system to solve
            if (TSparseSpace::Size(rDx) != 0) {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false) {
                    if (this->GetKeepSystemConstantDuringIterations() == false) {
                        // A = 0.00;
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
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (mStoreNonconvergedSolutionsFlag == true) {
                Vector ith;
                this->GetCurrentSolution(r_dof_set, ith);
                NonconvergedSolutions.push_back(ith);
            }

            residual_is_updated = false;

            if (is_converged == true) {
                if (mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }

            // ---- SEEPAGE SEAM 3 --------------------------------------------------------------
            if (any_switched) is_converged = false;
            // ----------------------------------------------------------------------------------
        }

        // plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            this->MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("GeoSeepageNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / " << mMaxIterationNumber
                << " iterations" << std::endl;
        }

        // calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            mNonconvergedSolutionsMatrix = Matrix(r_dof_set.size(), NonconvergedSolutions.size());
            for (std::size_t i = 0; i < NonconvergedSolutions.size(); ++i) {
                block_for_each(r_dof_set, [&](const auto& r_dof) {
                    mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) =
                        NonconvergedSolutions[i](r_dof.EquationId());
                });
            }
        }

        return is_converged;
    }

protected:
    // Decides whether one seepage node must change between a Dirichlet and a zero-flux Neumann
    // boundary, applies that change, and reports whether anything changed. At most one node
    // switches per non-linear iteration, across the whole model.
    virtual bool UpdateSeepageBoundaryConditions()
    {
        KRATOS_TRY

        if (mSeepageNodes.empty()) return false;

        const auto nodal_flows = Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(
            BaseType::GetModelPart(), BaseType::GetModelPart().GetProcessInfo());

        return Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(mSeepageNodes, nodal_flows);

        KRATOS_CATCH("")
    }


private:
    // Cached once in Initialize. The conditions of a model part do not change during a solve, so
    // there is no need to rediscover the seepage nodes every iteration.
    std::vector<Node*> mSeepageNodes;
}; // Class GeoSeepageNewtonRaphsonStrategy

} // namespace Kratos





