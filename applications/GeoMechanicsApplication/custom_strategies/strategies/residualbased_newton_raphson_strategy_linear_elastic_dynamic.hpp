//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// default builder and solver
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_linear_elastic_dynamic.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"

#include "custom_processes/deactivate_conditions_on_inactive_elements_process.hpp"
#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"

// Application includes

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
 * @class GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic
 * @ingroup KratosGeomechanicsApplication
 * @brief This is Newton Raphson strategy especially for Dynamic linear elastic systems.
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Aron Noordam
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic);

    using SolvingStrategyType = SolvingStrategy<TSparseSpace, TDenseSpace>;
    using BaseType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using ClassType =
        GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TDataType                = typename BaseType::TDataType;
    using SparseSpaceType          = TSparseSpace;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType    = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType    = typename BaseType::LocalSystemMatrixType;
    using TSystemMatrixPointerType = typename BaseType::TSystemMatrixPointerType;
    using TSystemVectorPointerType = typename BaseType::TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic(ModelPart& rModelPart,
                                                                  typename TSchemeType::Pointer pScheme,
                                                                  typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                                                  typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                                                  int  MaxIterations      = 30,
                                                                  bool CalculateReactions = false,
                                                                  bool MoveMeshFlag       = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, false, MoveMeshFlag)
    {
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY
        BaseType::Initialize();

        // Note that FindNeighbourElementsOfConditionsProcess and DeactivateConditionsOnInactiveElements are required to be performed before initializing the System and State
        // this means that these operations are done twice in the GeomechanicsSolver in python
        FindNeighbourElementsOfConditionsProcess{ BaseType::GetModelPart() }.Execute();

        DeactivateConditionsOnInactiveElements{ BaseType::GetModelPart() }.Execute();

        if (!BaseType::mStiffnessMatrixIsBuilt)
            // initialize the system matrices and the initial second derivative
            this->InititalizeSystemAndState();

        KRATOS_CATCH("")
    }

    void Predict() override
    {
        KRATOS_TRY

        BaseType::GetScheme()->Predict(BaseType::GetModelPart(), BaseType::GetBuilderAndSolver()->GetDofSet(),
                                       *BaseType::mpA, *BaseType::mpDx, *BaseType::mpb);

        // Note that constraints are not applied in this predict, nor is an update performed, constrains are added in the builder and solver

        if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep() override

    {
        KRATOS_TRY
        BaseType::InitializeSolutionStep();

        // it is required to initialize the mDxTot vector here, as SolveSolutionStep can be called
        // multiple times, in a single time step (from an overlaying strategy)
        mDxTot = TSystemVectorType(BaseType::GetBuilderAndSolver()->GetDofSet().size(), 0.0);
        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart&                    r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme     = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        auto&               r_dof_set = p_builder_and_solver->GetDofSet();
        std::vector<Vector> non_converged_solutions;

        if (BaseType::mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            BaseType::GetCurrentSolution(r_dof_set, initial);
            non_converged_solutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        // initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

        bool is_converged = this->PerformIterationCycle(rA, rDx, rb, mDxTot, non_converged_solutions, iteration_number);

        if (is_converged) {
			// here only the derivatives are updated in the scheme that is used, generally both derivatives and solution step are updated
            p_scheme->Update(r_model_part, r_dof_set, rA, mDxTot, rb);
        }

        // plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            this->MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyLinearElasticDynamic", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << BaseType::mMaxIterationNumber << " iterations" << std::endl;
        }

        // calculate reactions if required
        if (BaseType::mCalculateReactionsFlag)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, mDxTot, rb);

        if (BaseType::mStoreNonconvergedSolutionsFlag) {
            BaseType::mNonconvergedSolutionsMatrix =
                Matrix(r_dof_set.size(), non_converged_solutions.size());
            for (std::size_t i = 0; i < non_converged_solutions.size(); ++i) {
                block_for_each(r_dof_set, [non_converged_solutions, i, this](const auto& r_dof) {
                    BaseType::mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) =
                        non_converged_solutions[i](r_dof.EquationId());
                });
            }
        }

        return is_converged;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        auto default_parameters = Parameters(R"(
        {
            "name"                                : "newton_raphson_strategy_linear_elastic_dynamic"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name() { return "newton_raphson_strategy_linear_elastic_dynamic"; }

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
        return "ResidualBasedNewtonRaphsonStrategyLinearElasticDynamic";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@}
    ///@name Static Member Variables
    ///@{
    TSystemVectorType mDxTot;

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method prints information after reach the max number of iterations
     */
    void MaxIterationsExceeded() override
    {
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyLinearElasticDynamic", BaseType::GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << BaseType::mMaxIterationNumber << " ) exceeded!"
            << std::endl;
    }

    bool PerformIterationCycle(TSystemMatrixType&   rA,
                               TSystemVectorType&   rDx,
                               TSystemVectorType&   rb,
                               TSystemVectorType&   rDxTot,
                               std::vector<Vector>& rNonconvergedSolutions,
                               unsigned int&        rIterationNumber)
    {
        ModelPart&                    r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme     = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        auto&       r_dof_set              = p_builder_and_solver->GetDofSet();
        const auto& r_current_process_info = r_model_part.GetProcessInfo();

        bool is_converged = false;
        for (; rIterationNumber < BaseType::mMaxIterationNumber; rIterationNumber++) {
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = rIterationNumber;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            if (SparseSpaceType::Size(rDx) != 0) {
                TSparseSpace::SetToZero(rDx);
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);

            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(rIterationNumber);

            // Updating the results stored in the database
            this->UpdateSolutionStepValue(rDx, rDxTot);

            // only finalize condition non linear iteration
			p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);

            BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (BaseType::mStoreNonconvergedSolutionsFlag) {
                Vector ith;
                BaseType::GetCurrentSolution(r_dof_set, ith);
                rNonconvergedSolutions.push_back(ith);
            }

            if (is_converged) {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                }

                is_converged =
                    BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
            if (is_converged) {
                return true;
            }
        }

        return false;
    }

    void UpdateSolutionStepValue(TSystemVectorType& rDx, TSystemVectorType& rDx_tot)
    {
        // performs: rDx_tot += rDx;
        TSparseSpace::UnaliasedAdd(rDx_tot, 1.0, rDx);

        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        DofsArrayType& r_dof_set = p_builder_and_solver->GetDofSet();

        block_for_each(r_dof_set, [&rDx](auto& dof) {
            if (dof.IsFree()) {
                dof.GetSolutionStepValue() += TSparseSpace::GetValue(rDx, dof.EquationId());
            }
        });
    }

    /// <summary>
    /// Initializes the system matrices and the initial state
    /// </summary>
    void InititalizeSystemAndState()
    {
        KRATOS_TRY

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        ModelPart&                    r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme     = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        this->InitializeSolutionStep();

        // Initialize non linear iteration for elements here, as the scheme only initializes conditions.
        const auto& r_current_process_info = r_model_part.GetProcessInfo();
        block_for_each(r_model_part.Elements(), [&r_current_process_info](Element& r_element) {
            if (r_element.IsActive()) {
                r_element.InitializeNonLinearIteration(r_current_process_info);
            }
        });

        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        p_builder_and_solver->Build(p_scheme, r_model_part, rA, rb);

        this->FinalizeSolutionStep();
        BaseType::mStiffnessMatrixIsBuilt = true;

        KRATOS_CATCH("")
    }

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

}; /* Class GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */
