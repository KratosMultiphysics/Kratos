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
#include <iostream>

// External includes

// Project includes

#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// #include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "utilities/builtin_timer.h"

// default builder and solver
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_linear_elastic_dynamic.h"
// #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"

#include "custom_processes/deactivate_conditions_on_inactive_elements_process.hpp"
#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"

// Application includes
#include "geo_mechanics_application_variables.h"

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
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Riccardo Rossi
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    using BaseType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    typedef GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic(
        ModelPart&                                 rModelPart,
        typename TSchemeType::Pointer              pScheme,
        typename TLinearSolver::Pointer            pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer    pNewBuilderAndSolver,
        Parameters&                                rParameters,
        double                                     Beta,
        double                                     Gamma,
        int                                        MaxIterations      = 30,
        bool                                       CalculateReactions = false,
        bool                                       MoveMeshFlag       = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              rModelPart,
              pScheme,
              /*pNewLinearSolver,*/
              pNewConvergenceCriteria,
              pNewBuilderAndSolver,
              MaxIterations,
              CalculateReactions,
              false,
              MoveMeshFlag),
          mBeta(Beta),
          mGamma(Gamma)
    {
        // new constructor
    }

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used
    in which the values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        const DataCommunicator& r_comm = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator();
        // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        // if the operations needed were already performed this does nothing
        if (BaseType::mInitializeWasPerformed == false) Initialize();

        const auto equation_size = BaseType::GetBuilderAndSolver()->GetDofSet().size();

        ModelPart& model_part = BaseType::GetModelPart();

        double delta_time = model_part.GetProcessInfo()[DELTA_TIME];

        TSystemVectorType& rFirstDerivativeVector  = TSystemVectorType(equation_size, 0.0);
        TSystemVectorType& rSecondDerivativeVector = TSystemVectorType(equation_size, 0.0);

        TSystemVectorType& rUpdatedFirstDerivativeVector  = TSystemVectorType(equation_size, 0.0);
        TSystemVectorType& rUpdatedSecondDerivativeVector = TSystemVectorType(equation_size, 0.0);

        this->GetFirstAndSecondDerivativeVector(rFirstDerivativeVector, rSecondDerivativeVector, model_part, 0);

        // rUpdatedFirstDerivativeVector = rFirstDerivativeVector * (mGamma / mBeta) + rSecondDerivativeVector * (delta_time * (mGamma / (2 * mBeta) - 1));
        TSparseSpace::ScaleAndAdd(mGamma / mBeta, rFirstDerivativeVector,
                                  delta_time * (mGamma / (2 * mBeta) - 1), rSecondDerivativeVector,
                                  rUpdatedFirstDerivativeVector);

        // rUpdatedSecondDerivativeVector = rFirstDerivativeVector * (1.0 / (mBeta * delta_time)) + rSecondDerivativeVector * (1.0 / (2 * mBeta));
        TSparseSpace::ScaleAndAdd(1.0 / (mBeta * delta_time), rFirstDerivativeVector, 1.0 / (2 * mBeta),
                                  rSecondDerivativeVector, rUpdatedSecondDerivativeVector);

        this->SetFirstAndSecondDerivativeVector(rUpdatedFirstDerivativeVector,
                                                rUpdatedSecondDerivativeVector, model_part);

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true) BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (BaseType::mInitializeWasPerformed == false) {
            // pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
            typename TConvergenceCriteriaType::Pointer p_convergence_criteria = BaseType::mpConvergenceCriteria;
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            // Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            // initialisation of the convergence criteria
            if (p_convergence_criteria->IsInitialized() == false)
                p_convergence_criteria->Initialize(BaseType::GetModelPart());

            BaseType::mInitializeWasPerformed = true;
        }

        auto find_neighbour_elements_of_conditions_process =
            FindNeighbourElementsOfConditionsProcess(BaseType::GetModelPart());
        find_neighbour_elements_of_conditions_process.Execute();

        auto deactivate_conditions_on_inactive_elements_process =
            DeactivateConditionsOnInactiveElements(BaseType::GetModelPart());
        deactivate_conditions_on_inactive_elements_process.Execute();

        if (!BaseType::mStiffnessMatrixIsBuilt)
            // initialize the system matrices and the initial second derivative
            this->InititalizeSystemAndState();

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        // Pointers needed in the solution
        typename TSchemeType::Pointer           p_scheme             = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        ModelPart&                              r_model_part         = BaseType::GetModelPart();

        // Set up the system, operation performed just once unless it is required
        // to reform the dof set at each iteration
        BuiltinTimer system_construction_time;
        if (!p_builder_and_solver->GetDofSetIsInitializedFlag()) {
            // Setting up the list of the DOFs to be solved
            BuiltinTimer setup_dofs_time;
            p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0)
                << "Setup Dofs Time: " << setup_dofs_time << std::endl;

            // Shaping correctly the system
            BuiltinTimer setup_system_time;
            p_builder_and_solver->SetUpSystem(r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0)
                << "Setup System Time: " << setup_system_time << std::endl;

            // Setting up the Vectors involved to the correct size
            BuiltinTimer system_matrix_resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb, r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0)
                << "System Matrix Resize Time: " << system_matrix_resize_time << std::endl;
        }

        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0)
            << "System Construction Time: " << system_construction_time << std::endl;

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        // Initial operations ... things that are constant over the Solution Step
        p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

        // only initialize solution step of conditions
        const auto& r_current_process_info = r_model_part.GetProcessInfo();
        block_for_each(r_model_part.Conditions(), [&r_current_process_info](Condition& r_condition) {
            if (r_condition.IsActive()) {
                r_condition.InitializeSolutionStep(r_current_process_info);
            }
        });

        // Initialisation of the convergence criteria
        if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        }

        BaseType::mpConvergenceCriteria->InitializeSolutionStep(
            r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
            TSparseSpace::SetToZero(rb);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = BaseType::GetModelPart();

        typename TSchemeType::Pointer           p_scheme             = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        // Finalisation of the solution step,
        // operations to be done after achieving convergence, for example the
        // Final Residual Vector (mb) has to be saved in there
        // to avoid error accumulation

        // just finialize conditions and not elements as defined in the scheme
        const auto& r_current_process_info = r_model_part.GetProcessInfo();

        block_for_each(r_model_part.Conditions(), [&r_current_process_info](Condition& r_condition) {
            if (r_condition.IsActive()) {
                r_condition.FinalizeSolutionStep(r_current_process_info);
            }
        });

        p_builder_and_solver->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->FinalizeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(),
                                                    rA, rDx, rb);

        // Cleaning memory after the solution
        p_scheme->Clean();

        if (BaseType::mReformDofSetAtEachStep == true) // deallocate the systemvectors
        {
            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart&                              r_model_part         = BaseType::GetModelPart();
        typename TSchemeType::Pointer           p_scheme             = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        auto&                                   r_dof_set = p_builder_and_solver->GetDofSet();
        std::vector<Vector>                     NonconvergedSolutions;

        if (mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            BaseType::GetCurrentSolution(r_dof_set, initial);
            NonconvergedSolutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        TSystemVectorType rDx_tot = TSystemVectorType(r_dof_set.size(), 0.0);

        // initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number                      = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated                           = false;

        // only initialize conditions
        const auto& r_current_process_info = r_model_part.GetProcessInfo();
        block_for_each(r_model_part.Conditions(), [&r_current_process_info](Condition& r_condition) {
            if (r_condition.IsActive()) {
                r_condition.InitializeNonLinearIteration(r_current_process_info);
            }
        });

        BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        bool is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        TSparseSpace::SetToZero(rDx);
        TSparseSpace::SetToZero(rb);

        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);

        // Debugging info
        BaseType::EchoInfo(iteration_number);

        // Updating the results stored in the database

        this->UpdateSolutionStepValue(rDx, rDx_tot);

        // only finalize condition non linear iteration
        block_for_each(r_model_part.Conditions(), [&r_current_process_info](Condition& r_condition) {
            if (r_condition.IsActive()) {
                r_condition.FinalizeNonLinearIteration(r_current_process_info);
            }
        });

        // p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (BaseType::mStoreNonconvergedSolutionsFlag) {
            Vector first;
            BaseType::GetCurrentSolution(r_dof_set, first);
            NonconvergedSolutions.push_back(first);
        }

        if (is_converged) {
            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        // Iteration Cycle... performed only for non linear RHS
        while (is_converged == false && iteration_number++ < BaseType::mMaxIterationNumber) {
            // setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            // call the linear system solver to find the correction mDx for the
            // it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0) {
                TSparseSpace::SetToZero(rDx);
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);

            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(iteration_number);

            // Updating the results stored in the database
            this->UpdateSolutionStepValue(rDx, rDx_tot);

            // only finalize condition non linear iteration
            block_for_each(r_model_part.Conditions(), [&r_current_process_info](Condition& r_condition) {
                if (r_condition.IsActive()) {
                    r_condition.FinalizeNonLinearIteration(r_current_process_info);
                }
            });

            // p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (BaseType::mStoreNonconvergedSolutionsFlag == true) {
                Vector ith;
                BaseType::GetCurrentSolution(r_dof_set, ith);
                NonconvergedSolutions.push_back(ith);
            }

            residual_is_updated = false;

            if (is_converged == true) {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        if (is_converged) {
            this->UpdateSolutionStepDerivative(rDx_tot, r_model_part);
        }

        // plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            this->MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / " << BaseType::mMaxIterationNumber
                << " iterations" << std::endl;
        }

        // calculate reactions if required
        if (BaseType::mCalculateReactionsFlag)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx_tot, rb);

        if (BaseType::mStoreNonconvergedSolutionsFlag) {
            BaseType::mNonconvergedSolutionsMatrix = Matrix(r_dof_set.size(), NonconvergedSolutions.size());
            for (std::size_t i = 0; i < NonconvergedSolutions.size(); ++i) {
                block_for_each(r_dof_set, [&](const auto& r_dof) {
                    BaseType::mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) =
                        NonconvergedSolutions[i](r_dof.EquationId());
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
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                : "newton_raphson_strategy_linear_elastic_dynamic",
            "use_old_stiffness_in_first_iteration": false,
            "max_iteration"                       : 10,
            "reform_dofs_at_each_step"            : false,
            "compute_reactions"                   : false,
            "builder_and_solver_settings"         : {},
            "convergence_criteria_settings"       : {},
            "linear_solver_settings"              : {},
            "scheme_settings"                     : {}
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

    double mBeta;
    double mGamma;

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method prints information after reach the max number of iterations
     */
    void MaxIterationsExceeded() override
    {
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyLinearElasticDynamic", BaseType::GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << BaseType::mMaxIterationNumber << " ) exceeded!" << std::endl;
    }

    void GetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
                                           TSystemVectorType& rSecondDerivativeVector,
                                           ModelPart&         rModelPart,
                                           IndexType          i)
    {
        block_for_each(rModelPart.Nodes(),
                       [&rFirstDerivativeVector, &rSecondDerivativeVector, i, this](Node& rNode) {
            if (rNode.IsActive()) {
                GetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector,
                                          rSecondDerivativeVector, i);
                GetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector,
                                          rSecondDerivativeVector, i);

                const std::vector<const Variable<double>*> optional_variables = {
                    &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

                for (const auto p_variable : optional_variables) {
                    GetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                                                      rSecondDerivativeVector, i);
                }
            }
        });
    }

    void SetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
                                           TSystemVectorType& rSecondDerivativeVector,
                                           ModelPart&         rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&rFirstDerivativeVector, &rSecondDerivativeVector, this](Node& rNode) {
            if (rNode.IsActive()) {
                SetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
                SetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector, rSecondDerivativeVector);

                const std::vector<const Variable<double>*> optional_variables = {
                    &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

                for (const auto p_variable : optional_variables) {
                    SetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                                                      rSecondDerivativeVector);
                }
            }
        });
    }

    void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                           const Node&             rNode,
                                           TSystemVectorType&      rFirstDerivativeVector,
                                           TSystemVectorType&      rSecondDerivativeVector,
                                           IndexType               i) const
    {
        if (rNode.HasDofFor(rVariable)) {
            GetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector, i);
        }
    }

    void SetDerivativesForOptionalVariable(const Variable<double>&  rVariable,
                                           Node&                    rNode,
                                           const TSystemVectorType& rFirstDerivativeVector,
                                           const TSystemVectorType& rSecondDerivativeVector)
    {
        if (rNode.HasDofFor(rVariable)) {
            SetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
        }
    }

    void GetDerivativesForVariable(const Variable<double>& rVariable,
                                   const Node&             rNode,
                                   TSystemVectorType&      rFirstDerivativeVector,
                                   TSystemVectorType&      rSecondDerivativeVector,
                                   IndexType               i) const
    {
        const auto& r_first_derivative  = rVariable.GetTimeDerivative();
        const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

        const auto equation_id              = rNode.GetDof(rVariable).EquationId();
        rFirstDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_first_derivative, i);
        rSecondDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_second_derivative, i);
    }

    void SetDerivativesForVariable(const Variable<double>&  rVariable,
                                   Node&                    rNode,
                                   const TSystemVectorType& rFirstDerivativeVector,
                                   const TSystemVectorType& rSecondDerivativeVector)
    {
        const auto& r_first_derivative  = rVariable.GetTimeDerivative();
        const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

        const auto equation_id                              = rNode.GetDof(rVariable).EquationId();
        rNode.FastGetSolutionStepValue(r_first_derivative)  = rFirstDerivativeVector[equation_id];
        rNode.FastGetSolutionStepValue(r_second_derivative) = rSecondDerivativeVector[equation_id];
    }

    void UpdateSolutionStepValue(TSystemVectorType& rDx, TSystemVectorType& rDx_tot)
    {
        // rDx_tot += rDx;
        TSparseSpace::UnaliasedAdd(rDx_tot, 1.0, rDx);

        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        DofsArrayType&                          r_dof_set = p_builder_and_solver->GetDofSet();

        block_for_each(r_dof_set, [&rDx](auto& dof) {
            if (dof.IsFree()) {
                dof.GetSolutionStepValue() += TSparseSpace::GetValue(rDx, dof.EquationId());
            }
        });
    }

    void UpdateSolutionStepDerivative(TSystemVectorType& rDx_tot, ModelPart& rModelPart)
    {
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        const DofsArrayType&                    r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemVectorType& r_first_derivative_vector  = TSystemVectorType(r_dof_set.size(), 0.0);
        TSystemVectorType& r_second_derivative_vector = TSystemVectorType(r_dof_set.size(), 0.0);

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        // get values from previous time step as the derivatives are already updated in the Predict step
        this->GetFirstAndSecondDerivativeVector(r_first_derivative_vector,
                                                r_second_derivative_vector, rModelPart, 1);

        // const TSystemVectorType& r_delta_first_derivative_vector = rDx_tot * (mGamma / (mBeta * delta_time)) - r_first_derivative_vector * (mGamma/mBeta) + r_second_derivative_vector * (delta_time * (1-mGamma / (2 * mBeta)));
        TSystemVectorType& r_delta_first_derivative_vector = TSystemVectorType(r_dof_set.size(), 0.0);
        TSparseSpace::UnaliasedAdd(r_delta_first_derivative_vector, (mGamma / (mBeta * delta_time)), rDx_tot);
        TSparseSpace::UnaliasedAdd(r_delta_first_derivative_vector, -(mGamma / mBeta), r_first_derivative_vector);
        TSparseSpace::UnaliasedAdd(r_delta_first_derivative_vector,
                                   delta_time * (1 - mGamma / (2 * mBeta)), r_second_derivative_vector);

        // const TSystemVectorType& r_delta_second_derivative_vector = rDx_tot * (1 / (mBeta * delta_time * delta_time)) - r_first_derivative_vector * (1 / (mBeta * delta_time)) - r_second_derivative_vector * (1 / (2 * mBeta));
        TSystemVectorType& r_delta_second_derivative_vector = TSystemVectorType(r_dof_set.size(), 0.0);
        TSparseSpace::UnaliasedAdd(r_delta_second_derivative_vector,
                                   1 / (mBeta * delta_time * delta_time), rDx_tot);
        TSparseSpace::UnaliasedAdd(r_delta_second_derivative_vector, -1 / (mBeta * delta_time),
                                   r_first_derivative_vector);
        TSparseSpace::UnaliasedAdd(r_delta_second_derivative_vector, -1 / (2 * mBeta), r_second_derivative_vector);

        // r_first_derivative_vector += r_delta_first_derivative_vector;
        TSparseSpace::UnaliasedAdd(r_first_derivative_vector, 1.0, r_delta_first_derivative_vector);

        // r_second_derivative_vector += r_delta_second_derivative_vector;
        TSparseSpace::UnaliasedAdd(r_second_derivative_vector, 1.0, r_delta_second_derivative_vector);

        this->SetFirstAndSecondDerivativeVector(r_first_derivative_vector, r_second_derivative_vector, rModelPart);
    }

    /// <summary>
    /// Initializes the system matrices and the initial state
    /// </summary>
    void InititalizeSystemAndState()
    {
        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        ModelPart&                              r_model_part         = BaseType::GetModelPart();
        typename TSchemeType::Pointer           p_scheme             = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        this->InitializeSolutionStep();
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        p_builder_and_solver->Build(p_scheme, r_model_part, rA, rb);
        this->FinalizeSolutionStep();
        BaseType::mStiffnessMatrixIsBuilt = true;
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

    /**
     * Copy constructor.
     */

    GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic(const GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic& Other) {};

    ///@}

}; /* Class GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */
