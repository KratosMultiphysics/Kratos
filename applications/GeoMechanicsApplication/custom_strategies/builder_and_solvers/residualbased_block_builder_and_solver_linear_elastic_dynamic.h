// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         BSD License
//                   geo_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi, Aron Noordam
//  Collaborators:   Vicente Mataix
//
//

#pragma once

/* System includes */
#include <unordered_set>

/* External includes */
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "custom_utilities/sparse_system_utilities.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "includes/model_part.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/atomic_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/timer.h"
#include "utilities/variable_utils.h"

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
 * @class ResidualBasedBlockBuilderAndSolverLinearElasticDynamic
 * @ingroup GeoMechanicsApplication
 * @brief Current class provides an implementation for builder and solving operations, especially
 * for linear elastic dynamic systems
 * @details When the LHS is build, the global mass and damping matrices are built separately. When
 * building the RHS, the mass and damping matrices are multiplied with respectively the second and
 * first derivative vector to calculate the mass and damping contribution. The RHS is constituted by
 * the unbalanced loads (residual) Degrees of freedom are reordered putting the restrained degrees
 * of freedom at the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information. Calculation of the reactions involves a cost very similar to the calculation of
 * the total residual. This class is intended to be used when the stiffness, mass and damping matrices are
 * constant throughout the iterations. If the matrices are not constant, this class cannot be used.
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author Aron Noordam
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ResidualBasedBlockBuilderAndSolverLinearElasticDynamic
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG(SILENT_WARNINGS);

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverLinearElasticDynamic);

    /// Definition of the base class
    using BaseType = ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definition of the classes from the base class
    using TSchemeType           = typename BaseType::TSchemeType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using NodesArrayType        = typename BaseType::NodesArrayType;
    using ElementsArrayType     = typename BaseType::ElementsArrayType;
    using ConditionsArrayType   = typename BaseType::ConditionsArrayType;

    /// Additional definitions
    using ElementsContainerType = PointerVectorSet<Element, IndexedObject>;
    using EquationIdVectorType  = Element::EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverLinearElasticDynamic(typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                                                    double Beta,
                                                                    double Gamma,
                                                                    bool CalculateInitialSecondDerivative)
        : BaseType(pNewLinearSystemSolver),
          mBeta(Beta),
          mGamma(Gamma),
          mCalculateInitialSecondDerivative(CalculateInitialSecondDerivative)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverLinearElasticDynamic() override = default;

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        mPreviousOutOfBalanceVector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        mCurrentOutOfBalanceVector  = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);

        if (mPreviousExternalForceVector.size() == 0) {
            mPreviousExternalForceVector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        }
    }

    void Build(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rb) override
    {
        Timer::Start("Build");
        this->BuildLHS(pScheme, rModelPart, rA);
        this->BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        Timer::Stop("Build");

        TSystemVectorType dummy_b(rA.size1(), 0.0);
        TSystemVectorType dummy_rDx(rA.size1(), 0.0);

        // apply constraints
        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            const auto timer_constraints = BuiltinTimer();
            Timer::Start("ApplyConstraints");
            BaseType::ApplyConstraints(pScheme, rModelPart, rA, rb);
            BaseType::ApplyConstraints(pScheme, rModelPart, mMassMatrix, dummy_b);
            BaseType::ApplyConstraints(pScheme, rModelPart, mDampingMatrix, dummy_b);
            Timer::Stop("ApplyConstraints");
            KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 1)
                << "Constraints build time: " << timer_constraints << std::endl;
        }

        // apply dirichlet conditions
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, dummy_rDx, rb);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mMassMatrix, dummy_rDx, dummy_b);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mDampingMatrix, dummy_rDx, dummy_b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << dummy_rDx
            << "\nRHS vector = " << rb << std::endl;

        if (mCalculateInitialSecondDerivative) {
            this->CalculateInitialSecondDerivative(rModelPart, rA, rb, pScheme);
            mCopyExternalForceVector = true;
        }

        // only add dynamics to lhs after calculating initial second derivative
        this->AddDynamicsToLhs(rA, rModelPart);

        // Initialize the linear solver, such that the solver can factorize the matrices already. In
        // case the matrices can be pre-factorized, this step is not performed during calculation.
        BaseType::mpLinearSystemSolver->InitializeSolutionStep(rA, dummy_rDx, rb);

        // check if PerformSolutionStep can be performed instead of Solve, this is more efficient for solvers, which can be pre-factorized
        try {
            BaseType::mpLinearSystemSolver->PerformSolutionStep(rA, dummy_rDx, rb);
            mUsePerformSolutionStep = true;
        } catch (const Kratos::Exception& e) {
            std::string error_message = e.what();

            // if PerformSolutionStep is not implemented, the following error is thrown, in this case, use Solve
            if (error_message.find("Error: Calling linear solver base class") != std::string::npos) {
                mUsePerformSolutionStep = false;
            }
            // Re-throw the exception if it's not the specific error we're looking for
            else {
                throw;
            }
        }
    }

    /**
     * @brief Function to perform the build of the LHS, mass matrix and damping matrix
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(pScheme) << "No scheme provided!" << std::endl;

        this->InitializeDynamicMatrix(mMassMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);
        this->InitializeDynamicMatrix(mDampingMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);

        const auto timer = BuiltinTimer();

        const ElementsArrayType& r_elements = rModelPart.Elements();
        this->CalculateGlobalMatrices(r_elements, rA, rModelPart);

        const ConditionsArrayType& r_conditions = rModelPart.Conditions();
        this->CalculateGlobalMatrices(r_conditions, rA, rModelPart);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 1)
            << "Build time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic",
                       (BaseType::GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0))
            << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the RHS and solves the system with an already defined LHS
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(typename TSchemeType::Pointer pScheme,
                          ModelPart&                    rModelPart,
                          TSystemMatrixType&            rA,
                          TSystemVectorType&            rDx,
                          TSystemVectorType&            rb) override
    {
        this->BuildRHS(pScheme, rModelPart, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            TSystemVectorType Dxmodified(rb.size());

            // Initialize the vector
            TSparseSpace::SetToZero(Dxmodified);

            this->InternalSystemSolveWithPhysics(rA, Dxmodified, rb, rModelPart);

            // recover solution of the original problem
            TSparseSpace::Mult(BaseType::mT, Dxmodified, rDx);
        } else {
            this->InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }

        TSparseSpace::Copy(mCurrentOutOfBalanceVector, mPreviousOutOfBalanceVector);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 1)
            << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() >= 3)
            << "After the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void InternalSystemSolveWithPhysics(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, ModelPart& rModelPart)
    {
        double norm_b = 0.00;
        if (TSparseSpace::Size(rb) != 0) norm_b = TSparseSpace::TwoNorm(rb);

        if (norm_b != 0.00) {
            // if the system is already factorized, perform solution step. In case the solver does not support this, use Solve
            if (mUsePerformSolutionStep) {
                BaseType::mpLinearSystemSolver->PerformSolutionStep(rA, rDx, rb);
            } else {
                BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
            }
        } else {
            KRATOS_WARNING_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic",
                              BaseType::mOptions.IsNot(SILENT_WARNINGS))
                << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // Prints information about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", BaseType::GetEchoLevel() > 1)
            << *(BaseType::mpLinearSystemSolver) << std::endl;
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb) override
    {
        KRATOS_TRY

        Timer::Start("BuildRHS");

        this->BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        // add dirichlet conditions to RHS
        this->ApplyDirichletConditionsRhs(rb);

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }

    void CalculateReactions(typename TSchemeType::Pointer pScheme,
                            ModelPart&                    rModelPart,
                            TSystemMatrixType&            A,
                            TSystemVectorType&            Dx,
                            TSystemVectorType&            b) override
    {
        TSparseSpace::SetToZero(b);

        // refresh RHS to have the correct reactions
        this->BuildRHSNoDirichlet(pScheme, rModelPart, b);

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof) {
            const std::size_t i = rDof.EquationId();

            rDof.GetSolutionStepReactionValue() = -b[i];
        });
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // intitial copy should only happen if second derivative vector is calculated
        if (mCopyExternalForceVector) {
            TSparseSpace::Copy(mCurrentExternalForceVector, mPreviousExternalForceVector);
        }
        mCopyExternalForceVector = true;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        auto default_parameters = Parameters(R"(
        {
            "name"                                 : "block_builder_and_solver_linear_elastic_dynamic",
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
        })");

        // Getting base class default parameters
        const auto base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name() { return "block_builder_and_solver_linear_elastic_dynamic"; }

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
        return "ResidualBasedBlockBuilderAndSolverLinearElasticDynamic";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    TSystemMatrixType mMassMatrix;
    TSystemMatrixType mDampingMatrix;
    TSystemVectorType mPreviousExternalForceVector;
    TSystemVectorType mCurrentExternalForceVector;

    TSystemVectorType mPreviousOutOfBalanceVector;
    TSystemVectorType mCurrentOutOfBalanceVector;

    double mBeta;
    double mGamma;
    bool   mCalculateInitialSecondDerivative;
    bool   mCopyExternalForceVector = false;
    bool   mUsePerformSolutionStep  = false;

    void BuildRHSNoDirichlet(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb)
    {
        // getting the array of the conditions
        const ConditionsArrayType& r_conditions = rModelPart.Conditions();

        mCurrentExternalForceVector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);

        // assemble all conditions
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        block_for_each(r_conditions, [&r_current_process_info, this](Condition& r_condition) {
            LocalSystemVectorType           local_external_force = LocalSystemVectorType(0);
            Condition::EquationIdVectorType equation_ids;

            if (r_condition.IsActive()) {
                r_condition.CalculateRightHandSide(local_external_force, r_current_process_info);
                r_condition.EquationIdVector(equation_ids, r_current_process_info);

                // assemble the elemental contribution
                BaseType::AssembleRHS(mCurrentExternalForceVector, local_external_force, equation_ids);
            }
        });

        // Does: mCurrentOutOfBalanceVector = mCurrentExternalForceVector - mPreviousExternalForceVector;
        TSparseSpace::ScaleAndAdd(1.0, mCurrentExternalForceVector, -1.0,
                                  mPreviousExternalForceVector, mCurrentOutOfBalanceVector);

        // Add constraint to the out of balance force before mass and damping components are added, since the mass and damping components are
        // already constraint by the constraint mass and damping matrix.
        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            BaseType::ApplyRHSConstraints(pScheme, rModelPart, mCurrentOutOfBalanceVector);
            Timer::Stop("ApplyRHSConstraints");
        }

        this->AddMassAndDampingToRhs(rModelPart, mCurrentOutOfBalanceVector);

        // Does: rb = mCurrentOutOfBalanceVector - mPreviousOutOfBalanceVector;
        TSparseSpace::ScaleAndAdd(1.0, mCurrentOutOfBalanceVector, -1.0, mPreviousOutOfBalanceVector, rb);
    }

    void CalculateGlobalMatrices(const ElementsArrayType& rElements, TSystemMatrixType& rA, ModelPart& rModelPart)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rElements, [&r_current_process_info, &rA, this](Element& r_element) {
            LocalSystemMatrixType lhs_contribution(0, 0);
            LocalSystemMatrixType mass_contribution(0, 0);
            LocalSystemMatrixType damping_contribution(0, 0);

            std::vector<std::size_t> equation_ids;

            if (r_element.IsActive()) {
                r_element.CalculateLeftHandSide(lhs_contribution, r_current_process_info);
                r_element.CalculateMassMatrix(mass_contribution, r_current_process_info);
                r_element.CalculateDampingMatrix(damping_contribution, r_current_process_info);

                r_element.EquationIdVector(equation_ids, r_current_process_info);

                if (mass_contribution.size1() != 0) {
                    BaseType::AssembleLHS(mMassMatrix, mass_contribution, equation_ids);
                }
                if (damping_contribution.size1() != 0) {
                    BaseType::AssembleLHS(mDampingMatrix, damping_contribution, equation_ids);
                }

                // Assemble the elemental contribution
                BaseType::AssembleLHS(rA, lhs_contribution, equation_ids);
            }
        });
    }

    void CalculateGlobalMatrices(const ConditionsArrayType& rConditions, TSystemMatrixType& rA, ModelPart& rModelPart)
    {
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        block_for_each(rConditions, [&r_current_process_info, &rA, this](Condition& r_condition) {
            LocalSystemMatrixType lhs_contribution(0, 0);
            LocalSystemMatrixType mass_contribution(0, 0);
            LocalSystemMatrixType damping_contribution(0, 0);

            std::vector<std::size_t> equation_ids;

            if (r_condition.IsActive()) {
                r_condition.CalculateLeftHandSide(lhs_contribution, r_current_process_info);
                r_condition.CalculateMassMatrix(mass_contribution, r_current_process_info);
                r_condition.CalculateDampingMatrix(damping_contribution, r_current_process_info);

                r_condition.EquationIdVector(equation_ids, r_current_process_info);

                if (mass_contribution.size1() != 0) {
                    BaseType::AssembleLHS(mMassMatrix, mass_contribution, equation_ids);
                }
                if (damping_contribution.size1() != 0) {
                    BaseType::AssembleLHS(mDampingMatrix, damping_contribution, equation_ids);
                }

                // Assemble the elemental contribution
                BaseType::AssembleLHS(rA, lhs_contribution, equation_ids);
            }
        });
    }

    void AddDynamicsToLhs(TSystemMatrixType& rA, const ModelPart& rModelPart)
    {
        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        double*       a_values = rA.value_data().begin();
        const double* m_values = mMassMatrix.value_data().begin();
        const double* c_values = mDampingMatrix.value_data().begin();

        // add mass and damping contribution to LHS sparse matrix
        // mass contribution: 1.0 / (mBeta * delta_time * delta_time) * M
        // damping contribution: mGamma / (mBeta * delta_time) * C
        for (std::size_t i = 0; i < rA.size1(); i++) {
            const std::size_t col_begin = rA.index1_data()[i];
            const std::size_t col_end   = rA.index1_data()[i + 1];

            for (std::size_t j = col_begin; j < col_end; ++j) {
                a_values[j] += (1.0 / (mBeta * delta_time * delta_time)) * m_values[j];
                a_values[j] += (mGamma / (mBeta * delta_time)) * c_values[j];
            }
        }
    }

    void CalculateInitialSecondDerivative(ModelPart&                    rModelPart,
                                          TSystemMatrixType&            rStiffnessMatrix,
                                          TSystemVectorType&            rExternalForce,
                                          typename TSchemeType::Pointer pScheme)
    {
        TSystemVectorType solution_step_values = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSystemVectorType first_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSystemVectorType second_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);

        auto& r_dof_set = BaseType::GetDofSet();

        block_for_each(r_dof_set, [&solution_step_values](Dof<double>& r_dof) {
            solution_step_values[r_dof.EquationId()] = r_dof.GetSolutionStepValue(0);
        });

        Geo::SparseSystemUtilities::GetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, r_dof_set, rModelPart, 0);

        // calculate initial second derivative vector
        TSystemVectorType stiffness_contribution = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::Mult(rStiffnessMatrix, solution_step_values, stiffness_contribution);

        TSystemVectorType damping_contribution = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::Mult(mDampingMatrix, first_derivative_vector, damping_contribution);

        // performs: initial_force_vector = rExternalForce - stiffness_contribution - damping_contribution;
        TSystemVectorType initial_force_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::ScaleAndAdd(1.0, rExternalForce, -1.0, stiffness_contribution, initial_force_vector);
        TSparseSpace::UnaliasedAdd(initial_force_vector, -1.0, damping_contribution);

        // apply constraint to initial force vector, as the mMassmatrix is also constrained
        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            BaseType::ApplyRHSConstraints(pScheme, rModelPart, initial_force_vector);
            Timer::Stop("ApplyRHSConstraints");
        }

        // add dirichlet conditions to initial_force_vector
        this->ApplyDirichletConditionsRhs(initial_force_vector);

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            TSystemVectorType Dxmodified(initial_force_vector.size());

            // Initialize the vector
            TSparseSpace::SetToZero(Dxmodified);
            BaseType::mpLinearSystemSolver->Solve(mMassMatrix, second_derivative_vector, initial_force_vector);

            // recover solution of the original problem
            TSparseSpace::Mult(BaseType::mT, Dxmodified, second_derivative_vector);
        } else {
            BaseType::mpLinearSystemSolver->Solve(mMassMatrix, second_derivative_vector, initial_force_vector);
        }

        BaseType::mpLinearSystemSolver->Solve(mMassMatrix, second_derivative_vector, initial_force_vector);

        Geo::SparseSystemUtilities::SetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, rModelPart);
    }

    void InitializeDynamicMatrix(TSystemMatrixType&            rMatrix,
                                 std::size_t                   MatrixSize,
                                 typename TSchemeType::Pointer pScheme,
                                 ModelPart&                    rModelPart)
    {
        BaseType::ConstructMatrixStructure(pScheme, rMatrix, rModelPart);
        TSparseSpace::SetToZero(rMatrix);
    }

    void CalculateAndAddDynamicContributionToRhs(TSystemVectorType& rSolutionVector,
                                                 TSystemMatrixType& rGlobalMatrix,
                                                 TSystemVectorType& rb)
    {
        TSystemVectorType contribution = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::Mult(rGlobalMatrix, rSolutionVector, contribution);

        TSparseSpace::UnaliasedAdd(rb, 1.0, contribution);
    }

    /**
     * @brief Function to add the mass and damping contribution to the rhs.
     * @details Damping contribution is the dot product of the global damping matrix and the first
     * derivative vector, Mass contribution is the dot product of the global mass matrix and the
     * second derivative vector
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void AddMassAndDampingToRhs(ModelPart& rModelPart, TSystemVectorType& rb)
    {
        // Get first and second derivative vector
        TSystemVectorType first_derivative_vector;
        TSystemVectorType second_derivative_vector;

        Geo::SparseSystemUtilities::GetUFirstAndSecondDerivativeVector(
            first_derivative_vector, second_derivative_vector, BaseType::mDofSet, rModelPart, 0);

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        TSystemVectorType m_part_vector =
            first_derivative_vector * (1.0 / (mBeta * delta_time)) +
            second_derivative_vector * (1.0 / (2.0 * mBeta));

        TSystemVectorType c_part_vector =
            first_derivative_vector * (mGamma / mBeta) +
            second_derivative_vector *
            (delta_time * (mGamma / (2 * mBeta) - 1));

        // calculate and add mass and damping contribution to rhs
        this->CalculateAndAddDynamicContributionToRhs(m_part_vector, mMassMatrix, rb);
        this->CalculateAndAddDynamicContributionToRhs(c_part_vector, mDampingMatrix, rb);
    }

    void ApplyDirichletConditionsRhs(TSystemVectorType& rb)
    {
        // add dirichlet conditions to RHS
        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&rb](const Dof<double>& r_dof) {
            if (r_dof.IsFixed()) {
                const std::size_t i = r_dof.EquationId();
                rb[i]               = 0.0;
            }
        });
    }

}; /* Class ResidualBasedBlockBuilderAndSolverLinearElasticDynamic */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolverLinearElasticDynamic<TSparseSpace, TDenseSpace, TLinearSolver>::SILENT_WARNINGS(
    Kratos::Flags::Create(0));

///@}

} /* namespace Kratos.*/
