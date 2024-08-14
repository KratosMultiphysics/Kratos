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
 * @class ResidualBasedBlockBuilderAndSolverWithMassAndDamping
 * @ingroup GeoMechanicsApplication
 * @brief Current class provides an implementation for builder and solving operations, while the
 * global mass and damping matrices are stored.
 * @details When the LHS is build, the global mass and damping matrices are build separately. When
 * building the RHS, the mass and damping matrices are multiplied with respectively the second and
 * first derivative vector to calculate the mass and damping contribution. The RHS is constituted by
 * the unbalanced loads (residual) Degrees of freedom are reordered putting the restrained degrees
 * of freedom at the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information. Calculation of the reactions involves a cost very similar to the calculation of
 * the total residual. This class is intended to be used when the mass and damping matrices are
 * constant throughout the iterations, using this class, when rebuilding the LHS every iteration is
 * not efficient.
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

        // mOutOfBalanceVector = ZeroVector(BaseType::mEquationSystemSize);
        mPreviousOutOfBalanceVector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        mCurrentOutOfBalanceVector  = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);

        if (mPreviousExternalForceVector.size() == 0) {
            mPreviousExternalForceVector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        }
    }

    void Build(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rb) override
    {
        Timer::Start("Build");

        BuildLHS(pScheme, rModelPart, rA);
        BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        Timer::Stop("Build");

        // apply dirichlet conditions
        TSystemVectorType dummy_b(rA.size1(), 0.0);
        TSystemVectorType dummy_rDx(rA.size1(), 0.0);

        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, dummy_rDx, rb);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mMassMatrix, dummy_rDx, dummy_b);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mDampingMatrix, dummy_rDx, dummy_b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << dummy_rDx
            << "\nRHS vector = " << rb << std::endl;

        if (mCalculateInitialSecondDerivative) {
            CalculateInitialSecondDerivative(rModelPart, rA, rb);
            mCopyExternalForceVector = true;
        }

        // only add dynamics to lhs after calculating intial force vector
        this->AddDynamicsToLhs(rA, rModelPart);

        // this approach is not working for all solvers, this approach is meant for solvers which can be prefactorized.
        // For future reference, use BaseType::SystemSolveWithPhysics(rA, rDx, rb, rModelPart) instead of the following lines if a non compatible solver is required.
        BaseType::mpLinearSystemSolver->InitializeSolutionStep(rA, dummy_rDx, rb);
    }

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total
     * number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        InitializeDynamicMatrix(mMassMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);
        InitializeDynamicMatrix(mDampingMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);

        // Assemble all elements
        const auto timer = BuiltinTimer();

        // getting the array of the conditions
        const ElementsArrayType& r_elements = rModelPart.Elements();

        this->CalculateGlobalMatrices(r_elements, rA, rModelPart);

        const ConditionsArrayType& r_conditions = rModelPart.Conditions();
        this->CalculateGlobalMatrices(r_conditions, rA, rModelPart);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() >= 1)
            << "Build time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic",
                       (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0))
            << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and
     * only the RHS is built again
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
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        // this approach is not working for all solvers, this approach is meant for solvers which can be prefactorized.
        // For future reference, use BaseType::SystemSolveWithPhysics(rA, rDx, rb, rModelPart) instead of the following line if a non compatible solver is required.
        BaseType::mpLinearSystemSolver->PerformSolutionStep(rA, rDx, rb);

        TSparseSpace::Copy(mCurrentOutOfBalanceVector, mPreviousOutOfBalanceVector);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() >= 1)
            << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "After the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb) override
    {
        KRATOS_TRY

        Timer::Start("BuildRHS");

        BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        // add dirichlet conditions to RHS
        this->ApplyDirichletConditionsRhs(rb);

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb)
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
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "block_builder_and_solver_linear_elastic_dynamic",
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
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

protected:
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

    void GetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
                                           TSystemVectorType& rSecondDerivativeVector,
                                           ModelPart&         rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&rFirstDerivativeVector, &rSecondDerivativeVector, this](Node& r_node) {
            if (r_node.IsActive()) {
                GetDerivativesForVariable(DISPLACEMENT_X, r_node, rFirstDerivativeVector, rSecondDerivativeVector);
                GetDerivativesForVariable(DISPLACEMENT_Y, r_node, rFirstDerivativeVector, rSecondDerivativeVector);

                const std::vector<const Variable<double>*> optional_variables = {
                    &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

                for (const auto p_variable : optional_variables) {
                    GetDerivativesForOptionalVariable(*p_variable, r_node, rFirstDerivativeVector,
                                                      rSecondDerivativeVector);
                }
            }
        });
    }

    void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                           const Node&             rNode,
                                           TSystemVectorType&      rFirstDerivativeVector,
                                           TSystemVectorType&      rSecondDerivativeVector) const
    {
        if (rNode.HasDofFor(rVariable)) {
            GetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
        }
    }

    void GetDerivativesForVariable(const Variable<double>& rVariable,
                                   const Node&             rNode,
                                   TSystemVectorType&      rFirstDerivativeVector,
                                   TSystemVectorType&      rSecondDerivativeVector) const
    {
        const auto& first_derivative  = rVariable.GetTimeDerivative();
        const auto& second_derivative = first_derivative.GetTimeDerivative();

        const auto equation_id               = rNode.GetDof(rVariable).EquationId();
        rFirstDerivativeVector[equation_id]  = rNode.FastGetSolutionStepValue(first_derivative);
        rSecondDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(second_derivative);
    }

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

        AddMassAndDampingToRhs(rModelPart, mCurrentOutOfBalanceVector);

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

        double* a_values = rA.value_data().begin();
        double* m_values = mMassMatrix.value_data().begin();
        double* c_values = mDampingMatrix.value_data().begin();

        // add mass and damping contribution to LHS sparse matrix
        // mass contribution: 1.0 / (mBeta * delta_time * delta_time) * M
        // damping contribution: mGamma / (mBeta * delta_time) * C
        for (std::size_t i = 0; i < rA.size1(); i++) {
            const std::size_t col_begin = rA.index1_data()[i];
            const std::size_t col_end   = rA.index1_data()[i + 1];

            for (std::size_t j = col_begin; j < col_end; ++j) {
                const std::size_t col = rA.index2_data()[j];
                a_values[j] += (1.0 / (mBeta * delta_time * delta_time)) * m_values[j];
                a_values[j] += (mGamma / (mBeta * delta_time)) * c_values[j];
            }
        }
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

    void SetDerivativesForOptionalVariable(const Variable<double>&  rVariable,
                                           Node&                    rNode,
                                           const TSystemVectorType& rFirstDerivativeVector,
                                           const TSystemVectorType& rSecondDerivativeVector)
    {
        if (rNode.HasDofFor(rVariable)) {
            SetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
        }
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

    void CalculateInitialSecondDerivative(ModelPart&         rModelPart,
                                          TSystemMatrixType& rStiffnessMatrix,
                                          TSystemVectorType& rExternalForce)
    {
        TSystemVectorType& r_solution_step_values = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSystemVectorType& r_first_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSystemVectorType& r_second_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);

        auto& r_dof_set = this->GetDofSet();

        block_for_each(r_dof_set, [&r_solution_step_values](Dof<double>& r_dof) {
            r_solution_step_values[r_dof.EquationId()] = r_dof.GetSolutionStepValue(0);
        });

        this->GetFirstAndSecondDerivativeVector(r_first_derivative_vector, r_second_derivative_vector, rModelPart);

        // calculate initial second derivative vector
        TSystemVectorType& r_stiffness_contribution = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::Mult(rStiffnessMatrix, r_solution_step_values, r_stiffness_contribution);

        TSystemVectorType& r_damping_contribution = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::Mult(mDampingMatrix, r_first_derivative_vector, r_damping_contribution);

        // initial_force_vector = rExternalForce - r_stiffness_contribution - r_damping_contribution;
        TSystemVectorType& initial_force_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSparseSpace::ScaleAndAdd(1.0, rExternalForce, -1.0, r_stiffness_contribution, initial_force_vector);
        TSparseSpace::UnaliasedAdd(initial_force_vector, -1.0, r_damping_contribution);

        // add dirichlet conditions to initial_force_vector
        this->ApplyDirichletConditionsRhs(initial_force_vector);

        // solve for initial second derivative vector
        BaseType::mpLinearSystemSolver->Solve(mMassMatrix, r_second_derivative_vector, initial_force_vector);

        this->SetFirstAndSecondDerivativeVector(r_first_derivative_vector, r_second_derivative_vector, rModelPart);
    }

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
        TSystemVectorType first_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        TSystemVectorType second_derivative_vector = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        GetFirstAndSecondDerivativeVector(first_derivative_vector, second_derivative_vector, rModelPart);

        // calculate and add mass and damping contribution to rhs
        CalculateAndAddDynamicContributionToRhs(second_derivative_vector, mMassMatrix, rb);
        CalculateAndAddDynamicContributionToRhs(first_derivative_vector, mDampingMatrix, rb);
    }

    void ApplyDirichletConditionsRhs(TSystemVectorType& rb)
    {
        // add dirichlet conditions to RHS
        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& r_dof) {
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
