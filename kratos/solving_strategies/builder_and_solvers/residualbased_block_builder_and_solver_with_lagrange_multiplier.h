//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix
//
//

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "utilities/atomic_utilities.h"

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
 * @class ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similar to the calculation of the total residual
 * Additionally the constraints are solver considering Lagrange multiplier (or double Lagrange multiplier)
 * @note Based on https://www.code-aster.org/V2/doc/default/en/man_r/r3/r3.03.01.pdf
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG( DOUBLE_LAGRANGE_MULTIPLIER );
    KRATOS_DEFINE_LOCAL_FLAG( TRANSFORMATION_MATRIX_COMPUTED );

    // Constraint enum
    enum class CONSTRAINT_FACTOR {CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR = 0, CONSIDER_MEAN_DIAGONAL_CONSTRAINT_FACTOR = 1, CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR = 2};
    enum class AUXILIAR_CONSTRAINT_FACTOR {CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR = 0, CONSIDER_MEAN_DIAGONAL_CONSTRAINT_FACTOR = 1, CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR = 2};

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>   BaseBuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The definition of the current class
    typedef ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    /// DoF types definition
    typedef typename Node::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
        // Setting flags
        BaseType::mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING;
        BaseType::mOptions.Set(BaseType::SILENT_WARNINGS, false);
        mConstraintFactorConsidered = CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR;
        mAuxiliarConstraintFactorConsidered = AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR;
        BaseType::mOptions.Set(DOUBLE_LAGRANGE_MULTIPLIER, true);
        BaseType::mOptions.Set(TRANSFORMATION_MATRIX_COMPUTED, false);
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier() override
    {
    }

    /**
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    typename BaseBuilderAndSolverType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
    }

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Resize again the system to the original size
        if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        // Base build
        BaseType::Build(pScheme, rModelPart, rA, rb);

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void SystemSolve(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Compute the norm
        const double norm_b = (TSparseSpace::Size(rb) != 0) ? TSparseSpace::TwoNorm(rb) : 0.0;
        if (norm_b < std::numeric_limits<double>::epsilon()) {
            // Do solve
            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx);
        }

        // Prints information about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void SystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
        ) override
    {
        BaseType::InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Resize to the proper size
        const SizeType total_system_size = (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) ? BaseType::mEquationSystemSize + 2 * BaseType::mSlaveIds.size() : BaseType::mEquationSystemSize + BaseType::mSlaveIds.size();
        if (rDx.size() != total_system_size) {
            rDx.resize(total_system_size,  false);
            TSparseSpace::SetToZero(rDx);
        }

        // Base build and solve
        BaseType::BuildAndSolve(pScheme, rModelPart, rA, rDx, rb);

        // Update the Lagrange multiplier solution
        IndexPartition<std::size_t>(BaseType::mSlaveIds.size()).for_each([&, this](std::size_t Index){
            mLagrangeMultiplierVector[Index] += rDx[this->mEquationSystemSize + Index];
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Resize to the proper size
        const SizeType total_system_size = (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) ? BaseType::mEquationSystemSize + 2 * BaseType::mSlaveIds.size() : BaseType::mEquationSystemSize + BaseType::mSlaveIds.size();
        if (rDx.size() != total_system_size) {
            rDx.resize(total_system_size,  false);
            TSparseSpace::SetToZero(rDx);
        }

        // Base build and solve
        BaseType::BuildRHSAndSolve(pScheme, rModelPart, rA, rDx, rb);

        // Update the Lagrange multiplier solution
        IndexPartition<std::size_t>(BaseType::mSlaveIds.size()).for_each([&, this](std::size_t Index){
            mLagrangeMultiplierVector[Index] += rDx[this->mEquationSystemSize + Index];
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Resize again the system to the original size
        if (rb.size() != BaseType::mEquationSystemSize) {
            rb.resize(BaseType::mEquationSystemSize, false);
        }

        // Build the base RHS
        BaseType::BuildRHS(pScheme, rModelPart, rb);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method computes the reactions of the system
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        TSparseSpace::SetToZero(rb);

        // Refresh RHS to have the correct reactions
        BaseType::BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        // First iterator
        const auto it_dof_begin = BaseType::mDofSet.begin();

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof){
            if (rDof.IsFixed()) {
                rDof.GetSolutionStepReactionValue() = -rb[rDof.EquationId()];
            }
        });

        // NOTE: The constraints reactions are already computed when solving the dofs
        IndexPartition<std::size_t>(BaseType::mSlaveIds.size()).for_each([&, this](std::size_t Index){
            const IndexType equation_id = this->mSlaveIds[Index];
            auto it_dof = it_dof_begin + equation_id;
            it_dof->GetSolutionStepReactionValue() = mLagrangeMultiplierVector[mCorrespondanceDofsSlave[equation_id]];
        });
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation chosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver chosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        const std::size_t system_size = rA.size1();
        Vector scaling_factors (system_size);

        const auto it_dof_iterator_begin = BaseType::mDofSet.begin();
        const std::size_t ndofs = BaseType::mDofSet.size();

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        IndexPartition<std::size_t>(ndofs).for_each([&](std::size_t Index){
            auto it_dof_iterator = it_dof_iterator_begin + Index;
            if (it_dof_iterator->IsFixed()) {
                scaling_factors[Index] = 0.0;
            } else {
                scaling_factors[Index] = 1.0;
            }
        });

        // Filling with ones the LM dofs
        const std::size_t loop_size = system_size - ndofs;

        IndexPartition<std::size_t>(loop_size).for_each([&](std::size_t Index){
            scaling_factors[ndofs + Index] = 1.0;
        });


        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        BaseType::mScaleFactor = TSparseSpace::GetScaleNorm(rModelPart.GetProcessInfo(), rA, BaseType::mScalingDiagonal);

        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
            std::size_t col_begin = 0, col_end  = 0;
            bool empty = true;

            col_begin = Arow_indices[Index];
            col_end = Arow_indices[Index + 1];
            empty = true;
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if(std::abs(Avalues[j]) > zero_tolerance) {
                    empty = false;
                    break;
                }
            }

            if(empty) {
                rA(Index, Index) = this->mScaleFactor;
                rb[Index] = 0.0;
            }
        });

        IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index+1];
            const double k_factor = scaling_factors[Index];
            if (k_factor == 0.0) {
                // Zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (Acol_indices[j] != Index )
                        Avalues[j] = 0.0;

                // Zero out the RHS
                rb[Index] = 0.0;
            } else {
                // Zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        });
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix (RHS only)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void ApplyRHSConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            // First we check if CONSTRAINT_SCALE_FACTOR is defined
            if (mConstraintFactorConsidered != CONSTRAINT_FACTOR::CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR) {
                TSystemMatrixType A(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize);
                BaseType::ConstructMatrixStructure(pScheme, A, rModelPart);
                this->BuildLHS(pScheme, rModelPart, A);
                const double constraint_scale_factor = mConstraintFactorConsidered == CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR ? TSparseSpace::GetMaxDiagonal(A) : TSparseSpace::GetDiagonalNorm(A);
                mConstraintFactor = constraint_scale_factor;
            }

            // Check T has been computed
            if (TSparseSpace::Size1(BaseType::mT) != BaseType::mSlaveIds.size() || TSparseSpace::Size2(BaseType::mT) != BaseType::mEquationSystemSize) {
                BaseType::mT.resize(BaseType::mSlaveIds.size(), BaseType::mEquationSystemSize, false);
                ConstructMasterSlaveConstraintsStructure(rModelPart);
            }

            // If not previously computed we compute
            if (BaseType::mOptions.IsNot(TRANSFORMATION_MATRIX_COMPUTED)) {
                BuildMasterSlaveConstraints(rModelPart);
            }

            // Extend with the LM constribution
            const SizeType number_of_slave_dofs = TSparseSpace::Size1(BaseType::mT);

            // Definition of the total size of the system
            const SizeType total_size_of_system = BaseType::mEquationSystemSize + (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 2 * number_of_slave_dofs : number_of_slave_dofs);
            TSystemVectorType b_modified(total_size_of_system);

            // Copy the RHS
            IndexPartition<std::size_t>(this->mEquationSystemSize).for_each([&](std::size_t Index){
                b_modified[Index] = rb[Index];
            });

            const SizeType loop_size = total_size_of_system - BaseType::mEquationSystemSize;
            const SizeType start_index = BaseType::mEquationSystemSize;

            // Fill with zeros
            IndexPartition<std::size_t>(loop_size).for_each([&](std::size_t Index){
                b_modified[start_index + Index] = 0.0;
            });

            rb.resize(total_size_of_system, false);

            // Compute LM contributions
            TSystemVectorType b_lm(total_size_of_system);
            ComputeRHSLMContributions(b_lm, mConstraintFactor);

            // Fill auxiliar vector
            TSparseSpace::UnaliasedAdd(b_modified, 1.0, b_lm);

            // Finally reassign
            TSparseSpace::Copy(b_modified, rb);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            // Getting process info
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // First build the relation matrix
            BuildMasterSlaveConstraints(rModelPart);

            // Copy the LHS to avoid memory errors
            TSystemMatrixType copy_of_A;
            copy_of_A.swap(rA);

            TSystemMatrixType copy_of_T(BaseType::mT);
            TSystemMatrixType transpose_of_T(TSparseSpace::Size2(BaseType::mT), TSparseSpace::Size1(BaseType::mT));
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(transpose_of_T, BaseType::mT);

            // Some common values
            const SizeType number_of_slave_dofs = TSparseSpace::Size1(BaseType::mT);

            // Definition of the total size of the system
            const SizeType total_size_of_system = BaseType::mEquationSystemSize + (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 2 * number_of_slave_dofs : number_of_slave_dofs);
            TSystemVectorType b_modified(total_size_of_system);

            // Copy the RHS
            IndexPartition<std::size_t>(this->mEquationSystemSize).for_each([&](std::size_t Index){
                b_modified[Index] = rb[Index];
            });

            auto loop_size = static_cast<int>(total_size_of_system) - static_cast<int>(BaseType::mEquationSystemSize);
            auto start_index = BaseType::mEquationSystemSize;

            // Fill with zeros
            IndexPartition<std::size_t>(loop_size).for_each([&](std::size_t Index){
                b_modified[start_index + Index] = 0.0;
            });
            rb.resize(total_size_of_system, false);

            // Definition of the number of blocks
            const SizeType number_of_blocks = BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 3 : 2;

            // Create blocks
            DenseMatrix<TSystemMatrixType*> matrices_p_blocks(number_of_blocks, number_of_blocks);
            DenseMatrix<double> contribution_coefficients(number_of_blocks, number_of_blocks);
            DenseMatrix<bool> transpose_blocks(number_of_blocks, number_of_blocks);

            // Definition of the auxiliar values
            const bool has_constraint_scale_factor = mConstraintFactorConsidered == CONSTRAINT_FACTOR::CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR ? true : false;
            KRATOS_ERROR_IF(has_constraint_scale_factor && !r_current_process_info.Has(CONSTRAINT_SCALE_FACTOR)) << "Constraint scale factor not defined at process info" << std::endl;
            const double constraint_scale_factor = has_constraint_scale_factor ? r_current_process_info.GetValue(CONSTRAINT_SCALE_FACTOR) : mConstraintFactorConsidered == CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR ? TSparseSpace::GetDiagonalNorm(copy_of_A) : TSparseSpace::GetAveragevalueDiagonal(copy_of_A);
            mConstraintFactor = constraint_scale_factor;

            /* Fill common blocks */
            // Fill blocks
            matrices_p_blocks(0,0) = &copy_of_A;
            matrices_p_blocks(0,1) = &transpose_of_T;
            matrices_p_blocks(1,0) = &copy_of_T;

            // Fill coefficients
            contribution_coefficients(0, 0) = 1.0;
            contribution_coefficients(0, 1) = mConstraintFactor;
            contribution_coefficients(1, 0) = mConstraintFactor;

            // Fill transpose positions
            for (IndexType i = 0; i < number_of_blocks; ++i) {
                for (IndexType j = 0; j < number_of_blocks; ++j) {
                    transpose_blocks(i, j) = false;
                }
            }

            // Assemble the blocks
            if (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) {
                // Definition of the build scale factor auxiliar value
                const bool has_auxiliar_constraint_scale_factor = mAuxiliarConstraintFactorConsidered == AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR ? true : false;
                KRATOS_ERROR_IF(has_auxiliar_constraint_scale_factor && !r_current_process_info.Has(AUXILIAR_CONSTRAINT_SCALE_FACTOR)) << "Auxiliar constraint scale factor not defined at process info" << std::endl;
                const double auxiliar_constraint_scale_factor = has_auxiliar_constraint_scale_factor ? r_current_process_info.GetValue(AUXILIAR_CONSTRAINT_SCALE_FACTOR) : mAuxiliarConstraintFactorConsidered == AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR ? TSparseSpace::GetDiagonalNorm(copy_of_A) : TSparseSpace::GetAveragevalueDiagonal(copy_of_A);
                mAuxiliarConstraintFactor = auxiliar_constraint_scale_factor;

                // Create auxiliar identity matrix
                TSystemMatrixType identity_matrix(number_of_slave_dofs, number_of_slave_dofs);
                for (IndexType i = 0; i < number_of_slave_dofs; ++i) {
                    identity_matrix.push_back(i, i, 1.0);
                }

                KRATOS_ERROR_IF_NOT(identity_matrix.nnz() == number_of_slave_dofs) << "Inconsistent number of non-zero values in the identity matrix: " << number_of_slave_dofs << " vs " << identity_matrix.nnz() << std::endl;

                // Fill blocks
                matrices_p_blocks(0,2) = &transpose_of_T;
                matrices_p_blocks(2,0) = &copy_of_T;
                matrices_p_blocks(1,1) = &identity_matrix;
                matrices_p_blocks(1,2) = &identity_matrix;
                matrices_p_blocks(2,1) = &identity_matrix;
                matrices_p_blocks(2,2) = &identity_matrix;

                // Fill coefficients
                contribution_coefficients(0, 2) = mConstraintFactor;
                contribution_coefficients(2, 0) = mConstraintFactor;
                contribution_coefficients(1, 1) = -mAuxiliarConstraintFactor;
                contribution_coefficients(1, 2) = mAuxiliarConstraintFactor;
                contribution_coefficients(2, 1) = mAuxiliarConstraintFactor;
                contribution_coefficients(2, 2) = -mAuxiliarConstraintFactor;

                // Assemble the matrix (NOTE: Like the identity matrix is created inside the condition must be used meanwhile is alive, so inside the condition)
                SparseMatrixMultiplicationUtility::AssembleSparseMatrixByBlocks(rA, matrices_p_blocks, contribution_coefficients, transpose_blocks);
            } else {
                // Create auxiliar zero matrix
                TSystemMatrixType zero_matrix(number_of_slave_dofs, number_of_slave_dofs);

                // Fill blocks
                matrices_p_blocks(1,1) = &zero_matrix;

                // Fill coefficients
                contribution_coefficients(1, 1) = 0.0;

                // Assemble the matrix (NOTE: Like the zero matrix is created inside the condition must be used meanwhile is alive, so inside the condition)
                SparseMatrixMultiplicationUtility::AssembleSparseMatrixByBlocks(rA, matrices_p_blocks, contribution_coefficients, transpose_blocks);
            }

            // Compute LM contributions
            TSystemVectorType b_lm(total_size_of_system);
            ComputeRHSLMContributions(b_lm, constraint_scale_factor);

            // Fill auxiliar vector
            TSparseSpace::UnaliasedAdd(b_modified, 1.0, b_lm);

            // Finally reassign
            TSparseSpace::Copy(b_modified, rb);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        // Clear member variables
        mCorrespondanceDofsSlave.clear();
        mLagrangeMultiplierVector.resize(0,false);

        // Reset flag
        BaseType::mOptions.Set(TRANSFORMATION_MATRIX_COMPUTED, false);
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part of the problem to solve
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        return BaseType::Check(rModelPart);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                               : "ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier",
            "consider_lagrange_multiplier_constraint_resolution" : "double",
            "constraint_scale_factor"                            : "use_mean_diagonal",
            "auxiliar_constraint_scale_factor"                   : "use_mean_diagonal"
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
    static std::string Name()
    {
        return "block_builder_and_solver_with_lagrange_multiplier";
    }

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
        return "ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    std::unordered_map<IndexType, IndexType> mCorrespondanceDofsSlave; /// A map of the correspondence between the slave dofs
    TSystemVectorType mLagrangeMultiplierVector;                       /// This is vector containing the Lagrange multiplier solution
    double mConstraintFactor = 0.0;                                    /// The constraint scale factor
    double mAuxiliarConstraintFactor = 0.0;                            /// The auxiliary constraint scale factor

    CONSTRAINT_FACTOR mConstraintFactorConsidered;                  /// The value considered for the constraint factor
    AUXILIAR_CONSTRAINT_FACTOR mAuxiliarConstraintFactorConsidered; /// The value considered for the auxiliary constraint factor

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method constructs the master slaeve constraint structure
     * @param rModelPart The problem model part
     */
    void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() > 0) {
            Timer::Start("ConstraintsRelationMatrixStructure");
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Vector containing the localization in the system of the different terms
            DofsVectorType slave_dof_list, master_dof_list;

            // Constraint initial iterator
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();

            const std::size_t size_indices = BaseType::mDofSet.size();
            std::vector<std::unordered_set<IndexType>> indices(size_indices);

            std::vector<LockObject> lock_array(size_indices);

            #pragma omp parallel firstprivate(slave_dof_list, master_dof_list)
            {
                Element::EquationIdVectorType slave_ids(3);
                Element::EquationIdVectorType master_ids(3);
                std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

                #pragma omp for schedule(guided, 512) nowait
                for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                    auto it_const = it_const_begin + i_const;

                    // If the constraint is active
                    if(it_const->IsActive()) {
                        it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);

                        // Slave DoFs
                        for (auto &id_i : slave_ids) {
                            temp_indices[id_i].insert(id_i);
                            temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
                        }
                    }
                }

                // Merging all the temporal indexes
                for (int i = 0; i < static_cast<int>(temp_indices.size()); ++i) {
                    lock_array[i].lock();
                    indices[i].insert(temp_indices[i].begin(), temp_indices[i].end());
                    lock_array[i].unlock();
                }
            }

            IndexType counter = 0;
            mCorrespondanceDofsSlave.clear();
            BaseType::mSlaveIds.clear();
            BaseType::mMasterIds.clear();
            for (int i = 0; i < static_cast<int>(size_indices); ++i) {
                if (indices[i].size() == 0) { // Master dof!
                    BaseType::mMasterIds.push_back(i);
                } else { // Slave dof
                    BaseType::mSlaveIds.push_back(i);
                    mCorrespondanceDofsSlave.insert(std::pair<IndexType, IndexType>(i, counter));
                    ++counter;
                }
            }

            // The slave size
            const std::size_t slave_size = BaseType::mSlaveIds.size();

            // Count the row sizes
            std::size_t nnz = 0;
            nnz = IndexPartition<std::size_t>(slave_size).for_each<SumReduction<std::size_t>>([&, this](std::size_t Index){
                return indices[this->mSlaveIds[Index]].size();
            });

            BaseType::mT = TSystemMatrixType(slave_size, size_indices, nnz);
            BaseType::mConstantVector.resize(slave_size, false);
            mLagrangeMultiplierVector.resize(slave_size, false);
            TSparseSpace::SetToZero(mLagrangeMultiplierVector);

            double *Tvalues = BaseType::mT.value_data().begin();
            IndexType *Trow_indices = BaseType::mT.index1_data().begin();
            IndexType *Tcol_indices = BaseType::mT.index2_data().begin();

            // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
            Trow_indices[0] = 0;
            for (int i = 0; i < static_cast<int>(slave_size); ++i) {
                Trow_indices[i + 1] = Trow_indices[i] + indices[BaseType::mSlaveIds[i]].size();
            }

            IndexPartition<std::size_t>(slave_size).for_each([&, this](std::size_t Index){
                const IndexType row_begin = Trow_indices[Index];
                const IndexType row_end = Trow_indices[Index + 1];
                IndexType k = row_begin;
                const IndexType i_slave = this->mSlaveIds[Index];
                for (auto it = indices[i_slave].begin(); it != indices[i_slave].end(); ++it) {
                    Tcol_indices[k] = *it;
                    Tvalues[k] = 0.0;
                    ++k;
                }

                indices[i_slave].clear(); //deallocating the memory

                std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
            });

            BaseType::mT.set_filled(slave_size + 1, nnz);

            // Reset flag
            BaseType::mOptions.Set(TRANSFORMATION_MATRIX_COMPUTED, false);

            Timer::Stop("ConstraintsRelationMatrixStructure");
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method builds the master slave relation matrix and vector
     * @param rModelPart The problem model part
     */
    void BuildMasterSlaveConstraints(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        TSparseSpace::SetToZero(BaseType::mT);
        TSparseSpace::SetToZero(BaseType::mConstantVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Contributions to the system
        Matrix transformation_matrix = LocalSystemMatrixType(0, 0);
        Vector constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_equation_ids, master_equation_ids)
        {
            #pragma omp for schedule(guided, 512)
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                // If the constraint is active
                if (it_const->IsActive()) {
                    it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = mCorrespondanceDofsSlave[slave_equation_ids[i]];

                        // Assemble matrix row
                        BaseType::AssembleRowContribution(BaseType::mT, - transformation_matrix, i_global, i, master_equation_ids);

                        // Assemble constant vector
                        const double constant_value = constant_vector[i];
                        double& r_value = BaseType::mConstantVector[i_global];
                        AtomicAdd(r_value, constant_value);
                    }
                }
            }
        }

        // Setting the slave dofs into the T system
        for (auto eq_id : BaseType::mSlaveIds) {
            BaseType::mT(mCorrespondanceDofsSlave[eq_id], eq_id) = 1.0;
        }

        // Set flag
        BaseType::mOptions.Set(TRANSFORMATION_MATRIX_COMPUTED, true);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters
        ) const override
    {
        ThisParameters.RecursivelyValidateAndAssignDefaults(DefaultParameters);
        return ThisParameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Auxiliar set for constraints
        std::set<std::string> available_options_for_constraints_scale = {"use_mean_diagonal","use_diagonal_norm","defined_in_process_info"};

        // Definition of the constraint scale factor
        const std::string& r_constraint_scale_factor = ThisParameters["constraint_scale_factor"].GetString();

        // Check the values
        if (available_options_for_constraints_scale.find(r_constraint_scale_factor) == available_options_for_constraints_scale.end()) {
            std::stringstream msg;
            msg << "Currently prescribed constraint scale factor : " << r_constraint_scale_factor << "\n";
            msg << "Admissible values for the constraint scale factor are : use_mean_diagonal, use_diagonal_norm, or defined_in_process_info" << "\n";
            KRATOS_ERROR << msg.str() << std::endl;
        }

        // This case will consider the mean value in the diagonal as a scaling value
        if (r_constraint_scale_factor == "use_mean_diagonal") {
            mConstraintFactorConsidered = CONSTRAINT_FACTOR::CONSIDER_MEAN_DIAGONAL_CONSTRAINT_FACTOR;
        } else if (r_constraint_scale_factor == "use_diagonal_norm") { // On this case the norm of the diagonal will be considered
            mConstraintFactorConsidered = CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR;
        } else { // Otherwise we will assume we impose a numerical value
            mConstraintFactorConsidered = CONSTRAINT_FACTOR::CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR;
        }

        // Definition of the auxiliar constraint scale factor
        const std::string& r_auxiliar_constraint_scale_factor = ThisParameters["auxiliar_constraint_scale_factor"].GetString();

        // Check the values
        if (available_options_for_constraints_scale.find(r_auxiliar_constraint_scale_factor) == available_options_for_constraints_scale.end()) {
            std::stringstream msg;
            msg << "Currently prescribed constraint scale factor : " << r_auxiliar_constraint_scale_factor << "\n";
            msg << "Admissible values for the constraint scale factor are : use_mean_diagonal, use_diagonal_norm, or defined_in_process_info" << "\n";
            KRATOS_ERROR << msg.str() << std::endl;
        }

        // This case will consider the mean value in the diagonal as a scaling value
        if (r_auxiliar_constraint_scale_factor == "use_mean_diagonal") {
            mAuxiliarConstraintFactorConsidered = AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_MEAN_DIAGONAL_CONSTRAINT_FACTOR;
        } else if (r_auxiliar_constraint_scale_factor == "use_diagonal_norm") { // On this case the norm of the diagonal will be considered
            mAuxiliarConstraintFactorConsidered = AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_NORM_DIAGONAL_CONSTRAINT_FACTOR;
        } else { // Otherwise we will assume we impose a numerical value
            mAuxiliarConstraintFactorConsidered = AUXILIAR_CONSTRAINT_FACTOR::CONSIDER_PRESCRIBED_CONSTRAINT_FACTOR;
        }

        // Type of LM
        if (ThisParameters["consider_lagrange_multiplier_constraint_resolution"].GetString() == "double") {
            BaseType::mOptions.Set(DOUBLE_LAGRANGE_MULTIPLIER, true);
        } else {
            BaseType::mOptions.Set(DOUBLE_LAGRANGE_MULTIPLIER, false);
        }

        // Initialize flag
        BaseType::mOptions.Set(TRANSFORMATION_MATRIX_COMPUTED, false);
    }

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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies the constraints LM contribution to the RHS
     * @param rbLM The RHS vector
     * @param ScaleFactor The scale factor considered
     */
    void ComputeRHSLMContributions(
        TSystemVectorType& rbLM,
        const double ScaleFactor = 1.0
        )
    {
        KRATOS_TRY

        const auto it_dof_begin = BaseType::mDofSet.begin();
        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // Our auxiliar vector
        const SizeType number_of_slave_dofs = TSparseSpace::Size1(BaseType::mT);
        const SizeType total_size_of_system = BaseType::mEquationSystemSize + (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 2 * number_of_slave_dofs : number_of_slave_dofs);
        if (TSparseSpace::Size(rbLM) != total_size_of_system)
            rbLM.resize(total_size_of_system, false);
        TSystemVectorType aux_lm_rhs_contribution(number_of_slave_dofs);
        TSystemVectorType aux_whole_dof_vector(ndofs);

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        IndexPartition<std::size_t>(ndofs).for_each([&](std::size_t Index){
            auto it_dof = it_dof_begin + Index;
            aux_whole_dof_vector[Index] = it_dof->GetSolutionStepValue();
        });

        // Compute auxiliar contribution
        TSystemVectorType aux_slave_dof_vector(number_of_slave_dofs);
        TSparseSpace::Mult(BaseType::mT, aux_whole_dof_vector, aux_slave_dof_vector);

        // Finally compute the RHS LM contribution
        noalias(aux_lm_rhs_contribution) = ScaleFactor * (BaseType::mConstantVector -  aux_slave_dof_vector);

        if (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) {
            IndexPartition<std::size_t>(number_of_slave_dofs).for_each([&](std::size_t Index){
                rbLM[ndofs + Index] = aux_lm_rhs_contribution[Index];
                rbLM[ndofs + number_of_slave_dofs + Index] = aux_lm_rhs_contribution[Index];
            });
        } else {
            IndexPartition<std::size_t>(number_of_slave_dofs).for_each([&](std::size_t Index){
                rbLM[ndofs + Index] = aux_lm_rhs_contribution[Index];
            });
        }

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(ndofs, number_of_slave_dofs);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, BaseType::mT, -ScaleFactor);

        TSparseSpace::Mult(T_transpose_matrix, mLagrangeMultiplierVector, aux_whole_dof_vector);
        IndexPartition<std::size_t>(ndofs).for_each([&](std::size_t Index){
            rbLM[Index] = aux_whole_dof_vector[Index];
        });

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

    ///@}

}; /* Class ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier<TSparseSpace, TDenseSpace, TLinearSolver>::DOUBLE_LAGRANGE_MULTIPLIER(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier<TSparseSpace, TDenseSpace, TLinearSolver>::TRANSFORMATION_MATRIX_COMPUTED(Kratos::Flags::Create(3));

///@}

} /* namespace Kratos.*/
