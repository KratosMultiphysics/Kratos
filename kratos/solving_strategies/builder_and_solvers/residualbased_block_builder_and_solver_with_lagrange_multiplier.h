//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix
//
//
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_LAGRANGE_MULTIPLIER )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_LAGRANGE_MULTIPLIER

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
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

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>   BaseBuilderAndSolverType;
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                : "ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier",
            "scale_diagonal"                      : true,
            "silent_warnings"                     : false,
            "consider_double_lagrange_multiplier" : true
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Setting flags
        BaseType::mOptions.Set(BaseType::SCALE_DIAGONAL, ThisParameters["scale_diagonal"].GetBool());
        BaseType::mOptions.Set(BaseType::SILENT_WARNINGS, ThisParameters["silent_warnings"].GetBool());
        BaseType::mOptions.Set(DOUBLE_LAGRANGE_MULTIPLIER, ThisParameters["consider_double_lagrange_multiplier"].GetBool());
    }

    /**
     * @brief Default constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        const bool ScaleDiagonal = true,
        const bool SilentWarnings = false,
        const bool ConsiderDoubleLagrangeMultiplier = true
        ) : BaseType(pNewLinearSystemSolver, ScaleDiagonal, SilentWarnings)
    {
        // Setting flags
        BaseType::mOptions.Set(DOUBLE_LAGRANGE_MULTIPLIER, ConsiderDoubleLagrangeMultiplier);
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier() override
    {
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

        // Prints informations about the current time
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
        }

        // Base build and solve
        BaseType::BuildAndSolve(pScheme, rModelPart, rA, rDx, rb);

        // Update the Lagrange multiplier solution
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(BaseType::mSlaveIds.size()); ++i) {
            mLagrangeMultiplierVector[i] += rDx[BaseType::mEquationSystemSize + i];
        }
        
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
        }

        // Base build and solve
        BaseType::BuildRHSAndSolve(pScheme, rModelPart, rA, rDx, rb);

        // Update the Lagrange multiplier solution
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(BaseType::mSlaveIds.size()); ++i) {
            mLagrangeMultiplierVector[i] += rDx[BaseType::mEquationSystemSize + i];
        }
        
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
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        ) override
    {
        KRATOS_TRY;

        BaseType::SetUpDofSet(pScheme, rModelPart);

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPart& rModelPart
        ) override
    {
        BaseType::SetUpSystem(rModelPart);
    }

    /**
     * @brief We call this method in order to do an initial resize of the system of equations
     * @param pScheme The integration scheme considered
     * @param pA Pointer to the LHS
     * @param pDx Pointer to the unknowns increment vector
     * @param pb The pointer to the RHS
     * @param rModelPart The model part of the problem to solve
     */
    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
        ) override
    {
        // Base method
        BaseType::ResizeAndInitializeVectors(pScheme, pA, pDx, pb, rModelPart);
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

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // First iterator
        const auto it_dof_begin = BaseType::mDofSet.begin();

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for
        for (int k = 0; k<ndofs; k++) {
            auto it_dof =  it_dof_begin + k;

            if (it_dof->IsFixed()) {
                it_dof->GetSolutionStepReactionValue() = -rb[it_dof->EquationId()];
            }
        }

        // NOTE: The constraints reactions are already computed when solving the dofs
        const int number_slave_dofs = BaseType::mSlaveIds.size();
        #pragma omp parallel for
        for (int k = 0; k<number_slave_dofs; k++) {
            const IndexType equation_id = BaseType::mSlaveIds[k];
            auto it_dof = it_dof_begin + equation_id;
            it_dof->GetSolutionStepReactionValue() = mLagrangeMultiplierVector[mCorrespondanceDofsSlave[equation_id]];
        }

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
            auto& r_process_info = rModelPart.GetProcessInfo();
            if (!r_process_info.Has(CONSTRAINT_SCALE_FACTOR)) {
                TSystemMatrixType A(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize);
                BaseType::ConstructMatrixStructure(pScheme, A, rModelPart);
                this->BuildLHS(pScheme, rModelPart, A);
                r_process_info.SetValue(CONSTRAINT_SCALE_FACTOR, TSparseSpace::TwoNorm(A));
            }

            // Check T has been computed
            if (TSparseSpace::Size1(BaseType::mT) != BaseType::mSlaveIds.size() || TSparseSpace::Size2(BaseType::mT) != BaseType::mEquationSystemSize) {
                BaseType::mT.resize(BaseType::mSlaveIds.size(), BaseType::mEquationSystemSize, false);
                ConstructMasterSlaveConstraintsStructure(rModelPart);
            }

            // If not previously computed we compute
            if (TSparseSpace::TwoNorm(BaseType::mT) < std::numeric_limits<double>::epsilon()) {
                BuildMasterSlaveConstraints(rModelPart);
            }

            // Extend with the LM constribution
            const SizeType number_of_slave_dofs = TSparseSpace::Size1(BaseType::mT);

            // Auxiliar values
            const double constraint_scale_factor = r_process_info[CONSTRAINT_SCALE_FACTOR];

            // Definition of the total size of the system
            const SizeType total_size_of_system = BaseType::mEquationSystemSize + (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 2 * number_of_slave_dofs : number_of_slave_dofs);
            TSystemVectorType b_modified(total_size_of_system);

            // Copy the RHS
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(BaseType::mEquationSystemSize); ++i) {
                b_modified[i] = rb[i];
            }
            rb.resize(total_size_of_system, false);

            // Compute LM contributions
            TSystemVectorType b_lm(total_size_of_system);
            ComputeRHSLMContributions(b_lm, constraint_scale_factor);

            // Fill auxiliar vector
            TSparseSpace::UnaliasedAdd(b_modified, 1.0, b_lm);

            // Finally reassign
            TSparseSpace::Copy(b_modified, rb);
            b_modified.resize(0, false); // Free memory
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
            // First build the relation matrix
            BuildMasterSlaveConstraints(rModelPart);

            // Copy the LHS to avoid memory errors
            TSystemMatrixType copy_of_A(rA);
            TSystemMatrixType copy_of_T(BaseType::mT);
            TSystemMatrixType transpose_of_T(TSparseSpace::Size2(BaseType::mT), TSparseSpace::Size1(BaseType::mT));
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(transpose_of_T, BaseType::mT);

            // Some common values
            const SizeType number_of_slave_dofs = TSparseSpace::Size1(BaseType::mT);

            // Definition of the total size of the system
            const SizeType total_size_of_system = BaseType::mEquationSystemSize + (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 2 * number_of_slave_dofs : number_of_slave_dofs);
            TSystemVectorType b_modified(total_size_of_system);

            // Copy the RHS
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(BaseType::mEquationSystemSize); ++i) {
                b_modified[i] = rb[i];
            }
            rb.resize(total_size_of_system, false);

            // Definition of the number of blocks
            const SizeType number_of_blocks = BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER) ? 3 : 2;

            // Create blocks
            DenseMatrix<TSystemMatrixType*> matrices_p_blocks(number_of_blocks, number_of_blocks);
            DenseMatrix<double> contribution_coefficients(number_of_blocks, number_of_blocks);
            DenseMatrix<bool> transpose_blocks(number_of_blocks, number_of_blocks);

            // Definition of the auxiliar values
            auto& r_process_info = rModelPart.GetProcessInfo();
            const bool has_constraint_scale_factor = r_process_info.Has(CONSTRAINT_SCALE_FACTOR);
            const double constraint_scale_factor = has_constraint_scale_factor ? r_process_info[CONSTRAINT_SCALE_FACTOR] : TSparseSpace::TwoNorm(rA);
            if (!has_constraint_scale_factor) {
                r_process_info.SetValue(CONSTRAINT_SCALE_FACTOR, constraint_scale_factor);
            }

            /* Fill common blocks */
            // Fill blocks
            matrices_p_blocks(0,0) = &copy_of_A;
            matrices_p_blocks(0,1) = &transpose_of_T;
            matrices_p_blocks(1,0) = &copy_of_T;

            // Fill coefficients
            contribution_coefficients(0, 0) = 1.0;
            contribution_coefficients(0, 1) = constraint_scale_factor;
            contribution_coefficients(1, 0) = constraint_scale_factor;

            // Fill transpose positions
            for (IndexType i = 0; i < number_of_blocks; ++i) {
                for (IndexType j = 0; j < number_of_blocks; ++j) {
                    transpose_blocks(i, j) = false;
                }
            }

            // Assemble the blocks
            if (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) {
                // Definition of the build scale factor auxiliar value
                const bool has_build_scale_factor = r_process_info.Has(BUILD_SCALE_FACTOR);
                const double build_scale_factor = has_build_scale_factor ? r_process_info[BUILD_SCALE_FACTOR] : constraint_scale_factor;
                if (!has_build_scale_factor) {
                    r_process_info.SetValue(BUILD_SCALE_FACTOR, build_scale_factor);
                }

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
                contribution_coefficients(0, 2) = constraint_scale_factor;
                contribution_coefficients(2, 0) = constraint_scale_factor;
                contribution_coefficients(1, 1) = -build_scale_factor;
                contribution_coefficients(1, 2) = build_scale_factor;
                contribution_coefficients(2, 1) = build_scale_factor;
                contribution_coefficients(2, 2) = -build_scale_factor;

                // Assemble the matrix (NOTE: Like the identity matrix is created inside the condition must be used meanwhile is alive, so inside the condition)
                SparseMatrixMultiplicationUtility::AssembleSparseMatrixByBlocks<TSystemMatrixType>(rA, matrices_p_blocks, contribution_coefficients, transpose_blocks);
            } else {
                // Create auxiliar zero matrix
                TSystemMatrixType zero_matrix(number_of_slave_dofs, number_of_slave_dofs);

                // Fill blocks
                matrices_p_blocks(1,1) = &zero_matrix;

                // Fill coefficients
                contribution_coefficients(1, 1) = 0.0;

                // Assemble the matrix (NOTE: Like the zero matrix is created inside the condition must be used meanwhile is alive, so inside the condition)
                SparseMatrixMultiplicationUtility::AssembleSparseMatrixByBlocks<TSystemMatrixType>(rA, matrices_p_blocks, contribution_coefficients, transpose_blocks);
            }

            // Compute LM contributions
            TSystemVectorType b_lm(total_size_of_system);
            ComputeRHSLMContributions(b_lm, constraint_scale_factor);

            // Fill auxiliar vector
            TSparseSpace::UnaliasedAdd(b_modified, 1.0, b_lm);

            // Finally reassign
            TSparseSpace::Copy(b_modified, rb);
            b_modified.resize(0, false); // Free memory
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

    std::unordered_map<IndexType, IndexType> mCorrespondanceDofsSlave;
    
    TSystemVectorType mLagrangeMultiplierVector; /// This is vector containing the Lagrange multiplier solution

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

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if( it_const->IsDefined(ACTIVE) ) {
                        constraint_is_active = it_const->Is(ACTIVE);
                    }

                    if(constraint_is_active) {
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
                    lock_array[i].SetLock();
                    indices[i].insert(temp_indices[i].begin(), temp_indices[i].end());
                    lock_array[i].UnSetLock();
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
            #pragma omp parallel for reduction(+:nnz)
            for (int i = 0; i < static_cast<int>(slave_size); ++i) {
                nnz += indices[BaseType::mSlaveIds[i]].size();
            }

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

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(slave_size); ++i) {
                const IndexType row_begin = Trow_indices[i];
                const IndexType row_end = Trow_indices[i + 1];
                IndexType k = row_begin;
                const IndexType i_slave = BaseType::mSlaveIds[i];
                for (auto it = indices[i_slave].begin(); it != indices[i_slave].end(); ++it) {
                    Tcol_indices[k] = *it;
                    Tvalues[k] = 0.0;
                    ++k;
                }

                indices[i_slave].clear(); //deallocating the memory

                std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
            }

            BaseType::mT.set_filled(slave_size + 1, nnz);

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

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if (it_const->IsDefined(ACTIVE))
                    constraint_is_active = it_const->Is(ACTIVE);

                if (constraint_is_active) {
                    it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = mCorrespondanceDofsSlave[slave_equation_ids[i]];

                        // Assemble matrix row
                        BaseType::AssembleRowContribution(BaseType::mT, - transformation_matrix, i_global, i, master_equation_ids);

                        // Assemble constant vector
                        const double constant_value = constant_vector[i];
                        double& r_value = BaseType::mConstantVector[i_global];
                        #pragma omp atomic
                        r_value += constant_value;
                    }
                }
            }
        }

        // Setting the slave dofs into the T system
        for (auto eq_id : BaseType::mSlaveIds) {
            BaseType::mT(mCorrespondanceDofsSlave[eq_id], eq_id) = 1.0;
        }

        KRATOS_CATCH("")
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
        #pragma omp parallel for
        for (int i = 0; i < ndofs; ++i) {
            auto it_dof = it_dof_begin + i;
            aux_whole_dof_vector[i] = it_dof->GetSolutionStepValue();
        }

        // Compute auxiliar contribution
        TSystemVectorType aux_slave_dof_vector(number_of_slave_dofs);
        TSparseSpace::Mult(BaseType::mT, aux_whole_dof_vector, aux_slave_dof_vector);

        // Finally compute the RHS LM contribution
        noalias(aux_lm_rhs_contribution) = ScaleFactor * (BaseType::mConstantVector -  aux_slave_dof_vector);
        aux_slave_dof_vector.resize(0, false); // Free memory

        if (BaseType::mOptions.Is(DOUBLE_LAGRANGE_MULTIPLIER)) {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(number_of_slave_dofs); ++i) {
                rbLM[ndofs + i] = aux_lm_rhs_contribution[i];
                rbLM[ndofs + number_of_slave_dofs + i] = aux_lm_rhs_contribution[i];
            }
        } else {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(number_of_slave_dofs); ++i) {
                rbLM[ndofs + i] = aux_lm_rhs_contribution[i];
            }
        }
        aux_lm_rhs_contribution.resize(0, false); // Free memory

        // We compute the transposed matrix of the global relation matrix
        TSystemMatrixType T_transpose_matrix(ndofs, number_of_slave_dofs);
        SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, BaseType::mT, -ScaleFactor);

        TSparseSpace::Mult(T_transpose_matrix, mLagrangeMultiplierVector, aux_whole_dof_vector);
        #pragma omp parallel for
        for (int i = 0; i < ndofs; ++i) {
            rbLM[i] = aux_whole_dof_vector[i];
        }
        aux_whole_dof_vector.resize(0, false); // Free memory

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
const Kratos::Flags ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier<TSparseSpace, TDenseSpace, TLinearSolver>::DOUBLE_LAGRANGE_MULTIPLIER(Kratos::Flags::Create(3));

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
