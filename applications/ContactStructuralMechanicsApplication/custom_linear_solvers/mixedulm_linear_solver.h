// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MIXEDULM_SOLVER_H_INCLUDED )
#define  KRATOS_MIXEDULM_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/reorderer.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "custom_utilities/logging_settings.hpp"

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
 * @class MixedULMLinearSolver
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This solver is designed for the solution of mixed U-LM problems (this solver in particular is optimized for dual LM, to avoid the resolution).
 * @details It uses a block structure diving the matrix in UU LMLM ULM LMU blocks
 * and uses "standard" linear solvers for the different blocks as well as a GMRES for the outer part
 * @author Vicente Mataix Ferrandiz
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MixedULMLinearSolver :
    public IterativeSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>
{
public:
    ///@}
    ///@name Enums
    ///@{

    /// This enum is used to identify each index whick kind is
    enum class BlockType {
            OTHER,
            MASTER,
            SLAVE_INACTIVE,
            SLAVE_ACTIVE,
            LM_INACTIVE,
            LM_ACTIVE
            };

    ///@name Type Definitions
    ///@{

    /// Pointer definition of MixedULMLinearSolver
    KRATOS_CLASS_POINTER_DEFINITION (MixedULMLinearSolver);

    /// The base class corresponds to the an iterative solver
    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    /// The base class for the linear solver
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> LinearSolverType;

    /// The pointer to a linear solver
    typedef typename LinearSolverType::Pointer LinearSolverPointerType;

    /// The sparse matrix type
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// The vector type
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// The dense matrix type
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// The dense vector type
    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    /// The node type
    typedef Node<3> NodeType;

    /// The definition of the dof type
    typedef typename ModelPart::DofType DofType;

    /// The array containing the dofs
    typedef typename ModelPart::DofsArrayType DofsArrayType;

    /// An array of conditions
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// An array of nodes
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// The size type
    typedef std::size_t SizeType;

    /// The index type
    typedef std::size_t IndexType;

    /// A vector of indexes
    typedef DenseVector<IndexType> IndexVectorType;

    /// A vector of types
    typedef DenseVector<BlockType> BlockTypeVectorType;

    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param pSolverDispBlock The linear solver used for the displacement block
     * @param MaxTolerance The maximal tolrance considered
     * @param MaxIterationNumber The maximal number of iterations
     */
    MixedULMLinearSolver (
        LinearSolverPointerType pSolverDispBlock,
        const double MaxTolerance,
        const std::size_t MaxIterationNumber
        ) : BaseType (MaxTolerance, MaxIterationNumber),
            mpSolverDispBlock(pSolverDispBlock)
    {
        // Initializing the remaining variables
        mBlocksAreAllocated = false;
        mIsInitialized = false;
    }

    /**
     * @brief Second constructor, it uses a Kratos parameters as input instead of direct input
     * @param pSolverDispBlock The linear solver used for the displacement block
     * @param ThisParameters The configuration parameters considered
     */

    MixedULMLinearSolver(
        LinearSolverPointerType pSolverDispBlock,
        Parameters ThisParameters =  Parameters(R"({})")
        ): BaseType (),
            mpSolverDispBlock(pSolverDispBlock)

    {
        KRATOS_TRY

        // Now validate agains defaults -- this also ensures no type mismatch
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Initializing the remaining variables
        this->SetTolerance( ThisParameters["tolerance"].GetDouble() );
        this->SetMaxIterationsNumber( ThisParameters["max_iteration_number"].GetInt() );
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        mBlocksAreAllocated = false;
        mIsInitialized = false;

        KRATOS_CATCH("")
    }


    /// Copy constructor.
    MixedULMLinearSolver (const MixedULMLinearSolver& rOther)
        : BaseType(rOther),
          mpSolverDispBlock(rOther.mpSolverDispBlock),
          mBlocksAreAllocated(rOther.mBlocksAreAllocated),
          mIsInitialized(rOther.mIsInitialized),
          mMasterIndices(rOther.mMasterIndices),
          mSlaveInactiveIndices(rOther.mSlaveInactiveIndices),
          mSlaveActiveIndices(rOther.mSlaveActiveIndices),
          mLMInactiveIndices(rOther.mLMInactiveIndices),
          mLMActiveIndices(rOther.mLMActiveIndices),
          mOtherIndices(rOther.mOtherIndices),
          mGlobalToLocalIndexing(rOther.mGlobalToLocalIndexing),
          mWhichBlockType(rOther.mWhichBlockType),
          mKDispModified(rOther.mKDispModified),
          mKLMAModified(rOther.mKLMAModified),
          mKLMIModified(rOther.mKLMIModified),
          mKSAN(rOther.mKSAN),
          mKSAM(rOther.mKSAM),
          mKSASI(rOther.mKSASI),
          mKSASA(rOther.mKSASA),
          mPOperator(rOther.mPOperator),
          mCOperator(rOther.mCOperator),
          mResidualLMActive(rOther.mResidualLMActive),
          mResidualLMInactive(rOther.mResidualLMInactive),
          mResidualDisp(rOther.mResidualDisp),
          mLMActive(rOther.mLMActive),
          mLMInactive(rOther.mLMInactive),
          mDisp(rOther.mDisp),
          mEchoLevel(rOther.mEchoLevel),
          mFileCreated(rOther.mFileCreated)
    {
    }

    /// Destructor.
    ~MixedULMLinearSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MixedULMLinearSolver& operator= (const MixedULMLinearSolver& Other)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called as few times as possible. It creates the data structures
     * that only depend on the connectivity of the matrix (and not on its coefficients)
     * @details So that the memory can be allocated once and expensive operations can be done only when strictly
     * needed
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
    */
    void Initialize (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        if (mBlocksAreAllocated == true) {
            mpSolverDispBlock->Initialize(mKDispModified, mDisp, mResidualDisp);
            mIsInitialized = true;
        } else
            KRATOS_DETAIL("MixedULM Initialize") << "Linear solver intialization is deferred to the moment at which blocks are available" << std::endl;
    }

    /**
     * @brief This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
    */
    void InitializeSolutionStep (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        // Copy to local matrices
        if (mBlocksAreAllocated == false) {
            FillBlockMatrices (true, rA, rX, rB);
            mBlocksAreAllocated = true;
        } else {
            FillBlockMatrices (false, rA, rX, rB);
            mBlocksAreAllocated = true;
        }

        if(mIsInitialized == false)
            this->Initialize(rA,rX,rB);

        mpSolverDispBlock->InitializeSolutionStep(mKDispModified, mDisp, mResidualDisp);
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the
     * @details Initialize and InitializeSolutionStep functions.
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
    */
    void PerformSolutionStep (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        // Auxiliar size
        const SizeType lm_active_size = mLMActiveIndices.size();
        const SizeType lm_inactive_size = mLMInactiveIndices.size();
        const SizeType total_disp_size = mOtherIndices.size() + mMasterIndices.size() + mSlaveInactiveIndices.size() + mSlaveActiveIndices.size();

        // Get the u and lm residuals
        GetUPart (rB, mResidualDisp);

        // Solve u block
        if (mDisp.size() != total_disp_size)
            mDisp.resize(total_disp_size, false);
        mpSolverDispBlock->Solve (mKDispModified, mDisp, mResidualDisp);

        // Write back solution
        SetUPart(rX, mDisp);

        // Solve LM
        if (lm_active_size > 0) {
            // Now we compute the residual of the LM
            GetLMAPart (rB, mResidualLMActive);

            // LM = D⁻1*rLM
            if (mLMActive.size() != lm_active_size)
                mLMActive.resize(lm_active_size, false);
            TSparseSpaceType::Mult (mKLMAModified, mResidualLMActive, mLMActive);

            // Write back solution
            SetLMAPart(rX, mLMActive);
        }

        if (lm_inactive_size > 0) {
            // Now we compute the residual of the LM
            GetLMIPart (rB, mResidualLMInactive);

            // LM = D⁻1*rLM
            if (mLMInactive.size() != lm_inactive_size)
                mLMInactive.resize(lm_inactive_size, false);
            TSparseSpaceType::Mult (mKLMIModified, mResidualLMInactive, mLMInactive);

            // Write back solution
            SetLMIPart(rX, mLMInactive);
        }
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details For example this is the place to remove any data that we do not want to save for later
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
    */
    void FinalizeSolutionStep (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        mpSolverDispBlock->FinalizeSolutionStep(mKDispModified, mDisp, mResidualDisp);
    }

    /**
     * @brief This function is designed to clean up all internal data in the solver.
     * @details Clear is designed to leave the solver object as if newly created. After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mBlocksAreAllocated = false;
        mpSolverDispBlock->Clear();

        // We clear the matrixes and vectors
        mKDispModified.clear(); /// The modified displacement block
        mKLMAModified.clear();  /// The modified active LM block (diagonal)
        mKLMIModified.clear();  /// The modified inaactive LM block (diagonal)

        mKSAN.clear();  /// The slave active-displacement block
        mKSAM.clear();  /// The active slave-master block
        mKSASI.clear(); /// The active slave-inactive slave block
        mKSASA.clear(); /// The active slave-slave active block

        mPOperator.clear(); /// The operator used for the master blocks
        mCOperator.clear(); /// The operator used for the active slave block

        mResidualLMActive.clear();   /// The residual corresponding the active LM
        mResidualLMInactive.clear(); /// The residual corresponding the inactive LM
        mResidualDisp.clear();       /// The residual of the displacements

        mLMActive.clear();   /// The solution of the active LM
        mLMInactive.clear(); /// The solution of the inactive LM
        mDisp.clear();       /// The solution of the displacement

        mIsInitialized = false;
    }

    /**
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rVectorx is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        if (mIsInitialized == false)
            this->Initialize (rA,rX,rB);

        this->InitializeSolutionStep (rA,rX,rB);

        this->PerformSolutionStep (rA,rX,rB);

        this->FinalizeSolutionStep (rA,rX,rB);

        // We print the resulting system (if needed)
        if (mEchoLevel == 2) { //if it is needed to print the debug info
            KRATOS_INFO("Dx")  << "Solution obtained = " << mDisp << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << mResidualDisp << std::endl;
        } else if (mEchoLevel == 3) { //if it is needed to print the debug info
            KRATOS_INFO("LHS") << "SystemMatrix = " << mKDispModified << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << mDisp << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << mResidualDisp << std::endl;
        } else if (mEchoLevel >= 4) { //print to matrix market file
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << mFileCreated << ".mm";
            TSparseSpaceType::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), mKDispModified, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << mFileCreated << ".mm.rhs";
            TSparseSpaceType::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), mResidualDisp);
            mFileCreated++;
        }

        return false;
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rVectorx is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    bool Solve (
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB
        ) override
    {
        return false;
    }

    /**
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * @details Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /**
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * @details To make an example when solving a mixed u-p problem, it is important to identify the row associated to v and p. Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers which require knowledge on the spatial position of the nodes associated to a given dof. This function is the place to eventually provide such data
     * @param rA System matrix
     * @param rX Solution vector. It's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        DofsArrayType& rDofSet,
        ModelPart& rModelPart
        ) override
    {
        // Allocating auxiliar parameters
        IndexType node_id;

        // Count LM dofs
        SizeType n_lm_inactive_dofs = 0, n_lm_active_dofs = 0;
        SizeType n_master_dofs = 0;
        SizeType n_slave_inactive_dofs = 0, n_slave_active_dofs = 0;
        SizeType tot_active_dofs = 0;

        // We separate if we consider a block builder and solver or an elimination builder and solver
        if (rModelPart.IsNot(TO_SPLIT)) {
            // In case of block builder and solver
            for (auto& i_dof : rDofSet) {
                node_id = i_dof.Id();
                const NodeType& node = rModelPart.GetNode(node_id);
                if (i_dof.EquationId() < rA.size1()) {
                    tot_active_dofs++;
                    if (IsLMDof(i_dof)) {
                        if (node.Is(ACTIVE))
                            n_lm_active_dofs++;
                        else
                            n_lm_inactive_dofs++;
                    } else if (node.Is(INTERFACE) && IsDisplacementDof(i_dof)) {
                        if (node.Is(MASTER)) {
                            n_master_dofs++;
                        } else if (node.Is(SLAVE)) {
                            if (node.Is(ACTIVE))
                                n_slave_active_dofs++;
                            else
                                n_slave_inactive_dofs++;
                        }
                    }
                }
            }
        } else {
            // In case of elimination builder and solver
            for (auto& i_dof : rDofSet) {
                node_id = i_dof.Id();
                const NodeType& node = rModelPart.GetNode(node_id);
                tot_active_dofs++;
                if (IsLMDof(i_dof)) {
                    if (node.Is(ACTIVE))
                        n_lm_active_dofs++;
                    else
                        n_lm_inactive_dofs++;
                } else if (node.Is(INTERFACE) && IsDisplacementDof(i_dof)) {
                    if (node.Is(MASTER)) {
                        n_master_dofs++;
                    } else if (node.Is(SLAVE)) {
                        if (node.Is(ACTIVE))
                            n_slave_active_dofs++;
                        else
                            n_slave_inactive_dofs++;
                    }
                }
            }
        }

        KRATOS_ERROR_IF(tot_active_dofs != rA.size1()) << "Total system size does not coincide with the free dof map: " << tot_active_dofs << " vs " << rA.size1() << std::endl;

        // Resize arrays as needed
        if (mMasterIndices.size() != n_master_dofs)
            mMasterIndices.resize (n_master_dofs,false);
        if (mSlaveInactiveIndices.size() != n_slave_inactive_dofs)
            mSlaveInactiveIndices.resize (n_slave_inactive_dofs,false);
        if (mSlaveActiveIndices.size() != n_slave_active_dofs)
            mSlaveActiveIndices.resize (n_slave_active_dofs,false);
        if (mLMInactiveIndices.size() != n_lm_inactive_dofs)
            mLMInactiveIndices.resize (n_lm_inactive_dofs,false);
        if (mLMActiveIndices.size() != n_lm_active_dofs)
            mLMActiveIndices.resize (n_lm_active_dofs,false);

        const SizeType n_other_dofs = tot_active_dofs - n_lm_inactive_dofs - n_lm_active_dofs - n_master_dofs - n_slave_inactive_dofs - n_slave_active_dofs;
        if (mOtherIndices.size() != n_other_dofs)
            mOtherIndices.resize (n_other_dofs, false);
        if (mGlobalToLocalIndexing.size() != tot_active_dofs)
            mGlobalToLocalIndexing.resize (tot_active_dofs,false);
        if (mWhichBlockType.size() != tot_active_dofs)
            mWhichBlockType.resize(tot_active_dofs, false);

        // Size check
        KRATOS_ERROR_IF_NOT(n_lm_active_dofs == n_slave_active_dofs) << "The number of active LM dofs: " << n_lm_active_dofs << " and active slave nodes dofs: " << n_slave_active_dofs << " does not coincide" << std::endl;

        /**
         * Construct aux_lists as needed
         * "other_counter[i]" i will contain the position in the global system of the i-th NON-LM node
         * "lm_active_counter[i]" will contain the in the global system of the i-th NON-LM node
         * mGlobalToLocalIndexing[i] will contain the position in the local blocks of the
         */
        SizeType lm_inactive_counter = 0, lm_active_counter = 0;
        SizeType master_counter = 0;
        SizeType slave_inactive_counter = 0, slave_active_counter = 0;
        SizeType other_counter = 0;
        IndexType global_pos = 0;

        // We separate if we consider a block builder and solver or an elimination builder and solver
        if (rModelPart.IsNot(TO_SPLIT)) {
            // In case of block builder and solver
            for (auto& i_dof : rDofSet) {
                node_id = i_dof.Id();
                const NodeType& r_node = rModelPart.GetNode(node_id);
                if (i_dof.EquationId() < rA.size1()) {
                    if (IsLMDof(i_dof)) {
                        if (r_node.Is(ACTIVE)) {
                            mLMActiveIndices[lm_active_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = lm_active_counter;
                            mWhichBlockType[global_pos] = BlockType::LM_ACTIVE;
                            ++lm_active_counter;
                        } else {
                            mLMInactiveIndices[lm_inactive_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = lm_inactive_counter;
                            mWhichBlockType[global_pos] = BlockType::LM_INACTIVE;
                            ++lm_inactive_counter;
                        }
                    } else if ( r_node.Is(INTERFACE) && IsDisplacementDof(i_dof)) {
                        if (r_node.Is(MASTER)) {
                            mMasterIndices[master_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = master_counter;
                            mWhichBlockType[global_pos] = BlockType::MASTER;
                            ++master_counter;
                        } else if (r_node.Is(SLAVE)) {
                            if (r_node.Is(ACTIVE)) {
                                mSlaveActiveIndices[slave_active_counter] = global_pos;
                                mGlobalToLocalIndexing[global_pos] = slave_active_counter;
                                mWhichBlockType[global_pos] = BlockType::SLAVE_ACTIVE;
                                ++slave_active_counter;
                            } else {
                                mSlaveInactiveIndices[slave_inactive_counter] = global_pos;
                                mGlobalToLocalIndexing[global_pos] = slave_inactive_counter;
                                mWhichBlockType[global_pos] = BlockType::SLAVE_INACTIVE;
                                ++slave_inactive_counter;
                            }
                        } else { // We need to consider always an else to ensure that the system size is consistent
                            mOtherIndices[other_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = other_counter;
                            mWhichBlockType[global_pos] = BlockType::OTHER;
                            ++other_counter;
                        }
                    } else {
                        mOtherIndices[other_counter] = global_pos;
                        mGlobalToLocalIndexing[global_pos] = other_counter;
                        mWhichBlockType[global_pos] = BlockType::OTHER;
                        ++other_counter;
                    }
                    ++global_pos;
                }
            }
        } else {
            // In case of elimination builder and solver
            for (auto& i_dof : rDofSet) {
                node_id = i_dof.Id();
                const NodeType& r_node = rModelPart.GetNode(node_id);
                if (IsLMDof(i_dof)) {
                    if (r_node.Is(ACTIVE)) {
                        mLMActiveIndices[lm_active_counter] = global_pos;
                        mGlobalToLocalIndexing[global_pos] = lm_active_counter;
                        mWhichBlockType[global_pos] = BlockType::LM_ACTIVE;
                        ++lm_active_counter;
                    } else {
                        mLMInactiveIndices[lm_inactive_counter] = global_pos;
                        mGlobalToLocalIndexing[global_pos] = lm_inactive_counter;
                        mWhichBlockType[global_pos] = BlockType::LM_INACTIVE;
                        ++lm_inactive_counter;
                    }
                } else if ( r_node.Is(INTERFACE) && IsDisplacementDof(i_dof)) {
                    if (r_node.Is(MASTER)) {
                        mMasterIndices[master_counter] = global_pos;
                        mGlobalToLocalIndexing[global_pos] = master_counter;
                        mWhichBlockType[global_pos] = BlockType::MASTER;
                        ++master_counter;
                    } else if (r_node.Is(SLAVE)) {
                        if (r_node.Is(ACTIVE)) {
                            mSlaveActiveIndices[slave_active_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = slave_active_counter;
                            mWhichBlockType[global_pos] = BlockType::SLAVE_ACTIVE;
                            ++slave_active_counter;
                        } else {
                            mSlaveInactiveIndices[slave_inactive_counter] = global_pos;
                            mGlobalToLocalIndexing[global_pos] = slave_inactive_counter;
                            mWhichBlockType[global_pos] = BlockType::SLAVE_INACTIVE;
                            ++slave_inactive_counter;
                        }
                    } else { // We need to consider always an else to ensure that the system size is consistent
                        mOtherIndices[other_counter] = global_pos;
                        mGlobalToLocalIndexing[global_pos] = other_counter;
                        mWhichBlockType[global_pos] = BlockType::OTHER;
                        ++other_counter;
                    }
                } else {
                    mOtherIndices[other_counter] = global_pos;
                    mGlobalToLocalIndexing[global_pos] = other_counter;
                    mWhichBlockType[global_pos] = BlockType::OTHER;
                    ++other_counter;
                }
                ++global_pos;
            }
        }

        KRATOS_DEBUG_ERROR_IF(master_counter != n_master_dofs) << "The number of active slave dofs counter : " << master_counter << "is higher than the expected: " << n_master_dofs << std::endl;
        KRATOS_DEBUG_ERROR_IF(slave_active_counter != n_slave_active_dofs) << "The number of active slave dofs counter : " << slave_active_counter << "is higher than the expected: " << n_slave_active_dofs << std::endl;
        KRATOS_DEBUG_ERROR_IF(slave_inactive_counter != n_slave_inactive_dofs) << "The number of inactive slave dofs counter : " << slave_inactive_counter << "is higher than the expected: " << n_slave_inactive_dofs << std::endl;
        KRATOS_DEBUG_ERROR_IF(lm_active_counter != n_lm_active_dofs) << "The number of active LM dofs counter : " << lm_active_counter << "is higher than the expected: " << n_lm_active_dofs << std::endl;
        KRATOS_DEBUG_ERROR_IF(lm_inactive_counter != n_lm_inactive_dofs) << "The number of inactive LM dofs counter : " << lm_inactive_counter << "is higher than the expected: " << n_lm_inactive_dofs << std::endl;
        KRATOS_DEBUG_ERROR_IF(other_counter != n_other_dofs) << "The number of other dofs counter : " << other_counter << "is higher than the expected: " << n_other_dofs << std::endl;
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
        return "Mixed displacement LM linear solver";
    }

    /// Print information about this object.
    void PrintInfo (std::ostream& rOStream) const override
    {
        rOStream << "Mixed displacement LM linear solver";
    }

    /// Print object's data.
    void PrintData (std::ostream& rOStream) const override
    {
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

    /**
     * @brief T his function generates the subblocks of matrix A
     * @details as A = ( KNN  KNM    KNSI    KNSA     0        0    ) u
     *                 ( KMN  KMM    KMSI    KMSA    -MI^T    -MA^T ) u_master
     *                 ( KSIN KSIM   KSISI   KSISA   DII^T    DIA^T ) u_slave_inactive
     *                 ( KSAN KSAM   KSASI   KSASA   DAI^T    DAA^T ) u_slave_active
     *                 (  0    0      0      0       ALMI        0  ) LMInactive
     *                 (  0   KLMAM  KLMASI  KLMASA   0     KLMALMA ) LMActive
     * We will call as A = ( KNN  KNM    KNSI    KNSA     0      0      ) u
     *                     ( KMN  KMM    KMSI    KMSA   KMLMI   KMLMA   ) u_master
     *                     ( KSIN KSIM   KSISI   KSISA  KSILMI  KSILMA  ) u_slave_inactive
     *                     ( KSAN KSAM   KSASI   KSASA  KSALMI  KSALMA  ) u_slave_active
     *                     (  0    0      0      0      KLMILMI   0     ) LMInactive
     *                     (  0   KLMAM  KLMASI  KLMASA   0     KLMALMA ) LMActive
     * Subblocks are allocated or nor depending on the value of "NeedAllocation"
     * @param rA System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    void FillBlockMatrices (
        const bool NeedAllocation,
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        )
    {
        KRATOS_TRY

        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();
        const SizeType slave_active_size = mSlaveActiveIndices.size();
        const SizeType lm_active_size = mLMActiveIndices.size();
        const SizeType lm_inactive_size = mLMInactiveIndices.size();

        if (NeedAllocation)
            AllocateBlocks();

        // Get access to A data
        const IndexType* index1 = rA.index1_data().begin();
        const IndexType* index2 = rA.index2_data().begin();
        const double* values = rA.value_data().begin();

        // Allocate the auxiliar blocks by push_back
        SparseMatrixType KMLMA(master_size, lm_active_size);            /// The master-active LM block (this is the big block of M)
        SparseMatrixType KLMALMA(lm_active_size, lm_active_size);       /// The active LM-active LM block
        SparseMatrixType KSALMA(slave_active_size, lm_active_size);     /// The active slave-active LM block (this is the big block of D, diagonal)
        SparseMatrixType KLMILMI(lm_inactive_size, lm_inactive_size);   /// The inactive LM- inactive LM block (diagonal)

        IndexType* KMLMA_ptr = new IndexType[master_size + 1];
        IndexType* mKSAN_ptr = new IndexType[slave_active_size + 1];
        IndexType* mKSAM_ptr = new IndexType[slave_active_size + 1];
        IndexType* mKSASI_ptr = new IndexType[slave_active_size + 1];
        IndexType* mKSASA_ptr = new IndexType[slave_active_size + 1];
        IndexType* KSALMA_ptr = new IndexType[slave_active_size + 1];
        IndexType* KLMILMI_ptr = new IndexType[lm_inactive_size + 1];
        IndexType* KLMALMA_ptr = new IndexType[lm_active_size + 1];

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(master_size + 1); i++)
            KMLMA_ptr[i] = 0;
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(slave_active_size + 1); i++) {
            mKSAN_ptr[i] = 0;
            mKSAM_ptr[i] = 0;
            mKSASI_ptr[i] = 0;
            mKSASA_ptr[i] = 0;
            KSALMA_ptr[i] = 0;
        }
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(lm_inactive_size + 1); i++)
            KLMILMI_ptr[i] = 0;
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(lm_active_size + 1); i++)
            KLMALMA_ptr[i] = 0;

        #pragma omp parallel
        {
            // We iterate over original matrix
            #pragma omp for
            for (int i=0; i<static_cast<int>(rA.size1()); i++) {
                const IndexType row_begin = index1[i];
                const IndexType row_end   = index1[i+1];
                const IndexType local_row_id = mGlobalToLocalIndexing[i];

                IndexType KMLMA_cols = 0;
                IndexType mKSAN_cols = 0;
                IndexType mKSAM_cols = 0;
                IndexType mKSASI_cols = 0;
                IndexType mKSASA_cols = 0;
                IndexType KSALMA_cols = 0;
                IndexType KLMILMI_cols = 0;
                IndexType KLMALMA_cols = 0;

                if ( mWhichBlockType[i] == BlockType::MASTER) { // KMLMA
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if ( mWhichBlockType[col_index] == BlockType::LM_ACTIVE) { // KMLMA block
                            ++KMLMA_cols;
                        }
                    }
                    KRATOS_DEBUG_ERROR_IF(local_row_id > master_size) << "MASTER:: Local row ID: " << local_row_id <<" is greater than the number of rows " << master_size << std::endl;
                    KMLMA_ptr[local_row_id + 1] = KMLMA_cols;
                } else if ( mWhichBlockType[i] == BlockType::SLAVE_ACTIVE) { //either KSAN or KSAM or KSASA or KSASA or KSALM
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if (mWhichBlockType[col_index] == BlockType::OTHER) {                 // KSAN block
                            ++mKSAN_cols;
                        } else if (mWhichBlockType[col_index] == BlockType::MASTER) {         // KSAM block
                            ++mKSAM_cols;
                        } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) { // KSASI block
                            ++mKSASI_cols;
                        } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {   // KSASA block
                            ++mKSASA_cols;
                        } else if ( mWhichBlockType[col_index] == BlockType::LM_ACTIVE) {     // KSALMA block (diagonal)
                            ++KSALMA_cols;
                        }
                    }
                    KRATOS_DEBUG_ERROR_IF(local_row_id > slave_active_size) << "SLAVE_ACTIVE:: Local row ID: " << local_row_id <<" is greater than the number of rows " << slave_active_size << std::endl;
                    mKSAN_ptr[local_row_id + 1]  = mKSAN_cols;
                    mKSAM_ptr[local_row_id + 1]  = mKSAM_cols;
                    mKSASI_ptr[local_row_id + 1] = mKSASI_cols;
                    mKSASA_ptr[local_row_id + 1] = mKSASA_cols;
                    KSALMA_ptr[local_row_id + 1] = KSALMA_cols;
                } else if ( mWhichBlockType[i] == BlockType::LM_INACTIVE) { // KLMILMI
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if (mWhichBlockType[col_index] == BlockType::LM_INACTIVE) { // KLMILMI block (diagonal)
                            ++KLMILMI_cols;
                        }
                    }
                    KRATOS_DEBUG_ERROR_IF(local_row_id > lm_inactive_size) << "LM_INACTIVE:: Local row ID: " << local_row_id <<" is greater than the number of rows " << lm_inactive_size << std::endl;
                    KLMILMI_ptr[local_row_id + 1] = KLMILMI_cols;
                } else if ( mWhichBlockType[i] == BlockType::LM_ACTIVE) { // KLMALMA
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if (mWhichBlockType[col_index] == BlockType::LM_ACTIVE) { // KLMALMA block
                            ++KLMALMA_cols;
                        }
                    }
                    KRATOS_DEBUG_ERROR_IF(local_row_id > lm_active_size) << "LM_ACTIVE:: Local row ID: " << local_row_id <<" is greater than the number of rows " << lm_active_size << std::endl;
                    KLMALMA_ptr[local_row_id + 1] = KLMALMA_cols;
                }
            }
        }

        // We initialize the blocks sparse matrix
        std::partial_sum(KMLMA_ptr, KMLMA_ptr + master_size + 1, KMLMA_ptr);
        const std::size_t KMLMA_nonzero_values = KMLMA_ptr[master_size];
        IndexType* aux_index2_KMLMA= new IndexType[KMLMA_nonzero_values];
        double* aux_val_KMLMA= new double[KMLMA_nonzero_values];

        std::partial_sum(mKSAN_ptr, mKSAN_ptr + slave_active_size + 1, mKSAN_ptr);
        const std::size_t mKSAN_nonzero_values = mKSAN_ptr[slave_active_size];
        IndexType* aux_index2_mKSAN= new IndexType[mKSAN_nonzero_values];
        double* aux_val_mKSAN= new double[mKSAN_nonzero_values];

        std::partial_sum(mKSAM_ptr, mKSAM_ptr + slave_active_size + 1, mKSAM_ptr);
        const std::size_t mKSAM_nonzero_values = mKSAM_ptr[slave_active_size];
        IndexType* aux_index2_mKSAM= new IndexType[mKSAM_nonzero_values];
        double* aux_val_mKSAM= new double[mKSAM_nonzero_values];

        std::partial_sum(mKSASI_ptr, mKSASI_ptr + slave_active_size + 1, mKSASI_ptr);
        const std::size_t mKSASI_nonzero_values = mKSASI_ptr[slave_active_size];
        IndexType* aux_index2_mKSASI= new IndexType[mKSASI_nonzero_values];
        double* aux_val_mKSASI= new double[mKSASI_nonzero_values];

        std::partial_sum(mKSASA_ptr, mKSASA_ptr + slave_active_size + 1, mKSASA_ptr);
        const std::size_t mKSASA_nonzero_values = mKSASA_ptr[slave_active_size];
        IndexType* aux_index2_mKSASA= new IndexType[mKSASA_nonzero_values];
        double* aux_val_mKSASA = new double[mKSASA_nonzero_values];

        std::partial_sum(KSALMA_ptr, KSALMA_ptr + slave_active_size + 1, KSALMA_ptr);
        const std::size_t KSALMA_nonzero_values = KSALMA_ptr[slave_active_size];
        IndexType* aux_index2_KSALMA= new IndexType[KSALMA_nonzero_values];
        double* aux_val_KSALMA = new double[KSALMA_nonzero_values];

        std::partial_sum(KLMILMI_ptr, KLMILMI_ptr + lm_inactive_size + 1, KLMILMI_ptr);
        const std::size_t KLMILMI_nonzero_values = KLMILMI_ptr[lm_inactive_size];
        IndexType* aux_index2_KLMILMI= new IndexType[KLMILMI_nonzero_values];
        double* aux_val_KLMILMI = new double[KLMILMI_nonzero_values];

        std::partial_sum(KLMALMA_ptr, KLMALMA_ptr + lm_active_size + 1, KLMALMA_ptr);
        const std::size_t KLMALMA_nonzero_values = KLMALMA_ptr[lm_active_size];
        IndexType* aux_index2_KLMALMA = new IndexType[KLMALMA_nonzero_values];
        double* aux_val_KLMALMA = new double[KLMALMA_nonzero_values];

        #pragma omp parallel
        {
            // We iterate over original matrix
            #pragma omp for
            for (int i=0; i<static_cast<int>(rA.size1()); i++) {
                const IndexType row_begin = index1[i];
                const IndexType row_end   = index1[i+1];
                const IndexType local_row_id = mGlobalToLocalIndexing[i];

                if ( mWhichBlockType[i] == BlockType::MASTER) { // KMLMA
                    IndexType KMLMA_row_beg = KMLMA_ptr[local_row_id];
                    IndexType KMLMA_row_end = KMLMA_row_beg;
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if ( mWhichBlockType[col_index] == BlockType::LM_ACTIVE) { // KMLMA block
                            const double value = values[j];
                            const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
                            aux_index2_KMLMA[KMLMA_row_end] = local_col_id;
                            aux_val_KMLMA[KMLMA_row_end] = value;
                            ++KMLMA_row_end;
                        }
                    }
                } else if ( mWhichBlockType[i] == BlockType::SLAVE_ACTIVE) { //either KSAN or KSAM or KSASA or KSASA or KSALM
                    IndexType mKSAN_row_beg = mKSAN_ptr[local_row_id];
                    IndexType mKSAN_row_end = mKSAN_row_beg;
                    IndexType mKSAM_row_beg = mKSAM_ptr[local_row_id];
                    IndexType mKSAM_row_end = mKSAM_row_beg;
                    IndexType mKSASI_row_beg = mKSASI_ptr[local_row_id];
                    IndexType mKSASI_row_end = mKSASI_row_beg;
                    IndexType mKSASA_row_beg = mKSASA_ptr[local_row_id];
                    IndexType mKSASA_row_end = mKSASA_row_beg;
                    IndexType KSALMA_row_beg = KSALMA_ptr[local_row_id];
                    IndexType KSALMA_row_end = KSALMA_row_beg;
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        const double value = values[j];
                        const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
                        if (mWhichBlockType[col_index] == BlockType::OTHER) {                 // KSAN block
                            aux_index2_mKSAN[mKSAN_row_end] = local_col_id;
                            aux_val_mKSAN[mKSAN_row_end] = value;
                            ++mKSAN_row_end;
                        } else if (mWhichBlockType[col_index] == BlockType::MASTER) {         // KSAM block
                            aux_index2_mKSAM[mKSAM_row_end] = local_col_id;
                            aux_val_mKSAM[mKSAM_row_end] = value;
                            ++mKSAM_row_end;
                        } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) { // KSASI block
                            aux_index2_mKSASI[mKSASI_row_end] = local_col_id;
                            aux_val_mKSASI[mKSASI_row_end] = value;
                            ++mKSASI_row_end;
                        } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {   // KSASA block
                            aux_index2_mKSASA[mKSASA_row_end] = local_col_id;
                            aux_val_mKSASA[mKSASA_row_end] = value;
                            ++mKSASA_row_end;
                        } else if ( mWhichBlockType[col_index] == BlockType::LM_ACTIVE) {     // KSALMA block (diagonal)
                            aux_index2_KSALMA[KSALMA_row_end] = local_col_id;
                            aux_val_KSALMA[KSALMA_row_end] = value;
                            ++KSALMA_row_end;
                        }
                    }
                } else if ( mWhichBlockType[i] == BlockType::LM_INACTIVE) { // KLMILMI
                    IndexType KLMILMI_row_beg = KLMILMI_ptr[local_row_id];
                    IndexType KLMILMI_row_end = KLMILMI_row_beg;
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if (mWhichBlockType[col_index] == BlockType::LM_INACTIVE) { // KLMILMI block (diagonal)
                            const double value = values[j];
                            const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
                            aux_index2_KLMILMI[KLMILMI_row_end] = local_col_id;
                            aux_val_KLMILMI[KLMILMI_row_end] = value;
                            ++KLMILMI_row_end;
                        }
                    }
                } else if ( mWhichBlockType[i] == BlockType::LM_ACTIVE) { // KLMALMA
                    IndexType KLMALMA_row_beg = KLMALMA_ptr[local_row_id];
                    IndexType KLMALMA_row_end = KLMALMA_row_beg;
                    for (IndexType j=row_begin; j<row_end; j++) {
                        const IndexType col_index = index2[j];
                        if (mWhichBlockType[col_index] == BlockType::LM_ACTIVE) { // KLMALMA block
                            const double value = values[j];
                            const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
                            aux_index2_KLMALMA[KLMALMA_row_end] = local_col_id;
                            aux_val_KLMALMA[KLMALMA_row_end] = value;
                            ++KLMALMA_row_end;
                        }
                    }
                }
            }
        }

        CreateMatrix(KMLMA, master_size, lm_active_size, KMLMA_ptr, aux_index2_KMLMA, aux_val_KMLMA);
        CreateMatrix(mKSAN, slave_active_size, other_dof_size, mKSAN_ptr, aux_index2_mKSAN, aux_val_mKSAN);
        CreateMatrix(mKSAM, slave_active_size, master_size, mKSAM_ptr, aux_index2_mKSAM, aux_val_mKSAM);
        CreateMatrix(mKSASI, slave_active_size, slave_inactive_size, mKSASI_ptr, aux_index2_mKSASI, aux_val_mKSASI);
        CreateMatrix(mKSASA, slave_active_size, slave_active_size, mKSASA_ptr, aux_index2_mKSASA, aux_val_mKSASA);
        CreateMatrix(KSALMA, slave_active_size, lm_active_size, KSALMA_ptr, aux_index2_KSALMA, aux_val_KSALMA);
        CreateMatrix(KLMILMI, lm_inactive_size, lm_inactive_size, KLMILMI_ptr, aux_index2_KLMILMI, aux_val_KLMILMI);
        CreateMatrix(KLMALMA, lm_active_size, lm_active_size, KLMALMA_ptr, aux_index2_KLMALMA, aux_val_KLMALMA);

        // We compute directly the inverse of the KSALMA matrix
        // KSALMA it is supposed to be a diagonal matrix (in fact it is the key point of this formulation)
        // (NOTE: technically it is not a stiffness matrix, we give that name)
        if (lm_active_size > 0) {
            ComputeDiagonalByLumping(KSALMA, mKLMAModified, ZeroTolerance);
        }

        // We compute directly the inverse of the KLMILMI matrix
        // KLMILMI it is supposed to be a diagonal matrix (in fact it is the key point of this formulation)
        // (NOTE: technically it is not a stiffness matrix, we give that name)
        if (lm_inactive_size > 0) {
            ComputeDiagonalByLumping(KLMILMI, mKLMIModified, ZeroTolerance);
        }

        // Compute the P and C operators
        if (slave_active_size > 0) {
            SparseMatrixMultiplicationUtility::MatrixMultiplication(KMLMA,   mKLMAModified, mPOperator);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(KLMALMA, mKLMAModified, mCOperator);
        }

        // We proceed with the auxiliar products for the master blocks
        SparseMatrixType master_auxKSAN(master_size, other_dof_size);
        SparseMatrixType master_auxKSAM(master_size, master_size);
        SparseMatrixType master_auxKSASI(master_size, slave_inactive_size);
        SparseMatrixType master_auxKSASA(master_size, slave_active_size);

        if (slave_active_size > 0) {
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mPOperator, mKSAN, master_auxKSAN);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mPOperator, mKSAM, master_auxKSAM);
            if (slave_inactive_size > 0)
                SparseMatrixMultiplicationUtility::MatrixMultiplication(mPOperator, mKSASI, master_auxKSASI);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mPOperator, mKSASA, master_auxKSASA);
        }

        // We proceed with the auxiliar products for the active slave blocks
        SparseMatrixType aslave_auxKSAN(slave_active_size, other_dof_size);
        SparseMatrixType aslave_auxKSAM(slave_active_size, master_size);
        SparseMatrixType aslave_auxKSASI(slave_active_size, slave_inactive_size);
        SparseMatrixType aslave_auxKSASA(slave_active_size, slave_active_size);

        if (slave_active_size > 0) {
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mCOperator, mKSAN, aslave_auxKSAN);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mCOperator, mKSAM, aslave_auxKSAM);
            if (slave_inactive_size > 0)
                SparseMatrixMultiplicationUtility::MatrixMultiplication(mCOperator, mKSASI, aslave_auxKSASI);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(mCOperator, mKSASA, aslave_auxKSASA);
        }

        // Auxiliar indexes
        const SizeType other_dof_initial_index = 0;
        const SizeType master_dof_initial_index = other_dof_size;
        const SizeType slave_inactive_dof_initial_index = master_dof_initial_index + master_size;
        const SizeType assembling_slave_dof_initial_index = slave_inactive_dof_initial_index + slave_inactive_size;

        // The auxiliar index structure
        const SizeType nrows = mKDispModified.size1();
        const SizeType ncols = mKDispModified.size2();
        IndexType* K_disp_modified_ptr_aux1 = new IndexType[nrows + 1];
        K_disp_modified_ptr_aux1[0] = 0;

        #pragma omp parallel
        {
            #pragma omp for
            for (int i=0; i<static_cast<int>(rA.size1()); i++) {
                if ( mWhichBlockType[i] == BlockType::OTHER) { //either KNN or KNM or KNSI or KNSA
                    ComputeNonZeroColumnsDispDoFs( index1, index2, values,  i, other_dof_initial_index, K_disp_modified_ptr_aux1);
                } else if ( mWhichBlockType[i] == BlockType::MASTER) { //either KMN or KMM or KMSI or KMLM
                    ComputeNonZeroColumnsDispDoFs( index1, index2, values,  i, master_dof_initial_index, K_disp_modified_ptr_aux1);
                } else if ( mWhichBlockType[i] == BlockType::SLAVE_INACTIVE) { //either KSIN or KSIM or KSISI or KSISA
                    ComputeNonZeroColumnsDispDoFs( index1, index2, values,  i, slave_inactive_dof_initial_index, K_disp_modified_ptr_aux1);
                } else if ( mWhichBlockType[i] == BlockType::LM_ACTIVE) { //either KLMAM or KLMASI or KLMASA
                    ComputeNonZeroColumnsPartialDispDoFs( index1, index2, values,  i, assembling_slave_dof_initial_index, K_disp_modified_ptr_aux1);
                }
            }
        }

        // We initialize the final sparse matrix
        std::partial_sum(K_disp_modified_ptr_aux1, K_disp_modified_ptr_aux1 + nrows + 1, K_disp_modified_ptr_aux1);
        const SizeType nonzero_values_aux1 = K_disp_modified_ptr_aux1[nrows];
        IndexType* aux_index2_K_disp_modified_aux1 = new IndexType[nonzero_values_aux1];
        double* aux_val_K_disp_modified_aux1 = new double[nonzero_values_aux1];

        #pragma omp parallel
        {
            #pragma omp for
            for (int i=0; i<static_cast<int>(rA.size1()); i++) {
                if ( mWhichBlockType[i] == BlockType::OTHER) { //either KNN or KNM or KNSI or KNSA
                    ComputeAuxiliarValuesDispDoFs( index1, index2, values,  i, other_dof_initial_index, K_disp_modified_ptr_aux1, aux_index2_K_disp_modified_aux1, aux_val_K_disp_modified_aux1);
                } else if ( mWhichBlockType[i] == BlockType::MASTER) { //either KMN or KMM or KMSI or KMLM
                    ComputeAuxiliarValuesDispDoFs( index1, index2, values,  i, master_dof_initial_index, K_disp_modified_ptr_aux1, aux_index2_K_disp_modified_aux1, aux_val_K_disp_modified_aux1);
                } else if ( mWhichBlockType[i] == BlockType::SLAVE_INACTIVE) { //either KSIN or KSIM or KSISI or KSISA
                    ComputeAuxiliarValuesDispDoFs( index1, index2, values,  i, slave_inactive_dof_initial_index, K_disp_modified_ptr_aux1, aux_index2_K_disp_modified_aux1, aux_val_K_disp_modified_aux1);
                } else if ( mWhichBlockType[i] == BlockType::LM_ACTIVE) { //either KLMAM or KLMASI or KLMASA
                    ComputeAuxiliarValuesPartialDispDoFs( index1, index2, values,  i, assembling_slave_dof_initial_index, K_disp_modified_ptr_aux1, aux_index2_K_disp_modified_aux1, aux_val_K_disp_modified_aux1);
                }
            }
        }

        // Create the first auxiliar matrix
        CreateMatrix(mKDispModified, nrows, ncols, K_disp_modified_ptr_aux1, aux_index2_K_disp_modified_aux1, aux_val_K_disp_modified_aux1);

        // Now we create the second matrix block to sum
        IndexType* K_disp_modified_ptr_aux2 = new IndexType[nrows + 1];
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nrows + 1); i++)
            K_disp_modified_ptr_aux2[i] = 0;

        #pragma omp parallel
        {
            #pragma omp for
            for (int i=0; i<static_cast<int>(master_size); i++) {

                IndexType K_disp_modified_cols_aux2 = 0;

                // Get access to master_auxKSAN data
                if (master_auxKSAN.nnz() > 0 && other_dof_size > 0) {
                    ComputeNonZeroBlocks(master_auxKSAN, i, K_disp_modified_cols_aux2);
                }

                // Get access to master_auxKSAM data
                if (master_auxKSAM.nnz() > 0) {
                    ComputeNonZeroBlocks(master_auxKSAM, i, K_disp_modified_cols_aux2);
                }

                // Get access to master_auxKSASI data
                if (master_auxKSASI.nnz() > 0 && slave_inactive_size > 0) {
                    ComputeNonZeroBlocks(master_auxKSASI, i, K_disp_modified_cols_aux2);
                }

                // Get access to master_auxKSASA data
                if (master_auxKSASA.nnz() > 0 && slave_active_size > 0) {
                    ComputeNonZeroBlocks(master_auxKSASA, i, K_disp_modified_cols_aux2);
                }

                K_disp_modified_ptr_aux2[master_dof_initial_index + i + 1] = K_disp_modified_cols_aux2;
            }

            #pragma omp for
            for (int i=0; i<static_cast<int>(slave_active_size); i++) {

                IndexType K_disp_modified_cols_aux2 = 0;

                // Get access to aslave_auxKSAN data
                if (aslave_auxKSAN.nnz() > 0 && other_dof_size > 0) {
                    ComputeNonZeroBlocks(aslave_auxKSAN, i, K_disp_modified_cols_aux2);
                }

                // Get access to aslave_auxKSAM data
                if (aslave_auxKSAM.nnz() > 0 && master_size > 0) {
                    ComputeNonZeroBlocks(aslave_auxKSAM, i, K_disp_modified_cols_aux2);
                }

                // Get access to aslave_auxKSASI data
                if (aslave_auxKSASI.nnz() > 0 && slave_inactive_size > 0) {
                    ComputeNonZeroBlocks(aslave_auxKSASI, i, K_disp_modified_cols_aux2);
                }

                // Get access to aslave_auxKSASA data
                if (aslave_auxKSASA.nnz() > 0) {
                    ComputeNonZeroBlocks(aslave_auxKSASA, i, K_disp_modified_cols_aux2);
                }

                K_disp_modified_ptr_aux2[assembling_slave_dof_initial_index + i + 1] = K_disp_modified_cols_aux2;
            }
        }

        // We initialize the final sparse matrix
        std::partial_sum(K_disp_modified_ptr_aux2, K_disp_modified_ptr_aux2 + nrows + 1, K_disp_modified_ptr_aux2);
        const SizeType nonzero_values_aux2 = K_disp_modified_ptr_aux2[nrows];
        IndexType* aux_index2_K_disp_modified_aux2 = new IndexType[nonzero_values_aux2];
        double* aux_val_K_disp_modified_aux2 = new double[nonzero_values_aux2];

        #pragma omp parallel
        {
            #pragma omp for
            for (int i=0; i<static_cast<int>(master_size); i++) {
                const IndexType row_beg = K_disp_modified_ptr_aux2[master_dof_initial_index + i];
                IndexType row_end = row_beg;

                // Get access to master_auxKSAN data
                if (master_auxKSAN.nnz() > 0 && other_dof_size > 0) {
                    ComputeAuxiliarValuesBlocks(master_auxKSAN, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, other_dof_initial_index);
                }

                // Get access to master_auxKSAM data
                if (master_auxKSAM.nnz() > 0) {
                    ComputeAuxiliarValuesBlocks(master_auxKSAM, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, master_dof_initial_index);
                }

                // Get access to master_auxKSASI data
                if (master_auxKSASI.nnz() > 0 && slave_inactive_size > 0) {
                    ComputeAuxiliarValuesBlocks(master_auxKSASI, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, slave_inactive_dof_initial_index);
                }

                // Get access to master_auxKSASA data
                if (master_auxKSASA.nnz() > 0 && slave_active_size > 0) {
                    ComputeAuxiliarValuesBlocks(master_auxKSASA, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, assembling_slave_dof_initial_index);
                }
            }

            #pragma omp for
            for (int i=0; i<static_cast<int>(slave_active_size); i++) {
                const IndexType row_beg = K_disp_modified_ptr_aux2[assembling_slave_dof_initial_index + i];
                IndexType row_end = row_beg;

                // Get access to aslave_auxKSAN data
                if (aslave_auxKSAN.nnz() > 0 && other_dof_size > 0) {
                    ComputeAuxiliarValuesBlocks(aslave_auxKSAN, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, other_dof_initial_index);
                }

                // Get access to aslave_auxKSAM data
                if (aslave_auxKSAM.nnz() > 0 && master_size > 0) {
                    ComputeAuxiliarValuesBlocks(aslave_auxKSAM, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, master_dof_initial_index);
                }

                // Get access to aslave_auxKSASI data
                if (aslave_auxKSASI.nnz() > 0 && slave_inactive_size > 0) {
                    ComputeAuxiliarValuesBlocks(aslave_auxKSASI, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, slave_inactive_dof_initial_index);
                }

                // Get access to aslave_auxKSASA data
                if (aslave_auxKSASA.nnz() > 0) {
                    ComputeAuxiliarValuesBlocks(aslave_auxKSASA, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2, i, row_end, assembling_slave_dof_initial_index);
                }
            }
        }

        // Create the second auxiliar matrix
        SparseMatrixType K_disp_modified_aux2(nrows, ncols);
        CreateMatrix(K_disp_modified_aux2, nrows, ncols, K_disp_modified_ptr_aux2, aux_index2_K_disp_modified_aux2, aux_val_K_disp_modified_aux2);

        // We sum the auxiliar matrices
        SparseMatrixMultiplicationUtility::MatrixAdd<SparseMatrixType, SparseMatrixType>(mKDispModified, K_disp_modified_aux2, 1.0);

        // Finally we ensure that the matrix is structurally symmetric
        EnsureStructuralSymmetryMatrix(mKDispModified);

    #ifdef KRATOS_DEBUG
        CheckMatrix(mKDispModified);
    #endif

//         // DEBUG
//         LOG_MATRIX_PRETTY(rA)
//         LOG_MATRIX_PRETTY(mKDispModified)

        KRATOS_CATCH ("")
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

    LinearSolverPointerType mpSolverDispBlock; /// The pointer to the displacement linear solver

    bool mBlocksAreAllocated; /// The flag that indicates if the blocks are allocated
    bool mIsInitialized;      /// The flag that indicates if the solution is mIsInitialized

    IndexVectorType mMasterIndices;         /// The vector storing the indices of the master nodes in contact
    IndexVectorType mSlaveInactiveIndices;  /// The vector storing the indices of the slave nodes in contact (Inactive)
    IndexVectorType mSlaveActiveIndices;    /// The vector storing the indices of the slave nodes in contact (Active)
    IndexVectorType mLMInactiveIndices;     /// The vector storing the indices of the LM (Inactive)
    IndexVectorType mLMActiveIndices;       /// The vector storing the indices of the LM (Active)
    IndexVectorType mOtherIndices;          /// The vector containing the indices for other DoF
    IndexVectorType mGlobalToLocalIndexing; /// This vector stores the correspondance between the local and global
    BlockTypeVectorType mWhichBlockType; /// This vector stores the LM block belongings

    SparseMatrixType mKDispModified; /// The modified displacement block
    SparseMatrixType mKLMAModified;  /// The modified active LM block (inverted diagonal)
    SparseMatrixType mKLMIModified;  /// The modified inactive LM block (inverted diagonal)

    SparseMatrixType mKSAN;    /// The slave active-displacement block
    SparseMatrixType mKSAM;    /// The active slave-master block
    SparseMatrixType mKSASI;   /// The active slave-inactive slave block
    SparseMatrixType mKSASA;   /// The inactive slave-active slave block

    SparseMatrixType mPOperator; /// The operator used for the master blocks
    SparseMatrixType mCOperator; /// The operator used for the active slave block

    VectorType mResidualLMActive;   /// The residual of the active lagrange multipliers
    VectorType mResidualLMInactive; /// The residual of the inactive lagrange multipliers
    VectorType mResidualDisp;       /// The residual of the rest of displacements

    VectorType mLMActive;           /// The solution of the active lagrange multiplies
    VectorType mLMInactive;         /// The solution of the inactive lagrange multiplies
    VectorType mDisp;               /// The solution of the rest of displacements

    IndexType mEchoLevel = 0;       /// The echo level of the solver
    IndexType mFileCreated = 0;     /// The index used to identify the file created

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method is mean to avoid code duplication when computing the non zero terms in the Aux1 matrix
     * @param Index1 The indexes of nonzero rows
     * @param Index2 The indexes of nonzero columns
     * @param Values The array containing the values of the matrix
     * @param CurrentRow The current row computed
     * @param InitialIndex The index corresponding to the current row in the global contribution
     * @param Ptr The nonzero terms of each column
     */
    inline void ComputeNonZeroColumnsDispDoFs(
        const IndexType* Index1,
        const IndexType* Index2,
        const double* Values,
        const int CurrentRow,
        const IndexType InitialIndex,
        IndexType* Ptr
        )
    {
        const IndexType row_begin = Index1[CurrentRow];
        const IndexType row_end   = Index1[CurrentRow + 1];

        IndexType cols = 0;

        const IndexType local_row_id = mGlobalToLocalIndexing[CurrentRow] + InitialIndex;
        for (IndexType j=row_begin; j<row_end; j++) {
            const IndexType col_index = Index2[j];
            if (mWhichBlockType[col_index] == BlockType::OTHER) {
                ++cols;
            } else if (mWhichBlockType[col_index] == BlockType::MASTER) {
                ++cols;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) {
                ++cols;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {
                ++cols;
            }
        }
        Ptr[local_row_id + 1] = cols;
    }

    /**
     * @brief This method is mean to avoid code duplication when computing the non zero terms in the Aux1 matrix
     * @details The same as the previous one but not taking into account the contribution of the other dofs
     * @param Index1 The indexes of nonzero rows
     * @param Index2 The indexes of nonzero columns
     * @param Values The array containing the values of the matrix
     * @param CurrentRow The current row computed
     * @param InitialIndex The index corresponding to the current row in the global contribution
     * @param Ptr The nonzero terms of each column
     */
    inline void ComputeNonZeroColumnsPartialDispDoFs(
        const IndexType* Index1,
        const IndexType* Index2,
        const double* Values,
        const int CurrentRow,
        const IndexType InitialIndex,
        IndexType* Ptr
        )
    {
        const IndexType row_begin = Index1[CurrentRow];
        const IndexType row_end   = Index1[CurrentRow + 1];

        IndexType cols = 0;

        const IndexType local_row_id = mGlobalToLocalIndexing[CurrentRow] + InitialIndex;
        for (IndexType j=row_begin; j<row_end; j++) {
            const IndexType col_index = Index2[j];
            if (mWhichBlockType[col_index] == BlockType::MASTER) {
                ++cols;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) {
                ++cols;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {
                ++cols;
            }
        }
        Ptr[local_row_id + 1] = cols;
    }

    /**
     * @brief This method is mean to avoid code duplication when evaluate the terms of the Aux1 matrix
     * @param Index1 The indexes of nonzero rows
     * @param Index2 The indexes of nonzero columns
     * @param Values The array containing the values of the matrix
     * @param CurrentRow The current row computed
     * @param InitialIndex The index corresponding to the current row in the global contribution
     * @param Ptr The nonzero terms of each column
     * @param AuxIndex2 The indexes of the non zero columns
     * @param AuxVals The values of the final matrix
     */
    inline void ComputeAuxiliarValuesDispDoFs(
        const IndexType* Index1,
        const IndexType* Index2,
        const double* Values,
        const int CurrentRow,
        const IndexType InitialIndex,
        IndexType* Ptr,
        IndexType* AuxIndex2,
        double* AuxVals
        )
    {
        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();

        // Auxiliar indexes
        const SizeType other_dof_initial_index = 0;
        const SizeType master_dof_initial_index = other_dof_size;
        const SizeType slave_inactive_dof_initial_index = master_dof_initial_index + master_size;
        const SizeType assembling_slave_dof_initial_index = slave_inactive_dof_initial_index + slave_inactive_size;

        // Some indexes
        const IndexType local_row_id = mGlobalToLocalIndexing[CurrentRow] + InitialIndex;

        const IndexType row_begin_A = Index1[CurrentRow];
        const IndexType row_end_A   = Index1[CurrentRow + 1];

        const IndexType row_beg = Ptr[local_row_id];
        IndexType row_end = row_beg;

        for (IndexType j=row_begin_A; j<row_end_A; j++) {
            const IndexType col_index = Index2[j];
            const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
            const double value = Values[j];
            if (mWhichBlockType[col_index] == BlockType::OTHER) {
                AuxIndex2[row_end] = local_col_id + other_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            } else if (mWhichBlockType[col_index] == BlockType::MASTER) {
                AuxIndex2[row_end] = local_col_id + master_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) {
                AuxIndex2[row_end] = local_col_id + slave_inactive_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {
                AuxIndex2[row_end] = local_col_id + assembling_slave_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            }
        }
    }

    /**
     * @brief This method is mean to avoid code duplication when evaluate the terms of the Aux1 matrix
     * @details The same as the previous one but not taking into account the contribution of the other dofs
     * @param Index1 The indexes of nonzero rows
     * @param Index2 The indexes of nonzero columns
     * @param Values The array containing the values of the matrix
     * @param CurrentRow The current row computed
     * @param InitialIndex The index corresponding to the current row in the global contribution
     * @param Ptr The nonzero terms of each column
     * @param AuxIndex2 The indexes of the non zero columns
     * @param AuxVals The values of the final matrix
     */
    inline void ComputeAuxiliarValuesPartialDispDoFs(
        const IndexType* Index1,
        const IndexType* Index2,
        const double* Values,
        const int CurrentRow,
        const IndexType InitialIndex,
        IndexType* Ptr,
        IndexType* AuxIndex2,
        double* AuxVals
        )
    {
        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();

        // Auxiliar indexes
        const SizeType master_dof_initial_index = other_dof_size;
        const SizeType slave_inactive_dof_initial_index = master_dof_initial_index + master_size;
        const SizeType assembling_slave_dof_initial_index = slave_inactive_dof_initial_index + slave_inactive_size;

        // Some indexes
        const IndexType local_row_id = mGlobalToLocalIndexing[CurrentRow] + InitialIndex;

        const IndexType row_begin_A = Index1[CurrentRow];
        const IndexType row_end_A   = Index1[CurrentRow + 1];

        const IndexType row_beg = Ptr[local_row_id];
        IndexType row_end = row_beg;

        for (IndexType j=row_begin_A; j<row_end_A; j++) {
            const IndexType col_index = Index2[j];
            const IndexType local_col_id = mGlobalToLocalIndexing[col_index];
            const double value = Values[j];
            if (mWhichBlockType[col_index] == BlockType::MASTER) {
                AuxIndex2[row_end] = local_col_id + master_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_INACTIVE) {
                AuxIndex2[row_end] = local_col_id + slave_inactive_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            } else if (mWhichBlockType[col_index] == BlockType::SLAVE_ACTIVE) {
                AuxIndex2[row_end] = local_col_id + assembling_slave_dof_initial_index;
                AuxVals[row_end] = value;
                ++row_end;
            }
        }
    }

    /**
     * @brief This is a method to check the block containing nonzero values
     * @param AuxK The auxiliar block
     * @param CurrentRow The current row computed
     * @param KDispModifiedColsAux2 The nonzero rows array
     */
    inline void ComputeNonZeroBlocks(
        const SparseMatrixType& AuxK,
        const int CurrentRow,
        IndexType& KDispModifiedColsAux2
        )
    {
        // Get access to aux_K data
        const IndexType* aux_K_index1 = AuxK.index1_data().begin();

        const IndexType row_begin = aux_K_index1[CurrentRow];
        const IndexType row_end   = aux_K_index1[CurrentRow + 1];

        for (IndexType j=row_begin; j<row_end; j++) {
            ++KDispModifiedColsAux2;
        }
    }

    /**
     * @brief This is a method to compute the contribution of the auxiliar blocks
     * @param AuxK The auxiliar block
     * @param AuxIndex2 The indexes of the non zero columns
     * @param AuxVals The values of the final matrix
     * @param CurrentRow The current row computed
     * @param RowEnd The last column computed
     * @param InitialIndexColumn The initial column index of the auxiliar block in the final matrix
     */
    inline void ComputeAuxiliarValuesBlocks(
        const SparseMatrixType& AuxK,
        IndexType* AuxIndex2,
        double* AuxVals,
        const int CurrentRow,
        IndexType& RowEnd,
        const SizeType InitialIndexColumn
        )
    {
        // Get access to aux_K data
        const double* aux_values = AuxK.value_data().begin();
        const IndexType* aux_K_index1 = AuxK.index1_data().begin();
        const IndexType* aux_K_index2 = AuxK.index2_data().begin();

        const IndexType aux_K_row_begin = aux_K_index1[CurrentRow];
        const IndexType aux_K_row_end   = aux_K_index1[CurrentRow + 1];

        for (IndexType j=aux_K_row_begin; j<aux_K_row_end; j++) {
            const IndexType col_index = InitialIndexColumn + aux_K_index2[j];
            AuxIndex2[RowEnd] = col_index;
            AuxVals[RowEnd] = -aux_values[j];
            ++RowEnd;
        }
    }

    /**
     * @brief It allocates all the blocks and operators
     */
    inline void AllocateBlocks()
    {
        // We clear the matrixes
        mKDispModified.clear(); /// The modified displacement block
        mKLMAModified.clear();  /// The modified active LM block (diagonal)
        mKLMIModified.clear();  /// The modified inaactive LM block (diagonal)

        mKSAN.clear();  /// The slave active-displacement block
        mKSAM.clear();  /// The active slave-master block
        mKSASI.clear(); /// The active slave-inactive slave block
        mKSASA.clear(); /// The active slave-slave active block

        mPOperator.clear(); /// The operator used for the master blocks
        mCOperator.clear(); /// The operator used for the active slave block

        mResidualLMActive.clear();   /// The residual corresponding the active LM
        mResidualLMInactive.clear(); /// The residual corresponding the inactive LM
        mResidualDisp.clear();       /// The residual of the displacements

        mLMActive.clear();   /// The solution of the active LM
        mLMInactive.clear(); /// The solution of the inactive LM
        mDisp.clear();       /// The solution of the displacement

        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();
        const SizeType slave_active_size = mSlaveActiveIndices.size();
        const SizeType lm_active_size = mLMActiveIndices.size();
        const SizeType lm_inactive_size = mLMInactiveIndices.size();
        const SizeType total_size = other_dof_size + master_size + slave_inactive_size + slave_active_size;

        // We do the allocation
        mKDispModified.resize(total_size, total_size, false);            /// The modified displacement block
        mKLMAModified.resize(lm_active_size, lm_active_size, false);     /// The modified active LM block (diagonal)
        mKLMAModified.reserve(lm_active_size);
        mKLMIModified.resize(lm_inactive_size, lm_inactive_size, false); /// The modified inactve LM block (diagonal)
        mKLMIModified.reserve(lm_inactive_size);

        mKSAN.resize(slave_active_size, other_dof_size, false);       /// The slave active-displacement block
        mKSAM.resize(slave_active_size, master_size, false);          /// The active slave-master block
        mKSASI.resize(slave_active_size, slave_inactive_size, false); /// The active slave-inactive slave block
        mKSASA.resize(slave_active_size, slave_active_size, false);   /// The active slave-slave active block

        mPOperator.resize(master_size, slave_active_size, false);    /// The operator used for the master blocks
        mCOperator.resize(lm_active_size, slave_active_size, false); /// The operator used for the active slave block

        mResidualLMActive.resize(lm_active_size, false );     /// The residual corresponding the active LM
        mResidualLMInactive.resize(lm_inactive_size, false ); /// The residual corresponding the inactive LM
        mResidualDisp.resize(total_size );             /// The residual of the displacements

        mLMActive.resize(lm_active_size, false);     /// The solution of the active LM
        mLMInactive.resize(lm_inactive_size, false); /// The solution of the inactive LM
        mDisp.resize(total_size, false);             /// The solution of the displacement
    }

    /**
     * @brief This function extracts from a vector which has the size of the overall r, the part that corresponds to u-dofs
     * @param rTotalResidual The total residual of the problem
     * @param ResidualU The vector containing the residual relative to the displacements
     */
    inline void GetUPart (
        const VectorType& rTotalResidual,
        VectorType& ResidualU
        )
    {
        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();
        const SizeType slave_active_size = mSlaveActiveIndices.size();
        const SizeType lm_active_size = mLMActiveIndices.size();
        const SizeType total_size = other_dof_size + master_size + slave_inactive_size + slave_active_size;

        // Resize in case the size is not correct
        if (ResidualU.size() != total_size )
            ResidualU.resize (total_size, false);

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(other_dof_size); i++)
            ResidualU[i] = rTotalResidual[mOtherIndices[i]];

        // The corresponding residual for the active slave DoF's
        VectorType aux_res_active_slave(slave_active_size);
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(slave_active_size); i++)
            aux_res_active_slave[i] = rTotalResidual[mSlaveActiveIndices[i]];

        if (slave_active_size > 0) {
            // We compute the complementary residual for the master dofs
            VectorType aux_complement_master_residual(master_size);
            TSparseSpaceType::Mult(mPOperator, aux_res_active_slave, aux_complement_master_residual);

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(master_size); i++)
                ResidualU[other_dof_size + i] = rTotalResidual[mMasterIndices[i]] - aux_complement_master_residual[i];
        } else {
            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(master_size); i++)
                ResidualU[other_dof_size + i] = rTotalResidual[mMasterIndices[i]];
        }

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(slave_inactive_size); i++)
            ResidualU[other_dof_size + master_size + i] = rTotalResidual[mSlaveInactiveIndices[i]];

        if (slave_active_size > 0) {
            // We compute the complementary residual for the master dofs
            VectorType aux_complement_active_lm_residual(lm_active_size);
            TSparseSpaceType::Mult(mCOperator, aux_res_active_slave, aux_complement_active_lm_residual);

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(lm_active_size); i++)
                ResidualU[other_dof_size + master_size + slave_inactive_size + i] = rTotalResidual[mLMActiveIndices[i]] - aux_complement_active_lm_residual[i];
        } else {
            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(lm_active_size); i++)
                ResidualU[other_dof_size + master_size + slave_inactive_size + i] = rTotalResidual[mLMActiveIndices[i]];
        }
    }

    /**
     * @brief This function extracts from a vector which has the size of the overall r, the part that corresponds to active lm-dofs
     * @param rTotalResidual The total residual of the problem
     * @param rResidualLMA The vector containing the residual relative to the active LM
     */
    inline void GetLMAPart(
        const VectorType& rTotalResidual,
        VectorType& rResidualLMA
        )
    {
        // Auxiliar sizes
        const SizeType other_dof_size = mOtherIndices.size();
        const SizeType master_size = mMasterIndices.size();
        const SizeType slave_inactive_size = mSlaveInactiveIndices.size();
        const SizeType slave_active_size = mSlaveActiveIndices.size();

        // We add the other
        if (slave_active_size > 0) {

            // We get the displacement residual of the active slave nodes
            if (rResidualLMA.size() != slave_active_size )
                rResidualLMA.resize (slave_active_size, false);

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(rResidualLMA.size()); i++)
                rResidualLMA[i] = rTotalResidual[mSlaveActiveIndices[i]];

            // From the computed displacements we get the components of the displacements for each block
            VectorType disp_N(other_dof_size);
            VectorType disp_M(master_size);
            VectorType disp_SI(slave_inactive_size);
            VectorType disp_SA(slave_active_size);

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(other_dof_size); i++)
                disp_N[i] = mDisp[i];

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(master_size); i++)
                disp_M[i] = mDisp[other_dof_size + i];

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(slave_inactive_size); i++)
                disp_SI[i] = mDisp[other_dof_size + master_size + i];

            #pragma omp parallel for
            for (int i = 0; i<static_cast<int>(slave_active_size); i++)
                disp_SA[i] = mDisp[other_dof_size + master_size + slave_inactive_size + i];

            VectorType aux_mult(slave_active_size);
            TSparseSpaceType::Mult(mKSAN, disp_N, aux_mult);
            TSparseSpaceType::UnaliasedAdd (rResidualLMA, -1.0, aux_mult);
            TSparseSpaceType::Mult(mKSAM, disp_M, aux_mult);
            TSparseSpaceType::UnaliasedAdd (rResidualLMA, -1.0, aux_mult);
            if (slave_inactive_size > 0) {
                TSparseSpaceType::Mult(mKSASI, disp_SI, aux_mult);
                TSparseSpaceType::UnaliasedAdd (rResidualLMA, -1.0, aux_mult);
            }
            TSparseSpaceType::Mult(mKSASA, disp_SA, aux_mult);
            TSparseSpaceType::UnaliasedAdd (rResidualLMA, -1.0, aux_mult);
        }
    }

    /**
     * @brief This function extracts from a vector which has the size of the overall r, the part that corresponds to inactive lm-dofs
     * @param rTotalResidual The total residual of the problem
     * @param rResidualLMI The vector containing the residual relative to the inactive LM
     */
    inline void GetLMIPart (
        const VectorType& rTotalResidual,
        VectorType& rResidualLMI
        )
    {
        // Auxiliar size
        const SizeType lm_inactive_size = mLMInactiveIndices.size();

        // We get the displacement residual of the active slave nodes
        if (rResidualLMI.size() != lm_inactive_size )
            rResidualLMI.resize (lm_inactive_size, false);

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(lm_inactive_size); i++)
            rResidualLMI[i] = rTotalResidual[mLMInactiveIndices[i]];
    }

    /**
     * @brief This method writes the displacement part
     * @param rTotalResidual The total residual of the problem
     * @param ResidualU The vector containing the residual relative to the displacements
     */
    inline void SetUPart (
        VectorType& rTotalResidual,
        const VectorType& ResidualU
        )
    {
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(mOtherIndices.size()); i++)
            rTotalResidual[mOtherIndices[i]] = ResidualU[i];

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(mMasterIndices.size()); i++)
            rTotalResidual[mMasterIndices[i]] = ResidualU[mOtherIndices.size() + i];

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(mSlaveInactiveIndices.size()); i++)
            rTotalResidual[mSlaveInactiveIndices[i]] = ResidualU[mOtherIndices.size() + mMasterIndices.size() + i];

        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(mSlaveActiveIndices.size()); i++)
            rTotalResidual[mSlaveActiveIndices[i]] = ResidualU[mOtherIndices.size() + mMasterIndices.size() + mSlaveInactiveIndices.size() + i];
    }

    /**
     * @brief This method writes the active Lagrange Multiplier part
     * @param rTotalResidual The total residual of the problem
     * @param ResidualLMA The vector containing the residual relative to the active LM
     */
    inline void SetLMAPart (
        VectorType& rTotalResidual,
        const VectorType& ResidualLMA
        )
    {
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(ResidualLMA.size()); i++)
            rTotalResidual[mLMActiveIndices[i]] = ResidualLMA[i];
    }

    /**
     * @brief This method writes the inaactive Lagrange Multiplier part
     * @param rTotalResidual The total residual of the problem
     * @param ResidualLMI The vector containing the residual relative to the inactive LM
     */
    inline void SetLMIPart (
        VectorType& rTotalResidual,
        const VectorType& ResidualLMI
        )
    {
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int>(ResidualLMI.size()); i++)
            rTotalResidual[mLMInactiveIndices[i]] = ResidualLMI[i];
    }

    /**
     * @brief This method is intended to use to ensure the matrix is structurally symmetric
     * @param rA The matrix to be checked
     */
    void EnsureStructuralSymmetryMatrix (SparseMatrixType& rA)
    {
        // We compute the transposed matrix
        const SizeType size_system_1 = rA.size1();
        const SizeType size_system_2 = rA.size2();
        SparseMatrixType transpose(size_system_2, size_system_1);

        SparseMatrixMultiplicationUtility::TransposeMatrix<SparseMatrixType, SparseMatrixType>(transpose, rA, 0.0);

        // Finally we sum the auxiliar matrices
        SparseMatrixMultiplicationUtility::MatrixAdd<SparseMatrixType, SparseMatrixType>(rA, transpose, 1.0);
    }

    /**
     * @brief This method is intended to use to check the matrix
     * @param rA The matrix to be checked
     */
    double CheckMatrix (const SparseMatrixType& rA)
    {
        // Get access to A data
        const std::size_t* index1 = rA.index1_data().begin();
        const std::size_t* index2 = rA.index2_data().begin();
        const double* values = rA.value_data().begin();
        double norm = 0.0;
        for (std::size_t i=0; i<rA.size1(); ++i) {
            std::size_t row_begin = index1[i];
            std::size_t row_end   = index1[i+1];
            if (row_end - row_begin == 0)
                KRATOS_WARNING("Checking sparse matrix") << "Line " << i << " has no elements" << std::endl;

            for (std::size_t j=row_begin; j<row_end; j++) {
                KRATOS_ERROR_IF( index2[j] > rA.size2() ) << "Array above size of A" << std::endl;
                norm += values[j]*values[j];
            }
        }

        return std::sqrt (norm);
    }

    /**
     * @brief This method is designed to create the final solution sparse matrix from the auxiliar values
     * @detail Before create it reorder the columns. It deletes the auxiliar values after compute the matrix
     * @param AuxK The matrix solution
     * @param NRows The number of rows of the matrix
     * @param NCols The number of columns of the matrix
     * @param Ptr The indexes taht indicate the number of nonzero values in each column
     * @param AuxIndex2 The indexes of the nonzero columns
     * @param AuxVal The array containing the values of the sparse matrix
     */
    void CreateMatrix(
        SparseMatrixType& AuxK,
        const SizeType NRows,
        const SizeType NCols,
        IndexType* Ptr,
        IndexType* AuxIndex2,
        double* AuxVal
        )
    {
        // We reorder the rows
        SparseMatrixMultiplicationUtility::SortRows(Ptr, NRows, NCols, AuxIndex2, AuxVal);

        // Finally we build the final matrix
        SparseMatrixMultiplicationUtility::CreateSolutionMatrix(AuxK, NRows, NCols, Ptr, AuxIndex2, AuxVal);

        // Release memory
        delete[] Ptr;
        delete[] AuxIndex2;
        delete[] AuxVal;
    }

    /**
     * @brief This method is intended to lump an existing matrix
     * @param rA The matrix to be lumped
     * @param rdiagA The resulting matrix
     * @param Tolerance The tolerance considered to check if the values are almost 0
     * @todo Improve the lumping in case of not pure diagonal matrix
     */
    void ComputeDiagonalByLumping (
        const SparseMatrixType& rA,
        SparseMatrixType& rdiagA,
        const double Tolerance = ZeroTolerance
        )
    {
        // Aux values
        const std::size_t size_A = rA.size1();
//         VectorType diagA_vector(size_A);
//
//         // In case of not pure lumped matrix
//         if (rA.nnz() > size_A) {
//             // Get access to A data
//             const std::size_t* index1 = rA.index1_data().begin();
//             const double* values = rA.value_data().begin();
//
//             #pragma omp parallel for
//             for (int i=0; i< static_cast<int>(size_A); i++) {
//                 const std::size_t row_begin = index1[i];
//                 const std::size_t row_end   = index1[i+1];
//                 double temp = 0.0;
//                 for (std::size_t j=row_begin; j<row_end; j++)
//                     temp += values[j]*values[j];
//
//                 diagA_vector[i] = std::sqrt(temp);
//             }
//         } else { // Otherwise
//             #pragma omp parallel for
//             for (int i=0; i< static_cast<int>(size_A); i++) {
//                 diagA_vector[i] = rA(i, i);
//             }
//         }

        IndexType* ptr = new IndexType[size_A + 1];
        ptr[0] = 0;
        IndexType* aux_index2 = new IndexType[size_A];
        double* aux_val = new double[size_A];

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(size_A); i++) {
            ptr[i+1] = i+1;
            aux_index2[i] = i;
            const double value = rA(i, i);
//             const double value = diagA_vector[i];
            if (std::abs(value) > Tolerance)
                aux_val[i] = 1.0/value;
            else // Auxiliar value
                aux_val[i] = 1.0;
        }

        SparseMatrixMultiplicationUtility::CreateSolutionMatrix(rdiagA, size_A, size_A, ptr, aux_index2, aux_val);

        delete[] ptr;
        delete[] aux_index2;
        delete[] aux_val;
    }

    /**
     * @brief Checks if the degree of freedom belongs to a displacement DoF
     * @param rDoF The degree of freedom
     * @return True if the DoF corresponds with a displacement dof
     */
    static inline bool IsDisplacementDof(const DofType& rDoF)
    {
        const auto& r_variable = rDoF.GetVariable();
        if (r_variable == DISPLACEMENT_X ||
            r_variable == DISPLACEMENT_Y ||
            r_variable == DISPLACEMENT_Z) {
                return true;
        }

        return false;
    }

    /**
     * @brief Checks if the degree of freedom belongs to a LM DoF
     * @param rDoF The degree of freedom
     * @return True if the DoF corresponds with a LM dof
     */
    static inline bool IsLMDof(const DofType& rDoF)
    {
        const auto& r_variable = rDoF.GetVariable();
        if (r_variable == VECTOR_LAGRANGE_MULTIPLIER_X ||
            r_variable == VECTOR_LAGRANGE_MULTIPLIER_Y ||
            r_variable == VECTOR_LAGRANGE_MULTIPLIER_Z) {
                return true;
        }

        return false;
    }

    /**
     * @brief This method returns the defaulr parameters in order to avoid code duplication
     * @return Returns the default parameters
     */

    Parameters GetDefaultParameters()
    {
        Parameters default_parameters( R"(
        {
            "solver_type"          : "mixed_ulm_linear_solver",
            "tolerance"            : 1.0e-6,
            "max_iteration_number" : 200,
            "echo_level"           : 0
        }  )" );

        return default_parameters;
    }

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
}; // Class MixedULMLinearSolver
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  MixedULMLinearSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}
/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MixedULMLinearSolver<TSparseSpaceType, TDenseSpaceType,TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo (rOStream);
    rOStream << std::endl;
    rThis.PrintData (rOStream);
    return rOStream;
}
///@}
}  // namespace Kratos.
#endif // KRATOS_MIXEDULM_SOLVER_H_INCLUDED  defined
