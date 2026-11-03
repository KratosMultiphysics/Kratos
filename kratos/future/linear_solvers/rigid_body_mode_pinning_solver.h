//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "future/linear_solvers/linear_solver.h"

namespace Kratos::Future
{
///@name Kratos Classes
///@{

/**
 * @class RigidBodyModePinningSolver
 * @ingroup KratosCore
 * @brief Linear solver wrapper that removes the (effective) rigid body modes of the
 *        unconstrained components of the system by pinning a minimal set of DOFs.
 * @details Before each solve, the sparsity graph of A is split into connected components.
 *          For each component the residual of the rigid body modes Z is computed as R = A Z
 *          (accumulated row by row, never stored) and the unconstrained modes are the
 *          eigenvectors of the Gram matrix R_c^T R_c with (numerically) zero residual.
 *          A greedy pivoted Gram-Schmidt on the rows of Z_eff selects exactly rank(Z_eff)
 *          DOFs which are pinned to zero. Pinning is applied IN PLACE on A and b:
 *          the off-diagonal entries of pinned rows/columns are set to zero and the diagonal
 *          is kept, so the sparsity pattern of A is never changed and no dense matrix of
 *          the size of A is built.
 *          Number of rigid body modes: 3 in 2D (u_x, u_y, theta_z), 6 in 3D.
 * @note The rows of A are modified in place (as ScalingSolver does): the matrix must be
 *       re-assembled before the next call to PerformSolutionStep.
 * @note The small 3x3/6x6 eigen-decomposition is a cyclic Jacobi routine, since Eigen is not
 *       part of the Kratos external libraries.
 * @tparam TLinearAlgebra The version of the linear algebra to be used
 */
template<class TLinearAlgebra>
class RigidBodyModePinningSolver
    : public LinearSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidBodyModePinningSolver
    KRATOS_CLASS_POINTER_DEFINITION(RigidBodyModePinningSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TLinearAlgebra>;

    /// The definition of the linear system type
    using LinearSystemType = typename BaseType::LinearSystemType;

    /// The definition of the spaces (sparse matrix)
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TLinearAlgebra::VectorType;

    /// The index type definition to be consistent
    using IndexType = typename BaseType::IndexType;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor without parameters
     * @param pLinearSolver The linear solver to be wrapped
     * @param Tolerance Relative tolerance used to detect the unconstrained rigid body modes
     */
    RigidBodyModePinningSolver(
        typename BaseType::Pointer pLinearSolver,
        const double Tolerance = 1.0e-6
        ) : BaseType(),
            mpLinearSolver(pLinearSolver),
            mTolerance(Tolerance)
    {
    }

    /// Copy constructor.
    RigidBodyModePinningSolver(const RigidBodyModePinningSolver& Other) = delete;

    /// Destructor.
    ~RigidBodyModePinningSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RigidBodyModePinningSolver& operator=(const RigidBodyModePinningSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Delegates the initialization to the wrapped solver
    void Initialize(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->Initialize(rLinearSystem);
    }

    /// Delegates the initialization of the solution step to the wrapped solver
    void InitializeSolutionStep(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->InitializeSolutionStep(rLinearSystem);
    }

    /**
     * @brief Pins the rigid body modes of the unconstrained components in place and solves.
     * @param rLinearSystem The linear system to be solved (the matrix and the vectors are
     *                      modified in place: pinned rows/columns are zeroed off-diagonal,
     *                      pinned entries of the rhs and of the initial guess are set to zero)
     */
    bool PerformSolutionStep(LinearSystemType& rLinearSystem) override
    {
        if (!rLinearSystem.IsConsistent(SparseMatrixTag::LHS, DenseVectorTag::RHS, DenseVectorTag::Dx)) {
            return false;
        }

        auto& r_lin_op = *rLinearSystem.pGetLinearOperator(SparseMatrixTag::LHS);
        KRATOS_ERROR_IF(r_lin_op.IsMatrixFree()) << "RigidBodyModePinningSolver cannot be used with matrix-free linear operators." << std::endl;
        auto& rA = r_lin_op.GetMatrix();
        auto& rX = *rLinearSystem.pGetVector(DenseVectorTag::Dx);
        auto& rB = *rLinearSystem.pGetVector(DenseVectorTag::RHS);

        KRATOS_ERROR_IF(mRigidBodyModes.size1() != rA.size1() || mNumRigidBodyModes == 0)
            << "RigidBodyModePinningSolver: PrepareAdditionalData must be called before PerformSolutionStep (node data missing)" << std::endl;

        ComputePinningMask(rA);
        ApplyPinning(rA, rX, rB);

        return mpLinearSolver->PerformSolutionStep(rLinearSystem);
    }

    /// Delegates the finalization of the solution step to the wrapped solver
    void FinalizeSolutionStep(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->FinalizeSolutionStep(rLinearSystem);
    }

    /// Delegates the cleanup to the wrapped solver
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
     * @brief Builds the rigid body mode matrix Z (n x M) from the node coordinates and DOF types.
     * @details M = 3 in 2D and M = 6 in 3D (see FillRigidBodyModeRow for the column layout).
     *          The model part is retrieved from the additional data of the linear system.
     */
    void PrepareAdditionalData(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->PrepareAdditionalData(rLinearSystem);

        auto& r_model_part = *rLinearSystem.pGetModelPart();
        auto& r_A = *rLinearSystem.pGetMatrix(SparseMatrixTag::LHS);

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_ERROR_IF_NOT(domain_size == 2 || domain_size == 3)
            << "RigidBodyModePinningSolver: DOMAIN_SIZE must be 2 or 3 (found " << domain_size << ")" << std::endl;
        mNumRigidBodyModes = (domain_size == 2) ? 3 : 6;

        mRigidBodyModes.resize(r_A.size1(), mNumRigidBodyModes, false);
        std::fill(mRigidBodyModes.data().begin(), mRigidBodyModes.data().begin() + mRigidBodyModes.size1() * mRigidBodyModes.size2(), 0.0);

        for (const auto& r_node : r_model_part.Nodes()) {
            for (const auto& p_dof : r_node.GetDofs()) {
                const std::size_t eq_id = p_dof->EquationId();
                if (eq_id < r_A.size1()) {
                    FillRigidBodyModeRow(p_dof->GetVariable(), r_node.X0(), r_node.Y0(), r_node.Z0(), row(mRigidBodyModes, eq_id));
                }
            }
        }
    }

    ///@}
    ///@name Access
    ///@{

    /// The rigid body modes need the spatial position of the nodes
    bool RequiresAdditionalData() const override
    {
        return true;
    }

    ///@}
    ///@name Inquiry
    ///@{

    IndexType GetIterationsNumber() override
    {
        return mpLinearSolver->GetIterationsNumber();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Rigid body mode pinning solver. Uses internally the following linear solver " << mpLinearSolver->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    ///@}

private:
    ///@name Private Static Helpers
    ///@{

    static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

    /// Flat, contiguous view of row i of a row-major dense matrix
    static Kratos::span<double> RowSpan(Matrix& rMatrix, const std::size_t Row)
    {
        return Kratos::span<double>(rMatrix.data().begin() + Row * rMatrix.size2(), rMatrix.size2());
    }

    /// Flat view of the M x M block of component c inside a stacked (n_comp*M) x M matrix
    static Kratos::span<double> BlockSpan(Matrix& rStacked, const std::size_t Component, const std::size_t M)
    {
        return Kratos::span<double>(rStacked.data().begin() + Component * M * M, M * M);
    }

    static std::size_t Find(std::vector<std::size_t>& rParent, std::size_t i)
    {
        while (rParent[i] != i) {
            rParent[i] = rParent[rParent[i]]; // path halving
            i = rParent[i];
        }
        return i;
    }

    static void Unite(std::vector<std::size_t>& rParent, const std::size_t I, const std::size_t J)
    {
        const std::size_t ri = Find(rParent, I);
        const std::size_t rj = Find(rParent, J);
        if (ri != rj) {
            rParent[std::max(ri, rj)] = std::min(ri, rj);
        }
    }

    /**
     * @brief Fills the rigid body mode row of one DOF.
     * @details 3D layout: [t_x, t_y, t_z, r_x, r_y, r_z]; 2D layout: [t_x, t_y, r_z].
     *          Rows of DOFs that are not displacement/rotation components stay zero.
     */
    static void FillRigidBodyModeRow(const VariableData& rVariable, const double X, const double Y, const double Z, matrix_row<Matrix> Row)
    {
        if (Row.size() == 3) {
            if (&rVariable == &DISPLACEMENT_X) {
                Row[0] = 1.0; Row[2] = -Y;
            } else if (&rVariable == &DISPLACEMENT_Y) {
                Row[1] = 1.0; Row[2] = X;
            } else if (&rVariable == &ROTATION_Z) {
                Row[2] = 1.0;
            }
        } else {
            if (&rVariable == &DISPLACEMENT_X) {
                Row[0] = 1.0; Row[4] = Z; Row[5] = -Y;
            } else if (&rVariable == &DISPLACEMENT_Y) {
                Row[1] = 1.0; Row[3] = -Z; Row[5] = X;
            } else if (&rVariable == &DISPLACEMENT_Z) {
                Row[2] = 1.0; Row[3] = Y; Row[4] = -X;
            } else if (&rVariable == &ROTATION_X) {
                Row[3] = 1.0;
            } else if (&rVariable == &ROTATION_Y) {
                Row[4] = 1.0;
            } else if (&rVariable == &ROTATION_Z) {
                Row[5] = 1.0;
            }
        }
    }

    /// Adds the outer product x x^T to the M x M block out (row-major)
    static void AddOuterProduct(Kratos::span<const double> X, Kratos::span<double> rOut)
    {
        const std::size_t m = X.size();
        for (std::size_t p = 0; p < m; ++p) {
            for (std::size_t q = 0; q < m; ++q) {
                rOut[p * m + q] += X[p] * X[q];
            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computes the pinned DOFs (mPinned / mPinnedDofs) from the current matrix A.
     * @details Works on the unmodified matrix. All buffers are reused between calls.
     */
    void ComputePinningMask(const MatrixType& rA)
    {
        const std::size_t n = rA.size1();
        const double global_k_norm = ComputeConnectedComponents(rA);
        AccumulateComponentGrams(rA);
        GroupDofsByComponent();

        mPinned.assign(n, 0);
        mPinnedDofs.clear();

        for (std::size_t c = 0; c < mNumComponents; ++c) {
            if (mComponentSize[c] == 0) {
                continue;
            }
            const double est_kc_norm = global_k_norm * std::sqrt(static_cast<double>(mComponentSize[c]) / static_cast<double>(n));
            const std::size_t r = ComputeUnconstrainedModes(c, est_kc_norm);
            if (r > 0) {
                GreedyPinComponent(c, r);
            }
        }
    }

    /**
     * @brief Splits the sparsity graph of A in connected components (union-find on nonzero entries).
     * @return The Frobenius norm of A
     */
    double ComputeConnectedComponents(const MatrixType& rA)
    {
        const std::size_t n = rA.size1();
        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        const auto& r_vals = rA.value_data();

        mParent.resize(n);
        std::iota(mParent.begin(), mParent.end(), std::size_t(0));

        double k_norm2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                if (r_vals[k] != 0.0) {
                    k_norm2 += r_vals[k] * r_vals[k];
                    Unite(mParent, i, r_col_idx[k]);
                }
            }
        }

        std::vector<std::size_t> root_to_component(n, npos);
        mDofComponent.resize(n);
        mNumComponents = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t root = Find(mParent, i);
            if (root_to_component[root] == npos) {
                root_to_component[root] = mNumComponents++;
            }
            mDofComponent[i] = root_to_component[root];
        }

        mComponentSize.assign(mNumComponents, 0);
        for (std::size_t i = 0; i < n; ++i) {
            ++mComponentSize[mDofComponent[i]];
        }

        return std::sqrt(k_norm2);
    }

    /**
     * @brief Accumulates per component H_c = sum_i R_i R_i^T (R_i = sum_j A_ij Z_j) and G_c = sum_i Z_i Z_i^T.
     * @details R is never stored: each row of A is contracted with Z on the fly.
     */
    void AccumulateComponentGrams(const MatrixType& rA)
    {
        const std::size_t m = mNumRigidBodyModes;
        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        const auto& r_vals = rA.value_data();

        ResetMatrix(mComponentH, mNumComponents * m, m);
        ResetMatrix(mComponentG, mNumComponents * m, m);
        mRowResidual.resize(m, false);

        for (std::size_t i = 0; i < rA.size1(); ++i) {
            const std::size_t c = mDofComponent[i];

            std::fill(mRowResidual.begin(), mRowResidual.end(), 0.0);
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                noalias(mRowResidual) += r_vals[k] * row(mRigidBodyModes, r_col_idx[k]);
            }

            AddOuterProduct(Kratos::span<const double>(mRowResidual.data().begin(), m), BlockSpan(mComponentH, c, m));
            AddOuterProduct(RowSpan(mRigidBodyModes, i), BlockSpan(mComponentG, c, m));
        }
    }

    /// Groups the DOF indices by component (counting sort), so each component can be visited directly
    void GroupDofsByComponent()
    {
        const std::size_t n = mDofComponent.size();
        mComponentOffsets.assign(mNumComponents + 1, 0);
        for (std::size_t c = 0; c < mNumComponents; ++c) {
            mComponentOffsets[c + 1] = mComponentOffsets[c] + mComponentSize[c];
        }

        mComponentDofs.resize(n);
        std::vector<std::size_t> cursor(mComponentOffsets.begin(), mComponentOffsets.end() - 1);
        for (std::size_t i = 0; i < n; ++i) {
            mComponentDofs[cursor[mDofComponent[i]]++] = i;
        }
    }

    /**
     * @brief Selects the unconstrained rigid body modes of component c.
     * @details Eigen-decomposes H_c and keeps the eigenvectors whose residual is negligible
     *          relative to the (estimated) stiffness scale of the component.
     *          The kept vectors are stored as rows of mUnconstrained (first r rows).
     * @return The number r of unconstrained modes
     */
    std::size_t ComputeUnconstrainedModes(const std::size_t c, const double EstKcNorm)
    {
        const std::size_t m = mNumRigidBodyModes;
        const double eps = std::numeric_limits<double>::epsilon();

        // eigen-decomposition of H_c
        mWorkH = subrange(mComponentH, c * m, (c + 1) * m, 0, m);
        JacobiEigenDecomposition(mWorkH, mLambda, mWorkV);

        mOrder.resize(m);
        std::iota(mOrder.begin(), mOrder.end(), std::size_t(0));
        std::sort(mOrder.begin(), mOrder.end(), [&](const std::size_t a, const std::size_t b) { return mLambda[a] < mLambda[b]; });

        mWorkG = subrange(mComponentG, c * m, (c + 1) * m, 0, m);

        mSelectedModes.clear();
        mMode.resize(m, false);
        for (const std::size_t k : mOrder) {
            mMode = column(mWorkV, k);

            const double res_norm = std::sqrt(std::max(mLambda[k], 0.0));
            const double mode_z_norm = std::sqrt(std::max(inner_prod(mMode, prod(mWorkG, mMode)), 0.0));
            const double scale = std::max(EstKcNorm * mode_z_norm, eps);

            if (res_norm / scale < mTolerance) {
                mSelectedModes.push_back(k);
            }
        }

        const std::size_t r = mSelectedModes.size();
        mUnconstrained.resize(r, m, false);
        for (std::size_t a = 0; a < r; ++a) {
            row(mUnconstrained, a) = column(mWorkV, mSelectedModes[a]);
        }
        return r;
    }

    /**
     * @brief Greedy pivoted Gram-Schmidt on the rows w_i = U^T Z_i of Z_eff.
     * @details Pins, in component c, the DOFs whose effective rows add the largest new direction
     *          until the rank r of Z_eff is exhausted. Equivalent to a pivoted QR on Z_eff^T.
     */
    void GreedyPinComponent(const std::size_t c, const std::size_t r)
    {
        const double rank_rel_tol = 1.0e-8;

        mTrialRow.resize(r, false);
        mBestRow.resize(r, false);
        mQBasis.resize(r, r, false);

        std::size_t n_q = 0;
        double max_norm2_0 = 0.0;
        while (n_q < r) {
            double best = -1.0;
            std::size_t best_dof = npos;
            for (std::size_t idx = mComponentOffsets[c]; idx < mComponentOffsets[c + 1]; ++idx) {
                const std::size_t i = mComponentDofs[idx];
                if (mPinned[i]) {
                    continue;
                }
                ComputeEffectiveRow(i, mTrialRow);
                ProjectOutSelected(mTrialRow, n_q);
                const double norm2 = inner_prod(mTrialRow, mTrialRow);
                if (norm2 > best) {
                    best = norm2;
                    best_dof = i;
                    mBestRow = mTrialRow;
                }
            }

            if (best_dof == npos) {
                break;
            }
            if (n_q == 0) {
                max_norm2_0 = best;
            }
            if (std::sqrt(std::max(best, 0.0)) <= rank_rel_tol * std::sqrt(max_norm2_0)) {
                break; // remaining rows are linearly dependent: rank reached
            }

            mPinned[best_dof] = 1;
            mPinnedDofs.push_back(best_dof);

            const double inv_norm = 1.0 / std::sqrt(best);
            for (std::size_t a = 0; a < r; ++a) {
                mQBasis(n_q, a) = mBestRow[a] * inv_norm;
            }
            ++n_q;
        }
    }

    /// w[a] = u_a . Z_i for the first r unconstrained modes u_a
    void ComputeEffectiveRow(const std::size_t I, Vector& rW) const
    {
        noalias(rW) = prod(mUnconstrained, row(mRigidBodyModes, I));
    }

    /// Removes from w the components along the first n_q orthonormal directions (two passes)
    void ProjectOutSelected(Vector& rW, const std::size_t n_q) const
    {
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t qk = 0; qk < n_q; ++qk) {
                const auto q_row = row(mQBasis, qk);
                noalias(rW) -= inner_prod(q_row, rW) * q_row;
            }
        }
    }

    /**
     * @brief Applies the pinning in place: zero the off-diagonal entries of pinned rows and columns
     *        (diagonal kept, so the pattern is unchanged), zero the pinned rhs and initial guess.
     */
    void ApplyPinning(MatrixType& rA, VectorType& rX, VectorType& rB) const
    {
        if (mPinnedDofs.empty()) {
            return;
        }

        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        auto& r_vals = rA.value_data();

        for (std::size_t i = 0; i < rA.size1(); ++i) {
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                const std::size_t j = r_col_idx[k];
                if (j != i && (mPinned[i] || mPinned[j])) {
                    r_vals[k] = 0.0;
                }
            }
        }

        for (const std::size_t p : mPinnedDofs) {
            rB[p] = 0.0;
            rX[p] = 0.0;
        }
    }

    /// Resizes a matrix to rows x cols and sets all its entries to zero (no reallocation if sizes match)
    static void ResetMatrix(Matrix& rMatrix, const std::size_t Rows, const std::size_t Cols)
    {
        rMatrix.resize(Rows, Cols, false);
        std::fill(rMatrix.data().begin(), rMatrix.data().begin() + Rows * Cols, 0.0);
    }

    /**
     * @brief Eigen-decomposition of a small symmetric matrix by cyclic Jacobi rotations.
     * @param rA Symmetric matrix (destroyed on exit)
     * @param rLambda Eigenvalues (not sorted)
     * @param rV Eigenvectors stored as columns
     */
    static void JacobiEigenDecomposition(Matrix& rA, Vector& rLambda, Matrix& rV)
    {
        const std::size_t m = rA.size1();
        constexpr std::size_t max_sweeps = 50;

        rV.resize(m, m, false);
        rLambda.resize(m, false);
        for (std::size_t p = 0; p < m; ++p) {
            for (std::size_t q = 0; q < m; ++q) {
                rV(p, q) = (p == q) ? 1.0 : 0.0;
            }
        }

        for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
            double off = 0.0;
            double diag = 0.0;
            for (std::size_t p = 0; p < m; ++p) {
                diag += rA(p, p) * rA(p, p);
                for (std::size_t q = p + 1; q < m; ++q) {
                    off += rA(p, q) * rA(p, q);
                }
            }
            if (off <= 1.0e-30 * (diag + 1.0e-300)) {
                break;
            }

            for (std::size_t p = 0; p < m; ++p) {
                for (std::size_t q = p + 1; q < m; ++q) {
                    JacobiRotation(rA, rV, p, q);
                }
            }
        }

        for (std::size_t k = 0; k < m; ++k) {
            rLambda[k] = rA(k, k);
        }
    }

    /// One Jacobi rotation annihilating the (p,q) entry of the symmetric matrix rA; rV accumulates the rotations
    static void JacobiRotation(Matrix& rA, Matrix& rV, const std::size_t P, const std::size_t Q)
    {
        const std::size_t m = rA.size1();
        const double apq = rA(P, Q);
        if (std::abs(apq) < 1.0e-300) {
            return;
        }

        const double theta = (rA(Q, Q) - rA(P, P)) / (2.0 * apq);
        const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;

        // A <- A J (columns P and Q)
        for (std::size_t k = 0; k < m; ++k) {
            const double akp = rA(k, P);
            const double akq = rA(k, Q);
            rA(k, P) = c * akp - s * akq;
            rA(k, Q) = s * akp + c * akq;
        }
        // A <- J^T A (rows P and Q)
        for (std::size_t k = 0; k < m; ++k) {
            const double apk = rA(P, k);
            const double aqk = rA(Q, k);
            rA(P, k) = c * apk - s * aqk;
            rA(Q, k) = s * apk + c * aqk;
        }
        // V <- V J
        for (std::size_t k = 0; k < m; ++k) {
            const double vkp = rV(k, P);
            const double vkq = rV(k, Q);
            rV(k, P) = c * vkp - s * vkq;
            rV(k, Q) = s * vkp + c * vkq;
        }
    }

    ///@}
    ///@name Member Variables
    ///@{

    typename BaseType::Pointer mpLinearSolver; /// Pointer to the internal linear solver
    double mTolerance = 1.0e-6;                /// Relative tolerance to detect unconstrained rigid body modes
    std::size_t mNumRigidBodyModes = 0;        /// 3 in 2D, 6 in 3D

    Matrix mRigidBodyModes;        /// Z (n x M), built in PrepareAdditionalData
    Matrix mComponentH;            /// Stacked (n_comp*M x M) blocks R_c^T R_c
    Matrix mComponentG;            /// Stacked (n_comp*M x M) blocks Z_c^T Z_c
    Matrix mUnconstrained;         /// Unconstrained modes of the current component, as rows (r x M)
    Matrix mWorkG;                 /// Work matrix: Z_c^T Z_c of the current component
    std::vector<std::size_t> mSelectedModes; /// Indices of the unconstrained eigenvectors of the current component
    Matrix mWorkH;                 /// Work matrix for the eigen-decomposition
    Matrix mWorkV;                 /// Eigenvectors of the work matrix
    Matrix mQBasis;                /// Orthonormal basis of the selected rows (greedy selection)
    Vector mLambda;                /// Eigenvalues
    Vector mRowResidual;           /// Work vector: R_i of the current row
    Vector mTrialRow;              /// Work vector: candidate effective row
    Vector mBestRow;               /// Work vector: best effective row of the current step
    Vector mMode;                  /// Work vector: eigenvector under test

    std::vector<std::size_t> mParent;           /// Union-find buffer (connected components)
    std::vector<std::size_t> mDofComponent;     /// Component id of each DOF
    std::vector<std::size_t> mComponentSize;    /// Number of DOFs per component
    std::vector<std::size_t> mComponentOffsets; /// CSR-like offsets of the DOFs grouped by component
    std::vector<std::size_t> mComponentDofs;    /// DOFs grouped by component
    std::vector<std::size_t> mOrder;            /// Sorting permutation of the eigenvalues
    std::vector<char> mPinned;                  /// Flag of pinned DOFs
    std::vector<std::size_t> mPinnedDofs;       /// List of pinned DOFs
    std::size_t mNumComponents = 0;             /// Number of connected components

    ///@}

}; // Class RigidBodyModePinningSolver

///@}

}  // namespace Kratos::Future.
