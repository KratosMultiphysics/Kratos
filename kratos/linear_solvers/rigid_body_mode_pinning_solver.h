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

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "input_output/logger.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
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
 *          DOFs which are pinned to zero. Pinning is applied IN PLACE on A and b in
 *          InitializeSolutionStep, before the wrapped solver factorizes: the off-diagonal
 *          entries of pinned rows/columns are set to zero and the diagonal is kept, so the
 *          sparsity pattern of A is never changed and no dense matrix of the size of A
 *          is built.
 *          After the wrapped solver performs the solve, the rigid body content of the
 *          unconstrained components is projected out of the solution
 *          (x <- x - Z_eff (Z_eff^T Z_eff)^-1 Z_eff^T x, with Z_eff = Z_c U_c the part of Z
 *          corresponding to the unconstrained modes of the component), so the returned
 *          displacement is the unique minimal-norm one, independent of the DOFs that
 *          were pinned.
 *          Number of rigid body modes: 3 in 2D (u_x, u_y, theta_z), 6 in 3D.
 * @details The Parameters constructor builds the wrapped solver from the nested
 *          "inner_solver_settings" parameter through the LinearSolverFactory
 *          (default {"solver_type": "skyline_lu_factorization"}); any solver registered in
 *          the factory may be used, and any inner-solver parameter may be set inside it.
 *          With "echo_level" > 0 the number of connected components and the number of pinned
 *          DOFs is printed; with "echo_level" > 1 the pinned DOFs (dof id and node id) are
 *          listed.
 * @note The rows of A are modified in place (as ScalingSolver does): the matrix must be
 *       re-assembled before the next call to InitializeSolutionStep.
 * @note The small 3x3/6x6 eigen-decomposition is a cyclic Jacobi routine, since Eigen is not
 *       part of the Kratos external libraries.
 * @tparam TSparseSpaceType The sparse space definition
 * @tparam TDenseSpaceType The dense space definition
 * @tparam TReordererType The reorder considered
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class RigidBodyModePinningSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidBodyModePinningSolver
    KRATOS_CLASS_POINTER_DEFINITION(RigidBodyModePinningSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// The definition of the spaces (sparse matrix)
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TSparseSpaceType::VectorType;

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType, TDenseSpaceType>;

    /// The index type definition to be consistent
    using IndexType = typename BaseType::IndexType;

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
        mEchoLevel = 0;
    }

    /**
     * @brief Constructor from parameters
     * @details The wrapped solver is selected through the nested "inner_solver_settings"
     *          parameter (default {"solver_type": "skyline_lu_factorization"}) and built with
     *          the LinearSolverFactory; any inner-solver parameter may be set inside it.
     * @param ThisParameters The configuration parameters
     */
    RigidBodyModePinningSolver(Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        KRATOS_TRY

        // Validate and assign default parameters
        Parameters default_parameters(R"(
        {
            "solver_type"           : "rigid_body_mode_pinning",
            "inner_solver_settings" : {
                "solver_type" : "skyline_lu_factorization"
            },
            "tolerance"             : 1.0e-6,
            "echo_level"            : 0
        }
        )");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mTolerance = ThisParameters["tolerance"].GetDouble();
        mEchoLevel = ThisParameters["echo_level"].GetInt();

        // Build the wrapped solver from the nested settings, as for the other factory-based
        // wrapper solvers (e.g. MonotonicityPreservingSolver): "inner_solver_settings" already
        // carries the inner "solver_type" and any inner-solver parameter, so it is delegated
        // to the factory as is
        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters["inner_solver_settings"]);

        KRATOS_CATCH("")
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
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->Initialize(rA, rX, rB);
    }

    /**
     * @brief Pins the rigid body modes of the unconstrained components in place and
     *        delegates the initialization of the solution step to the wrapped solver.
     * @details The pinning is applied before delegating, so a wrapped direct solver
     *          factorizes the already-pinned matrix in its own InitializeSolutionStep.
     * @param rA. System matrix (modified in place: pinned rows/columns are zeroed off-diagonal)
     * @param rX. Solution vector (pinned entries are set to zero)
     * @param rB. Right hand side vector (pinned entries are set to zero)
     */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (this->IsConsistent(rA, rX, rB)) {
            KRATOS_ERROR_IF(mRigidBodyModes.size1() != rA.size1() || mNumRigidBodyModes == 0)
                << "RigidBodyModePinningSolver: ProvideAdditionalData must be called before InitializeSolutionStep (node data missing)" << std::endl;

            // The work arrays are transient: they live only for the duration of this call
            PinningWork work;
            ComputePinningMask(rA, work);
            ApplyPinning(rA, rX, rB, work);

            // Echo the pinning summary
            KRATOS_INFO_IF("RigidBodyModePinningSolver", mEchoLevel > 0)
                << "Pinned " << work.num_pinned << " dofs of " << work.component_size.size()
                << " connected components (the unconstrained rigid body modes were removed)." << std::endl;
            if (mEchoLevel > 1) {
                std::stringstream buffer;
                for (std::size_t p = 0; p < work.pinned.size(); ++p) {
                    if (work.pinned[p]) {
                        buffer << "  dof " << p << " (node " << mDofNodeId[p] << ")" << '\n';
                    }
                }
                KRATOS_INFO("RigidBodyModePinningSolver") << "Pinned dofs (dof id, node id):\n" << buffer.str() << std::endl;
            }
        }

        mpLinearSolver->InitializeSolutionStep(rA, rX, rB);
    }

    /**
     * @brief Delegates the solution of the system to the wrapped solver and removes the
     *        rigid body content of the solution.
     * @details After a successful solve, the rigid body content of the unconstrained
     *          components is projected out of the solution, so the returned solution is
     *          the unique minimal-norm one, independent of the DOFs that were pinned.
     */
    bool PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        const bool success = mpLinearSolver->PerformSolutionStep(rA, rX, rB);
        if (success && mNumRigidBodyModes != 0) {
            FilterRigidBodyContent(rX);
        }
        return success;
    }

    /// Delegates the finalization of the solution step to the wrapped solver
    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA, rX, rB);
    }

    /// Delegates the cleanup to the wrapped solver
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
     * @brief Builds the rigid body mode matrix Z (n x M) from the node coordinates and DOF types.
     * @details M = 3 in 2D and M = 6 in 3D (see FillRigidBodyModeRow for the column layout).
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        mpLinearSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_ERROR_IF_NOT(domain_size == 2 || domain_size == 3)
            << "RigidBodyModePinningSolver: DOMAIN_SIZE must be 2 or 3 (found " << domain_size << ")" << std::endl;
        mNumRigidBodyModes = (domain_size == 2) ? 3 : 6;

        mRigidBodyModes.resize(rA.size1(), mNumRigidBodyModes, false);
        mRigidBodyModes.clear();

        mDofNodeId.resize(rA.size1(), false);
        std::fill(mDofNodeId.data().begin(), mDofNodeId.data().begin() + rA.size1(), std::size_t(-1));

        for (const auto& r_node : r_model_part.Nodes()) {
            for (const auto& p_dof : r_node.GetDofs()) {
                const std::size_t eq_id = p_dof->EquationId();
                if (eq_id < rA.size1()) {
                    mDofNodeId[eq_id] = r_node.Id();
                    FillRigidBodyModeRow(p_dof->GetVariable(), r_node.X0(), r_node.Y0(), r_node.Z0(), row(mRigidBodyModes, eq_id));
                }
            }
        }
    }

    /// The rigid body modes need the spatial position of the nodes
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{

    void SetTolerance(double NewTolerance) override
    {
        mTolerance = NewTolerance;
    }

    double GetTolerance() override
    {
        return mTolerance;
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
    ///@name Type Definitions
    ///@{

    /// Transient work arrays for one pinning computation: allocated in InitializeSolutionStep
    /// and living only for the duration of that call
    struct PinningWork
    {
        DenseVector<std::size_t> dof_component;     /// Component id of each DOF
        DenseVector<std::size_t> component_size;    /// Number of DOFs per component
        Matrix component_h;                         /// Stacked (n_comp*M x M) blocks R_c^T R_c
        DenseVector<char> pinned;                   /// Flag of pinned DOFs
        std::size_t num_pinned = 0;                 /// Number of pinned DOFs
    };

    ///@}
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

    static std::size_t Find(DenseVector<std::size_t>& rParent, std::size_t i)
    {
        while (rParent[i] != i) {
            rParent[i] = rParent[rParent[i]]; // path halving
            i = rParent[i];
        }
        return i;
    }

    static void Unite(DenseVector<std::size_t>& rParent, const std::size_t I, const std::size_t J)
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
            if (rVariable == DISPLACEMENT_X) {
                Row[0] = 1.0; Row[2] = -Y;
            } else if (rVariable == DISPLACEMENT_Y) {
                Row[1] = 1.0; Row[2] = X;
            } else if (rVariable == ROTATION_Z) {
                Row[2] = 1.0;
            }
        } else {
            if (rVariable == DISPLACEMENT_X) {
                Row[0] = 1.0; Row[4] = Z; Row[5] = -Y;
            } else if (rVariable == DISPLACEMENT_Y) {
                Row[1] = 1.0; Row[3] = -Z; Row[5] = X;
            } else if (rVariable == DISPLACEMENT_Z) {
                Row[2] = 1.0; Row[3] = Y; Row[4] = -X;
            } else if (rVariable == ROTATION_X) {
                Row[3] = 1.0;
            } else if (rVariable == ROTATION_Y) {
                Row[4] = 1.0;
            } else if (rVariable == ROTATION_Z) {
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
     * @brief Computes the pinned DOFs (rWork.pinned) from the current matrix A.
     * @details Works on the unmodified matrix; all work arrays are provided by the caller
     *          through rWork.
     */
    void ComputePinningMask(const SparseMatrixType& rA, PinningWork& rWork)
    {
        const std::size_t n = rA.size1();
        const double global_k_norm = ComputeConnectedComponents(rA, rWork);
        AccumulateComponentGrams(rA, rWork);
        GroupDofsByComponent(rWork);

        rWork.pinned.resize(n, false);
        rWork.pinned.clear();
        rWork.num_pinned = 0;

        // The unconstrained modes of all components are stacked compactly in mComponentModes
        // (at most m rows per component); the offsets record where each component's rows start
        const std::size_t num_components = rWork.component_size.size();
        const std::size_t m = mNumRigidBodyModes;
        mComponentModes.resize(num_components * m, m, false);
        mComponentModeOffsets.resize(num_components + 1, false);
        mComponentModeOffsets[0] = 0;
        ResetMatrix(mComponentModeGrams, num_components * m, m);

        Matrix unconstrained; /// Unconstrained modes of the current component, as rows (r x M)
        for (std::size_t c = 0; c < num_components; ++c) {
            std::size_t r = 0;
            if (rWork.component_size[c] > 0) {
                const double est_kc_norm = global_k_norm * std::sqrt(static_cast<double>(rWork.component_size[c]) / static_cast<double>(n));
                r = ComputeUnconstrainedModes(c, est_kc_norm, rWork, unconstrained);
                const std::size_t base = mComponentModeOffsets[c];
                for (std::size_t a = 0; a < r; ++a) {
                    row(mComponentModes, base + a) = row(unconstrained, a);
                }
                if (r > 0) {
                    GreedyPinComponent(c, r, unconstrained, rWork);
                    // Precompute the Gram Z_eff^T Z_eff = U^T G_c U of the unconstrained modes
                    // (symmetric: the upper triangle is formed and mirrored) and store it in the
                    // M x M block of the component, so the rigid body filter in
                    // PerformSolutionStep does not have to re-form it
                    for (std::size_t a = 0; a < r; ++a) {
                        for (std::size_t b = a; b < r; ++b) {
                            double s = 0.0;
                            for (std::size_t p = 0; p < m; ++p) {
                                for (std::size_t q = 0; q < m; ++q) {
                                    s += unconstrained(a, p) * mComponentG(c * m + p, q) * unconstrained(b, q);
                                }
                            }
                            mComponentModeGrams(c * m + a, b) = s;
                            if (a != b) {
                                mComponentModeGrams(c * m + b, a) = s;
                            }
                        }
                    }
                }
            }
            mComponentModeOffsets[c + 1] = mComponentModeOffsets[c] + r;
        }
    }

    /**
     * @brief Splits the sparsity graph of A in connected components (union-find on nonzero entries).
     * @return The Frobenius norm of A
     */
    double ComputeConnectedComponents(const SparseMatrixType& rA, PinningWork& rWork)
    {
        const std::size_t n = rA.size1();
        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        const auto& r_vals = rA.value_data();

        DenseVector<std::size_t> parent(n);
        std::iota(parent.begin(), parent.end(), std::size_t(0));

        double k_norm2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                if (r_vals[k] != 0.0) {
                    k_norm2 += r_vals[k] * r_vals[k];
                    Unite(parent, i, r_col_idx[k]);
                }
            }
        }

        DenseVector<std::size_t> root_to_component(n, npos);
        rWork.dof_component.resize(n, false);
        std::size_t num_components = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t root = Find(parent, i);
            if (root_to_component[root] == npos) {
                root_to_component[root] = num_components++;
            }
            rWork.dof_component[i] = root_to_component[root];
        }

        rWork.component_size.resize(num_components, false);
        rWork.component_size.clear();
        for (std::size_t i = 0; i < n; ++i) {
            ++rWork.component_size[rWork.dof_component[i]];
        }

        return std::sqrt(k_norm2);
    }

    /**
     * @brief Accumulates per component H_c = sum_i R_i R_i^T (R_i = sum_j A_ij Z_j) and G_c = sum_i Z_i Z_i^T.
     * @details R is never stored: each row of A is contracted with Z on the fly. The G_c blocks
     *          are kept in mComponentG, since they are needed for the selection of the
     *          unconstrained modes and for the precomputation of their Gram in ComputePinningMask.
     */
    void AccumulateComponentGrams(const SparseMatrixType& rA, PinningWork& rWork)
    {
        const std::size_t m = mNumRigidBodyModes;
        const std::size_t num_components = rWork.component_size.size();
        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        const auto& r_vals = rA.value_data();

        ResetMatrix(rWork.component_h, num_components * m, m);
        ResetMatrix(mComponentG, num_components * m, m);
        Vector row_residual(m);

        for (std::size_t i = 0; i < rA.size1(); ++i) {
            const std::size_t c = rWork.dof_component[i];

            row_residual.clear();
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                noalias(row_residual) += r_vals[k] * row(mRigidBodyModes, r_col_idx[k]);
            }

            AddOuterProduct(Kratos::span<const double>(row_residual.data().begin(), m), BlockSpan(rWork.component_h, c, m));
            AddOuterProduct(RowSpan(mRigidBodyModes, i), BlockSpan(mComponentG, c, m));
        }
    }

    /// Groups the DOF indices by component (counting sort), so each component can be visited directly
    void GroupDofsByComponent(const PinningWork& rWork)
    {
        const std::size_t n = rWork.dof_component.size();
        const std::size_t num_components = rWork.component_size.size();

        mComponentOffsets.resize(num_components + 1, false);
        mComponentOffsets[0] = 0;
        for (std::size_t c = 0; c < num_components; ++c) {
            mComponentOffsets[c + 1] = mComponentOffsets[c] + rWork.component_size[c];
        }

        mComponentDofs.resize(n, false);
        DenseVector<std::size_t> cursor(num_components);
        for (std::size_t c = 0; c < num_components; ++c) {
            cursor[c] = mComponentOffsets[c];
        }
        for (std::size_t i = 0; i < n; ++i) {
            mComponentDofs[cursor[rWork.dof_component[i]]++] = i;
        }
    }

    /**
     * @brief Selects the unconstrained rigid body modes of component c.
     * @details Eigen-decomposes H_c and keeps the eigenvectors whose residual is negligible
     *          relative to the (estimated) stiffness scale of the component.
     *          The kept vectors are stored as rows of rUnconstrained (r rows).
     * @return The number r of unconstrained modes
     */
    std::size_t ComputeUnconstrainedModes(
        const std::size_t c,
        const double EstKcNorm,
        const PinningWork& rWork,
        Matrix& rUnconstrained)
    {
        const std::size_t m = mNumRigidBodyModes;
        const double eps = std::numeric_limits<double>::epsilon();

        // eigen-decomposition of H_c
        Matrix work_h;
        work_h = subrange(rWork.component_h, c * m, (c + 1) * m, 0, m);
        Vector lambda;
        Matrix work_v;
        JacobiEigenDecomposition(work_h, lambda, work_v);

        DenseVector<std::size_t> order(m);
        std::iota(order.begin(), order.end(), std::size_t(0));
        std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) { return lambda[a] < lambda[b]; });

        Matrix work_g;
        work_g = subrange(mComponentG, c * m, (c + 1) * m, 0, m);

        // true for the eigenvectors with (numerically) zero residual
        const auto is_unconstrained = [&](const std::size_t k)
        {
            Vector mode;
            mode = column(work_v, k);
            const double res_norm = std::sqrt(std::max(lambda[k], 0.0));
            const double mode_z_norm = std::sqrt(std::max(inner_prod(mode, prod(work_g, mode)), 0.0));
            const double scale = std::max(EstKcNorm * mode_z_norm, eps);
            return res_norm / scale < mTolerance;
        };

        std::size_t r = 0;
        for (const std::size_t k : order) {
            if (is_unconstrained(k)) {
                ++r;
            }
        }

        rUnconstrained.resize(r, m, false);
        std::size_t a = 0;
        for (std::size_t k = 0; k < m; ++k) {
            if (is_unconstrained(k)) {
                row(rUnconstrained, a) = column(work_v, k);
                ++a;
            }
        }
        return r;
    }

    /**
     * @brief Greedy pivoted Gram-Schmidt on the rows w_i = U^T Z_i of Z_eff.
     * @details Pins, in component c, the DOFs whose effective rows add the largest new direction
     *          until the rank r of Z_eff is exhausted. Equivalent to a pivoted QR on Z_eff^T.
     */
    void GreedyPinComponent(const std::size_t c, const std::size_t r, const Matrix& rUnconstrained, PinningWork& rWork)
    {
        const double rank_rel_tol = 1.0e-8;

        Vector trial_row(r);
        Vector best_row(r);
        Matrix q_basis(r, r, false);

        std::size_t n_q = 0;
        double max_norm2_0 = 0.0;
        while (n_q < r) {
            double best = -1.0;
            std::size_t best_dof = npos;
            for (std::size_t idx = mComponentOffsets[c]; idx < mComponentOffsets[c + 1]; ++idx) {
                const std::size_t i = mComponentDofs[idx];
                if (rWork.pinned[i]) {
                    continue;
                }
                ComputeEffectiveRow(i, trial_row, rUnconstrained);
                ProjectOutSelected(trial_row, n_q, q_basis);
                const double norm2 = inner_prod(trial_row, trial_row);
                if (norm2 > best) {
                    best = norm2;
                    best_dof = i;
                    best_row = trial_row;
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

            rWork.pinned[best_dof] = 1;
            ++rWork.num_pinned;

            const double inv_norm = 1.0 / std::sqrt(best);
            for (std::size_t a = 0; a < r; ++a) {
                q_basis(n_q, a) = best_row[a] * inv_norm;
            }
            ++n_q;
        }
    }

    /// w[a] = u_a . Z_i for the first r unconstrained modes u_a
    void ComputeEffectiveRow(const std::size_t I, Vector& rW, const Matrix& rUnconstrained) const
    {
        noalias(rW) = prod(rUnconstrained, row(mRigidBodyModes, I));
    }

    /// Removes from w the components along the first n_q orthonormal directions (two passes)
    void ProjectOutSelected(Vector& rW, const std::size_t n_q, const Matrix& rQBasis) const
    {
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t qk = 0; qk < n_q; ++qk) {
                const auto q_row = row(rQBasis, qk);
                noalias(rW) -= inner_prod(q_row, rW) * q_row;
            }
        }
    }

    /**
     * @brief Applies the pinning in place: zero the off-diagonal entries of pinned rows and columns
     *        (diagonal kept, so the pattern is unchanged), zero the pinned rhs and initial guess.
     */
    void ApplyPinning(SparseMatrixType& rA, VectorType& rX, VectorType& rB, const PinningWork& rWork) const
    {
        if (rWork.num_pinned == 0) {
            return;
        }

        const auto& r_row_ptr = rA.index1_data();
        const auto& r_col_idx = rA.index2_data();
        auto& r_vals = rA.value_data();

        for (std::size_t i = 0; i < rA.size1(); ++i) {
            for (std::size_t k = r_row_ptr[i]; k < r_row_ptr[i + 1]; ++k) {
                const std::size_t j = r_col_idx[k];
                if (j != i && (rWork.pinned[i] || rWork.pinned[j])) {
                    r_vals[k] = 0.0;
                }
            }
        }

        for (std::size_t p = 0; p < rWork.pinned.size(); ++p) {
            if (rWork.pinned[p]) {
                rB[p] = 0.0;
                rX[p] = 0.0;
            }
        }
    }

    /**
     * @brief Removes the rigid body content of the unconstrained components from rX.
     * @details For each component c with r unconstrained modes u_a (rows of mComponentModes)
     *          the effective mode matrix Z_eff = Z_c U_c (n_c x r) spans the null space of the
     *          pinned system on that component, so the projection
     *              x <- x - Z_eff (Z_eff^T Z_eff)^-1 Z_eff^T x
     *          removes the rigid body content and the result no longer depends on the DOFs
     *          that were pinned. Z_eff is never stored: the Gram Z_eff^T Z_eff = U^T G_c U is
     *          precomputed at pinning time (from the per-component Gram G_c = Z_c^T Z_c and
     *          the stored modes) and kept in mComponentModeGrams, so this filter only copies
     *          it into a scratch matrix (the pivoted elimination overwrites it); the right
     *          hand side Z_eff^T x = U^T (Z_c^T x) is accumulated row by row over the DOFs
     *          of the component, and the r x r system (r <= M) is solved by pivoted
     *          elimination before the correction is applied in a second pass. The small work
     *          arrays are allocated once and reused across components.
     */
    void FilterRigidBodyContent(VectorType& rX) const
    {
        const std::size_t M = mNumRigidBodyModes;
        if (mComponentModeOffsets.size() < 2 || mRigidBodyModes.size1() != rX.size()) {
            return; // no pinning was computed for this system
        }

        Matrix gram(M, M, false);     /// Z_eff^T Z_eff of the current component (destroyed by the elimination)
        Vector z_cx(M, false);        /// Z_c^T x of the current component
        Vector rhs(M, false);         /// Z_eff^T x
        Vector gamma(M, false);       /// (Z_eff^T Z_eff)^-1 Z_eff^T x
        Vector mode_combo(M, false);  /// M-vector combination of the removed rigid displacement

        for (std::size_t c = 0; c + 1 < mComponentModeOffsets.size(); ++c) {
            const std::size_t r = mComponentModeOffsets[c + 1] - mComponentModeOffsets[c];
            if (r == 0) {
                continue;
            }
            const std::size_t mode_base = mComponentModeOffsets[c];

            // Copy the gram = Z_eff^T Z_eff precomputed at pinning time into the scratch matrix:
            // the pivoted elimination overwrites its matrix argument
            for (std::size_t a = 0; a < r; ++a) {
                for (std::size_t b = 0; b < r; ++b) {
                    gram(a, b) = mComponentModeGrams(c * M + a, b);
                }
            }

            // pass 1: rhs_a = u_a . (Z_c^T x), with Z_c^T x accumulated row by row
            z_cx.clear();
            for (std::size_t idx = mComponentOffsets[c]; idx < mComponentOffsets[c + 1]; ++idx) {
                const std::size_t i = mComponentDofs[idx];
                const double x_i = rX[i];
                for (std::size_t b = 0; b < M; ++b) {
                    z_cx[b] += x_i * mRigidBodyModes(i, b);
                }
            }
            for (std::size_t a = 0; a < r; ++a) {
                double s = 0.0;
                for (std::size_t p = 0; p < M; ++p) {
                    s += mComponentModes(mode_base + a, p) * z_cx[p];
                }
                rhs[a] = s;
            }

            // solve the small system gram * gamma = rhs
            const std::size_t rank = SolveSmallSystem(gram, rhs, gamma, r);

            // mode_combo = U^T gamma: the mode combination of the removed rigid displacement
            mode_combo.clear();
            for (std::size_t a = 0; a < rank; ++a) {
                for (std::size_t b = 0; b < M; ++b) {
                    mode_combo[b] += mComponentModes(mode_base + a, b) * gamma[a];
                }
            }

            // pass 2: remove Z_i . mode_combo from the solution of the DOFs of c
            for (std::size_t idx = mComponentOffsets[c]; idx < mComponentOffsets[c + 1]; ++idx) {
                const std::size_t i = mComponentDofs[idx];
                double s = 0.0;
                for (std::size_t b = 0; b < M; ++b) {
                    s += mRigidBodyModes(i, b) * mode_combo[b];
                }
                rX[i] -= s;
            }
        }
    }

    /**
     * @brief Solves the small dense system G * x = b by partial-pivoting Gaussian elimination.
     * @details G (r x r) and b (r x 1) are overwritten. If G is rank deficient the
     *          elimination stops and the remaining components of x are set to zero.
     * @return The number of pivots found
     */
    static std::size_t SolveSmallSystem(Matrix& rG, Vector& rB, Vector& rX, const std::size_t r)
    {
        double scale = 0.0;
        for (std::size_t p = 0; p < r; ++p) {
            for (std::size_t q = 0; q < r; ++q) {
                scale = std::max(scale, std::abs(rG(p, q)));
            }
        }
        const double rank_tol = 1.0e-12 * scale;

        std::size_t rank = 0;
        for (std::size_t col = 0; col < r; ++col) {
            std::size_t pivot = col;
            double best = std::abs(rG(col, col));
            for (std::size_t k = col + 1; k < r; ++k) {
                const double candidate = std::abs(rG(k, col));
                if (candidate > best) {
                    best = candidate;
                    pivot = k;
                }
            }
            if (best <= rank_tol) {
                break; // G is rank deficient from this column on
            }

            if (pivot != col) {
                for (std::size_t q = 0; q < r; ++q) {
                    std::swap(rG(col, q), rG(pivot, q));
                }
                std::swap(rB[col], rB[pivot]);
            }

            const double inv_pivot = 1.0 / rG(col, col);
            for (std::size_t k = col + 1; k < r; ++k) {
                const double factor = rG(k, col) * inv_pivot;
                if (factor != 0.0) {
                    for (std::size_t q = col + 1; q < r; ++q) {
                        rG(k, q) -= factor * rG(col, q);
                    }
                    rB[k] -= factor * rB[col];
                }
                rG(k, col) = 0.0;
            }
            ++rank;
        }

        for (std::size_t k = r; k > 0; --k) {
            const std::size_t i = k - 1;
            if (i < rank) {
                double s = rB[i];
                for (std::size_t q = i + 1; q < rank; ++q) {
                    s -= rG(i, q) * rX[q];
                }
                rX[i] = s / rG(i, i);
            } else {
                rX[i] = 0.0;
            }
        }
        return rank;
    }

    /// Resizes a matrix to rows x cols and sets all its entries to zero (no reallocation if sizes match)
    static void ResetMatrix(Matrix& rMatrix, const std::size_t Rows, const std::size_t Cols)
    {
        rMatrix.resize(Rows, Cols, false);
        rMatrix.clear();
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
    int mEchoLevel = 0;                       /// Echo level (0: silent, 1: pinning summary, 2: pinned dof list)
    std::size_t mNumRigidBodyModes = 0;        /// 3 in 2D, 6 in 3D
    Matrix mRigidBodyModes;                    /// Z (n x M), built in ProvideAdditionalData
    DenseVector<std::size_t> mDofNodeId;       /// Node id of each DOF (std::size_t(-1) if unknown), built in ProvideAdditionalData

    /// Pinning state of the current system: filled in InitializeSolutionStep, consumed by the
    /// rigid body filter in PerformSolutionStep
    DenseVector<std::size_t> mComponentDofs;        /// DOFs grouped by component
    DenseVector<std::size_t> mComponentOffsets;     /// CSR-like offsets of the DOFs grouped by component
    Matrix mComponentModes;                         /// Unconstrained modes of all components, stacked compactly as rows
    DenseVector<std::size_t> mComponentModeOffsets; /// Row offsets of the per-component modes inside mComponentModes
    Matrix mComponentG;                             /// Stacked (n_comp*M x M) blocks Z_c^T Z_c
    Matrix mComponentModeGrams;                     /// Stacked (n_comp*M x M) blocks: the top-left r x r holds the Gram Z_eff^T Z_eff of the component's unconstrained modes, precomputed at pinning time

    ///@}

}; // Class RigidBodyModePinningSolver

///@}

}  // namespace Kratos.
