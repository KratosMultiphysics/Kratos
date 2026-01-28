/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Authors: Armin Geiser
//           Máté Kelemen
*/

#pragma once

// External includes
#include <Eigen/Core>
#include <Spectra/SymGEigsShiftSolver.h>

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_utilities/ublas_wrapper.h"
#include "utilities/builtin_timer.h"
#include "factories/linear_solver_factory.h"
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp"

namespace Kratos {

template<
    class TSparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>,
    class TDenseSpaceType = UblasSpace<double, Matrix, Vector>,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SpectraSymGEigsShiftSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

    typename LinearSolver<
        TSparseSpaceType,
        TDenseSpaceType
    >::Pointer mpSolver;

public:
    KRATOS_CLASS_POINTER_DEFINITION(SpectraSymGEigsShiftSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    SpectraSymGEigsShiftSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "spectra_sym_g_eigs_shift",
            "number_of_eigenvalues": 1,
            "normalize_eigenvectors": false,
            "max_iteration": 1000,
            "shift": 0.0,
            "echo_level": 1,
            "linear_solver_settings" : {
                "solver_type" : "sparse_lu"
            }
        })");

        // In descending order of preference, the following
        // linear solvers are defaulted:
        // - Intel Pardiso (requires MKL)
        // - CHOLMOD (requires SuiteSparse)
        // - Eigen LU (always available in the LinearSolversApplication)
        using SolverRegistry = KratosComponents<
            LinearSolverFactory<
                TSparseSpaceType,
                TDenseSpaceType
            >
        >;
        if (SolverRegistry::Has("pardiso_llt")) {
            default_params["linear_solver_settings"]["solver_type"].SetString("pardiso_llt");
        } else if (SolverRegistry::Has("cholmod")) {
            default_params["linear_solver_settings"]["solver_type"].SetString("cholmod");
        }

        mParam.ValidateAndAssignDefaults(default_params);
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~SpectraSymGEigsShiftSolver() override {}

    /**
     * Solve the generalized eigenvalue problem of symmetric matrices K and M
     */
    void Solve(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        VectorType& rEigenvalues,
        DenseMatrixType& rEigenvectors) override
    {
        using scalar_t = double;
        using vector_t = Kratos::EigenDynamicVector<scalar_t>;
        using matrix_t = Kratos::EigenDynamicMatrix<scalar_t>;

        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        const int echo_level = mParam["echo_level"].GetInt();
        const double shift = mParam["shift"].GetDouble();

        // --- wrap ublas matrices

        UblasWrapper<scalar_t> b_wrapper(rM);
        const auto& b = b_wrapper.matrix();

        // --- timer
        const auto timer = BuiltinTimer();

        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        using OpType = OwnSymShiftInvert<scalar_t>;
        using BOpType = OwnSparseSymMatProd<scalar_t>;

        OpType op(rK, rM, mParam["linear_solver_settings"]);
        BOpType Bop(b);
        const int ncv = 3 * nroot;  // TODO find a good value
        Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert> eigs(
            op,
            Bop,
            nroot,
            ncv,
            shift);

        eigs.init();
        const int nconv = eigs.compute(Spectra::SortRule::LargestAlge, max_iteration);
        const int niter = eigs.num_iterations();
        const int nops = eigs.num_operations();

        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "nconv = " << nconv << std::endl;
        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "niter = " << niter << std::endl;
        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "nops  = " << nops << std::endl;

        KRATOS_ERROR_IF(eigs.info() != Spectra::CompInfo::Successful) << "SpectraSymGEigsShiftSolver: failed" << std::endl;

        if (static_cast<int>(rEigenvalues.size()) != nroot) {
            rEigenvalues.resize(nroot);
        }
        if (static_cast<int>(rEigenvectors.size1()) != nroot || rEigenvectors.size2() != rK.size1()) {
            rEigenvectors.resize(nroot, rK.size1());
        }

        Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());

        // Spectra::SortRule::LargestAlge results in values being in descending order
        eigvals = eigs.eigenvalues().reverse();
        // conversion back to RowMajor ordering of Kratos and transpose because Spectra returns eigvecs as columns
        eigvecs = eigs.eigenvectors().transpose().colwise().reverse();

        // --- normalization
        // TODO seems to be normalized by Spectra already -> confirm!
        // Given generalized eigenvalue problem (A - eigenvalue * B) * eigenvector = 0,
        // eigenvector is normalized such that eigenvector^T * B * eigenvector = 1
        if(mParam["normalize_eigenvectors"].GetBool())
        {
            for (int i = 0; i != nroot; ++i)
            {
                const double tmp = eigvecs.row(i) * b_wrapper.matrix() * eigvecs.row(i).transpose();
                const double factor = 1.0 / std::sqrt(tmp);
                eigvecs.row(i) *=  factor;
                KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
            }
        }

        // --- output
        if (echo_level > 0) {
            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO("SpectraSymGEigsShiftSolver:") << "Completed in " << timer << std::endl
                      << "                   Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "SpectraSymGEigsShiftSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

private:

    ///
    /// \ingroup MatOp
    ///
    /// This class defines the matrix-vector multiplication operation on a
    /// sparse real symmetric matrix \f$A\f$, i.e., calculating \f$y=Ax\f$ for any vector
    /// \f$x\f$. This is adapted from Spectra::OwnSparseSymMatProd
    ///
    template <typename Scalar_>
    class OwnSparseSymMatProd
    {
    public:
        ///
        /// Element type of the matrix.
        ///
        using Scalar = Scalar_;

    private:
        using Index = Eigen::Index;
        using Vector = Kratos::EigenDynamicVector<Scalar>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;
        using Matrix = Kratos::EigenDynamicMatrix<Scalar>;
        using SparseMatrix = Kratos::EigenSparseMatrix<Scalar>;
        using ConstGenericSparseMatrix = const Eigen::Ref<const SparseMatrix>;

        ConstGenericSparseMatrix m_mat;

    public:
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param mat An **Eigen** sparse matrix object, whose type can be
        /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
        /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
        ///
        OwnSparseSymMatProd(ConstGenericSparseMatrix& mat) :
            m_mat(mat)
        {}

        ///
        /// Return the number of rows of the underlying matrix.
        ///
        Index rows() const { return m_mat.rows(); }
        ///
        /// Return the number of columns of the underlying matrix.
        ///
        Index cols() const { return m_mat.cols(); }

        /// @brief Perform the matrix-vector multiplication operation \f$y=Ax\f$.
        /// @param x_in  Pointer to the \f$x\f$ vector.
        /// @param y_out Pointer to the \f$y\f$ vector.
        // y_out = A * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            block_for_each(y_out, y_out + m_mat.rows(), [](Scalar& rOut) {rOut = static_cast<Scalar>(0);});
            BalancedProduct(
                m_mat.outerIndexPtr(),
                m_mat.outerIndexPtr() + m_mat.outerSize() + 1,
                m_mat.innerIndexPtr(),
                m_mat.valuePtr(),
                x_in,
                y_out);
        }

        ///
        /// Perform the matrix-matrix multiplication operation \f$y=Ax\f$.
        ///
        Matrix operator*(const Eigen::Ref<const Matrix>& mat_in) const
        {
            return m_mat * mat_in;
        }

        ///
        /// Extract (i,j) element of the underlying matrix.
        ///
        Scalar operator()(Index i, Index j) const
        {
            return m_mat.coeff(i, j);
        }
    };



    ///
    /// \ingroup MatOp
    ///
    /// This class defines matrix operations required by the generalized eigen solver
    /// in the shift-and-invert mode. Given two symmetric matrices \f$A\f$ and \f$B\f$,
    /// it solves the linear equation \f$y=(A-\sigma B)^{-1}x\f$, where \f$\sigma\f$ is a real shift.
    /// This is adapted from Spectra::OwnSymShiftInvert
    ///
    /// This class is intended to be used with the SymGEigsShiftSolver generalized eigen solver.
    ///
    /// \tparam TScalar        The element type of the matrices.
    ///                        Currently supported types are `float`, `double`, and `long double`.
    template <typename TScalar>
    class OwnSymShiftInvert
    {
    public:
        ///
        /// Element type of the matrix.
        ///
        using Scalar = TScalar;

    private:
        using TSparseMatrix = typename TUblasSparseSpace<TScalar>::MatrixType;

        using TVector = typename TUblasSparseSpace<TScalar>::VectorType;

        TSparseMatrix* mpTarget;

        const TSparseMatrix* mpB;

        mutable TVector mDummyRhs;

        mutable TVector mDummySolution;

        TScalar mLastShift;

        typename LinearSolver<
            TUblasSparseSpace<TScalar>,
            TUblasDenseSpace<TScalar>
        >::Pointer mpSolver;

    public:
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param rA A sparse matrix object, whose type can be
        ///           `Eigen::SparseMatrix<...>`,
        ///           `Eigen::Map<Eigen::SparseMatrix<...>>`,
        ///           `Eigen::Ref<Eigen::SparseMatrix<...>>`, etc.
        /// \param rB A sparse matrix object.
        ///
        OwnSymShiftInvert(TSparseMatrix& rA,
                          const TSparseMatrix& rB,
                          Parameters Settings)
            : mpTarget(&rA),
              mpB(&rB),
              mDummyRhs(rA.size1()),
              mDummySolution(rA.size1()),
              mLastShift(static_cast<TScalar>(0)),
              mpSolver(nullptr)
        {
            const std::size_t system_size = rA.size1();
            KRATOS_ERROR_IF(system_size != rA.size2() || system_size != rB.size1() || system_size != rB.size2())
                << "SymShiftInvert: A and B must be square matrices of the same size";

            MergeMatrices<TScalar>(*mpTarget, *mpB);
            TUblasSparseSpace<TScalar>::SetToZero(mDummyRhs);
            //TUblasSparseSpace<TScalar>::SetToZero(mDummySolution);

            using SolverRegistry = KratosComponents<LinearSolverFactory<
                TUblasSparseSpace<TScalar>,
                TUblasDenseSpace<TScalar>
            >>;
            const std::string solver_name = Settings["solver_type"].GetString();
            KRATOS_ERROR_IF_NOT(SolverRegistry::Has(solver_name))
                << "'" << solver_name << "' does not name a registered linear solver. Options are:\n"
                << []() -> std::string {
                    std::stringstream message;
                    for (const auto& r_pair : SolverRegistry::GetComponents())
                        message << "\t" << r_pair.first << "\n";
                    return message.str();
                }();

            KRATOS_TRY
            mpSolver = SolverRegistry::Get(solver_name).Create(Settings);
            KRATOS_CATCH("")
        }

        /// Return the number of rows of the underlying matrix.
        std::ptrdiff_t rows() const {return mpTarget->size1();}

        /// Return the number of columns of the underlying matrix.
        std::ptrdiff_t cols() const {return mpTarget->size2();}

        /// Set the real shift \f$\sigma\f$.
        void set_shift(const Scalar& sigma)
        {
            const TScalar relative_shift = sigma - mLastShift;
            mLastShift = sigma;

            if (relative_shift == static_cast<TScalar>(0)) {
                mpSolver->InitializeSolutionStep(*mpTarget, mDummySolution, mDummyRhs);
            } else {
                InPlaceMatrixAdd(*mpTarget, *mpB, -relative_shift);
                mpSolver->InitializeSolutionStep(*mpTarget, mDummySolution, mDummyRhs);
            }
        }

        ///
        /// Perform the shift-invert operation \f$y=(A-\sigma B)^{-1}x\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = inv(A - sigma * B) * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            //MapConstVec x(x_in, m_n);
            //MapVec y(y_out, m_n);
            //y.noalias() = m_solver.solve(x);
            TUblasSparseSpace<TScalar>::SetToZero(mDummySolution);
            std::copy(x_in,
                      x_in + mpTarget->size1(),
                      mDummyRhs.begin());
            mpSolver->PerformSolutionStep(*mpTarget, mDummySolution, mDummyRhs);
            std::copy(mDummySolution.begin(),
                      mDummySolution.end(),
                      y_out);
        }
    };


}; // class SpectraSymGEigsShiftSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    SpectraSymGEigsShiftSolver<TSparseSpaceType,
    TDenseSpaceType,
    TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const SpectraSymGEigsShiftSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

