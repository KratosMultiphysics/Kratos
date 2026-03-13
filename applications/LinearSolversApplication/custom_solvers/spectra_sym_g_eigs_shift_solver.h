/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Authors: Armin Geiser
*/

#if !defined(KRATOS_SPECTRA_SYM_G_EIGS_SHIFT_SOLVER_H_INCLUDED)
#define KRATOS_SPECTRA_SYM_G_EIGS_SHIFT_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Spectra/SymGEigsShiftSolver.h>
#if defined EIGEN_USE_MKL_ALL
#include "eigen_pardiso_lu_solver.h"
#endif // defined EIGEN_USE_MKL_ALL

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_utilities/ublas_wrapper.h"
#include "utilities/builtin_timer.h"

namespace Kratos
{

template<
    class TSparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>,
    class TDenseSpaceType = UblasSpace<double, Matrix, Vector>,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SpectraSymGEigsShiftSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

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
            "echo_level": 1
        })");

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
        
       
        // --- timer
        const auto timer = BuiltinTimer();
        // If the eigenvalue solution is not sough after return trivial solution of zeros

        Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());
if (nroot !=0){
        // --- wrap ublas matrices
        UblasWrapper<scalar_t> a_wrapper(rK);
        UblasWrapper<scalar_t> b_wrapper(rM);

        const auto& a = a_wrapper.matrix();
        const auto& b = b_wrapper.matrix();


        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        using OpType = OwnSymShiftInvert<scalar_t>;
        using BOpType = OwnSparseSymMatProd<scalar_t>;

        OpType op(a, b);
        BOpType Bop(b);
        const int ncv = 3 * nroot;  // TODO find a good value
        Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert> eigs(op, Bop, nroot, ncv, shift);

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
    }
    else{
        if (static_cast<int>(rEigenvalues.size()) != nroot) {
            rEigenvalues.resize(nroot);
        }
        if (static_cast<int>(rEigenvectors.size1()) != nroot || rEigenvectors.size2() != rK.size1()) {
            rEigenvectors.resize(nroot, rK.size1());
        }
        eigvals.setZero();
        eigvecs.setZero();
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

        ///
        /// Perform the matrix-vector multiplication operation \f$y=Ax\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = A * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in, m_mat.cols());
            MapVec y(y_out, m_mat.rows());
            y.noalias() = m_mat * x;
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
    /// \tparam Scalar_        The element type of the matrices.
    ///                        Currently supported types are `float`, `double`, and `long double`.
    template <typename Scalar_>
    class OwnSymShiftInvert
    {
    public:
        ///
        /// Element type of the matrix.
        ///
        using Scalar = Scalar_;

    private:
        using Index = Eigen::Index;

        // type of the A matrix
        using MatrixA = Kratos::EigenSparseMatrix<Scalar>;

        // type of the B matrix
        using MatrixB = Kratos::EigenSparseMatrix<Scalar>;

        using Vector = Kratos::EigenDynamicVector<Scalar>;
        using MapConstVec = Eigen::Map<const Vector>;
        using MapVec = Eigen::Map<Vector>;

        using ResType = MatrixA;

        #if defined EIGEN_USE_MKL_ALL
        using FacType = Eigen::PardisoLU<ResType>;
        #else
        using FacType = Eigen::SparseLU<ResType>;
        #endif // defined EIGEN_USE_MKL_ALL

        using ConstGenericMatrixA = const Eigen::Ref<const MatrixA>;
        using ConstGenericMatrixB = const Eigen::Ref<const MatrixB>;

        ConstGenericMatrixA m_matA;
        ConstGenericMatrixB m_matB;
        const Index m_n;
        FacType m_solver;

    public:
        ///
        /// Constructor to create the matrix operation object.
        ///
        /// \param A A sparse matrix object, whose type can be
        ///          `Eigen::SparseMatrix<...>`,
        ///          `Eigen::Map<Eigen::SparseMatrix<...>>`,
        ///          `Eigen::Ref<Eigen::SparseMatrix<...>>`, etc.
        /// \param B A sparse matrix object.
        ///
        OwnSymShiftInvert(ConstGenericMatrixA& A, ConstGenericMatrixB& B) :
            m_matA(A), m_matB(B), m_n(A.rows())
        {
            KRATOS_ERROR_IF(m_n != A.cols() || m_n != B.rows() || m_n != B.cols()) << "SymShiftInvert: A and B must be square matrices of the same size" << std::endl;
        }

        ///
        /// Return the number of rows of the underlying matrix.
        ///
        Index rows() const { return m_n; }
        ///
        /// Return the number of columns of the underlying matrix.
        ///
        Index cols() const { return m_n; }

        ///
        /// Set the real shift \f$\sigma\f$.
        ///
        void set_shift(const Scalar& sigma)
        {
            if (sigma == 0.0) {
                m_solver.compute(m_matA);
            } else {
                m_solver.compute(m_matA - sigma * m_matB);
            }
            const bool success = m_solver.info() == Eigen::Success;
            KRATOS_ERROR_IF_NOT(success) << "SymShiftInvert: factorization failed with the given shift" << std::endl;
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
            MapConstVec x(x_in, m_n);
            MapVec y(y_out, m_n);
            y.noalias() = m_solver.solve(x);
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

#endif // defined(KRATOS_SPECTRA_SYM_G_EIGS_SHIFT_SOLVER_H_INCLUDED)
