/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Authors: Armin Geiser
*/

#if !defined(KRATOS_SPECTRA_G_EIGS_SHIFT_SOLVER_H_INCLUDED)
#define KRATOS_SPECTRA_G_EIGS_SHIFT_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Spectra/GenEigsRealShiftSolver.h>
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
#include "utilities/math_utils.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "eigen_sparse_lu_solver.h"

namespace Kratos
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType  = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;
using LU_SolverType   = SkylineLUFactorizationSolver<SparseSpaceType, DenseSpaceType>;
using AMGCLSolverType = AMGCLSolver< SparseSpaceType,DenseSpaceType>;
using BICGSTABSolverType = BICGSTABSolver< SparseSpaceType,DenseSpaceType>;


template<
    class TSparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>,
    class TDenseSpaceType = UblasSpace<double, Matrix, Vector>,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SpectraGEigsShiftSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(SpectraGEigsShiftSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    SpectraGEigsShiftSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "spectra_g_eigs_shift",
            "number_of_eigenvalues": 1,
            "normalize_eigenvectors": false,
            "max_iteration": 5000,
            "tol": 1e-6,
            "shift": 0.0,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);
        
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~SpectraGEigsShiftSolver() override {}

    /**
     * Solve the generalized eigenvalue problem of matrices K and M
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
        double shift = mParam["shift"].GetDouble();
        shift = 1.0 / (shift * shift);
        scalar_t tol = mParam["tol"].GetDouble();

        // --- calculate intermediate general matrix
        // --- timer
        const auto timer = BuiltinTimer();
        // If the eigenvalue solution is not sough after return trivial solution of zeros

        Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());
if (nroot !=0){
        SparseMatrixType rA;
        auto amgcl_settings = Parameters(R"({
        "solver_type" : "amgcl",
        "max_iteration" : 1000,
        "tolerance" : 1e-10
        })");

        auto p_mass_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);
        BuildGeneralizedSystemMatrix(rK, rM, rA, amgcl_settings);

        // --- wrap ublas matrices
       

        UblasWrapper<scalar_t> a_wrapper(rA);
        UblasWrapper<scalar_t> b_wrapper(rM);
        

        const auto& A = a_wrapper.matrix();
        const auto& B = b_wrapper.matrix();

        KRATOS_INFO_IF("SpectraGEigsShiftSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        using OpType = OwnShiftInvert<scalar_t>;

        OpType op(A);
        int ncv = 2 * nroot +2;  // TODO find a good value
        if (ncv < 30){ncv = 30;}
        
        Spectra::GenEigsRealShiftSolver<OpType> eigs(op, nroot, ncv, shift);

        eigs.init();
        // const int nconv = eigs.compute(Spectra::SortRule::SmallestMagn, 1e-16, max_iteration);
        const int nconv = eigs.compute(Spectra::SortRule::SmallestMagn , max_iteration, tol);
        const int niter = eigs.num_iterations();
        const int nops = eigs.num_operations();

        KRATOS_INFO_IF("SpectraGEigsShiftSolver:", echo_level > 0) << "nconv = " << nconv << std::endl;
        KRATOS_INFO_IF("SpectraGEigsShiftSolver:", echo_level > 0) << "niter = " << niter << std::endl;
        KRATOS_INFO_IF("SpectraGEigsShiftSolver:", echo_level > 0) << "nops  = " << nops << std::endl;
        KRATOS_WATCH(nconv)
        KRATOS_WATCH(niter)
        KRATOS_WATCH(nops)

        KRATOS_ERROR_IF(eigs.info() != Spectra::CompInfo::Successful) << "SpectraGEigsShiftSolver: failed" << std::endl;

        if (static_cast<int>(rEigenvalues.size()) != nroot) {
            rEigenvalues.resize(nroot);
        }
        if (static_cast<int>(rEigenvectors.size1()) != nroot || rEigenvectors.size2() != rK.size1()) {
            rEigenvectors.resize(nroot, rK.size1());
        }

        // Spectra::SortRule::LargestAlge results in values being in descending order
        eigvals = eigs.eigenvalues().real();

        const std::size_t n = eigvals.size();
        for (std::size_t i = 0; i < n; ++i){
            if(eigvals[i]!=0.0){eigvals[i]=1/eigvals[i];}
        }
        // conversion back to RowMajor ordering of Kratos and transpose because Spectra returns eigvecs as columns
        // eigvecs = eigs.eigenvectors().real().transpose().colwise().reverse();
        eigvecs = eigs.eigenvectors().real().transpose();

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
                KRATOS_INFO_IF("SpectraGEigsShiftSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
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

            KRATOS_INFO("SpectraGEigsShiftSolver:") << "Completed in " << timer << std::endl
                      << "                   Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "SpectraGEigsShiftSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

private:

    void BuildGeneralizedSystemMatrix(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        SparseMatrixType& rA, Parameters param = Parameters(R"({})"))
    {
        const std::size_t n = rK.size1();
        KRATOS_ERROR_IF(n != rK.size2() || n != rM.size1() || n != rM.size2())
            << "ShiftInvert: A and B must be square matrices of the same size" << std::endl;


        // LU_SolverType lu_solver;
        // AMGCLSolverType amgcl_solver(param);
        // BICGSTABSolverType bicg_solver;
        Eigen::SparseLU<Kratos::EigenSparseMatrix<double>> splu_solver;

        DenseMatrixType Kinv = ZeroMatrix(n,n);
        DenseMatrixType I = ZeroMatrix(n,n);
        for (std::size_t i = 0; i < n; ++i){I(i,i)=1.0;}
        
        // Initialize A as zero sparse matrix
        rA.clear();
        rA.resize(n, n, false);


        UblasWrapper<double> a_wrapper(rK);
        // UblasWrapper<double> x_wrapper(Minv);
        // UblasWrapper<double> b_wrapper(I);
        

        const auto& A = a_wrapper.matrix();
        // const auto& X = x_wrapper.matrix();
        // const auto& B = b_wrapper.matrix();

        OwnShiftInvert<double> opxd(A);
        opxd.set_shift(0);
        


        Kratos::Vector x(n,0);
        Kratos::Vector b(n,0);
        // UblasWrapper<double> b_wrapper(b);
        // UblasWrapper<double> x_wrapper(x);
        // const auto& B = b_wrapper.matrix();
        // const auto& X = x_wrapper.matrix();

        // opxd.perform_op(B.data(),X.data())
        

        // BuiltinTimer system_construction_time;
        
        // KRATOS_INFO_IF("Fill Time", true)
        //     << system_construction_time << std::endl;
        // BuiltinTimer system_construction;
        // lu_solver.Solve(rM, Minv, I);
        
        // KRATOS_INFO_IF("Solve Time", true)
        //     << system_construction << std::endl;

        for (std::size_t j = 0; j < n; ++j)
        {
            // Create a column vector of zeros and the corresponding row being 1.0
            std::fill(b.begin(), b.end(), 0.0);
            std::fill(x.begin(), x.end(), 0.0);
            b[j] = 1.0;

            // lu_solver.Solve(rM, x, b);
            opxd.perform_op(&b[0],&x[0]);
            // amgcl_solver.Solve(rM, x, b);
            // bicg_solver.Solve(rM, x, b);
            // Write x into column j of Minv
            double threshold = std::abs(x[j]) * 0;
            for (std::size_t i = 0; i < n; i++) {
                if(std::abs(x[i])>threshold){Kinv(i,j) = x[i];}
			} 
            // TDenseSpaceType::SetColumn(j,Minv, x);
        }        
        rA =  prod(Kinv,rM);
        // KRATOS_WATCH("Stiffness")
        // SparseSpaceType::WriteMatrixMarketMatrix("K.mm", rK, false);
        // KRATOS_WATCH("MinvK")
        // SparseSpaceType::WriteMatrixMarketMatrix("MinvK.mm", rA, false);
        // KRATOS_WATCH("Exit")
        rA.complete_index1_data();
        // SparseSpaceType::WriteMatrixMarketMatrix("MinvK.mm", rA, false);
    }

        ///
    /// \ingroup MatOp
    ///
    /// This class defines matrix operations required by the generalized eigen solver
    /// in the shift-and-invert mode. Given a matrix \f$A\f$,
    /// it solves the linear equation \f$y=(A-\sigma I)^{-1}x\f$, where \f$\sigma\f$ is a real shift.
    /// This is adapted from Spectra::SparseGenRealShiftSolve.h
    ///
    /// This class is intended to be used with the GEigsRealShiftSolver generalized eigen solver.
    ///
    /// \tparam Scalar_        The element type of the matrices.
    ///                        Currently supported types are `float`, `double`, and `long double`.
    template <typename Scalar_>
    class OwnShiftInvert
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
        using SparseMatrix = Kratos::EigenSparseMatrix<Scalar>;
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

        ConstGenericMatrixA m_matA;
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
        ///
        OwnShiftInvert(ConstGenericMatrixA& A) :
            m_matA(A), m_n(A.rows())
        {
            KRATOS_ERROR_IF(m_n != A.cols()) << "ShiftInvert: A must be square matrix" << std::endl;
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
                SparseMatrix I(m_n, m_n);
                I.setIdentity();
                m_solver.compute(m_matA - sigma * I);
            }
            const bool success = m_solver.info() == Eigen::Success;
            
            KRATOS_ERROR_IF_NOT(success) << "ShiftInvert: factorization failed with the given shift" << std::endl;
        }

        ///
        /// Perform the shift-invert operation \f$y=(A-\sigma I)^{-1}x\f$.
        ///
        /// \param x_in  Pointer to the \f$x\f$ vector.
        /// \param y_out Pointer to the \f$y\f$ vector.
        ///
        // y_out = inv(A - sigma * I) * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in, m_n);
            MapVec y(y_out, m_n);
            y.noalias() = m_solver.solve(x);
        }

    };


}; // class SpectraGEigsShiftSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    SpectraGEigsShiftSolver<TSparseSpaceType,
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
    const SpectraGEigsShiftSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_SPECTRA_G_EIGS_SHIFT_SOLVER_H_INCLUDED)
