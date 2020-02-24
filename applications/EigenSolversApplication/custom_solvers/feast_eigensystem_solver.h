/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author:  Quirin Aumann
*/

#if !defined(KRATOS_FEAST_EIGENSYSTEM_SOLVER_H_INCLUDED)
#define KRATOS_FEAST_EIGENSYSTEM_SOLVER_H_INCLUDED

// External includes
// #include <Eigen/Core>
// #include <Eigen/Eigenvalues>

// Project includes
#include "includes/define.h"
// #if defined EIGEN_USE_MKL_ALL
// #include "eigen_pardiso_ldlt_solver.h"
// #else // defined EIGEN_USE_MKL_ALL
// #include "eigen_sparse_lu_solver.h"
// #endif // defined EIGEN_USE_MKL_ALL
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/ublas_wrapper.h"

extern "C" {
    #include <feast.h>
    #include <feast_sparse.h>
}

namespace Kratos
{

template<
    class TSparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>,
    class TDenseSpaceType = UblasSpace<double, Matrix, Vector>,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class FEASTEigensystemSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(FEASTEigensystemSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    FEASTEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "feast_eigensystem",
            "number_of_eigenvalues": 1,
            "normalize_eigenvectors": false,
            "emin" : 0,
            "emax" : 0,
            "M0" : 0,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~FEASTEigensystemSolver() override {}

    /**
     * Solve the generalized eigenvalue problem using an eigen subspace iteration method
     * The implementation follows the code from
     * K. J. Bathe, Finite Element Procedures second Edition, ISBN-13: 978-0979004957
     * page 954 and following
     * The naming of the variables is chose according to the reference.
     *
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        VectorType& rEigenvalues,
        DenseMatrixType& rEigenvectors) override
    {
        // using scalar_t = double;
        // using vector_t = Eigen::VectorXd;
        // using matrix_t = Eigen::MatrixXd;

        std::cout << "FANCY FEAST EIGENSOLVER!!\n";
        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();

        if( rEigenvalues.size() != mParam["M0"].GetInt() )
            rEigenvalues.resize(mParam["M0"].GetInt(), false);

        if( rEigenvectors.size1() != rK.size1() || rEigenvectors.size2() != mParam["M0"].GetInt() )
            rEigenvectors.resize(rK.size1(), mParam["M0"].GetInt(), false);

        VectorType Residual(mParam["M0"].GetInt());
        KRATOS_WATCH(Residual)

        int fpm[64] = {};
        feastinit(fpm);
        KRATOS_WATCH(fpm)
        echo_level > 0 ? fpm[0] = 1 : fpm[0] = 0;
        // fpm[2] = 8;
        fpm[39] = -1; //M0/2 lowest eigenvalues in interval

        char UPLO = 'F';
        int N = static_cast<int>(rK.size1());
        double* A = rK.value_data().begin();
        int IA[N+1] = {};
        KRATOS_WATCH(N+1)
        KRATOS_WATCH(rK.index1_data().size())
        for( int i=0; i<N+1; ++i )
        {
            IA[i] = static_cast<int>(rK.index1_data()[i]) + 1;
        }
        int JA[IA[N]-1] = {};
        KRATOS_WATCH(IA[N])
        KRATOS_WATCH(rK.index2_data().size())
        for( int i=0; i<IA[N]-1; ++i )
        {
            JA[i] = static_cast<int>(rK.index2_data()[i]) + 1;
        }

        double* B = rM.value_data().begin();
        int IB[N+1] = {};
        KRATOS_WATCH(N+1)
        KRATOS_WATCH(rM.index1_data().size())
        for( int i=0; i<N+1; ++i )
        {
            IB[i] = static_cast<int>(rM.index1_data()[i]) + 1;
        }
        int JB[IB[N]-1] = {};
        KRATOS_WATCH(IB[N])
        KRATOS_WATCH(rM.index2_data().size())
        for( int i=0; i<IB[N]-1; ++i )
        {
            JB[i] = static_cast<int>(rM.index2_data()[i]) + 1;
        }
        double epsout;
        int loop;
        double Emin = mParam["emin"].GetDouble();
        double Emax = mParam["emax"].GetDouble();
        // int* M0 = (int*) mParam["M0"].GetInt();
        int M0 = mParam["M0"].GetInt();
        double* E = rEigenvalues.data().begin();
        double* X = rEigenvectors.data().begin();
        int M;
        double* res = Residual.data().begin();
        int info;
        std::cout << "geschafft\n";
        /**
        // KRATOS_WATCH(rM.index1_data())
        std::for_each(rM.index1_data().begin(), rM.index1_data().end(), [](size_t i) { std::cout << i << ","; });
        std::cout << "\n";
        KRATOS_WATCH(IB[0])
        KRATOS_WATCH(IB[1])
        KRATOS_WATCH(IB[2])
        KRATOS_WATCH(IB[3])
        std::for_each(rM.index2_data().begin(), rM.index2_data().end(), [](size_t i) { std::cout << i << ","; });
        std::cout << "\n";
        KRATOS_WATCH(JB[0])
        KRATOS_WATCH(JB[1])
        KRATOS_WATCH(JB[2])
        std::for_each(rM.value_data().begin(), rM.value_data().end(), [](double i) { std::cout << i << ","; });
        std::cout << "\n";
        std::for_each(rK.index1_data().begin(), rK.index1_data().end(), [](size_t i) { std::cout << i << ","; });
        std::cout << "\n";
        KRATOS_WATCH(IA[0])
        KRATOS_WATCH(IA[1])
        KRATOS_WATCH(IA[2])
        KRATOS_WATCH(IA[3])
        std::for_each(rK.index2_data().begin(), rK.index2_data().end(), [](size_t i) { std::cout << i << ","; });
        std::cout << "\n";
        KRATOS_WATCH(JA[0])
        KRATOS_WATCH(JA[1])
        KRATOS_WATCH(JA[2])
        KRATOS_WATCH(JA[3])
        std::for_each(rK.value_data().begin(), rK.value_data().end(), [](double i) { std::cout << i << ","; });
        std::cout << "\n";
        **/
        dfeast_scsrgv(&UPLO, &N, A, IA, JA, B, IB, JB, fpm, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);
        
        std::cout << "yeah\n";

        // // --- wrap ublas matrices

        // UblasWrapper<scalar_t> a_wrapper(rK);
        // UblasWrapper<scalar_t> b_wrapper(rM);

        // const auto& a = a_wrapper.matrix();
        // const auto& b = b_wrapper.matrix();


        // // --- timer

        // double start_time = OpenMPUtils::GetCurrentTime();

        // KRATOS_INFO_IF("FEASTEigensystemSolver:", echo_level > 0) << "Start"  << std::endl;

        // // --- calculation

        // int nn = a.rows();
        // int nc = std::min(2 * nroot, nroot + 8);

        // // projections
        // matrix_t ar(nc, nc);
        // matrix_t br(nc, nc);

        // // storage for eigenvalues
        // vector_t prev_eigv = vector_t::Zero(nc);

        // // storage for eigenvectors
        // matrix_t r = matrix_t::Zero(nn, nc);
        // for (int i = 0; i != nn; ++i) {
        //     r(i, 0) = b.coeff(i, i);
        // }

        // vector_t tmp(nn);

        // // working vector
        // vector_t w(nn);
        // for (int i = 0; i != nn; ++i) {
        //     w(i) = r(i, 0) / a.coeff(i, i);
        // }

        // int nd = nn / nc;
        // int l = nn - nd;

        // vector_t tt(nn);
        // int ij = 0;

        // tt(0) = 0.0;

        // for (int j = 1; j != nc; ++j) {
        //     double rt = 0.0;

        //     for (int i = 0; i != l; ++i) {
        //         if (w(i) >= rt) {
        //             rt = w(i);
        //             ij = i;
        //         }
        //     }

        //     for (int i = l - 1; i != nn; ++i) {
        //         if (w(i) > rt) {
        //             rt = w(i);
        //             ij = i;
        //         }
        //     }

        //     tt(j) = ij;
        //     w(ij) = 0.0;

        //     l -= nd;

        //     r(ij, j) = 1.0;
        // }

        // #if defined USE_EIGEN_MKL
        // EigenPardisoLDLTSolver<double> solver;
        // #else  // defined USE_EIGEN_MKL
        // EigenSparseLUSolver<double> solver;
        // #endif // defined USE_EIGEN_MKL

        // solver.Compute(a);

        // int iteration = 0;

        // Eigen::GeneralizedSelfAdjointEigenSolver<matrix_t> eig;

        // do {
        //     iteration++;

        //     KRATOS_INFO_IF("FEASTEigensystemSolver:", echo_level > 1) << "Iteration " << iteration <<std::endl;

        //     for (int j = 0; j != nc; ++j) {
        //         tmp = r.col(j);
        //         solver.Solve(tmp, tt);

        //         for (int i = j; i != nc; ++i) {
        //             ar(i, j) = r.col(i).dot(tt);
        //         }

        //         r.col(j) = tt;
        //     }

        //     for (int j = 0; j != nc; ++j) {
        //         tt = b * r.col(j);

        //         for (int i = j; i != nc; ++i) {
        //             br(i, j) = r.col(i).dot(tt);
        //         }

        //         r.col(j) = tt;
        //     }

        //     eig.compute(ar, br);

        //     if(eig.info() != Eigen::Success) {
        //         KRATOS_WARNING("FEASTEigensystemSolver:") << "Eigen solution was not successful!" << std::endl;
        //         break;
        //     }

        //     r *= eig.eigenvectors();

        //     bool is_converged = true;
        //     for (int i = 0; i != nc; i++) {
        //         double eigv = eig.eigenvalues()(i);
        //         double dif = eigv - prev_eigv(i);
        //         double rtolv = std::abs(dif / eigv);

        //         if (rtolv > tolerance) {
        //             is_converged = false;
        //             KRATOS_WARNING_IF("FEASTEigensystemSolver:", echo_level > 1) << "Convergence not reached for eigenvalue #"<<i+1<<": " << rtolv <<"." << std::endl;
        //             break;
        //         }
        //     }

        //     if (is_converged) {
        //         KRATOS_INFO_IF("FEASTEigensystemSolver:", echo_level > 0) << "Convergence reached after " << iteration << " iterations within a relative tolerance: " << tolerance << std::endl;
        //         break;
        //     } else if (iteration >= max_iteration) {
        //         KRATOS_INFO_IF("FEASTEigensystemSolver:", echo_level > 0) << "Convergence not reached in " << max_iteration << " iterations." << std::endl;
        //         break;
        //     }

        //     prev_eigv = eig.eigenvalues();
        // } while (true);


        // if (static_cast<int>(rEigenvalues.size()) != nroot) {
        //     rEigenvalues.resize(nroot);
        // }
        // if (static_cast<int>(rEigenvectors.size1()) != nroot || static_cast<int>(rEigenvectors.size2()) != nn) {
        //     rEigenvectors.resize(nroot, nn);
        // }

        // Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        // Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());

        // eigvals = eig.eigenvalues().head(nroot);

        // for (int i = 0; i != nroot; ++i) {
        //     tmp = r.col(i);
        //     solver.Solve(tmp, eigvecs.row(i));
        //     eigvecs.row(i).normalize();
        // }

        // // --- normalization
        // // Given generalized eigenvalue problem (A - eigenvalue * B) * eigenvector = 0,
        // // eigenvector is normalized such that eigenvector^T * B * eigenvector = 1
        // if(mParam["normalize_eigenvectors"].GetBool())
        // {
        //     for (int i = 0; i != nroot; ++i)
        //     {
        //         const double tmp = eigvecs.row(i) * b * eigvecs.row(i).transpose();
        //         const double factor = 1.0 / std::sqrt(tmp);
        //         eigvecs.row(i) *=  factor;
        //         KRATOS_INFO_IF("FEASTEigensystemSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
        //     }
        // }

        // // --- output
        // if (echo_level > 0) {
        //     double end_time = OpenMPUtils::GetCurrentTime();
        //     double duration = end_time - start_time;

        //     Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

        //     KRATOS_INFO("FEASTEigensystemSolver:") << "Completed in " << duration << " seconds" << std::endl
        //               << "                   Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        // }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "FEASTEigensystemSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

}; // class FEASTEigensystemSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    FEASTEigensystemSolver<TSparseSpaceType,
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
    const FEASTEigensystemSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_FEAST_EIGENSYSTEM_SOLVER_H_INCLUDED)
