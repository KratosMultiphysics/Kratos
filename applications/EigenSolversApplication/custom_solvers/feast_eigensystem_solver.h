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
        // noalias(rEigenvectors) = ZeroMatrix(rEigenvectors.size1(), rEigenvectors.size2());

        //create column based matrix for the fortran routine
        //TODO: change data type
        matrix<double,column_major> tmp_eigenvectors(rEigenvectors.size1(), rEigenvectors.size2());

        VectorType Residual(mParam["M0"].GetInt());
        KRATOS_WATCH(Residual)

        int fpm[64] = {};
        feastinit(fpm);
        KRATOS_WATCH(fpm)
        echo_level > 0 ? fpm[0] = 1 : fpm[0] = 0;
        // fpm[2] = 8;

        //TODO if this is used, the other half of eigenvalues/vectors should be excluded from the solution
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
        double* X = tmp_eigenvectors.data().begin();
        // double* X = rEigenvectors.data().begin();
        int M;
        double* res = Residual.data().begin();
        int info;

        // std::cout << "geschafft\n";
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
        
        // copy eigenvectors back to the provided row based matrix
        noalias(rEigenvectors) = tmp_eigenvectors;
        // std::cout << "yeah\n";

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
