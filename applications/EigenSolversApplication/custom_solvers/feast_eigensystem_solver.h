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

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
// #include "linear_solvers/iterative_solver.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
// #include "custom_utilities/ublas_wrapper.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

extern "C" {
    #include <feast.h>
    #include <feast_sparse.h>
}

namespace Kratos
{

template<
    typename TScalar = double,
    class TSparseSpaceType = TUblasSparseSpace<TScalar>,
    class TDenseSpaceType = TUblasDenseSpace<TScalar>>
// template <typename TScalar = double>
class FEASTEigensystemSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(FEASTEigensystemSolver);

    // typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef TScalar ValueType;

    FEASTEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "feast_eigensystem",
            "number_of_eigenvalues": 0,
            "search_lowest_eigenvalues": false,
            "search_highest_eigenvalues": false,
            "e_min" : 0.0,
            "e_max" : 1.0,
            "subspace_size" : 0,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "echo_level": 1
        })");

        //TODO: e_mid, r for complex!!!!!!!!!!

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        // BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());

        KRATOS_ERROR_IF( mParam["search_lowest_eigenvalues"].GetBool() && mParam["search_highest_eigenvalues"].GetBool() ) <<
            "Cannot search for highest and lowest eigenvalues at the same time\n";

        KRATOS_ERROR_IF( mParam["e_min"].GetDouble() > mParam["e_max"].GetDouble() ) <<
            "Invalid eigenvalue limits provided\n";

        KRATOS_INFO_IF( "FEASTEigensystemSolver", 
            (mParam["search_lowest_eigenvalues"].GetBool() || mParam["search_highest_eigenvalues"].GetBool()) 
            && (mParam["subspace_size"].GetInt() > 0) ) <<
            "Manually defined subspace size will be overwritten because extremal eigenvalues are sought\n";
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
        // settings
        const size_t system_size = rK.size1();
        size_t subspace_size;

        if( mParam["search_lowest_eigenvalues"].GetBool() || mParam["search_highest_eigenvalues"].GetBool() )
        {
            subspace_size = 2 * static_cast<size_t>(mParam["number_of_eigenvalues"].GetInt());
        }
        else
        {
            subspace_size = static_cast<size_t>(mParam["subspace_size"].GetInt());
        }

        const int echo_level = mParam["echo_level"].GetInt();

        if( rEigenvalues.size() != subspace_size )
            rEigenvalues.resize(subspace_size, false);

        if( rEigenvectors.size1() != system_size || rEigenvectors.size2() != subspace_size )
            rEigenvectors.resize(system_size, subspace_size, false);

        // create column based matrix for the fortran routine
        matrix<ValueType, column_major> tmp_eigenvectors(rEigenvectors.size1(), rEigenvectors.size2());
        // matrix<double, column_major> tmp_eigenvectors = FeastEigenvectorMatrix(rEigenvectors);
        VectorType residual(subspace_size);

        // set FEAST settings
        int fpm[64] = {};
        feastinit(fpm);

        echo_level > 0 ? fpm[0] = 1 : fpm[0] = 0;

        if( mParam["search_lowest_eigenvalues"].GetBool() )
        {
            fpm[39] = -1;
        }
        if( mParam["search_highest_eigenvalues"].GetBool() )
        {
            fpm[39] = 1;
        }

        char UPLO = 'F';
        int N = static_cast<int>(system_size);

        double* A = reinterpret_cast<double*>(rK.value_data().begin());
        int IA[N+1] = {};
        // KRATOS_WATCH(N+1)
        // KRATOS_WATCH(rK.index1_data().size())
        for( int i=0; i<N+1; ++i )
        {
            IA[i] = static_cast<int>(rK.index1_data()[i]) + 1;
        }
        int JA[IA[N]-1] = {};
        // KRATOS_WATCH(IA[N])
        // KRATOS_WATCH(rK.index2_data().size())
        for( int i=0; i<IA[N]-1; ++i )
        {
            JA[i] = static_cast<int>(rK.index2_data()[i]) + 1;
        }

        double* B = reinterpret_cast<double*>(rM.value_data().begin());
        int IB[N+1] = {};
        // KRATOS_WATCH(N+1)
        // KRATOS_WATCH(rM.index1_data().size())
        for( int i=0; i<N+1; ++i )
        {
            IB[i] = static_cast<int>(rM.index1_data()[i]) + 1;
        }
        int JB[IB[N]-1] = {};
        // KRATOS_WATCH(IB[N])
        // KRATOS_WATCH(rM.index2_data().size())
        for( int i=0; i<IB[N]-1; ++i )
        {
            JB[i] = static_cast<int>(rM.index2_data()[i]) + 1;
        }

        double epsout;
        int loop;
        double Emin = mParam["e_min"].GetDouble();
        double Emax = mParam["e_max"].GetDouble();
        int M0 = static_cast<int>(subspace_size);
        double* E = reinterpret_cast<double*>(rEigenvalues.data().begin());
        double* X = reinterpret_cast<double*>(tmp_eigenvectors.data().begin());
        int M;
        double* res = reinterpret_cast<double*>(residual.data().begin());
        int info;

        ValueType T;
        // dfeast_scsrgv(&UPLO, &N, A, IA, JA, B, IB, JB, fpm, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);
        fptr feast = CallFeast(T);

        feast(&UPLO, &N, A, IA, JA, B, IB, JB, fpm, &epsout, &loop, &Emin, &Emax, &M0, E, X, &M, res, &info);
        
        // copy eigenvectors back to the provided row based matrix
        noalias(rEigenvectors) = tmp_eigenvectors;

        // discard the spurious eigenvalues
        if( mParam["search_lowest_eigenvalues"].GetBool() || mParam["search_highest_eigenvalues"].GetBool() )
        {
            // rEigenvectors.resize(system_size, subspace_size/2, true);
            rEigenvalues.resize(subspace_size/2, true);
        }

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

  private:

    typedef void (*fptr)(char*, int*, double*, int*, int*, double*, int*, int*, int*, double*, int*, double*, double*, int*, double*, double*, int*, double*, int*);

    fptr CallFeast(double T)
    {
        return dfeast_scsrgv;
    }

    fptr CallFeast(std::complex<double> T)
    {
        return zfeast_scsrgv;
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
