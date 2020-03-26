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
#include "linear_solvers/linear_solver.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"

extern "C" {
    #include <feast.h>
    #include <feast_sparse.h>
}

namespace Kratos
{

template<
    bool TSymmetric,
    typename TScalarIn,
    typename TScalarOut,
    class TSparseSpaceTypeIn = TUblasSparseSpace<TScalarIn>,
    class TDenseSpaceTypeIn = TUblasDenseSpace<TScalarIn>,
    class TSparseSpaceTypeOut = TUblasSparseSpace<TScalarOut>,
    class TDenseSpaceTypeOut = TUblasDenseSpace<TScalarOut>>
class FEASTEigensystemSolver
    : public LinearSolver<TSparseSpaceTypeIn, TDenseSpaceTypeIn>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(FEASTEigensystemSolver);

    typedef LinearSolver<TSparseSpaceTypeIn, TDenseSpaceTypeIn> BaseType;

    typedef typename TSparseSpaceTypeIn::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceTypeIn::VectorType VectorType;

    typedef typename TDenseSpaceTypeIn::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceTypeOut::MatrixType DenseOutputMatrixType;

    typedef matrix<TScalarOut, column_major> FEASTMatrixType;

    typedef typename TSparseSpaceTypeOut::VectorType OutputVectorType;

    typedef TScalarIn ValueType;

    FEASTEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type" : "eigen_feast",
            "symmetric" : true,
            "keep_real_solution" : false,
            "number_of_eigenvalues" : 0,
            "search_lowest_eigenvalues" : false,
            "search_highest_eigenvalues" : false,
            "e_min" : 0.0,
            "e_max" : 0.0,
            "e_mid_re" : 0.0,
            "e_mid_im" : 0.0,
            "e_r" : 0.0,
            "subspace_size" : 0,
            "max_iteration" : 20,
            "tolerance" : 1e-12,
            "echo_level" : 0
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        KRATOS_ERROR_IF( mParam["number_of_eigenvalues"].GetInt() < 0 ) <<
            "Invalid number of eigenvalues provided" << std::endl;

        KRATOS_ERROR_IF( mParam["subspace_size"].GetInt() < 0 ) <<
            "Invalid subspace size provided" << std::endl;

        KRATOS_ERROR_IF( mParam["max_iteration"].GetInt() < 1 ) <<
            "Invalid maximal number of iterations provided" << std::endl;

        KRATOS_INFO_IF( "FEASTEigensystemSolver",
            mParam["number_of_eigenvalues"].GetInt() > 0  && mParam["subspace_size"].GetInt() > 0 ) <<
            "Manually defined number of eigenvalues will be overwritten to match the defined subspace size" << std::endl;

        const TScalarOut T = {};
        CheckParameters(T);
    }

    ~FEASTEigensystemSolver() override {}

    /**
     * Solve the generalized eigenvalue problem using FEAST
     * @param rK first input matrix
     * @param rM second input matrix
     * @param rEigenvalues eigenvalues
     * @param rEigenvectors row-aligned eigenvectors [n_evs,n_dofs]
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
        const TScalarIn Ti = {};
        const TScalarOut T = {};

        if( mParam["search_lowest_eigenvalues"].GetBool() || mParam["search_highest_eigenvalues"].GetBool() )
        {
            subspace_size = 2 * static_cast<size_t>(mParam["number_of_eigenvalues"].GetInt());
        }
        else if( mParam["subspace_size"].GetInt() == 0 )
        {
            subspace_size = 1.5 * static_cast<size_t>(mParam["number_of_eigenvalues"].GetInt());
        }
        else
        {
            subspace_size = static_cast<size_t>(mParam["subspace_size"].GetInt());
        }

        // create column based matrix for the fortran routine
        FEASTMatrixType tmp_eigenvectors(system_size, subspace_size);
        OutputVectorType tmp_eigenvalues(subspace_size);
        VectorType residual(subspace_size);

        // set FEAST settings
        int fpm[64] = {};
        feastinit(fpm);

        mParam["echo_level"].GetInt() > 0 ? fpm[0] = 1 : fpm[0] = 0;

        fpm[2] = -std::log10(mParam["tolerance"].GetDouble());
        fpm[3] = mParam["max_iteration"].GetInt();

        // compute only right eigenvectors
        if( !TSymmetric )
            fpm[14] = 1;

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

        // provide matrices in array form. fortran indices start with 1, must be int
        double* A = reinterpret_cast<double*>(rK.value_data().begin());
        std::vector<int> IA;
        IA.reserve(N+1);
        for( int i=0; i<N+1; ++i )
        {
            IA[i] = static_cast<int>(rK.index1_data()[i]) + 1;
        }
        std::vector<int> JA;
        JA.reserve(IA[N]-1);
        for( int i=0; i<IA[N]-1; ++i )
        {
            JA[i] = static_cast<int>(rK.index2_data()[i]) + 1;
        }

        double* B = reinterpret_cast<double*>(rM.value_data().begin());
        std::vector<int> IB;
        IB.reserve(N+1);
        for( int i=0; i<N+1; ++i )
        {
            IB[i] = static_cast<int>(rM.index1_data()[i]) + 1;
        }
        std::vector<int> JB;
        JB.reserve(IB[N]-1);
        for( int i=0; i<IB[N]-1; ++i )
        {
            JB[i] = static_cast<int>(rM.index2_data()[i]) + 1;
        }

        double epsout;
        int loop;
        TScalarOut E1 = GetE1(T);
        double E2 = GetE2(T);
        double* Emin = reinterpret_cast<double*>(&E1);
        double* Emax = reinterpret_cast<double*>(&E2);
        int M0 = static_cast<int>(subspace_size);
        double* E = reinterpret_cast<double*>(tmp_eigenvalues.data().begin());
        double* X = reinterpret_cast<double*>(tmp_eigenvectors.data().begin());
        int M;
        double* res = reinterpret_cast<double*>(residual.data().begin());
        int info;

        // call feast
        auto feast = CallFeast(Ti, TSymmetric);
        feast(&UPLO, &N, A, IA.data(), JA.data(), B, IB.data(), JB.data(), fpm, &epsout, &loop, Emin, Emax, &M0, E, X, &M, res, &info);

        KRATOS_ERROR_IF(info < 0 || info > 99) << "FEAST encounterd error " << info << ". Please check FEAST output." << std::endl;
        KRATOS_INFO_IF("FeastEigensystemSolver", info > 1 && info < 6) << "FEAST finished with warning " << info << ". Please check FEAST output." << std::endl;
        KRATOS_ERROR_IF(info == 1) << "FEAST finished with warning " << info << ", no eigenvalues could be found in the given search interval." << std::endl;
        KRATOS_ERROR_IF(info == 7) << "FEAST finished with warning " << info << ", no extremal eigenvalues could be found. Please check FEAST output." << std::endl;

        // truncate non-converged results
        tmp_eigenvalues.resize(M, true);
        tmp_eigenvectors.resize(system_size, M, true);

        RetrieveEigensolution<OutputVectorType, VectorType>(tmp_eigenvalues, rEigenvalues, tmp_eigenvectors, rEigenvectors);
    }

    DenseOutputMatrixType& GetEigenvectorSolution()
    {
        KRATOS_ERROR_IF( mpEigenvectorMatrix == nullptr ) << "The eigenvector matrix is not initialized. This means, that the real general generalized " <<
            "eigenvalue problem did not yield complex results or the solution has not yet been performed" << std::endl;
        return *mpEigenvectorMatrix;
    }

    OutputVectorType& GetEigenvalueSolution()
    {
        KRATOS_ERROR_IF( mpEigenvectorMatrix == nullptr ) << "The eigenvector matrix is not initialized. This means, that the real general generalized " <<
            "eigenvalue problem did not yield complex results or the solution has not yet been performed" << std::endl;
        return *mpEigenvalueVector;
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

    typename TDenseSpaceTypeOut::MatrixPointerType mpEigenvectorMatrix;

    typename TSparseSpaceTypeOut::VectorPointerType mpEigenvalueVector;

    typedef void (feast_ptr)(char*, int*, double*, int*, int*, double*, int*, int*, int*, double*, int*, double*, double*, int*, double*, double*, int*, double*, int*);

    std::function<feast_ptr> CallFeast(double T, bool symmetric)
    {
        using namespace std::placeholders;

        if( symmetric ) {
            return std::bind(dfeast_scsrgv, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19);
        } else {
            return std::bind(dfeast_gcsrgv, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19);
        }
    }

    std::function<feast_ptr> CallFeast(std::complex<double> T, bool symmetric)
    {
        using namespace std::placeholders;

        if( symmetric ) {
            return std::bind(zfeast_scsrgv, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19);
        } else {
            return std::bind(zfeast_gcsrgv, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19);
        }
    }

    double GetE1(double T)
    {
        return mParam["e_min"].GetDouble();
    }

    double GetE2(double T)
    {
        return mParam["e_max"].GetDouble();
    }

    std::complex<double> GetE1(std::complex<double> T)
    {
        return std::complex<double>(mParam["e_mid_re"].GetDouble(), mParam["e_mid_im"].GetDouble());
    }

    double GetE2(std::complex<double> T)
    {
        return mParam["e_r"].GetDouble();
    }

    template<class TVectorType, class TOtherVectorType>
    void RetrieveEigensolution(TVectorType &rFeastEigenvalues, TOtherVectorType &rEigenvalues, FEASTMatrixType& rFeastEigenvectors, DenseMatrixType& rEigenvectors)
    {
        bool complex_result = false;

        // copy real parts of eigenvalues to result vector
        if( rEigenvalues.size() != rFeastEigenvalues.size() )
            rEigenvalues.resize(rFeastEigenvalues.size(), false);

        #pragma omp parallel for
        for( size_t i=0; i<rEigenvalues.size(); ++i ) {
            rEigenvalues[i] = std::real(rFeastEigenvalues[i]);
            if( std::abs(std::imag(rFeastEigenvalues[i])) > 1000*std::numeric_limits<double>::epsilon() ) {
                KRATOS_WARNING("FeastEigensystemSolver") << "Complex eigenvalue detected! eigenvalue[" << i << "] = "
                    << rFeastEigenvalues[i] << std::endl;
                complex_result = true;
            }
        }

        if( !complex_result || mParam["keep_real_solution"].GetBool() ) {
            // copy real parts of eigenvectors to result matrix
            // the eigensolver strategy expects an eigenvector matrix of shape [n_eigenvalues, n_dofs], so FEAST's eigenvector matrix has to be transposed
            if( rEigenvectors.size1() != rFeastEigenvectors.size2() || rEigenvectors.size2() != rFeastEigenvectors.size1() )
                rEigenvectors.resize(rFeastEigenvectors.size2(), rFeastEigenvectors.size1(), false);

            #pragma omp parallel for
            for( size_t i=0; i<rEigenvectors.size1(); ++i ) {
                for( size_t j=0; j<rEigenvectors.size2(); ++j ) {
                    rEigenvectors(i,j) = std::real(rFeastEigenvectors(j,i));
                }
            }
        }

        if( complex_result ) {
            KRATOS_WARNING("FeastEigensystemSolver") << "The computed solution is complex. " <<
                "Please use GetEigenvectorSolution() and GetEigenvalueSolution() to retrieve the complex result." << std::endl;

            mpEigenvectorMatrix = TDenseSpaceTypeOut::CreateEmptyMatrixPointer();
            mpEigenvalueVector = TSparseSpaceTypeOut::CreateEmptyVectorPointer();

            // provide the complex solution
            DenseOutputMatrixType& ev = *mpEigenvectorMatrix;
            if( ev.size1() != rFeastEigenvectors.size2() || ev.size2() != rFeastEigenvectors.size1() )
                ev.resize(rFeastEigenvectors.size2(), rFeastEigenvectors.size1(), false);
            noalias(ev) = trans(rFeastEigenvectors);

            OutputVectorType& e = *mpEigenvalueVector;
            e.swap(rFeastEigenvalues);

            if( !mParam["keep_real_solution"].GetBool() ) {
                // invalidate the double solution
                rEigenvalues.resize(1, false);
                rEigenvalues[0] = std::nan("");
                rEigenvectors.resize(1, 1, false);
                rEigenvectors(0,0) = std::nan("");
            }
        }
    }

    template<class TVectorType, class TOtherVectorType>
    void RetrieveEigensolution(TVectorType &rFeastEigenvalues, TVectorType &rEigenvalues, FEASTMatrixType& rFeastEigenvectors, DenseMatrixType& rEigenvectors)
    {
        // copy eigenvalues to result vector
        rEigenvalues.swap(rFeastEigenvalues);

        // copy eigenvectors to result matrix
        if( rEigenvectors.size1() != rFeastEigenvectors.size1() || rEigenvectors.size2() != rFeastEigenvectors.size2() )
            rEigenvectors.resize(rFeastEigenvectors.size1(), rFeastEigenvectors.size2(), false);

        noalias(rEigenvectors) = rFeastEigenvectors;

        // the eigensolver strategy expects an eigenvector matrix of shape [n_eigenvalues, n_dofs], so FEAST's eigenvector matrix has to be transposed
        rEigenvectors = trans(rEigenvectors);
    }

    void CheckParameters(double T)
    {
        KRATOS_ERROR_IF( mParam["search_lowest_eigenvalues"].GetBool() && mParam["search_highest_eigenvalues"].GetBool() ) <<
            "Cannot search for highest and lowest eigenvalues at the same time" << std::endl;

        KRATOS_ERROR_IF( mParam["e_max"].GetDouble() <= mParam["e_min"].GetDouble() ) <<
            "Invalid eigenvalue limits provided" << std::endl;

        KRATOS_INFO_IF( "FEASTEigensystemSolver",
            mParam["e_mid_re"].GetDouble() != 0.0 || mParam["e_mid_im"].GetDouble() != 0.0 || mParam["e_r"].GetDouble() != 0.0 ) <<
            "Manually defined e_mid_re, e_mid_im, e_r are not used for real symmetric matrices" << std::endl;
    }

    void CheckParameters(std::complex<double> T)
    {
        KRATOS_ERROR_IF( mParam["e_r"].GetDouble() <= 0.0 ) <<
            "Invalid search radius provided" << std::endl;

        KRATOS_INFO_IF( "FEASTEigensystemSolver",
            mParam["e_min"].GetDouble() != 0.0 || mParam["e_max"].GetDouble() != 0.0 ) <<
            "Manually defined e_min, e_max are not used for complex symmetric matrices" << std::endl;

        KRATOS_INFO_IF( "FEASTEigensystemSolver",
            mParam["search_lowest_eigenvalues"].GetBool() || mParam["search_highest_eigenvalues"].GetBool() ) <<
            "Search for extremal eigenvalues is only available for Hermitian problems" << std::endl;

    }

}; // class FEASTEigensystemSolver


/**
 * input stream function
 */
template<bool TSymmetric, typename TScalarIn, typename TScalarOut>
inline std::istream& operator >>(
    std::istream& rIStream,
    FEASTEigensystemSolver<TSymmetric, TScalarIn, TScalarOut>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<bool TSymmetric, typename TScalarIn, typename TScalarOut>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const FEASTEigensystemSolver<TSymmetric, TScalarIn, TScalarOut>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_FEAST_EIGENSYSTEM_SOLVER_H_INCLUDED)
