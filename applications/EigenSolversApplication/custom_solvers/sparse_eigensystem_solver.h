/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Authors: Thomas Oberbichler
//           Armin Geiser
*/

#if !defined(KRATOS_SPARSE_EIGENSYSTEM_SOLVER_H_INCLUDED)
#define KRATOS_SPARSE_EIGENSYSTEM_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#if defined EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SparseEigensystemSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(SparseEigensystemSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    SparseEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "SparseEigensystemSolver",
            "number_of_eigenvalues": 1,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~SparseEigensystemSolver() override {}

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
        using scalar_t = double;
        using vector_t = Eigen::VectorXd;
        using matrix_t = Eigen::MatrixXd;
        using sparse_t = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
	    #if defined EIGEN_USE_MKL_ALL
        using ldlt_solver_t = Eigen::PardisoLDLT<sparse_t>;
        #else
        using ldlt_solver_t = Eigen::SparseLU<sparse_t>;
        #endif
        
        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt(); // number of eigenvalues requested
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();


        // --- wrap ublas matrices
        
        std::vector<int> index1_vector_a(rK.index1_data().size());
        std::vector<int> index2_vector_a(rK.index2_data().size());

        for (size_t i = 0; i < rK.index1_data().size(); i++) {
            index1_vector_a[i] = (int)rK.index1_data()[i];
        }

        for (size_t i = 0; i < rK.index2_data().size(); i++) {
            index2_vector_a[i] = (int)rK.index2_data()[i];
        }
        
        Eigen::Map<sparse_t> a(rK.size1(), rK.size2(), rK.nnz(), index1_vector_a.data(), index2_vector_a.data(), rK.value_data().begin());


        std::vector<int> index1_vector_b(rM.index1_data().size());
        std::vector<int> index2_vector_b(rM.index2_data().size());

        for (size_t i = 0; i < rM.index1_data().size(); i++) {
            index1_vector_b[i] = (int)rM.index1_data()[i];
        }

        for (size_t i = 0; i < rM.index2_data().size(); i++) {
            index2_vector_b[i] = (int)rM.index2_data()[i];
        }
        
        Eigen::Map<sparse_t> b(rM.size1(), rM.size2(), rM.nnz(), index1_vector_b.data(), index2_vector_b.data(), rM.value_data().begin());


        // --- timer      

        double start_time = OpenMPUtils::GetCurrentTime();

        if (echo_level > 0) {
            std::cout << "SparseEigensystemSolver: Start"  << std::endl;
        }


        // --- calculation

        int nn = a.rows();
        int nc = std::min(2 * nroot, nroot + 8);
        
        bool eigen_solver_successful = true;

        // projections
        matrix_t ar(nc, nc);
        matrix_t br(nc, nc);

        // storage for eigenvalues
        vector_t prev_eigv(nc);

        // storage for eigenvectors
        matrix_t r = matrix_t::Zero(nn, nc);
        for (int i = 0; i != nn; ++i) {
            r(i, 0) = b.coeff(i, i);
        }

        // woaing vector
        vector_t w(nn);
        for (int i = 0; i != nn; ++i) {
            w(i) = r(i, 0) / a.coeff(i, i);
        }

        int nd = nn / nc;
        int l = nn - nd;

        vector_t tt(nn);
        int ij = 0;

        tt(0) = 0.0;

        for (int j = 1; j != nc; ++j) {
            double rt = 0.0;

            for (int i = 0; i != l; ++i) {
                if (w(i) >= rt) {
                    rt = w(i);
                    ij = i;
                }
            }

            for (int i = l - 1; i != nn; ++i) {
                if (w(i) > rt) {
                    rt = w(i);
                    ij = i;
                }
            }

            tt(j) = ij;
            w(ij) = 0.0;

            l -= nd;

            r(ij, j) = 1.0;
        }
       
        ldlt_solver_t solver;
        solver.compute(a);

        int iteration = 0;

        Eigen::GeneralizedSelfAdjointEigenSolver<matrix_t> eig;

        do {   
            iteration++;

            if (echo_level > 1) {
                std::cout << "SparseEigensystem: Iteration " << iteration <<std::endl;
            }

            for (int j = 0; j != nc; ++j) {
                vector_t tmp(nn);
                tmp = r.col(j);
                tt = solver.solve(tmp);

                for (int i = j; i != nc; ++i) {                    
                    ar(i, j) = r.col(i).dot(tt);
                }

                r.col(j) = tt;
            }

            for (int j = 0; j != nc; ++j) {
                tt = b * r.col(j);

                for (int i = j; i != nc; ++i) {
                    br(i, j) = r.col(i).dot(tt);
                }

                r.col(j) = tt;
            }

            eig.compute(ar, br);

            if(eig.info() != Eigen::Success) {
                std::cout << "SparseEigensystem: Eigen solution was not successful!" << std::endl;
                break;
            }

            r *= eig.eigenvectors();

            bool is_converged = true;

            for (int i = 0; i != nc; i++) {
                double eigv = eig.eigenvalues()(i);
                double dif = eigv - prev_eigv(i);
                double rtolv = std::abs(dif / eigv);

                if (rtolv > tolerance) {
                    is_converged = false;
                    break;
                }
            }

            if (is_converged) {
                if (echo_level > 0) {
                    std::cout << "SparseEigensystem: Convergence reached after " << iteration << " iterations within a relative tolerance: " << tolerance << std::endl;
                }
                break;
            } else if (iteration >= max_iteration) {
                if (echo_level > 0) {
                    std::cout << "SparseEigensystem: Convergence not reached in " << max_iteration << " iterations." << std::endl;
                }
                break;
            }

            prev_eigv = eig.eigenvectors();
        } while (true);

        Eigen::Map<vector_t> (rEigenvalues.data().begin(), rEigenvalues.size()) = eig.eigenvalues().head(nroot);
        Eigen::Map<matrix_t> (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2()) = solver.solve(r.leftCols(nroot));
        

        // --- timer     

        if (echo_level > 0) {
            double end_time = OpenMPUtils::GetCurrentTime();
            double duration = end_time - start_time;
            std::cout << "SparseEigensystemSolver: Completed in " << duration << " seconds" << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "SparseEigensystemSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }
}; // class SparseEigensystemSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    SparseEigensystemSolver<TSparseSpaceType,
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
    const SparseEigensystemSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_SPARSE_EIGENSYSTEM_SOLVER_H_INCLUDED)