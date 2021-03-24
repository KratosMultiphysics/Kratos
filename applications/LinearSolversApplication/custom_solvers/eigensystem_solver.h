/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Authors: Thomas Oberbichler
//           Armin Geiser
*/

#if !defined(KRATOS_EIGENSYSTEM_SOLVER_H_INCLUDED)
#define KRATOS_EIGENSYSTEM_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

// Project includes
#include "includes/define.h"
#if defined EIGEN_USE_MKL_ALL
#include "eigen_pardiso_lu_solver.h"
#endif // defined EIGEN_USE_MKL_ALL
#include "eigen_sparse_lu_solver.h"
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
class EigensystemSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigensystemSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 1,
            "normalize_eigenvectors": false,
            "use_mkl_if_available": true,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~EigensystemSolver() override {}

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

        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();


        // --- wrap ublas matrices

        UblasWrapper<scalar_t> a_wrapper(rK);
        UblasWrapper<scalar_t> b_wrapper(rM);

        const auto& a = a_wrapper.matrix();
        const auto& b = b_wrapper.matrix();


        // --- timer
        const auto timer = BuiltinTimer();

        KRATOS_INFO_IF("EigensystemSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        int nn = a.rows();
        int nc = std::min(2 * nroot, nroot + 8);

        // projections
        matrix_t ar(nc, nc);
        matrix_t br(nc, nc);

        // storage for eigenvalues
        vector_t prev_eigv = vector_t::Zero(nc);

        // storage for eigenvectors
        matrix_t r = matrix_t::Zero(nn, nc);
        for (int i = 0; i != nn; ++i) {
            r(i, 0) = b.coeff(i, i);
        }

        vector_t tmp(nn);

        // working vector
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

        Kratos::unique_ptr<DirectSolverWrapperBase> p_solver;

        #if defined USE_EIGEN_MKL
        if (mParam["use_mkl_if_available"].GetBool()) {
            p_solver = Kratos::make_unique<DirectSolverWrapper<EigenPardisoLUSolver<double>>>();
        } else {
            p_solver = Kratos::make_unique<DirectSolverWrapper<EigenSparseLUSolver<double>>>();
        }
        #else  // defined USE_EIGEN_MKL
        p_solver = Kratos::make_unique<DirectSolverWrapper<EigenSparseLUSolver<double>>>();
        #endif // defined USE_EIGEN_MKL

        p_solver->Compute(a);

        int iteration = 0;

        Eigen::GeneralizedSelfAdjointEigenSolver<matrix_t> eig;

        do {
            iteration++;

            KRATOS_INFO_IF("EigensystemSolver:", echo_level > 1) << "Iteration " << iteration <<std::endl;

            for (int j = 0; j != nc; ++j) {
                tmp = r.col(j);
                p_solver->Solve(tmp, tt);

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
                KRATOS_WARNING("EigensystemSolver:") << "Eigen solution was not successful!" << std::endl;
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
                    KRATOS_WARNING_IF("EigensystemSolver:", echo_level > 1) << "Convergence not reached for eigenvalue #"<<i+1<<": " << rtolv <<"." << std::endl;
                    break;
                }
            }

            if (is_converged) {
                KRATOS_INFO_IF("EigensystemSolver:", echo_level > 0) << "Convergence reached after " << iteration << " iterations within a relative tolerance: " << tolerance << std::endl;
                break;
            } else if (iteration >= max_iteration) {
                KRATOS_INFO_IF("EigensystemSolver:", echo_level > 0) << "Convergence not reached in " << max_iteration << " iterations." << std::endl;
                break;
            }

            prev_eigv = eig.eigenvalues();
        } while (true);


        if (static_cast<int>(rEigenvalues.size()) != nroot) {
            rEigenvalues.resize(nroot);
        }
        if (static_cast<int>(rEigenvectors.size1()) != nroot || static_cast<int>(rEigenvectors.size2()) != nn) {
            rEigenvectors.resize(nroot, nn);
        }

        Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());

        eigvals = eig.eigenvalues().head(nroot);

        for (int i = 0; i != nroot; ++i) {
            tmp = r.col(i);
            p_solver->Solve(tmp, eigvecs.row(i));
            eigvecs.row(i).normalize();
        }

        // --- normalization
        // Given generalized eigenvalue problem (A - eigenvalue * B) * eigenvector = 0,
        // eigenvector is normalized such that eigenvector^T * B * eigenvector = 1
        if(mParam["normalize_eigenvectors"].GetBool())
        {
            for (int i = 0; i != nroot; ++i)
            {
                const double tmp = eigvecs.row(i) * b * eigvecs.row(i).transpose();
                const double factor = 1.0 / std::sqrt(tmp);
                eigvecs.row(i) *=  factor;
                KRATOS_INFO_IF("EigensystemSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
            }
        }

        // --- output
        if (echo_level > 0) {
            double duration = timer.ElapsedSeconds();

            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO("EigensystemSolver:") << "Completed in " << duration << " seconds" << std::endl
                      << "                   Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "EigensystemSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

private:

    struct DirectSolverWrapperBase
    {
        typedef Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, int>> MatrixMapType;
        typedef Eigen::Matrix<double, Eigen::Dynamic, 1> EigenVectorType;
        typedef Eigen::Ref<const EigenVectorType> ConstVectorRefType;
        typedef Eigen::Ref<EigenVectorType> VectorRefType;

        virtual ~DirectSolverWrapperBase() = default;
        virtual void Compute(MatrixMapType a) = 0;
        virtual void Solve(ConstVectorRefType b, VectorRefType x) = 0;
    };

    template<class TSolver>
    struct DirectSolverWrapper : DirectSolverWrapperBase
    {
        typedef DirectSolverWrapperBase BaseType;
        typedef typename BaseType::MatrixMapType MatrixMapType;
        typedef typename BaseType::VectorRefType VectorRefType;
        typedef typename BaseType::ConstVectorRefType ConstVectorRefType;

        void Compute(MatrixMapType a) override {mSolver.Compute(a);}
        void Solve(ConstVectorRefType b, VectorRefType x) override {mSolver.Solve(b, x);}
        private:
            TSolver mSolver;
    };

}; // class EigensystemSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    EigensystemSolver<TSparseSpaceType,
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
    const EigensystemSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_EIGENSYSTEM_SOLVER_H_INCLUDED)
