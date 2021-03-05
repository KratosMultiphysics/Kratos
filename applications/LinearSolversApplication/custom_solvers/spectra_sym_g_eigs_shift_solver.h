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
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>

// Project includes
#include "includes/define.h"
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

        // BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~SpectraSymGEigsShiftSolver() override {}

    /**
     * Solve the generalized eigenvalue problem
     */
    void Solve(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        VectorType& rEigenvalues,
        DenseMatrixType& rEigenvectors) override
    {
        using scalar_t = double;
        using vector_t = Eigen::VectorXd;
        using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        // const int max_iteration = BaseType::GetMaxIterationsNumber();
        // const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();
        const double shift = mParam["shift"].GetDouble();

        // --- wrap ublas matrices

        UblasWrapper<scalar_t> a_wrapper(rK);
        UblasWrapper<scalar_t> b_wrapper(rM);

        // Spectra needs ColMajor ordering
        Eigen::SparseMatrix<scalar_t, Eigen::ColMajor, int> a = a_wrapper.matrix();
        Eigen::SparseMatrix<scalar_t, Eigen::ColMajor, int> b = b_wrapper.matrix();

        // --- timer
        const auto timer = BuiltinTimer();

        KRATOS_INFO_IF("SpectraSymGEigsShiftSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        using OpType = Spectra::SymShiftInvert<scalar_t, Eigen::Sparse, Eigen::Sparse>;
        using BOpType = Spectra::SparseSymMatProd<scalar_t>;

        OpType op(a, b);
        BOpType Bop(b);
        const int ncv = 3 * nroot;  // TODO find a good value
        Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert> eigs(op, Bop, nroot, ncv, shift);

        eigs.init();
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        int nconv = eigs.compute(Spectra::SortRule::LargestAlge, max_iteration);
        int niter = eigs.num_iterations();
        int nops = eigs.num_operations();

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
            double duration = timer.ElapsedSeconds();

            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO("SpectraSymGEigsShiftSolver:") << "Completed in " << duration << " seconds" << std::endl
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
