/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Manuel Messmer
*/

#if !defined(KRATOS_EIGEN_DENSE_EIGENVALUE_SOLVER_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_EIGENVALUE_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "spaces/ublas_space.h"
#include "utilities/builtin_timer.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos {

template <typename TScalar = double,
          class TSparseSpaceType = TUblasDenseSpace<TScalar>,
          class TDenseSpaceType = TUblasDenseSpace<TScalar>>
class DenseEigenvalueSolver :
    public LinearSolver<TSparseSpaceType, TDenseSpaceType>
{
    Parameters mParam;

    int mEchoLevel;

public:
    KRATOS_CLASS_POINTER_DEFINITION(DenseEigenvalueSolver);

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    DenseEigenvalueSolver( Parameters param) : mParam(param){

        Parameters default_params(R"(
        {
            "solver_type"       : "dense_eigensolver",
            "ascending_order"   : true,
            "echo_level"        : 0
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        mEchoLevel = mParam["echo_level"].GetInt();
    }

    ~DenseEigenvalueSolver() override {}

    /**
     * @brief Dense eigenvalue solver
     * @details Computes eigenvalues and eigenvectors of selfadjoint matrices.
     * @param rA System matrix
     * @param rDummy Dummy matrix
     * @param rEigenvalues Eigenvalue vector
     * @param rEigenvectors Eigenvector matrix
     */
    void Solve( DenseMatrixType& rA,
                DenseMatrixType& rDummy,
                DenseVectorType& rEigenvalues,
                DenseMatrixType& rEigenvectors) override
    {

        BuiltinTimer eigensolver_timer;

        KRATOS_INFO_IF("DenseEigenvalueSolver", mEchoLevel > 0) << "Start"  << std::endl;

        using vector_t = Kratos::EigenDynamicVector<TScalar>;
        using matrix_t = Kratos::EigenDynamicMatrix<TScalar>;

        Eigen::Map<matrix_t> A(rA.data().begin(), rA.size1(), rA.size2());

        Eigen::SelfAdjointEigenSolver<matrix_t> solver;

        solver.compute(A, Eigen::ComputeEigenvectors);

        rEigenvalues.resize(rA.size1());
        rEigenvectors.resize(rA.size1(), rA.size1());

        Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());

        if( mParam["ascending_order"].GetBool() ){
            eigvals = solver.eigenvalues();
            eigvecs = solver.eigenvectors();
        }
        else{
            eigvals = solver.eigenvalues().reverse();
            eigvecs = solver.eigenvectors().rowwise().reverse();
        }

        const bool success = solver.info() == Eigen::Success;

        if( success ){
            KRATOS_INFO_IF("DenseEigenvalueSolver", mEchoLevel > 0) << "Completed in "
                << eigensolver_timer << std::endl;

            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO_IF("", mEchoLevel > 1)
                << "                  Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        }
        else{
            KRATOS_WARNING("DenseEigenvalueSolver") << "Decomposition failed!" << std::endl;
        }
    }

    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DenseEigenvalueSolver";
    }

};

} // namespace Kratos

#endif // defined(KRATOS_EIGEN_DENSE_EIGENVALUE_SOLVER_H_INCLUDED)