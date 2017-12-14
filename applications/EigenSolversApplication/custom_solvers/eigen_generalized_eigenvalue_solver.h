/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Armin Geiser
*/

#if !defined(EIGEN_GENERALIZED_EIGEN_SOLVER_H_INCLUDED)
#define EIGEN_GENERALIZED_EIGEN_SOLVER_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/direct_solver.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// This solver solves the generalized eigenvalue problem using Eigen.
/// it solves for all eigenvalues
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class EigenGeneralizedEigenvalueSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenGeneralizedEigenvalueSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenGeneralizedEigenvalueSolver(Parameters param) : mParam(param)
    {

        Parameters default_params(R"(
        {
            "solver_type": "GeneralizedEigenSolver",
            "verbosity": 1
        })");

        // don't validate linear_solver_settings here
        mParam.ValidateAndAssignDefaults(default_params);
    }

    ~EigenGeneralizedEigenvalueSolver() override {}


    /**
     * Solve for all eigenvalues and vectors of the generalized eigenvalue problem.
     * for self adjoint version:
     * rA is self adjoint matrix; rB positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rA,
            SparseMatrixType& rB,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        const int verbosity = mParam["verbosity"].GetInt();

        if (rA.size1() > 100)
        {
            std::cout << "WARNING: EigenGeneralizedEigenvalueSolver solves for all "<<
                        "eigen values of a dense matrix. "<<
                        "This might take long for large matrices!" << std::endl;
        }

        if (verbosity>1) std::cout << "Start solving for eigenvalues" << evalues << std::endl;
        // create Eigen matrix A
        Eigen::MatrixXd A(rA.size1(),rA.size2());
        Eigen::MatrixXd B(rA.size1(),rA.size2());

        for (size_t i=0; i<rA.size1(); i++)
        {
            for (size_t j=0; j<rA.size2(); j++)
            {
                A(i,j) = rA(i,j);
                B(i,j) = rB(i,j);
            }
        }

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);

        Eigen::VectorXd evalues = ges.eigenvalues();
        Eigen::MatrixXd evecs = ges.eigenvectors();

        if(ges.info() != Eigen::Success)
        {
            KRATOS_ERROR << "Eigen solution was not successful!" << std::endl;
        }

        if (verbosity>1) std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;

        rEigenvalues.resize(evalues.rows());
        rEigenvectors.resize(evecs.rows(), evecs.cols());

        for (int i=0; i<evalues.rows(); i++){
            rEigenvalues(i) = evalues(i);
            for (int j=0; j<evecs.cols(); j++){
                rEigenvectors(i,j) = evecs(i,j);
            }
        }

        return;
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "EigenGeneralizedEigenvalueSolver.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///@}

  private:
    ///@name Member Variables
    ///@{

    Parameters mParam;

    ///@}

}; // class EigenGeneralizedEigenvalueSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        EigenGeneralizedEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const EigenGeneralizedEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(EIGEN_GENERALIZED_EIGEN_SOLVER_H_INCLUDED)