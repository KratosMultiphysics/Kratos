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
#include <Eigen/Sparse>
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

/// Adapter Eigen generalized eigenvalue problem solvers.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class EigenGeneralizedEigenSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

  public:
    KRATOS_CLASS_POINTER_DEFINITION(EigenGeneralizedEigenSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    EigenGeneralizedEigenSolver(Parameters::Pointer pParam) : mpParam(pParam)
    {

        Parameters default_params(R"(
        {
            "solver_type": "GeneralizedEigenSolver",
            "verbosity": 1
        })");

        // don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);
    }

    ~EigenGeneralizedEigenSolver() override {}


    /**
     * Solve the generalized eigenvalue problem.
     * for self adjoint version:
     * rA is self adjoint matrix; rB positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rA,
            SparseMatrixType& rB,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        // settings
        //Parameters& settings = *mpParam;
        const int verbosity = 2;

        if (rA.size1() > 100)
        {
            std::cout << "WARNING: EigenGeneralizedEigenSolver solves for all "<<
                        "eigen values of a dense matrix. "<<
                        "This might take long for large matrices!" << std::endl;
        }

        // create Eigen matrix A
        Eigen::MatrixXd A(rA.size1(),rA.size2());
        Eigen::MatrixXd B(rA.size1(),rA.size2());

        for (size_t i=0; i<rA.size1(); i++)
        {
            for (size_t j=0; j<rA.size2(); j++)
            {
                A(i,j) = double(rA(i,j));
                B(i,j) = double(rB(i,j));
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
        rOStream << "EigenGeneralizedEigenSolver.";
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

    Parameters::Pointer mpParam;

    ///@}

}; // class EigenGeneralizedEigenSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        EigenGeneralizedEigenSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const EigenGeneralizedEigenSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(EIGEN_GENERALIZED_EIGEN_SOLVER_H_INCLUDED)