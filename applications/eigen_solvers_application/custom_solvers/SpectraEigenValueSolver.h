/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#if !defined(KRATOS_SPECTRA_EIGEN_VALUE_SOLVER_H_INCLUDED)
#define KRATOS_SPECTRA_EIGEN_VALUE_SOLVER_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <SymGEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <MatOp/SparseCholesky.h>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/direct_solver.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Adapter to Spectra eigenvalue problem solver using Eigen direct solvers.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SpectraEigenValueSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

  public:
    KRATOS_CLASS_POINTER_DEFINITION(SpectraEigenValueSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> EigenSparseMatrixTypeRow;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSparseMatrixTypeCol;

    SpectraEigenValueSolver(Parameters::Pointer pParam) : mpParam(pParam)
    {

        Parameters default_params(R"(
        {
            "solver_type": "SpectraEigenValueSolver",
            "number_of_eigenvalues": 1,
            "ncv": -1,
            "linear_solver_settings": {}
        })");

        // don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);
    }

    ~SpectraEigenValueSolver() override {}


    /**
     * Solve the generalized eigenvalue problem.
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rK,
            SparseMatrixType& rM,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        // settings
        Parameters& settings = *mpParam;

        const int n_rows = rK.size1();
        const int nev = settings["number_of_eigenvalues"].GetInt(); // number of eigenvalues requested
        int ncv;
        if (settings["ncv"].GetInt() == -1)
            ncv = std::min(nev*3, n_rows-1);
        else
            ncv = settings["ncv"].GetInt();
        KRATOS_ERROR_IF(nev>n_rows-1) << "Number of Eigenvalues has to be smaller then matrix size" << std::endl;
        KRATOS_ERROR_IF(ncv>n_rows-1) << "ncv has to be smaller then matrix size" << std::endl;

        // create Eigen matrix A
        std::vector<int> index1_vector_A(rK.index1_data().size());
        std::vector<int> index2_vector_A(rK.index2_data().size());

        for (size_t i = 0; i < rK.index1_data().size(); i++) {
            index1_vector_A[i] = (int)rK.index1_data()[i];
        }

        for (size_t i = 0; i < rK.index2_data().size(); i++) {
            index2_vector_A[i] = (int)rK.index2_data()[i];
        }
        Eigen::Map< EigenSparseMatrixTypeRow > Arow(rK.size1(), rK.size2(), rK.nnz(), index1_vector_A.data(), index2_vector_A.data(), rK.value_data().begin());
        EigenSparseMatrixTypeCol A = Arow;

        // create Eigen matrix B
        std::vector<int> index1_vector_B(rM.index1_data().size());
        std::vector<int> index2_vector_B(rM.index2_data().size());

        for (size_t i = 0; i < rM.index1_data().size(); i++) {
            index1_vector_B[i] = (int)rM.index1_data()[i];
        }

        for (size_t i = 0; i < rM.index2_data().size(); i++) {
            index2_vector_B[i] = (int)rM.index2_data()[i];
        }
        Eigen::Map< EigenSparseMatrixTypeRow > Brow(rM.size1(), rM.size2(), rM.nnz(), index1_vector_B.data(), index2_vector_B.data(), rM.value_data().begin());
        EigenSparseMatrixTypeCol B = Brow;

        typedef Spectra::SparseSymMatProd<double> OpType;
        typedef Spectra::SparseCholesky<double> BOpType; // ->Spectra::GEIGS_CHOLESKY
        //typename Spectra::SparseRegularInverse<double> BOpType; // -> Spectra::GEIGS_REGULAR_INVERSE

        OpType op(A);
        BOpType Bop(B);
        Spectra::SymGEigsSolver<double, Spectra::SMALLEST_MAGN, OpType, BOpType, Spectra::GEIGS_CHOLESKY>
            geigs(&op, &Bop, nev, ncv);

        // Initialize and compute
        geigs.init();
        geigs.compute();

        // Retrieve results
        Eigen::VectorXd evalues;
        Eigen::MatrixXd evecs;
        if(geigs.info() != Spectra::SUCCESSFUL)
        {
            KRATOS_ERROR << "Solution of Eigenvalues using Spectra failed!"<< std::endl;
        }

        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();

        std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;

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
        rOStream << "SpectraEigenValueSolver.";
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

}; // class SpectraEigenValueSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        SpectraEigenValueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const SpectraEigenValueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_SPECTRA_EIGEN_VALUE_SOLVER_H_INCLUDED)