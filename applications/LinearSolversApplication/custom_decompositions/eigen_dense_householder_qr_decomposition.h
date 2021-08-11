//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EIGEN_DENSE_HOUSEHOLDER_QR_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_HOUSEHOLDER_QR_H_INCLUDED

// External includes
#include <Eigen/QR>

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
#include "utilities/dense_qr_decomposition.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template<class TDenseSpace>
class EigenDenseHouseholderQRDecomposition : public DenseQRDecomposition<TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(EigenDenseHouseholderQRDecomposition);

    typedef typename TDenseSpace::DataType DataType;
    typedef typename TDenseSpace::VectorType VectorType;
    typedef typename TDenseSpace::MatrixType MatrixType;

    using EigenVector = Kratos::EigenDynamicVector<DataType>;
    using EigenMatrix = Kratos::EigenDynamicMatrix<DataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    EigenDenseHouseholderQRDecomposition() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static std::string Name()
    {
        return "dense_householder_qr_decomposition";
    }

    void Compute(MatrixType& rInputMatrix) override
    {
        // Householder QR requires m >= n
        const std::size_t m = rInputMatrix.size1();
        const std::size_t n = rInputMatrix.size2();
        KRATOS_ERROR_IF(m < n) << "Householder QR decomposition requires m >= n. Input matrix size is (" << m << "," << n << ")." << std::endl;

        // Compute the Householder QR decomposition
        Eigen::Map<EigenMatrix> eigen_input_matrix_map(rInputMatrix.data().begin(), m, n);
        mHouseholderQR.compute(eigen_input_matrix_map);
    }

    void Compute(
        MatrixType& rInputMatrix,
        MatrixType& rMatrixQ,
        MatrixType& rMatrixR) override
    {
        KRATOS_ERROR << "Householder QR decomposition does not implement the R matrix return. Call the \'Compute()\' and the \'MatrixQ\' methods sequentially." << std::endl;
    }

    void Solve(
        MatrixType& rB,
        MatrixType& rX) const override
    {
        // Solve the problem Ax = b
        Eigen::Map<EigenMatrix> eigen_rhs_map(rB.data().begin(), rB.size1(), rB.size2());
        const auto& r_x = mHouseholderQR.solve(eigen_rhs_map);

        // Set the output matrix
        std::size_t m = r_x.rows();
        std::size_t n = r_x.cols();
        if (rX.size1() != m || rX.size2() != n) {
            rX.resize(m,n);
        }
        Eigen::Map<EigenMatrix> eigen_x_map(rX.data().begin(), m, n);
        eigen_x_map = r_x;
    }

    void Solve(
        const VectorType& rB,
        VectorType& rX) const override
    {
        // Convert the vector data to matrix type
        std::size_t n_rows = rB.size();
        MatrixType aux_input(n_rows, 1);
        MatrixType aux_output(n_rows, 1);
        IndexPartition<std::size_t>(n_rows).for_each([&](std::size_t i){aux_input(i,0) = rB(i);});

        // Call the matrix type method
        Solve(aux_input, aux_output);

        // Convert the matrix output data to vector type
        IndexPartition<std::size_t>(n_rows).for_each([&](std::size_t i){rX(i) = aux_output(i,0);});
    }

    void MatrixQ(MatrixType& rMatrixQ) const override
    {
        // Get the complete unitary matrix
        const auto& r_Q = mHouseholderQR.householderQ();
        const std::size_t Q_rows = r_Q.rows();
        const std::size_t Q_cols = r_Q.cols();
        MatrixType complete_Q(Q_rows, Q_cols);
        Eigen::Map<EigenMatrix> matrix_Q_map(complete_Q.data().begin(), Q_rows, Q_cols);
        matrix_Q_map = r_Q;

        // Calculate the thin Q to be returned
        const std::size_t n = mHouseholderQR.matrixQR().cols();
        const MatrixType aux_identity = IdentityMatrix(Q_cols,n);
        if (rMatrixQ.size1() != Q_rows || rMatrixQ.size2() != n) {
            rMatrixQ.resize(Q_rows,n);
        }
        noalias(rMatrixQ) = prod(complete_Q, aux_identity);
    }

    void MatrixR(MatrixType& rMatrixR) const override
    {
        KRATOS_ERROR << "Householder QR decomposition does not implement the R matrix return" << std::endl;
    }

    void MatrixP(MatrixType& rMatrixP) const override
    {
        KRATOS_ERROR << "Householder QR decomposition does not implement the P matrix return" << std::endl;
    }

    std::size_t Rank() const override
    {
        KRATOS_ERROR << "Householder QR decomposition is not rank revealing." << std::endl;
    }

    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "EigenDecomposition <" << Name() << "> finished.";
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Private static Member Variables
    ///@{


    ///@}
    ///@name Private member Variables
    ///@{

    Eigen::HouseholderQR<EigenMatrix> mHouseholderQR;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private LifeCycle
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{


    ///@}
};

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
} // namespace Kratos

#endif // defined(KRATOS_EIGEN_DENSE_HOUSEHOLDER_QR_H_INCLUDED)
