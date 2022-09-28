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

#if !defined(KRATOS_EIGEN_DENSE_COLUMN_PIVOTING_HOUSEHOLDER_QR_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_COLUMN_PIVOTING_HOUSEHOLDER_QR_H_INCLUDED

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
class EigenDenseColumnPivotingHouseholderQRDecomposition : public DenseQRDecomposition<TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(EigenDenseColumnPivotingHouseholderQRDecomposition);

    typedef typename TDenseSpace::DataType DataType;
    typedef typename TDenseSpace::VectorType VectorType;
    typedef typename TDenseSpace::MatrixType MatrixType;

    using EigenVector = Kratos::EigenDynamicVector<DataType>;
    using EigenMatrix = Kratos::EigenDynamicMatrix<DataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    EigenDenseColumnPivotingHouseholderQRDecomposition() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static std::string Name()
    {
        return "dense_column_pivoting_householder_qr_decomposition";
    }

    void Compute(MatrixType& rInputMatrix) override
    {
        // Householder QR requires m >= n
        const std::size_t m = rInputMatrix.size1();
        const std::size_t n = rInputMatrix.size2();
        KRATOS_ERROR_IF(m < n) << "Householder QR decomposition requires m >= n. Input matrix size is (" << m << "," << n << ")." << std::endl;

        // Compute the Householder QR decomposition
        // Note that the QR is computed when constructing the pointer
        Eigen::Map<EigenMatrix> eigen_input_matrix_map(rInputMatrix.data().begin(), m, n);
        mpColPivHouseholderQR = Kratos::make_unique<Eigen::ColPivHouseholderQR<Eigen::Ref<EigenMatrix>>>(eigen_input_matrix_map);
    }

    void Compute(
        MatrixType& rInputMatrix,
        MatrixType& rMatrixQ,
        MatrixType& rMatrixR) override
    {
        Compute(rInputMatrix);
        MatrixQ(rMatrixQ);
        MatrixR(rMatrixR);
    }

    void Solve(
        MatrixType& rB,
        MatrixType& rX) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Check output matrix size
        std::size_t n = rB.size2();
        const std::size_t rank = Rank();
        if (rX.size1() != rank || rX.size2() != n) {
            rX.resize(rank,n,false);
        }

        // Solve the problem Ax = b
        Eigen::Map<EigenMatrix> eigen_x_map(rX.data().begin(), rank, n);
        Eigen::Map<EigenMatrix> eigen_rhs_map(rB.data().begin(), rB.size1(), n);
        eigen_x_map = mpColPivHouseholderQR->solve(eigen_rhs_map);
    }

    void Solve(
        const VectorType& rB,
        VectorType& rX) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Check output matrix size
        const std::size_t rank = Rank();
        if (rX.size() != rank) {
            rX.resize(rank,false);
        }

        // Solve the problem Ax = b
        Eigen::Map<EigenMatrix> eigen_x_map(rX.data().begin(), rank, 1);
        Eigen::Map<EigenMatrix> eigen_rhs_map(const_cast<VectorType&>(rB).data().begin(), rB.size(), 1);
        eigen_x_map = mpColPivHouseholderQR->solve(eigen_rhs_map);
    }

    void MatrixQ(MatrixType& rMatrixQ) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixQ'." << std::endl;

        // Set the thin Q to be returned
        const std::size_t Q_rows = mpColPivHouseholderQR->householderQ().rows();
        const std::size_t rank = Rank();
        if (rMatrixQ.size1() != Q_rows || rMatrixQ.size2() != rank) {
            rMatrixQ.resize(Q_rows,rank,false);
        }

        // Get the thin unitary matrix Q from the complete one
        // Note that Eigen stores it not as matrix type but as a sequence of Householder transformations (householderQ())
        Eigen::Map<EigenMatrix> thin_Q(rMatrixQ.data().begin(), Q_rows, rank);
        thin_Q = EigenMatrix::Identity(Q_rows,rank);
        thin_Q = mpColPivHouseholderQR->householderQ() * thin_Q;
    }

    void MatrixR(MatrixType& rMatrixR) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixR'." << std::endl;

        // Set the matrix R to be returned
        const std::size_t rank = Rank();
        if (rMatrixR.size1() != rank || rMatrixR.size2() != rank) {
            rMatrixR.resize(rank,rank,false);
        }

        // Get the upper triangular matrix
        // Note that we specify Eigen to return the upper triangular part as the bottom part are auxiliary internal values
        Eigen::Map<EigenMatrix> matrix_R_map(rMatrixR.data().begin(), rank, rank);
        matrix_R_map = mpColPivHouseholderQR->matrixR().topLeftCorner(rank, rank).template triangularView<Eigen::Upper>();
    }

    void MatrixP(MatrixType& rMatrixP) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixP'." << std::endl;

        // Get the permutation matrix
        const auto& r_P = mpColPivHouseholderQR->colsPermutation();
        const std::size_t m = r_P.rows();
        const std::size_t n = r_P.cols();

        // Output the permutation matrix
        if (rMatrixP.size1() != m || rMatrixP.size2() != n) {
            rMatrixP.resize(m,n,false);
        }
        Eigen::Map<EigenMatrix> matrix_P_map(rMatrixP.data().begin(), m, n);
        matrix_P_map = r_P;
    }

    std::size_t Rank() const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpColPivHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'Rank'." << std::endl;

        return mpColPivHouseholderQR->rank();
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

    std::unique_ptr<Eigen::ColPivHouseholderQR<Eigen::Ref<EigenMatrix>>> mpColPivHouseholderQR = std::unique_ptr<Eigen::ColPivHouseholderQR<Eigen::Ref<EigenMatrix>>>(nullptr);

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

#endif // defined(KRATOS_EIGEN_DENSE_COLUMN_PIVOTING_HOUSEHOLDER_QR_H_INCLUDED)
