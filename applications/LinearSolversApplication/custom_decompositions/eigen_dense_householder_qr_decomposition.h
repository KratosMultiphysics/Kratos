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
        // Note that the QR is computed when constructing the pointer
        Eigen::Map<EigenMatrix> eigen_input_matrix_map(rInputMatrix.data().begin(), m, n);
        mpHouseholderQR = Kratos::make_unique<Eigen::HouseholderQR<Eigen::Ref<EigenMatrix>>>(eigen_input_matrix_map);
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
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Check output matrix size
        std::size_t n = rB.size2();
        const std::size_t k = mpHouseholderQR->matrixQR().cols();
        if (rX.size1() != k || rX.size2() != n) {
            rX.resize(k,n,false);
        }

        // Solve the problem Ax = b
        Eigen::Map<EigenMatrix> eigen_x_map(rX.data().begin(), k, n);
        Eigen::Map<EigenMatrix> eigen_rhs_map(rB.data().begin(), rB.size1(), n);
        eigen_x_map = mpHouseholderQR->solve(eigen_rhs_map);
    }

    void Solve(
        const VectorType& rB,
        VectorType& rX) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Check output matrix size
        const std::size_t k = mpHouseholderQR->matrixQR().cols();
        if (rX.size() != k) {
            rX.resize(k,false);
        }

        // Solve the problem Ax = b
        Eigen::Map<EigenMatrix> eigen_x_map(rX.data().begin(), k, 1);
        Eigen::Map<EigenMatrix> eigen_rhs_map(const_cast<VectorType&>(rB).data().begin(), rB.size(), 1);
        eigen_x_map = mpHouseholderQR->solve(eigen_rhs_map);
    }

    void MatrixQ(MatrixType& rMatrixQ) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixQ'." << std::endl;

        // Set the thin Q to be returned
        const std::size_t Q_rows = mpHouseholderQR->householderQ().rows();
        const std::size_t n = mpHouseholderQR->matrixQR().cols();
        if (rMatrixQ.size1() != Q_rows || rMatrixQ.size2() != n) {
            rMatrixQ.resize(Q_rows,n,false);
        }

        // Get the thin unitary matrix Q from the complete one
        // Note that Eigen stores it not as matrix type but as a sequence of Householder transformations (householderQ())
        Eigen::Map<EigenMatrix> thin_Q(rMatrixQ.data().begin(), Q_rows, n);
        thin_Q = EigenMatrix::Identity(Q_rows,n);
        thin_Q = mpHouseholderQR->householderQ() * thin_Q;
    }

    void MatrixR(MatrixType& rMatrixR) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpHouseholderQR) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixR'." << std::endl;

        // Set the matrix R to be returned
        const std::size_t n = mpHouseholderQR->matrixQR().cols();
        if (rMatrixR.size1() != n || rMatrixR.size2() != n) {
            rMatrixR.resize(n,n,false);
        }

        // Get the upper triangular matrix
        // Note that we specify Eigen to return the upper triangular part as the bottom part are auxiliary internal values
        Eigen::Map<EigenMatrix> matrix_R_map(rMatrixR.data().begin(), n, n);
        matrix_R_map = mpHouseholderQR->matrixQR().topLeftCorner(n, n).template triangularView<Eigen::Upper>();
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

    std::unique_ptr<Eigen::HouseholderQR<Eigen::Ref<EigenMatrix>>> mpHouseholderQR = std::unique_ptr<Eigen::HouseholderQR<Eigen::Ref<EigenMatrix>>>(nullptr);

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
