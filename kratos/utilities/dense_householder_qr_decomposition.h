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

#if !defined(KRATOS_DENSE_HOUSEHOLDER_QR_H_INCLUDED)
#define KRATOS_DENSE_HOUSEHOLDER_QR_H_INCLUDED

// External includes
#include "amgcl/detail/qr.hpp"

// Project includes
#include "dense_qr_decomposition.h"

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

template <class TDenseSpaceType>
class DenseHouseholderQRDecomposition : public DenseQRDecomposition<TDenseSpaceType>
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(DenseHouseholderQRDecomposition);

    using DataType = typename TDenseSpaceType::DataType;
    using VectorType = typename TDenseSpaceType::VectorType;
    using MatrixType = typename TDenseSpaceType::MatrixType;
    using AMGCLQRType = amgcl::detail::QR<DataType>;

    ///@}
    ///@name Life Cycle
    ///@{

    DenseHouseholderQRDecomposition() = default;

    virtual ~DenseHouseholderQRDecomposition() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Name of the QR
     * Returns a string containing the name of the current QR decomposition
     * @return std::string Name of the QR decomposition
     */
    static std::string Name()
    {
        return "dense_householder_qr_decomposition";
    }

    /**
     * @brief Compute the QR
     * Computes the QR Decomposition (QR) of the given imput matrix
     * Note that the input matrix is modidifed
     * @param rInputMatrix Matrix to compute the QR decomposition
     */
    void Compute(MatrixType& rInputMatrix) override
    {
        // Set input data
        // Note that we need a copy as QR decomposition is modifying the input
        mpA = &rInputMatrix;
        DataType *p_0_0 = &((*mpA)(0,0));
        const std::size_t m = rInputMatrix.size1();
        const std::size_t n = rInputMatrix.size2();

        // Compute the Householder QR decomposition
        mHouseholderQR.compute(m, n, p_0_0);
    }

    /**
     * @brief Compute the QR
     * Computes the QR (QR) of the given input matrix
     * Note that the input matrix is modidifed
     * @param rInputMatrix Matrix to compute the QR decomposition
     * @param rMatrixQ Unitary matrix
     * @param rMatrixR Upper triangular matrix
     */
    void Compute(
        MatrixType& rInputMatrix,
        MatrixType& rMatrixQ,
        MatrixType& rMatrixR) override
    {
        // Set input data
        // Note that we need a copy as QR decomposition is modifying the input
        mpA = &rInputMatrix;
        DataType *p_0_0 = &((*mpA)(0,0));
        const std::size_t m = rInputMatrix.size1();
        const std::size_t n = rInputMatrix.size2();

        // Compute the Householder QR decomposition
        mHouseholderQR.factorize(m, n, p_0_0);

        // Check sizes and fill Q values
        if (rMatrixQ.size1() != m || rMatrixQ.size2() != n) {
            rMatrixQ.resize(m,n,false);
        }
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                rMatrixQ(i,j) = mHouseholderQR.Q(i,j);
            }
        }

        // Check sizes and fill R values
        MatrixR(rMatrixR);
    }

    /**
     * @brief Solves the problem Ax=b
     * Being A the input matrix, this method solves the problem Ax = b
     * @param rB The Right Hand Side (RHS) matrix
     * @param rX The solution matrix
     */
    void Solve(
        MatrixType& rB,
        MatrixType& rX) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpA) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Set input data
        DataType* p_A_0_0 = &((*mpA)(0,0));
        const std::size_t m = mpA->size1();
        const std::size_t n = mpA->size2();

        // Check output matrix size
        const std::size_t l = rB.size2();
        if (rX.size1() != n || rX.size2() != l) {
            rX.resize(n, l, false);
        }

        // Call the QR solve method for each column
        VectorType aux_B_col(l);
        VectorType aux_X_col(l);
        const bool qr_computed = true;
        auto storage_order = amgcl::detail::storage_order::row_major;
        for (std::size_t i_col = 0; i_col < l; ++i_col) {
            TDenseSpaceType::GetColumn(i_col, rB, aux_B_col);
            TDenseSpaceType::GetColumn(i_col, rX, aux_X_col);
            DataType* p_b_0 = &(aux_B_col(0));
            DataType* p_x_0 = &(aux_X_col(0));
            const_cast<AMGCLQRType&>(mHouseholderQR).solve(m, n, p_A_0_0, p_b_0, p_x_0, storage_order, qr_computed);
            TDenseSpaceType::SetColumn(i_col, rX, aux_X_col);
        }
    }

    /**
     * @brief Solves the problem Ax=b
     * Being A the input matrix, this method solves the problem Ax = b
     * @param rB The Right Hand Side (RHS) vector
     * @param rX The solution vector
     */
    void Solve(
        const VectorType& rB,
        VectorType& rX) const override
    {
        // Check that QR decomposition has been already computed
        KRATOS_ERROR_IF(!mpA) << "QR decomposition not computed yet. Please call 'Compute' before 'Solve'." << std::endl;

        // Set input data
        DataType* p_b_0 = &((const_cast<VectorType&>(rB))(0));
        DataType* p_x_0 = &(rX(0));
        DataType *p_A_0_0 = &((*mpA)(0,0));
        const std::size_t m = mpA->size1();
        const std::size_t n = mpA->size2();

        // Check output vector size
        if (rX.size() != n) {
            rX.resize(n, false);
        }

        // Call the QR solve method
        const bool qr_computed = true;
        auto storage_order = amgcl::detail::storage_order::row_major;
        const_cast<AMGCLQRType&>(mHouseholderQR).solve(m, n, p_A_0_0, p_b_0, p_x_0, storage_order, qr_computed);
    }

    /**
     * @brief Unitary matrix getter
     * If computed, this method sets the unitary matrix in the provided array
     * @param rMatrixQ Unitary matrix
     */
    void MatrixQ(MatrixType& rMatrixQ) const override
    {
        // We intentionally avoid to implement this method as it requires to call the factorize of QR include
        // Otherwise we would need to always do the factorize call, which is more expensive than the simpler compute one
        KRATOS_ERROR << "This method is not implemented. Please use the 'Compute' interface with Q and R return." << std::endl;
    }

    /**
     * @brief Upper triangular matrix getter
     * If computed, this method sets the upper triangular matrix in the provided array
     * @param rMatrixR Upper triangular matrix
     */
    void MatrixR(MatrixType& rMatrixR) const override
    {
        // Check that QR has been already computed
        KRATOS_ERROR_IF(!mpA) << "QR decomposition not computed yet. Please call 'Compute' before 'MatrixR'." << std::endl;

        // Check input size
        const std::size_t n = mpA->size2();
        if (rMatrixR.size1() != n || rMatrixR.size2() != n) {
            rMatrixR.resize(n,n,false);
        }

        // Get R values from QR util
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                rMatrixR(i,j) = mHouseholderQR.R(i,j);
            }
        }
    }

    /**
     * @brief Pivoting matrix getter
     * If computed, this method sets the pivoting matrix
     * @param rMatrixP Pivoting matrix
     */
    void MatrixP(MatrixType& rMatrixP) const override
    {
        KRATOS_ERROR << "Householder QR decomposition does not implement the P matrix return" << std::endl;
    }

    /**
     * @brief Rank of the provided array
     * Calculates and returns the rank of the array decomposed with the QR
     * @return std::size_t Rank of the provided array
     */
    std::size_t Rank() const override
    {
        KRATOS_ERROR << "Householder QR decomposition is not rank revealing." << std::endl;
    }

    /**
     * @brief QR information
     * Outputs the QR class information
     * @param rOStream Information output
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "Decomposition <" << Name() << "> finished.";
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

    AMGCLQRType mHouseholderQR;

    MatrixType* mpA = nullptr;

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

#endif // defined(KRATOS_DENSE_HOUSEHOLDER_QR_H_INCLUDED)
