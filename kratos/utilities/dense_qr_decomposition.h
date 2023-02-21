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

#if !defined(KRATOS_DENSE_QR_H_INCLUDED)
#define KRATOS_DENSE_QR_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"

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

template<class TDenseSpaceType>
class DenseQRDecomposition
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(DenseQRDecomposition);

    typedef typename TDenseSpaceType::DataType DataType;
    typedef typename TDenseSpaceType::VectorType VectorType;
    typedef typename TDenseSpaceType::MatrixType MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    DenseQRDecomposition() = default;

    virtual ~DenseQRDecomposition() = default;

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
        return "dense_qr_decomposition";
    }

    /**
     * @brief Compute the QR
     * Computes the QR Decomposition (QR) of the given imput matrix
     * @param rInputMatrix Matrix to compute the QR decomposition
     */
    virtual void Compute(MatrixType& rInputMatrix) = 0;

    /**
     * @brief Compute the QR
     * Computes the QR (QR) of the given input matrix
     * @param rInputMatrix Matrix to compute the QR decomposition
     * @param rMatrixQ Unitary matrix
     * @param rMatrixR Upper triangular matrix
     */
    virtual void Compute(
        MatrixType& rInputMatrix,
        MatrixType& rMatrixQ,
        MatrixType& rMatrixR) = 0;

    /**
     * @brief Solves the problem Ax=b
     * Being A the input matrix, this method solves the problem Ax = b
     * @param rB The Right Hand Side (RHS) matrix
     * @param rX The solution matrix
     */
    virtual void Solve(
        MatrixType& rB,
        MatrixType& rX) const = 0;

    /**
     * @brief Solves the problem Ax=b
     * Being A the input matrix, this method solves the problem Ax = b
     * @param rB The Right Hand Side (RHS) vector
     * @param rX The solution vector
     */
    virtual void Solve(
        const VectorType& rB,
        VectorType& rX) const = 0;

    /**
     * @brief Unitary matrix getter
     * If computed, this method sets the unitary matrix in the provided array
     * @param rMatrixU Unitary matrix
     */
    virtual void MatrixQ(MatrixType& rMatrixQ) const = 0;

    /**
     * @brief Upper triangular matrix getter
     * If computed, this method sets the upper triangular matrix in the provided array
     * @param rMatrixV Upper triangular matrix
     */
    virtual void MatrixR(MatrixType& rMatrixR) const = 0;

    /**
     * @brief Pivoting matrix getter
     * If computed, this method sets the pivoting matrix
     * @param rMatrixP Pivoting matrix
     */
    virtual void MatrixP(MatrixType& rMatrixP) const = 0;

    /**
     * @brief Rank of the provided array
     * Calculates and returns the rank of the array decomposed with the QR
     * @return std::size_t Rank of the provided array
     */
    virtual std::size_t Rank() const = 0;

    /**
     * @brief QR information
     * Outputs the QR class information
     * @param rOStream Information output
     */
    virtual void PrintInfo(std::ostream &rOStream) const = 0;

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
};

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
} // namespace Kratos

#endif // defined(KRATOS_DENSE_QR_H_INCLUDED)
