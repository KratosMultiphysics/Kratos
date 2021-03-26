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

#if !defined(KRATOS_DENSE_SVD_H_INCLUDED)
#define KRATOS_DENSE_SVD_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

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
class DenseSingularValueDecomposition
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(DenseSingularValueDecomposition);

    typedef typename TDenseSpaceType::DataType DataType;
    typedef typename TDenseSpaceType::VectorType VectorType;
    typedef typename TDenseSpaceType::MatrixType MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    DenseSingularValueDecomposition() = default;

    virtual ~DenseSingularValueDecomposition() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Name of the SVD
     * Returns a string containing the name of the current SVD decomposition
     * @return std::string Name of the SVD decomposition
     */
    static std::string Name()
    {
        return "dense_singular_value_decomposition";
    }

    /**
     * @brief Compute the SVD
     * Computes the Singular Value Decomposition (SVD) of the given imput matrix
     * @param rInputMatrix Matrix to compute the SVD decomposition
     * @param Settings Settings for the SVD decomposition
     */
    virtual void Compute(
        MatrixType& rInputMatrix,
        Parameters Settings) = 0;

    /**
     * @brief Compute the SVD
     * Computes the Singular Value Decomposition (SVD) of the given input matrix
     * @param rInputMatrix Matrix to compute the SVD decomposition
     * @param rVectorS Vector containing the singular values (sorted from the largest one to smallest one)
     * @param rMatrixU Left singular vectors matrix
     * @param rMatrixV Right singular vectors matrix
     * @param Settings Settings for the SVD decomposition
     */
    virtual void Compute(
        MatrixType& rInputMatrix,
        VectorType& rVectorS,
        MatrixType& rMatrixU,
        MatrixType& rMatrixV,
        Parameters Settings) = 0;

    /**
     * @brief Left singular vectors matrix getter
     * If computed, this method sets the left singular vectors matrix in the provided array
     * @param rMatrixU Left singular vectors matrix
     */
    virtual void MatrixU(MatrixType& rMatrixU) = 0;

    /**
     * @brief Right singular vectors matrix getter
     * If computed, this method sets the right singular vectors matrix in the provided array
     * Note that this method is understood to return V (not its transpose). This means that
     * this matrix needs to be transposed in order to reconstruct the input matrix
     * @param rMatrixV
     */
    virtual void MatrixV(MatrixType& rMatrixV) = 0;

    /**
     * @brief Singular values vector getter
     * This method sets the singular values vector in the provided array
     * @param rVectorS
     */
    virtual void SingularValues(VectorType& rVectorS) = 0;

    /**
     * @brief Number of non-zero singular values
     * This method returns the number of non-zero singular values
     * @return std::size_t Number of non-zero singular values
     */
    virtual std::size_t NonZeroSingularValues() = 0;

    /**
     * @brief Rank of the provided array
     * Calculates and returns the rank of the array decomposed with the SVD
     * @return std::size_t Rank of the provided array
     */
    virtual std::size_t Rank() = 0;

    /**
     * @brief SVD information
     * Outputs the SVD class information
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

#endif // defined(KRATOS_DENSE_SVD_H_INCLUDED)
