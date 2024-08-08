//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//  Collaborators:   Vicente Mataix Ferrandiz
//                   Pablo Becker
//

#pragma once

// System includes
#include <cmath>
#include <type_traits>

// External includes

// Project includes
#include "input_output/logger.h"
#include "includes/ublas_interface.h"
#include "includes/global_variables.h"
#include "containers/array_1d.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

/**
 * @brief Enum class that defines the possible scenarios for the clamp function
 * @details This enum class defines the possible scenarios for the clamp function
 * WITHIN: The value is within the limits
 * MINIMUM: The value is below the minimum
 * MAXIMUM: The value is above the maximum
 */
enum class ClampScenario
{
    WITHIN_BOUNDS,
    BELOW_MINIMUM,
    ABOVE_MAXIMUM
};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MathUtils
 * @ingroup KratosCore
 * @brief Various mathematical utilitiy functions
 * @details Various mathematical utilitiy functions. Defines several utility functions.
 * @author Riccardo Rossi
 * @author Pooyan Dadvand
 */
template<class TDataType = double>
class KRATOS_API(KRATOS_CORE) MathUtils
{
public:

    ///@name Type Definitions
    ///@{

    /// The matrix type
    using MatrixType = Matrix;

    /// The vector type
    using VectorType = Vector;

    /// The size type
    using SizeType = std::size_t;

    /// The index type
    using IndexType = std::size_t;

    /// The indirect array type
    using IndirectArrayType = boost::numeric::ublas::indirect_array<DenseVector<std::size_t>>;

    /// The machine precision
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function returns the machine precision
     * @return The corresponding epsilon for the double
     */
    static inline double GetZeroTolerance()
    {
        return ZeroTolerance;
    }

    /**
     * @brief In geometry, Heron's formula (sometimes called Hero's formula), named after Hero of Alexandria, gives the area of a triangle by requiring no arbitrary choice of side as base or vertex as origin, contrary to other formulas for the area of a triangle, such as half the base times the height or half the norm of a cross product of two sides.
     * @param a First length
     * @param b Second length
     * @param c Third length
     * @return Heron solution: Heron's formula states that the area of a triangle whose sides have lengths a, b, and c
     */
    template<bool TCheck>// = false>
    static inline double Heron(
        double a,
        double b,
        double c
        )
    {
        const double s = 0.5 * (a + b + c);
        const double A2 = s * (s - a) * (s - b) * (s - c);
        if constexpr(TCheck) {
            if(A2 < 0.0) {
                KRATOS_ERROR << "The square of area is negative, probably the triangle is in bad shape:" << A2 << std::endl;
            } else {
                return std::sqrt(A2);
            }
        } else {
            return std::sqrt(std::abs(A2));
        }
    }

    /**
     * @brief Calculates the cofactor
     * @param rMat The matrix to calculate
     * @param i The index i
     * @param j The index j
     * @return The cofactor of the matrix
     */
    template<class TMatrixType>
    static double Cofactor(const TMatrixType& rMat, IndexType i, IndexType j)
    {
        static_assert(std::is_same<typename TMatrixType::value_type, double>::value, "Bad value type.");

        KRATOS_ERROR_IF(rMat.size1() != rMat.size2() || rMat.size1() == 0) << "Bad matrix dimensions." << std::endl;

        if (rMat.size1() == 1)
            return 1;

        IndirectArrayType ia1(rMat.size1() - 1), ia2(rMat.size2() - 1);

        // Construct the submatrix structure for the first minor.
        unsigned i_sub = 0;
        for (unsigned k = 0; k < rMat.size1(); ++k)
            if (k != i)
                ia1(i_sub++) = k;

        unsigned j_sub = 0;
        for (unsigned k = 0; k < rMat.size2(); ++k)
            if (k != j)
                ia2(j_sub++) = k;

        boost::numeric::ublas::matrix_indirect<const TMatrixType, IndirectArrayType> sub_mat(rMat, ia1, ia2);
        const double first_minor = Det(sub_mat);
        return ((i + j) % 2) ? -first_minor : first_minor;
    }

    /**
     * @brief Calculates the cofactor matrix
     * @param rMat The matrix to calculate
     * @return The cofactor matrix
     */
    template<class TMatrixType>
    static MatrixType CofactorMatrix(const TMatrixType& rMat)
    {
        static_assert(std::is_same<typename TMatrixType::value_type, double>::value, "Bad value type.");

        MatrixType cofactor_matrix(rMat.size1(), rMat.size2());

        for (IndexType i = 0; i < rMat.size1(); ++i)
            for (IndexType j = 0; j < rMat.size2(); ++j)
                cofactor_matrix(i, j) = Cofactor(rMat, i, j);

        return cofactor_matrix;
    }

    /**
     * @brief Calculates the inverse of a 2x2, 3x3 and 4x4 matrices (using bounded matrix for performance)
     * @param rInputMatrix The matrix to invert
     * @param rInputMatrixDet The determinant of the matrix
     * @param Tolerance The maximum tolerance considered
     * @return InvertMatrix: The inverted matrix
     */
    template<SizeType TDim>
    KRATOS_DEPRECATED_MESSAGE("Please use InvertMatrix() instead")
    static inline BoundedMatrix<double, TDim, TDim> InvertMatrix(
        const BoundedMatrix<double, TDim, TDim>& rInputMatrix,
        double& rInputMatrixDet,
        const double Tolerance = ZeroTolerance
        )
    {
        BoundedMatrix<double, TDim, TDim> inverted_matrix;

        /* Compute Determinant of the matrix */
        rInputMatrixDet = Det(rInputMatrix);

        if constexpr (TDim == 1) {
            inverted_matrix(0,0) = 1.0/rInputMatrix(0,0);
            rInputMatrixDet = rInputMatrix(0,0);
        } else if constexpr (TDim == 2) {
            InvertMatrix2(rInputMatrix, inverted_matrix, rInputMatrixDet);
        } else if constexpr (TDim == 3) {
            InvertMatrix3(rInputMatrix, inverted_matrix, rInputMatrixDet);
        } else if constexpr (TDim == 4) {
            InvertMatrix4(rInputMatrix, inverted_matrix, rInputMatrixDet);
        } else {
            KRATOS_ERROR << "Size not implemented. Size: " << TDim << std::endl;
        }

        // Checking condition number
        if (Tolerance > 0.0) { // Check is skipped for negative tolerances
            CheckConditionNumber(rInputMatrix, inverted_matrix, Tolerance);
        }

        return inverted_matrix;
    }

    /**
     * @brief This method checks the condition number of matrix
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param Tolerance The maximum tolerance considered
     */
    template<class TMatrix1, class TMatrix2>
    static inline bool CheckConditionNumber(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        const double Tolerance = ZeroTolerance,
        const bool ThrowError = true
        )
    {
        // We want at least 4 significant digits
        const double max_condition_number = (1.0/Tolerance) * 1.0e-4;

        // Find the condition number to define is inverse is OK
        const double input_matrix_norm = norm_frobenius(rInputMatrix);
        const double inverted_matrix_norm = norm_frobenius(rInvertedMatrix);

        // Now the condition number is the product of both norms
        const double cond_number = input_matrix_norm * inverted_matrix_norm ;
        // Finally check if the condition number is low enough
        if (cond_number > max_condition_number) {
            if (ThrowError) {
                KRATOS_WATCH(rInputMatrix);
                KRATOS_ERROR << " Condition number of the matrix is too high!, cond_number = " << cond_number << std::endl;
            }
            return false;
        }

        return true;
    }

    /**
     * @brief It inverts non square matrices (https://en.wikipedia.org/wiki/Inverse_element#Matrices)
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param rInputMatrixDet Is the determinant of the input matrix
     * @param Tolerance The maximum tolerance considered
     * @tparam TMatrix1 The type of the input matrix
     * @tparam TMatrix2 Is the type of the output matrix
     */
    template<class TMatrix1, class TMatrix2>
    static void GeneralizedInvertMatrix(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        double& rInputMatrixDet,
        const double Tolerance = ZeroTolerance
        )
    {
        const SizeType size_1 = rInputMatrix.size1();
        const SizeType size_2 = rInputMatrix.size2();

        if (size_1 == size_2) {
            InvertMatrix(rInputMatrix, rInvertedMatrix, rInputMatrixDet, Tolerance);
        } else if (size_1 < size_2) { // Right inverse
            if (rInvertedMatrix.size1() != size_2 || rInvertedMatrix.size2() != size_1) {
                rInvertedMatrix.resize(size_2, size_1, false);
            }
            const Matrix aux = prod(rInputMatrix, trans(rInputMatrix));
            Matrix auxInv;
            InvertMatrix(aux, auxInv, rInputMatrixDet, Tolerance);
            rInputMatrixDet = std::sqrt(rInputMatrixDet);
            noalias(rInvertedMatrix) = prod(trans(rInputMatrix), auxInv);
        } else { // Left inverse
            if (rInvertedMatrix.size1() != size_2 || rInvertedMatrix.size2() != size_1) {
                rInvertedMatrix.resize(size_2, size_1, false);
            }
            const Matrix aux = prod(trans(rInputMatrix), rInputMatrix);
            Matrix auxInv;
            InvertMatrix(aux, auxInv, rInputMatrixDet, Tolerance);
            rInputMatrixDet = std::sqrt(rInputMatrixDet);
            noalias(rInvertedMatrix) = prod(auxInv, trans(rInputMatrix));
        }
    }

    /**
     * @brief This function is designed to be called when a dense linear system is needed to be solved
     * @param A System matrix
     * @param rX Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB Right hand side vector.
     */
    static void Solve(
        MatrixType A,
        VectorType& rX,
        const VectorType& rB
        );

    /**
     * @brief It inverts matrices of order 2, 3 and 4
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param rInputMatrixDet Is the determinant of the input matrix
     * @tparam TMatrix1 The type of the input matrix
     * @tparam TMatrix2 Is the type of the output matrix
     */
    template<class TMatrix1, class TMatrix2>
    static void InvertMatrix(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        double& rInputMatrixDet,
        const double Tolerance = ZeroTolerance
        )
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rInputMatrix.size1() == rInputMatrix.size2()) << "Matrix provided is non-square" << std::endl;

        const SizeType size = rInputMatrix.size2();

        if(size == 1) {
            if(rInvertedMatrix.size1() != 1 || rInvertedMatrix.size2() != 1) {
                rInvertedMatrix.resize(1,1,false);
            }
            rInvertedMatrix(0,0) = 1.0/rInputMatrix(0,0);
            rInputMatrixDet = rInputMatrix(0,0);
        } else if (size == 2) {
            InvertMatrix2(rInputMatrix, rInvertedMatrix, rInputMatrixDet);
        } else if (size == 3) {
            InvertMatrix3(rInputMatrix, rInvertedMatrix, rInputMatrixDet);
        } else if (size == 4) {
            InvertMatrix4(rInputMatrix, rInvertedMatrix, rInputMatrixDet);
        } else if (std::is_same<TMatrix1, Matrix>::value) {

            const SizeType size1 = rInputMatrix.size1();
            const SizeType size2 = rInputMatrix.size2();
            if(rInvertedMatrix.size1() != size1 || rInvertedMatrix.size2() != size2) {
                rInvertedMatrix.resize(size1, size2,false);
            }

            Matrix A(rInputMatrix);
            typedef permutation_matrix<SizeType> pmatrix;
            pmatrix pm(A.size1());
            const int singular = lu_factorize(A,pm);
            rInvertedMatrix.assign( IdentityMatrix(A.size1()));
            KRATOS_ERROR_IF(singular == 1) << "Matrix is singular: " << rInputMatrix << std::endl;
            lu_substitute(A, pm, rInvertedMatrix);

            // Calculating determinant
            rInputMatrixDet = 1.0;

            for (IndexType i = 0; i < size1;++i) {
                IndexType ki = pm[i] == i ? 0 : 1;
                rInputMatrixDet *= (ki == 0) ? A(i,i) : -A(i,i);
            }
       } else { // Bounded-matrix case
            const SizeType size1 = rInputMatrix.size1();
            const SizeType size2 = rInputMatrix.size2();

            Matrix A(rInputMatrix);
            Matrix invA(rInvertedMatrix);

            typedef permutation_matrix<SizeType> pmatrix;
            pmatrix pm(size1);
            const int singular = lu_factorize(A,pm);
            invA.assign( IdentityMatrix(size1));
            KRATOS_ERROR_IF(singular == 1) << "Matrix is singular: " << rInputMatrix << std::endl;
            lu_substitute(A, pm, invA);

            // Calculating determinant
            rInputMatrixDet = 1.0;

            for (IndexType i = 0; i < size1;++i) {
                IndexType ki = pm[i] == i ? 0 : 1;
                rInputMatrixDet *= (ki == 0) ? A(i,i) : -A(i,i);
            }

            for (IndexType i = 0; i < size1;++i) {
                for (IndexType j = 0; j < size2;++j) {
                    rInvertedMatrix(i,j) = invA(i,j);
                }
            }
       }

       // Checking condition number
       if (Tolerance > 0.0) { // Check is skipped for negative tolerances
            CheckConditionNumber(rInputMatrix, rInvertedMatrix, Tolerance);
       }
    }

    /**
     * @brief It inverts matrices of order 2
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param rInputMatrixDet Is the determinant of the input matrix
     */
    template<class TMatrix1, class TMatrix2>
    static void InvertMatrix2(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        double& rInputMatrixDet
        )
    {
        KRATOS_TRY;

        KRATOS_DEBUG_ERROR_IF_NOT(rInputMatrix.size1() == rInputMatrix.size2()) << "Matrix provided is non-square" << std::endl;

        if(rInvertedMatrix.size1() != 2 || rInvertedMatrix.size2() != 2) {
            rInvertedMatrix.resize(2,2,false);
        }

        rInputMatrixDet = rInputMatrix(0,0)*rInputMatrix(1,1)-rInputMatrix(0,1)*rInputMatrix(1,0);

        rInvertedMatrix(0,0) =  rInputMatrix(1,1);
        rInvertedMatrix(0,1) = -rInputMatrix(0,1);
        rInvertedMatrix(1,0) = -rInputMatrix(1,0);
        rInvertedMatrix(1,1) =  rInputMatrix(0,0);

        rInvertedMatrix/=rInputMatrixDet;

        KRATOS_CATCH("");
    }

    /**
     * @brief It inverts matrices of order 3
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param rInputMatrixDet Is the determinant of the input matrix
     */
    template<class TMatrix1, class TMatrix2>
    static void InvertMatrix3(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        double& rInputMatrixDet
        )
    {
        KRATOS_TRY;

        KRATOS_DEBUG_ERROR_IF_NOT(rInputMatrix.size1() == rInputMatrix.size2()) << "Matrix provided is non-square" << std::endl;

        if(rInvertedMatrix.size1() != 3 || rInvertedMatrix.size2() != 3) {
            rInvertedMatrix.resize(3,3,false);
        }

        // Filling the inverted matrix with the algebraic complements
        // First column
        rInvertedMatrix(0,0) = rInputMatrix(1,1)*rInputMatrix(2,2) - rInputMatrix(1,2)*rInputMatrix(2,1);
        rInvertedMatrix(1,0) = -rInputMatrix(1,0)*rInputMatrix(2,2) + rInputMatrix(1,2)*rInputMatrix(2,0);
        rInvertedMatrix(2,0) = rInputMatrix(1,0)*rInputMatrix(2,1) - rInputMatrix(1,1)*rInputMatrix(2,0);

        // Second column
        rInvertedMatrix(0,1) = -rInputMatrix(0,1)*rInputMatrix(2,2) + rInputMatrix(0,2)*rInputMatrix(2,1);
        rInvertedMatrix(1,1) = rInputMatrix(0,0)*rInputMatrix(2,2) - rInputMatrix(0,2)*rInputMatrix(2,0);
        rInvertedMatrix(2,1) = -rInputMatrix(0,0)*rInputMatrix(2,1) + rInputMatrix(0,1)*rInputMatrix(2,0);

        // Third column
        rInvertedMatrix(0,2) = rInputMatrix(0,1)*rInputMatrix(1,2) - rInputMatrix(0,2)*rInputMatrix(1,1);
        rInvertedMatrix(1,2) = -rInputMatrix(0,0)*rInputMatrix(1,2) + rInputMatrix(0,2)*rInputMatrix(1,0);
        rInvertedMatrix(2,2) = rInputMatrix(0,0)*rInputMatrix(1,1) - rInputMatrix(0,1)*rInputMatrix(1,0);

        // Calculation of determinant (of the input matrix)
        rInputMatrixDet = rInputMatrix(0,0)*rInvertedMatrix(0,0) + rInputMatrix(0,1)*rInvertedMatrix(1,0) + rInputMatrix(0,2)*rInvertedMatrix(2,0);

        // Finalizing the calculation of the inverted matrix
        rInvertedMatrix /= rInputMatrixDet;

        KRATOS_CATCH("")
    }

    /**
     * @brief It inverts matrices of order 4
     * @param rInputMatrix Is the input matrix (unchanged at output)
     * @param rInvertedMatrix Is the inverse of the input matrix
     * @param rInputMatrixDet Is the determinant of the input matrix
     */
    template<class TMatrix1, class TMatrix2>
    static void InvertMatrix4(
        const TMatrix1& rInputMatrix,
        TMatrix2& rInvertedMatrix,
        double& rInputMatrixDet
        )
    {
        KRATOS_TRY;

        KRATOS_DEBUG_ERROR_IF_NOT(rInputMatrix.size1() == rInputMatrix.size2()) << "Matrix provided is non-square" << std::endl;

        if (rInvertedMatrix.size1() != 4 || rInvertedMatrix.size2() != 4) {
            rInvertedMatrix.resize(4, 4, false);
        }

        /* Compute inverse of the Matrix */
        // First column
        rInvertedMatrix(0, 0) = -(rInputMatrix(1, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 1)) + rInputMatrix(1, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 1) + rInputMatrix(1, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 2) - rInputMatrix(1, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 2) - rInputMatrix(1, 2) * rInputMatrix(2, 1) * rInputMatrix(3, 3) + rInputMatrix(1, 1) * rInputMatrix(2, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(0, 1) = rInputMatrix(0, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 1) - rInputMatrix(0, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 1) - rInputMatrix(0, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 2) + rInputMatrix(0, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 2) + rInputMatrix(0, 2) * rInputMatrix(2, 1) * rInputMatrix(3, 3) - rInputMatrix(0, 1) * rInputMatrix(2, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(0, 2) = -(rInputMatrix(0, 3) * rInputMatrix(1, 2) * rInputMatrix(3, 1)) + rInputMatrix(0, 2) * rInputMatrix(1, 3) * rInputMatrix(3, 1) + rInputMatrix(0, 3) * rInputMatrix(1, 1) * rInputMatrix(3, 2) - rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(3, 2) - rInputMatrix(0, 2) * rInputMatrix(1, 1) * rInputMatrix(3, 3) + rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(0, 3) = rInputMatrix(0, 3) * rInputMatrix(1, 2) * rInputMatrix(2, 1) - rInputMatrix(0, 2) * rInputMatrix(1, 3) * rInputMatrix(2, 1) - rInputMatrix(0, 3) * rInputMatrix(1, 1) * rInputMatrix(2, 2) + rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(2, 2) + rInputMatrix(0, 2) * rInputMatrix(1, 1) * rInputMatrix(2, 3) - rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(2, 3);

        // Second column
        rInvertedMatrix(1, 0) = rInputMatrix(1, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 0) - rInputMatrix(1, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 0) - rInputMatrix(1, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 2) + rInputMatrix(1, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 2) + rInputMatrix(1, 2) * rInputMatrix(2, 0) * rInputMatrix(3, 3) - rInputMatrix(1, 0) * rInputMatrix(2, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(1, 1) = -(rInputMatrix(0, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 0)) + rInputMatrix(0, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 0) + rInputMatrix(0, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 2) - rInputMatrix(0, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 2) - rInputMatrix(0, 2) * rInputMatrix(2, 0) * rInputMatrix(3, 3) + rInputMatrix(0, 0) * rInputMatrix(2, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(1, 2) = rInputMatrix(0, 3) * rInputMatrix(1, 2) * rInputMatrix(3, 0) - rInputMatrix(0, 2) * rInputMatrix(1, 3) * rInputMatrix(3, 0) - rInputMatrix(0, 3) * rInputMatrix(1, 0) * rInputMatrix(3, 2) + rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(3, 2) + rInputMatrix(0, 2) * rInputMatrix(1, 0) * rInputMatrix(3, 3) - rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(3, 3);
        rInvertedMatrix(1, 3) = -(rInputMatrix(0, 3) * rInputMatrix(1, 2) * rInputMatrix(2, 0)) + rInputMatrix(0, 2) * rInputMatrix(1, 3) * rInputMatrix(2, 0) + rInputMatrix(0, 3) * rInputMatrix(1, 0) * rInputMatrix(2, 2) - rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(2, 2) - rInputMatrix(0, 2) * rInputMatrix(1, 0) * rInputMatrix(2, 3) + rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(2, 3);

        // Third column
        rInvertedMatrix(2, 0) = -(rInputMatrix(1, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 0)) + rInputMatrix(1, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 0) + rInputMatrix(1, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 1) - rInputMatrix(1, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 1) - rInputMatrix(1, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 3) + rInputMatrix(1, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 3);
        rInvertedMatrix(2, 1) = rInputMatrix(0, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 0) - rInputMatrix(0, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 0) - rInputMatrix(0, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 1) + rInputMatrix(0, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 1) + rInputMatrix(0, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 3) - rInputMatrix(0, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 3);
        rInvertedMatrix(2, 2) = -(rInputMatrix(0, 3) * rInputMatrix(1, 1) * rInputMatrix(3, 0)) + rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(3, 0) + rInputMatrix(0, 3) * rInputMatrix(1, 0) * rInputMatrix(3, 1) - rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(3, 1) - rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(3, 3) + rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(3, 3);
        rInvertedMatrix(2, 3) = rInputMatrix(0, 3) * rInputMatrix(1, 1) * rInputMatrix(2, 0) - rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(2, 0) - rInputMatrix(0, 3) * rInputMatrix(1, 0) * rInputMatrix(2, 1) + rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(2, 1) + rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(2, 3) - rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(2, 3);

        // Fourth column
        rInvertedMatrix(3, 0) = rInputMatrix(1, 2) * rInputMatrix(2, 1) * rInputMatrix(3, 0) - rInputMatrix(1, 1) * rInputMatrix(2, 2) * rInputMatrix(3, 0) - rInputMatrix(1, 2) * rInputMatrix(2, 0) * rInputMatrix(3, 1) + rInputMatrix(1, 0) * rInputMatrix(2, 2) * rInputMatrix(3, 1) + rInputMatrix(1, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 2) - rInputMatrix(1, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 2);
        rInvertedMatrix(3, 1) = -(rInputMatrix(0, 2) * rInputMatrix(2, 1) * rInputMatrix(3, 0)) + rInputMatrix(0, 1) * rInputMatrix(2, 2) * rInputMatrix(3, 0) + rInputMatrix(0, 2) * rInputMatrix(2, 0) * rInputMatrix(3, 1) - rInputMatrix(0, 0) * rInputMatrix(2, 2) * rInputMatrix(3, 1) - rInputMatrix(0, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 2) + rInputMatrix(0, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 2);
        rInvertedMatrix(3, 2) = rInputMatrix(0, 2) * rInputMatrix(1, 1) * rInputMatrix(3, 0) - rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(3, 0) - rInputMatrix(0, 2) * rInputMatrix(1, 0) * rInputMatrix(3, 1) + rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(3, 1) + rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(3, 2) - rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(3, 2);
        rInvertedMatrix(3, 3) = -(rInputMatrix(0, 2) * rInputMatrix(1, 1) * rInputMatrix(2, 0)) + rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(2, 0) + rInputMatrix(0, 2) * rInputMatrix(1, 0) * rInputMatrix(2, 1) - rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(2, 1) - rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(2, 2) + rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(2, 2);

        // Calculation of determinant (of the input matrix)
        rInputMatrixDet = rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 0) - rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 0) - rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(2, 2) * rInputMatrix(3, 1) + rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(2, 3) * rInputMatrix(3, 1) - rInputMatrix(0, 1) * rInputMatrix(1, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 2) + rInputMatrix(0, 0) * rInputMatrix(1, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 2) + rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 2) - rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 2) + rInputMatrix(0, 3) * (rInputMatrix(1, 2) * rInputMatrix(2, 1) * rInputMatrix(3, 0) - rInputMatrix(1, 1) * rInputMatrix(2, 2) * rInputMatrix(3, 0) - rInputMatrix(1, 2) * rInputMatrix(2, 0) * rInputMatrix(3, 1) + rInputMatrix(1, 0) * rInputMatrix(2, 2) * rInputMatrix(3, 1) + rInputMatrix(1, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 2) - rInputMatrix(1, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 2)) + (rInputMatrix(0, 1) * rInputMatrix(1, 2) * rInputMatrix(2, 0) - rInputMatrix(0, 0) * rInputMatrix(1, 2) * rInputMatrix(2, 1) - rInputMatrix(0, 1) * rInputMatrix(1, 0) * rInputMatrix(2, 2) + rInputMatrix(0, 0) * rInputMatrix(1, 1) * rInputMatrix(2, 2)) * rInputMatrix(3, 3) + rInputMatrix(0, 2) * (-(rInputMatrix(1, 3) * rInputMatrix(2, 1) * rInputMatrix(3, 0)) + rInputMatrix(1, 1) * rInputMatrix(2, 3) * rInputMatrix(3, 0) + rInputMatrix(1, 3) * rInputMatrix(2, 0) * rInputMatrix(3, 1) - rInputMatrix(1, 0) * rInputMatrix(2, 3) * rInputMatrix(3, 1) - rInputMatrix(1, 1) * rInputMatrix(2, 0) * rInputMatrix(3, 3) + rInputMatrix(1, 0) * rInputMatrix(2, 1) * rInputMatrix(3, 3));

        // Finalizing the calculation of the inverted matrix
        rInvertedMatrix /= rInputMatrixDet;

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates the determinant of a matrix of dimension 2x2 (no check performed on matrix size)
     * @param rA Is the input matrix
     * @return The determinant of the 2x2 matrix
     */
    template<class TMatrixType>
    static inline double Det2(const TMatrixType& rA)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rA.size1() == rA.size2()) << "Matrix provided is non-square" << std::endl;

        return (rA(0,0)*rA(1,1)-rA(0,1)*rA(1,0));
    }

    /**
     * @brief Calculates the determinant of a matrix of dimension 3*3 (no check performed on matrix size)
     * @param rA Is the input matrix
     * @return The determinant of the 3x3 matrix
     */
    template<class TMatrixType>
    static inline double Det3(const TMatrixType& rA)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rA.size1() == rA.size2()) << "Matrix provided is non-square" << std::endl;

        // Calculating the algebraic complements to the first line
        const double a = rA(1,1)*rA(2,2) - rA(1,2)*rA(2,1);
        const double b = rA(1,0)*rA(2,2) - rA(1,2)*rA(2,0);
        const double c = rA(1,0)*rA(2,1) - rA(1,1)*rA(2,0);

        return rA(0,0)*a - rA(0,1)*b + rA(0,2)*c;
    }

    /**
     * @brief Calculates the determinant of a matrix of dimension 4*4 (no check performed on matrix size)
     * @param rA Is the input matrix
     * @return The determinant of the 4x4 matrix
     */
    template<class TMatrixType>
    static inline double Det4(const TMatrixType& rA)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rA.size1() == rA.size2()) << "Matrix provided is non-square" << std::endl;

        const double det = rA(0,1)*rA(1,3)*rA(2,2)*rA(3,0)-rA(0,1)*rA(1,2)*rA(2,3)*rA(3,0)-rA(0,0)*rA(1,3)*rA(2,2)*rA(3,1)+rA(0,0)*rA(1,2)*rA(2,3)*rA(3,1)
                          -rA(0,1)*rA(1,3)*rA(2,0)*rA(3,2)+rA(0,0)*rA(1,3)*rA(2,1)*rA(3,2)+rA(0,1)*rA(1,0)*rA(2,3)*rA(3,2)-rA(0,0)*rA(1,1)*rA(2,3)*rA(3,2)+rA(0,3)*(rA(1,2)*rA(2,1)*rA(3,0)-rA(1,1)*rA(2,2)*rA(3,0)-rA(1,2)*rA(2,0)*rA(3,1)+rA(1,0)*rA(2,2)*rA(3,1)+rA(1,1)*rA(2,0)*rA(3,2)
                          -rA(1,0)*rA(2,1)*rA(3,2))+(rA(0,1)*rA(1,2)*rA(2,0)-rA(0,0)*rA(1,2)*rA(2,1)-rA(0,1)*rA(1,0)*rA(2,2)+rA(0,0)*rA(1,1)*rA(2,2))*rA(3,3)+rA(0,2)*(-(rA(1,3)*rA(2,1)*rA(3,0))+rA(1,1)*rA(2,3)*rA(3,0)+rA(1,3)*rA(2,0)*rA(3,1)-rA(1,0)*rA(2,3)*rA(3,1)-rA(1,1)*rA(2,0)*rA(3,3)+rA(1,0)*rA(2,1)*rA(3,3));
        return det;
    }

public:
    /**
     * @brief Calculates the determinant of a matrix of a square matrix of any size (no check performed on release mode)
     * @param rA Is the input matrix
     * @return The determinant of any size matrix
     */
    template<class TMatrixType>
    static inline double Det(const TMatrixType& rA)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rA.size1() == rA.size2()) << "Matrix provided is non-square" << std::endl;

        switch (rA.size1()) {
            case 2:
                return Det2(rA);
            case 3:
                return Det3(rA);
            case 4:
                return Det4(rA);
            default:
                double det = 1.0;
                using namespace boost::numeric::ublas;
                typedef permutation_matrix<SizeType> pmatrix;
                Matrix Aux(rA);
                pmatrix pm(Aux.size1());
                bool singular = lu_factorize(Aux,pm);

                if (singular) {
                    return 0.0;
                }

                for (IndexType i = 0; i < Aux.size1();++i) {
                    IndexType ki = pm[i] == i ? 0 : 1;
                    det *= std::pow(-1.0, ki) * Aux(i,i);
                }
                return det;
       }
    }

    /**
     * @brief Calculates the determinant of any matrix (no check performed on matrix size)
     * @param rA Is the input matrix
     * @return The determinant of the 2x2 matrix
     */
    template<class TMatrixType>
    static inline double GeneralizedDet(const TMatrixType& rA)
    {
        if (rA.size1() == rA.size2()) {
            return Det(rA);
        } else if (rA.size1() < rA.size2()) { // Right determinant
            const Matrix AAT = prod( rA, trans(rA) );
            return std::sqrt(Det(AAT));
        } else { // Left determinant
            const Matrix ATA = prod( trans(rA), rA );
            return std::sqrt(Det(ATA));
        }
    }

    /**
     * @brief Performs the dot product of two vectors of dimension 3
     * @details No check performed on vector sizes
     * @param a First input vector
     * @param b Second input vector
     * @return The resulting norm
     */
    static inline double Dot3(
        const Vector& a,
        const Vector& b
        )
    {
        KRATOS_DEBUG_ERROR_IF_NOT(a.size() == 3) << "The size of the first vector is not 3. Size: " << a.size() << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(b.size() == 3) << "The size of the second vector is not 3. Size: " << b.size() << std::endl;
        return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }

    /**
     * @brief Computes the dot product of two 1D arrays of doubles with size 3.
     * @param rA The first array.
     * @param rB The second array.
     * @return double The dot product of the two arrays.
     */
    static inline double Dot3(
        const array_1d<double, 3>& rA,
        const array_1d<double, 3>& rB
        )
    {
        return (rA[0]*rB[0] + rA[1]*rB[1] + rA[2]*rB[2]);
    }

    /**
     * @brief Performs the dot product of two vectors of arbitrary size
     * @details No check performed on vector sizes
     * @param rFirstVector First input vector
     * @param rSecondVector Second input vector
     * @return The resulting norm
     */
    static inline double Dot(
        const Vector& rFirstVector,
        const Vector& rSecondVector
        )
    {
        Vector::const_iterator i = rFirstVector.begin();
        Vector::const_iterator j = rSecondVector.begin();
        double temp = 0.0;
        while(i != rFirstVector.end()) {
            temp += *i++ * *j++;
        }
        return temp;
        //return std::inner_product(rFirstVector.begin(), rFirstVector.end(), rSecondVector.begin(), 0.0);
    }

    /**
    * @brief Clamps a value between a minimum and maximum range.
    * @details This function ensures that the given value `rX` lies within the specified
    * range `[rMinimum, rMaximum]`. If `rX` is less than `rMinimum`, it returns
    * `rMinimum`. If `rX` is greater than `rMaximum`, it returns `rMaximum`.
    * Otherwise, it returns `rX`.
    * @tparam T The type of the value and bounds, typically a numeric type.
    * @param rX The value to clamp.
    * @param rMinimum The minimum bound.
    * @param rMaximum The maximum bound.
    * @param Tolerance The tolerance value to consider for clamping, defaults to 0.
    * @return The clamped value.
    */
    template <typename T>
    static T Clamp(
        const T& rX,
        const T& rMinimum,
        const T& rMaximum,
        const T Tolerance = 0.0
        )
    {
        if (rX < rMinimum - Tolerance) {
            return rMinimum;
        } else if (rX > rMaximum + Tolerance) {
            return rMaximum;
        } else {
            return rX;
        }
    }

    /**
    * @brief Clamps a value between a minimum and maximum range with a scenario.
    * @details This function ensures that the given value `rX` lies within the specified
    * range `[rMinimum, rMaximum]`. If `rX` is less than `rMinimum`, it returns
    * `rMinimum`. If `rX` is greater than `rMaximum`, it returns `rMaximum`.
    * Otherwise, it returns `rX`. It also returns a scenario indicating if the value is below the minimum, above the maximum or within the bounds.
    * @tparam T The type of the value and bounds, typically a numeric type.
    * @param rX The value to clamp.
    * @param rMinimum The minimum bound.
    * @param rMaximum The maximum bound.
    * @param rScenario The scenario where the value is clamped. It can be BELOW_MINIMUM, ABOVE_MAXIMUM or WITHIN_BOUNDS.
    * @param Tolerance The tolerance value to consider for clamping, defaults to 0.
    * @return The clamped value.
    */
    template <typename T>
    static T Clamp(
        const T& rX,
        const T& rMinimum,
        const T& rMaximum,
        ClampScenario& rScenario,
        const T Tolerance = 0.0
        )
    {
        if (rX < rMinimum - Tolerance) {
            rScenario = ClampScenario::BELOW_MINIMUM;
            return rMinimum;
        } else if (rX > rMaximum + Tolerance) {
            rScenario = ClampScenario::ABOVE_MAXIMUM;
            return rMaximum;
        } else {
            rScenario = ClampScenario::WITHIN_BOUNDS;
            return rX;
        }
    }

    /**
     * @brief Calculates the norm of vector "a" which is assumed to be of size 3
     * @details No check is performed on the vector's size
     * @param a Input vector
     * @return The resulting norm
     */
    template<class TVectorType>
    static inline double Norm3(const TVectorType& a)
    {
        double temp = std::pow(a[0],2) + std::pow(a[1],2) + std::pow(a[2],2);
        temp = std::sqrt(temp);
        return temp;
    }

    /**
     * @brief Calculates the norm of vector "a"
     * @param a Input vector
     * @return The resulting norm
     */
    static inline double Norm(const Vector& a)
    {
        Vector::const_iterator i = a.begin();
        double temp = 0.0;
        while(i != a.end()) {
            temp += (*i) * (*i);
            i++;
        }
        return std::sqrt(temp);
    }

    /**
     * @brief Calculates the norm of vector "a" while avoiding underflow and overflow.
     * @param a Input vector
     * @return The resulting norm
     * @see http://www.netlib.org/lapack/explore-html/da/d7f/dnrm2_8f_source.html
     */
    static inline double StableNorm(const Vector& a)
    {
        if (a.size() == 0) {
            return 0;
        }

        if (a.size() == 1) {
            return a[0];
        }

        double scale {0};

        double sqr_sum_scaled {1};

        for (auto it = a.begin(); it != a.end(); ++it) {
            double x = *it;

            if (x != 0) {
                const double abs_x = std::abs(x);

                if (scale < abs_x) {
                    const double f = scale / abs_x;
                    sqr_sum_scaled = sqr_sum_scaled * (f * f) + 1.0;
                    scale = abs_x;
                } else {
                    x = abs_x / scale;
                    sqr_sum_scaled += x * x;
                }
            }
        }

        return scale * std::sqrt(sqr_sum_scaled);
    }

    /**
     * @brief Performs the vector product of the two input vectors a,b
     * @details a,b are assumed to be of size 3 (no check is performed on vector sizes)
     * @param a First input vector
     * @param b Second input vector
     * @return The resulting vector
     */
    template<class T>
    static inline T CrossProduct(
        const T& a,
        const T& b
        )
    {
        T c(a);

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        return c;
    }

    /**
    * @brief Checks there is aliasing
    * @param value1 The first value
    * @param value2 The second value
    */
    template< class T1, class T2>
    static inline typename std::enable_if<std::is_same<T1, T2>::value, bool>::type CheckIsAlias(T1& value1, T2& value2)
    {
        return value1 == value2;
    }

    /**
    * @brief Checks there is aliasing
    * @param value1 The first value
    * @param value2 The second value
    */
    template< class T1, class T2>
    static inline typename std::enable_if<!std::is_same<T1, T2>::value, bool>::type CheckIsAlias(T1& value1, T2& value2)
    {
        return false;
    }

    /**
     * @brief Performs the cross product of the two input vectors a,b
     * @details a,b are assumed to be of size 3 (check is only performed on vector sizes in debug mode)
     * @param a First input vector
     * @param b Second input vector
     * @param c The resulting vector
     */
    template< class T1, class T2 , class T3>
    static inline void CrossProduct(T1& c, const T2& a, const T3& b ){
        if (c.size() != 3) c.resize(3);

        KRATOS_DEBUG_ERROR_IF(a.size() != 3 || b.size() != 3 || c.size() != 3)
            << "The size of the vectors is different of 3: "
            << a << ", " << b << " and " << c << std::endl;
        KRATOS_DEBUG_ERROR_IF(CheckIsAlias(c, a))
            << "Aliasing between the output parameter and the first "
            << "input parameter" << std::endl;
        KRATOS_DEBUG_ERROR_IF(CheckIsAlias(c, b))  << "Aliasing between "
            << "the output parameter and the second input parameter"  << std::endl;

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    /**
     * @brief Performs the unitary cross product of the two input vectors a,b
     * @details a,b are assumed to be of size 3 (no check is performed on vector sizes)
     * @param a First input vector
     * @param b Second input vector
     * @param c The resulting vector
     */
    template< class T1, class T2 , class T3>
    static inline void UnitCrossProduct(T1& c, const T2& a, const T3& b ){
        CrossProduct(c,a,b);
        const double norm = norm_2(c);
        KRATOS_DEBUG_ERROR_IF(norm < 1000.0*ZeroTolerance)
            << "norm is 0 when making the UnitCrossProduct of the vectors "
            << a << " and " << b << std::endl;
        c/=norm;
    }

    /**
     * @brief This computes a orthonormal basis from a given vector (Frisvad method)
     * @param c The input vector
     * @param a First resulting vector
     * @param b Second resulting vector
     * @param Type The type of method employed, 0 is HughesMoeller, 1 is Frisvad and otherwise Naive
     */
    template< class T1, class T2 , class T3>
    static inline void OrthonormalBasis(const T1& c,T2& a,T3& b, const IndexType Type = 0 ){
        if (Type == 0)
            OrthonormalBasisHughesMoeller(c,a,b);
        else if (Type == 1)
            OrthonormalBasisFrisvad(c,a,b);
        else
            OrthonormalBasisNaive(c,a,b);
    }

    /**
     * @brief This computes a orthonormal basis from a given vector (Hughes Moeller method)
     * @param c The input vector
     * @param a First resulting vector
     * @param b Second resulting vector
     * @note Orthonormal basis taken from: http://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
     */
    template< class T1, class T2 , class T3>
    static inline void OrthonormalBasisHughesMoeller(const T1& c,T2& a,T3& b ){
        KRATOS_DEBUG_ERROR_IF(norm_2(c) < (1.0 - 1.0e-6) || norm_2(c) > (1.0 + 1.0e-6)) << "Input should be a normal vector" << std::endl;
        //  Choose a vector  orthogonal  to n as the  direction  of b2.
        if(std::abs(c[0]) > std::abs(c[2])) {
            b[0] =  c[1];
            b[1] = -c[0];
            b[2] =  0.0;
        } else {
            b[0] =   0.0;
            b[1] =   c[2];
            b[2]  = -c[1];
        }
        b /=  norm_2(b); //  Normalize  b
        UnitCrossProduct(a, b , c); //  Construct  a  using a cross  product
    }

    /**
     * @brief This computes a orthonormal basis from a given vector (Frisvad method)
     * @param c The input vector
     * @param a First resulting vector
     * @param b Second resulting vector
     * @note Orthonormal basis taken from: http://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
     */
    template< class T1, class T2 , class T3>
    static inline void OrthonormalBasisFrisvad(const T1& c,T2& a,T3& b ){
        KRATOS_DEBUG_ERROR_IF(norm_2(c) < (1.0 - 1.0e-3) || norm_2(c) > (1.0 + 1.0e-3)) << "Input should be a normal vector" << std::endl;
        if ((c[2] + 1.0) > 1.0e4 * ZeroTolerance) {
            a[0] = 1.0 - std::pow(c[0], 2)/(1.0 + c[2]);
            a[1] = - (c[0] * c[1])/(1.0 + c[2]);
            a[2] = - c[0];
            const double norm_a = norm_2(a);
            a /= norm_a;
            b[0] = - (c[0] * c[1])/(1.0 + c[2]);
            b[1] = 1.0 - std::pow(c[1], 2)/(1.0 + c[2]);
            b[2] = -c[1];
            const double norm_b = norm_2(b);
            b /= norm_b;
        } else { // In case that the vector is in negative Z direction
            a[0] = 1.0;
            a[1] = 0.0;
            a[2] = 0.0;
            b[0] = 0.0;
            b[1] = -1.0;
            b[2] = 0.0;
        }
    }

    /**
     * @brief This computes a orthonormal basis from a given vector (Naive method)
     * @param c The input vector
     * @param a First resulting vector
     * @param b Second resulting vector
     * @note Orthonormal basis taken from: http://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
     */
    template< class T1, class T2 , class T3>
    static inline void OrthonormalBasisNaive(const T1& c,T2& a,T3& b ){
        KRATOS_DEBUG_ERROR_IF(norm_2(c) < (1.0 - 1.0e-3) || norm_2(c) > (1.0 + 1.0e-3)) << "Input should be a normal vector" << std::endl;
        // If c is near  the x-axis , use  the y-axis. Otherwise  use  the x-axis.
        if(c[0] > 0.9f) {
            a[0] = 0.0;
            a[1] = 1.0;
            a[2] = 0.0;
        } else {
            a[0] = 1.0;
            a[1] = 0.0;
            a[2] = 0.0;
        }
        a  -= c * inner_prod(a, c); // Make a  orthogonal  to c
        a /=  norm_2(a);            //  Normalize  a
        UnitCrossProduct(b, c, a);  //  Construct  b  using a cross  product
    }

    /**
     * @brief Computes the angle between two vectors in 3D
     * @param rV1 First input vector
     * @param rV2 Second input vector
     */
    template< class T1, class T2>
    static inline double VectorsAngle(const T1& rV1, const T2& rV2 ){
        const T1 aux_1 = rV1 * norm_2(rV2);
        const T2 aux_2 = norm_2(rV1) * rV2;
        const double num = norm_2(aux_1 - aux_2);
        const double denom = norm_2(aux_1 + aux_2);
        return 2.0 * std::atan2( num , denom);
    }

    /**
     * @brief Returns a matrix :
     * A = a.tensorproduct.b
     * @details a,b are assumed to be of order 3, no check is performed on the size of the vectors
     * @param a First input vector
     * @param b Second input vector
     * @return Returns A = a.tensorproduct.b
     */
    static inline MatrixType TensorProduct3(
        const Vector& a,
        const Vector& b
        )
    {
        MatrixType A(3,3);

        A(0,0) = a[0]*b[0];
        A(0,1) = a[0]*b[1];
        A(0,2) = a[0]*b[2];
        A(1,0) = a[1]*b[0];
        A(1,1) = a[1]*b[1];
        A(1,2) = a[1]*b[2];
        A(2,0) = a[2]*b[0];
        A(2,1) = a[2]*b[1];
        A(2,2) = a[2]*b[2];

        return A;
    }

    /**
     * @brief "rInputMatrix" is ADDED to "Destination" matrix starting from InitialRow and InitialCol of the destination matrix
     * @details "Destination" is assumed to be able to contain the "input matrix" (no check is performed on the bounds)
     * @param rDestination The matrix destination
     * @param rInputMatrix The input matrix to be added
     * @param InitialRow The initial row
     * @param InitialCol The initial column
     */
    template<class TMatrixType1, class TMatrixType2>
    static inline void AddMatrix(
        TMatrixType1& rDestination,
        const TMatrixType2& rInputMatrix,
        const IndexType InitialRow,
        const IndexType InitialCol
        )
    {
        KRATOS_TRY

        for(IndexType i = 0; i < rInputMatrix.size1(); ++i) {
            for(IndexType j = 0; j < rInputMatrix.size2(); ++j) {
                rDestination(InitialRow+i, InitialCol+j) += rInputMatrix(i,j);
            }
        }
        KRATOS_CATCH("")
    }

    /**
     * @brief "rInputVector" is ADDED to "Destination" vector starting from InitialIndex of the destination matrix
     * @details "Destination" is assumed to be able to contain the "input vector" (no check is performed on the bounds)
     * @param rDestination The vector destination
     * @param rInputVector The input vector to be added
     * @param InitialIndex The initial index
     */
    template<class TVectorType1, class TVectorType2>
    static inline void AddVector(
        TVectorType1& rDestination,
        const TVectorType2& rInputVector,
        const IndexType InitialIndex
        )
    {
        KRATOS_TRY

        for(IndexType i = 0; i < rInputVector.size(); ++i) {
            rDestination[InitialIndex+i] += rInputVector[i];
        }
        KRATOS_CATCH("")
    }

    /**
     * @brief "rInputMatrix" is SUBTRACTED to "rDestination" matrix starting from InitialRow and InitialCol of the destination matrix
     * @details "rDestination" is assumed to be able to contain the "input matrix" (no check is performed on the bounds)
     * @param rDestination The matric destination
     * @param rInputMatrix The input matrix to be computed
     * @param InitialRow The initial row to compute
     * @param InitialCol The initial column to compute
     */
    static inline void  SubtractMatrix(
        MatrixType& rDestination,
        const MatrixType& rInputMatrix,
        const IndexType InitialRow,
        const IndexType InitialCol
        )
    {
        KRATOS_TRY;

        for(IndexType i = 0; i<rInputMatrix.size1(); ++i) {
            for(IndexType j = 0; j<rInputMatrix.size2(); ++j) {
                rDestination(InitialRow+i, InitialCol+j) -= rInputMatrix(i,j);
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief "rInputMatrix" is WRITTEN on "Destination" matrix starting from InitialRow and InitialCol of the destination matrix
     * @details "Destination" is assumed to be able to contain the "input matrix" (no check is performed on the bounds)
     * @warning Destination is overwritten!!
     * @param rDestination The matric destination
     * @param rrInputMatrix The input matrix to be computed
     * @param InitialRow The initial row to compute
     * @param InitialCol The initial column to compute
     */
    static inline void  WriteMatrix(
        MatrixType& rDestination,
        const MatrixType& rInputMatrix,
        const IndexType InitialRow,
        const IndexType InitialCol
        )
    {
        KRATOS_TRY;

        for(IndexType i = 0; i < rInputMatrix.size1(); ++i) {
            for(IndexType j = 0; j < rInputMatrix.size2(); ++j) {
                rDestination(InitialRow+i, InitialCol+j) = rInputMatrix(i,j);
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs the Kroneker product of the Reduced Matrix with the identity matrix of size "dimension"
     * @param rDestination The matric destination
     * @param rReducedMatrix The reduced matrix to be computed
     * @param Dimension The dimension where we work
     */
    static inline void ExpandReducedMatrix(
        MatrixType& rDestination,
        const MatrixType& rReducedMatrix,
        const SizeType Dimension
        )
    {
        KRATOS_TRY;

        const SizeType size = rReducedMatrix.size2();
        IndexType rowindex = 0;
        IndexType colindex = 0;

        for (IndexType i = 0; i < size; ++i) {
            rowindex = i * Dimension;
            for (IndexType j = 0; j < size; ++j) {
                colindex = j * Dimension;
                for(IndexType ii = 0; ii < Dimension; ++ii) {
                    rDestination(rowindex+ii, colindex+ii) = rReducedMatrix(i, j);
                }
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs the Kroneker product of the Reduced Matrix with the identity matrix of size "dimension" ADDING to the destination matrix
     * @param rDestination The matric destination
     * @param rReducedMatrix The reduced matrix to be added
     * @param Dimension The dimension where we work
     */
    static inline void  ExpandAndAddReducedMatrix(
        MatrixType& rDestination,
        const MatrixType& rReducedMatrix,
        const SizeType Dimension
        )
    {
        KRATOS_TRY;

        const SizeType size = rReducedMatrix.size2();
        IndexType rowindex = 0;
        IndexType colindex = 0;

        for (IndexType i = 0; i < size; ++i) {
            rowindex = i * Dimension;
            for (IndexType j = 0; j < size; ++j) {
                colindex = j * Dimension;
                for(IndexType ii = 0; ii < Dimension; ++ii) {
                    rDestination(rowindex+ii, colindex+ii) += rReducedMatrix(i, j);
                }
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs rX += coeff*rY. no check on bounds is performed
     * @param rX The vector destination
     * @param rY The vector to be added
     * @param coeff The proportion to be added
     */
    static inline void  VecAdd(
        Vector& rX,
        const double coeff,
        Vector& rY
        )
    {
        KRATOS_TRY
        SizeType size=rX.size();

        for (IndexType i=0; i<size; ++i) {
            rX[i] += coeff * rY[i];
        }
        KRATOS_CATCH("")
    }

   /**
     * @brief Transforms a stess vector into a matrix. Stresses are assumed to be stored in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rStressVector the given stress vector
     * @return the corresponding stress tensor in matrix form
     * @tparam TVector The vector type considered
     * @tparam TMatrixType The matrix returning type
     */
    template<class TVector, class TMatrixType = MatrixType>
    static inline TMatrixType StressVectorToTensor(const TVector& rStressVector)
    {
        KRATOS_TRY;

        const SizeType matrix_size = rStressVector.size() == 3 ? 2 : 3;
        TMatrixType stress_tensor(matrix_size, matrix_size);

        if (rStressVector.size()==3) {
            stress_tensor(0,0) = rStressVector[0];
            stress_tensor(0,1) = rStressVector[2];
            stress_tensor(1,0) = rStressVector[2];
            stress_tensor(1,1) = rStressVector[1];
        } else if (rStressVector.size()==4) {
            stress_tensor(0,0) = rStressVector[0];
            stress_tensor(0,1) = rStressVector[3];
            stress_tensor(0,2) = 0.0;
            stress_tensor(1,0) = rStressVector[3];
            stress_tensor(1,1) = rStressVector[1];
            stress_tensor(1,2) = 0.0;
            stress_tensor(2,0) = 0.0;
            stress_tensor(2,1) = 0.0;
            stress_tensor(2,2) = rStressVector[2];
        } else if (rStressVector.size()==6) {
            stress_tensor(0,0) = rStressVector[0];
            stress_tensor(0,1) = rStressVector[3];
            stress_tensor(0,2) = rStressVector[5];
            stress_tensor(1,0) = rStressVector[3];
            stress_tensor(1,1) = rStressVector[1];
            stress_tensor(1,2) = rStressVector[4];
            stress_tensor(2,0) = rStressVector[5];
            stress_tensor(2,1) = rStressVector[4];
            stress_tensor(2,2) = rStressVector[2];
        }

        return stress_tensor;

        KRATOS_CATCH("");
    }

   /**
     * @brief Transforms a  vector into a symmetric matrix.
     * @details Components are assumed to be stored in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rVector the given stress vector
     * @return The corresponding Tensor in matrix form
     * @tparam TVector The vector type considered
     * @tparam TMatrixType The matrix returning type
     */
    template<class TVector, class TMatrixType = MatrixType>
    static inline TMatrixType VectorToSymmetricTensor(const TVector& rVector)
    {
        KRATOS_TRY;

        const SizeType matrix_size = rVector.size() == 3 ? 2 : 3;
        TMatrixType tensor(matrix_size, matrix_size);

        if (rVector.size() == 3) {
            tensor(0,0) = rVector[0];
            tensor(0,1) = rVector[2];
            tensor(1,0) = rVector[2];
            tensor(1,1) = rVector[1];
        } else if (rVector.size() == 4) {
            tensor(0,0) = rVector[0];
            tensor(0,1) = rVector[3];
            tensor(0,2) = 0.0;
            tensor(1,0) = rVector[3];
            tensor(1,1) = rVector[1];
            tensor(1,2) = 0.0;
            tensor(2,0) = 0.0;
            tensor(2,1) = 0.0;
            tensor(2,2) = rVector[2];
        } else if (rVector.size() == 6) {
            tensor(0,0) = rVector[0];
            tensor(0,1) = rVector[3];
            tensor(0,2) = rVector[5];
            tensor(1,0) = rVector[3];
            tensor(1,1) = rVector[1];
            tensor(1,2) = rVector[4];
            tensor(2,0) = rVector[5];
            tensor(2,1) = rVector[4];
            tensor(2,2) = rVector[2];
        }

        return tensor;

        KRATOS_CATCH("");
    }

    /**
     * @brief Sign function
     * @param ThisDataType The value to extract the sign
     * @return The sign of the value
     */
    static inline int Sign(const double& ThisDataType)
    {
        KRATOS_TRY;
        const double& x = ThisDataType;
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
        KRATOS_CATCH("");
    }


    /**
     * @brief Transforms a strain vector into a matrix. Strains are assumed to be stored in the following way:
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * \f$ [ e11, e22, e33, 2*e12 ] \f$ for 2D case.
     * \f$ [ e11, e22, 2*e12 ] \f$ for 2D case.
     * @details Hence the deviatoric components of the strain vector are divided by 2 while they are stored into the matrix
     * @param rStrainVector the given strain vector
     * @return the corresponding strain tensor in matrix form
     * @tparam TVector The vector type considered
     * @tparam TMatrixType The matrix returning type
     */
    template<class TVector, class TMatrixType = MatrixType>
    static inline TMatrixType StrainVectorToTensor( const TVector& rStrainVector)
    {
        KRATOS_TRY

        const SizeType matrix_size = rStrainVector.size() == 3 ? 2 : 3;
        TMatrixType strain_tensor(matrix_size, matrix_size);

        if (rStrainVector.size()==3) {
            strain_tensor(0,0) = rStrainVector[0];
            strain_tensor(0,1) = 0.5*rStrainVector[2];
            strain_tensor(1,0) = 0.5*rStrainVector[2];
            strain_tensor(1,1) = rStrainVector[1];
        } else if (rStrainVector.size()==4) {
            strain_tensor(0,0) = rStrainVector[0];
            strain_tensor(0,1) = 0.5*rStrainVector[3];
            strain_tensor(0,2) = 0;
            strain_tensor(1,0) = 0.5*rStrainVector[3];
            strain_tensor(1,1) = rStrainVector[1];
            strain_tensor(1,2) = 0;
            strain_tensor(2,0) = 0;
            strain_tensor(2,1) = 0;
            strain_tensor(2,2) = rStrainVector[2];
        } else if (rStrainVector.size()==6) {
            strain_tensor(0,0) = rStrainVector[0];
            strain_tensor(0,1) = 0.5*rStrainVector[3];
            strain_tensor(0,2) = 0.5*rStrainVector[5];
            strain_tensor(1,0) = 0.5*rStrainVector[3];
            strain_tensor(1,1) = rStrainVector[1];
            strain_tensor(1,2) = 0.5*rStrainVector[4];
            strain_tensor(2,0) = 0.5*rStrainVector[5];
            strain_tensor(2,1) = 0.5*rStrainVector[4];
            strain_tensor(2,2) = rStrainVector[2];
        }

        return strain_tensor;

        KRATOS_CATCH("");
    }

    /**
     * @brief Transforms a given symmetric Strain Tensor to Voigt Notation:
     * @details The following cases:
     *  - In the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     *    \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     *  - In the 2D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     *    \f$ [ e11, e22, e33, 2*e12 ] \f$ fir 2D case.
     *  - In the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     *    \f$ [ e11, e22, 2*e12 ] \f$ fir 2D case.
     * @param rStrainTensor the given symmetric second order strain tensor
     * @return the corresponding strain tensor in vector form
     * @tparam TMatrixType The matrix type considered
     * @tparam TVector The vector returning type
     */
    template<class TMatrixType, class TVector = Vector>
    static inline Vector StrainTensorToVector(
        const TMatrixType& rStrainTensor,
        SizeType rSize = 0
        )
    {
        KRATOS_TRY;

        if(rSize == 0) {
            if(rStrainTensor.size1() == 2) {
                rSize = 3;
            } else if(rStrainTensor.size1() == 3) {
                rSize = 6;
            }
        }

        Vector strain_vector(rSize);

        if (rSize == 3) {
            strain_vector[0] = rStrainTensor(0,0);
            strain_vector[1] = rStrainTensor(1,1);
            strain_vector[2] = 2.0*rStrainTensor(0,1);
        } else if (rSize == 4) {
            strain_vector[0] = rStrainTensor(0,0);
            strain_vector[1] = rStrainTensor(1,1);
            strain_vector[2] = rStrainTensor(2,2);
            strain_vector[3] = 2.0*rStrainTensor(0,1);
        } else if (rSize == 6) {
            strain_vector[0] = rStrainTensor(0,0);
            strain_vector[1] = rStrainTensor(1,1);
            strain_vector[2] = rStrainTensor(2,2);
            strain_vector[3] = 2.0*rStrainTensor(0,1);
            strain_vector[4] = 2.0*rStrainTensor(1,2);
            strain_vector[5] = 2.0*rStrainTensor(0,2);
        }

        return strain_vector;

        KRATOS_CATCH("");
     }

    /**
     * @brief Transforms a given symmetric Stress Tensor to Voigt Notation:
     * @details Components are assumed to be stored in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * In the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * In the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * In the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     * @param rStressTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     * @tparam TMatrixType The matrix type considered
     * @tparam TVector The vector returning type
     */
    template<class TMatrixType, class TVector = Vector>
    static inline TVector StressTensorToVector(
        const TMatrixType& rStressTensor,
        SizeType rSize = 0
        )
    {
        KRATOS_TRY;

        if(rSize == 0) {
            if(rStressTensor.size1() == 2) {
                rSize = 3;
            } else if(rStressTensor.size1() == 3) {
                rSize = 6;
            }
        }

        TVector stress_vector(rSize);

        if (rSize == 3) {
            stress_vector[0] = rStressTensor(0,0);
            stress_vector[1] = rStressTensor(1,1);
            stress_vector[2] = rStressTensor(0,1);
        } else if (rSize == 4) {
            stress_vector[0] = rStressTensor(0,0);
            stress_vector[1] = rStressTensor(1,1);
            stress_vector[2] = rStressTensor(2,2);
            stress_vector[3] = rStressTensor(0,1);
        } else if (rSize == 6) {
            stress_vector[0] = rStressTensor(0,0);
            stress_vector[1] = rStressTensor(1,1);
            stress_vector[2] = rStressTensor(2,2);
            stress_vector[3] = rStressTensor(0,1);
            stress_vector[4] = rStressTensor(1,2);
            stress_vector[5] = rStressTensor(0,2);
        }

        return stress_vector;

        KRATOS_CATCH("");
     }

    /**
     * @brief Transforms a given symmetric Tensor to Voigt Notation:
     * @details The following cases:
     *  - In the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     *  - In the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     *  - In the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     * @param rTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     * @tparam TMatrixType The matrix type considered
     * @tparam TVector The vector returning type
     */
    template<class TMatrixType, class TVector = Vector>
    static inline TVector SymmetricTensorToVector(
        const TMatrixType& rTensor,
        SizeType rSize = 0
        )
    {
        KRATOS_TRY;

        if(rSize == 0) {
            if(rTensor.size1() == 2) {
                rSize = 3;
            } else if(rTensor.size1() == 3) {
                rSize = 6;
            }
        }

        Vector vector(rSize);

        if (rSize == 3) {
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(0,1);

        } else if (rSize==4) {
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(2,2);
            vector[3]= rTensor(0,1);
        } else if (rSize==6) {
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(2,2);
            vector[3]= rTensor(0,1);
            vector[4]= rTensor(1,2);
            vector[5]= rTensor(0,2);
        }

        return vector;

        KRATOS_CATCH("");
     }

    /**
     * @brief Calculates the product operation B'DB
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @tparam TMatrixType1 The type of matrix considered (1)
     * @tparam TMatrixType2 The type of matrix considered (2)
     * @tparam TMatrixType3 The type of matrix considered (3)
     */
    template<class TMatrixType1, class TMatrixType2, class TMatrixType3>
    static inline void BtDBProductOperation(
        TMatrixType1& rA,
        const TMatrixType2& rD,
        const TMatrixType3& rB
        )
    {
        // The sizes
        const SizeType size1 = rB.size2();
        const SizeType size2 = rB.size2();

        if (rA.size1() != size1 || rA.size2() != size2)
            rA.resize(size1, size2, false);

        // Direct multiplication
        // noalias(rA) = prod( trans( rB ), MatrixType(prod(rD, rB)));

        // Manual multiplication
        rA.clear();
        for(IndexType k = 0; k< rD.size1(); ++k) {
            for(IndexType l = 0; l < rD.size2(); ++l) {
                const double Dkl = rD(k, l);
                for(IndexType j = 0; j < rB.size2(); ++j) {
                    const double DklBlj = Dkl * rB(l, j);
                    for(IndexType i = 0; i< rB.size2(); ++i) {
                        rA(i, j) += rB(k, i) * DklBlj;
                    }
                }
            }
        }
    }

    /**
     * @brief Calculates the product operation BDB'
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @tparam TMatrixType1 The type of matrix considered (1)
     * @tparam TMatrixType2 The type of matrix considered (2)
     * @tparam TMatrixType3 The type of matrix considered (3)
     */
    template<class TMatrixType1, class TMatrixType2, class TMatrixType3>
    static inline void BDBtProductOperation(
        TMatrixType1& rA,
        const TMatrixType2& rD,
        const TMatrixType3& rB
        )
    {
        // The sizes
        const SizeType size1 = rB.size1();
        const SizeType size2 = rB.size1();

        if (rA.size1() != size1 || rA.size2() != size2)
            rA.resize(size1, size2, false);

        // Direct multiplication
        // noalias(rA) = prod(rB, MatrixType(prod(rD, trans(rB))));

        // Manual multiplication
        rA.clear();
        for(IndexType k = 0; k< rD.size1(); ++k) {
            for(IndexType l = 0; l < rD.size2(); ++l) {
                const double Dkl = rD(k,l);
                for(IndexType j = 0; j < rB.size1(); ++j) {
                    const double DklBjl = Dkl * rB(j,l);
                    for(IndexType i = 0; i< rB.size1(); ++i) {
                        rA(i, j) += rB(i, k) * DklBjl;
                    }
                }
            }
        }
    }

    /**
     * @brief Calculates the eigenvectors and eigenvalues of given symmetric matrix
     * @details The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method. The resulting decomposition is LDL'
     * @note See https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
     * @param rA The given symmetric matrix the eigenvectors are to be calculated.
     * @param rEigenVectorsMatrix The result matrix (will be overwritten with the eigenvectors)
     * @param rEigenValuesMatrix The result diagonal matrix with the eigenvalues
     * @param Tolerance The largest value considered to be zero
     * @param MaxIterations Maximum number of iterations
     * @tparam TMatrixType1 The type of matrix considered (1)
     * @tparam TMatrixType2 The type of matrix considered (2)
     */
    template<class TMatrixType1, class TMatrixType2>
    static inline bool GaussSeidelEigenSystem(
        const TMatrixType1& rA,
        TMatrixType2& rEigenVectorsMatrix,
        TMatrixType2& rEigenValuesMatrix,
        const double Tolerance = 1.0e-18,
        const SizeType MaxIterations = 20
        )
    {
        bool is_converged = false;

        const SizeType size = rA.size1();

        if (rEigenVectorsMatrix.size1() != size || rEigenVectorsMatrix.size2() != size)
            rEigenVectorsMatrix.resize(size, size, false);
        if (rEigenValuesMatrix.size1() != size || rEigenValuesMatrix.size2() != size)
            rEigenValuesMatrix.resize(size, size, false);

        const TMatrixType2 identity_matrix = IdentityMatrix(size);
        noalias(rEigenVectorsMatrix) = identity_matrix;
        noalias(rEigenValuesMatrix) = rA;

        // Auxiliar values
        TMatrixType2 aux_A, aux_V_matrix, rotation_matrix;
        double a, u, c, s, gamma, teta;
        IndexType index1, index2;

        aux_A.resize(size,size,false);
        aux_V_matrix.resize(size,size,false);
        rotation_matrix.resize(size,size,false);

        for(IndexType iterations = 0; iterations < MaxIterations; ++iterations) {
            is_converged = true;

            a = 0.0;
            index1 = 0;
            index2 = 1;

            for(IndexType i = 0; i < size; ++i) {
                for(IndexType j = (i + 1); j < size; ++j) {
                    if((std::abs(rEigenValuesMatrix(i, j)) > a ) && (std::abs(rEigenValuesMatrix(i, j)) > Tolerance)) {
                        a = std::abs(rEigenValuesMatrix(i,j));
                        index1 = i;
                        index2 = j;
                        is_converged = false;
                    }
                }
            }

            if(is_converged) {
                break;
            }

            // Calculation of Rotation angle
            gamma = (rEigenValuesMatrix(index2, index2)-rEigenValuesMatrix(index1, index1)) / (2 * rEigenValuesMatrix(index1, index2));
            u = 1.0;

            if(std::abs(gamma) > Tolerance && std::abs(gamma)< (1.0/Tolerance)) {
                u = gamma / std::abs(gamma) * 1.0 / (std::abs(gamma) + std::sqrt(1.0 + gamma * gamma));
            } else {
                if  (std::abs(gamma) >= (1.0/Tolerance)) {
                    u = 0.5 / gamma;
                }
            }

            c = 1.0 / (std::sqrt(1.0 + u * u));
            s = c * u;
            teta = s / (1.0 + c);

            // Rotation of the Matrix
            noalias(aux_A) = rEigenValuesMatrix;
            aux_A(index2, index2) = rEigenValuesMatrix(index2,index2) + u * rEigenValuesMatrix(index1, index2);
            aux_A(index1, index1) = rEigenValuesMatrix(index1,index1) - u * rEigenValuesMatrix(index1, index2);
            aux_A(index1, index2) = 0.0;
            aux_A(index2, index1) = 0.0;

            for(IndexType i = 0; i < size; ++i) {
                if((i!= index1) && (i!= index2)) {
                    aux_A(index2, i) = rEigenValuesMatrix(index2, i) + s * (rEigenValuesMatrix(index1, i)- teta * rEigenValuesMatrix(index2, i));
                    aux_A(i, index2) = rEigenValuesMatrix(index2, i) + s * (rEigenValuesMatrix(index1, i)- teta * rEigenValuesMatrix(index2, i));
                    aux_A(index1, i) = rEigenValuesMatrix(index1, i) - s * (rEigenValuesMatrix(index2, i) + teta * rEigenValuesMatrix(index1, i));
                    aux_A(i, index1) = rEigenValuesMatrix(index1, i) - s * (rEigenValuesMatrix(index2, i) + teta * rEigenValuesMatrix(index1, i));
                }
            }

            noalias(rEigenValuesMatrix) = aux_A;

            // Calculation of the eigeneigen vector matrix V
            noalias(rotation_matrix) = identity_matrix;
            rotation_matrix(index2, index1) = -s;
            rotation_matrix(index1, index2) =  s;
            rotation_matrix(index1, index1) =  c;
            rotation_matrix(index2, index2) =  c;

            noalias(aux_V_matrix) = ZeroMatrix(size, size);

            for(IndexType i = 0; i < size; ++i) {
                for(IndexType j = 0; j < size; ++j) {
                    for(IndexType k = 0; k < size; ++k) {
                        aux_V_matrix(i, j) += rEigenVectorsMatrix(i, k) * rotation_matrix(k, j);
                    }
                }
            }
            noalias(rEigenVectorsMatrix) = aux_V_matrix;
        }

        KRATOS_WARNING_IF("MathUtils::EigenSystem", !is_converged) << "Spectral decomposition not converged " << std::endl;

        return is_converged;
    }

    /**
     * @brief Calculates the eigenvectors and eigenvalues of given symmetric TDimxTDim matrix
     * @details The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method. The resulting decomposition is L'DL
     * @note See https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
     * @param A The given symmetric matrix the eigenvectors are to be calculated.
     * @param rEigenVectorsMatrix The result matrix (will be overwritten with the eigenvectors)
     * @param rEigenValuesMatrix The result diagonal matrix with the eigenvalues
     * @param Tolerance The largest value considered to be zero
     * @param MaxIterations Maximum number of iterations
     * @tparam TDim The size of the bounded matrix
     * @warning This method is deprecated. Will be removed soon
     */
    template<SizeType TDim>
    KRATOS_DEPRECATED_MESSAGE("Please use GaussSeidelEigenSystem() instead. Note the resulting EigenVectors matrix is transposed respect GaussSeidelEigenSystem()")
    static inline bool EigenSystem(
        const BoundedMatrix<double, TDim, TDim>& rA,
        BoundedMatrix<double, TDim, TDim>& rEigenVectorsMatrix,
        BoundedMatrix<double, TDim, TDim>& rEigenValuesMatrix,
        const double Tolerance = 1.0e-18,
        const SizeType MaxIterations = 20
        )
    {
        const bool is_converged = GaussSeidelEigenSystem(rA, rEigenVectorsMatrix, rEigenValuesMatrix, Tolerance, MaxIterations);

        const BoundedMatrix<double, TDim, TDim> V_matrix = rEigenVectorsMatrix;
        for(IndexType i = 0; i < TDim; ++i) {
            for(IndexType j = 0; j < TDim; ++j) {
                rEigenVectorsMatrix(i, j) = V_matrix(j, i);
            }
        }

        return is_converged;
    }

    /**
     * @brief Calculates the square root of a matrix
     * @details This function calculates the square root of a matrix by doing an eigenvalue decomposition
     * The square root of a matrix A is defined as A = V*S*inv(V) where A is the eigenvectors matrix
     * and S the diagonal matrix containing the square root of the eigenvalues. Note that the previous
     * expression can be rewritten as A = V*S*trans(V) since V is orthogonal.
     * @tparam TMatrixType1 Input matrix type
     * @tparam TMatrixType2 Output matrix type
     * @param rA Input matrix
     * @param rMatrixSquareRoot Square root output matrix
     * @param Tolerance Tolerance of the eigenvalue decomposition
     * @param MaxIterations Maximum iterations of the eigenvalue decomposition
     * @return true The eigenvalue decomposition problem converged
     * @return false The eigenvalue decomposition problem did not converge
     */
    template<class TMatrixType1, class TMatrixType2>
    static inline bool MatrixSquareRoot(
        const TMatrixType1 &rA,
        TMatrixType2 &rMatrixSquareRoot,
        const double Tolerance = 1.0e-18,
        const SizeType MaxIterations = 20
        )
    {
        // Do an eigenvalue decomposition of the input matrix
        TMatrixType2 eigenvectors_matrix, eigenvalues_matrix;
        const bool is_converged = GaussSeidelEigenSystem(rA, eigenvectors_matrix, eigenvalues_matrix, Tolerance, MaxIterations);
        KRATOS_WARNING_IF("MatrixSquareRoot", !is_converged) << "GaussSeidelEigenSystem did not converge.\n";

        // Get the square root of the eigenvalues
        SizeType size = eigenvalues_matrix.size1();
        for (SizeType i = 0; i < size; ++i) {
            KRATOS_ERROR_IF(eigenvalues_matrix(i,i) < 0) << "Eigenvalue " << i << " is negative. Square root matrix cannot be computed" << std::endl;
            eigenvalues_matrix(i,i) = std::sqrt(eigenvalues_matrix(i,i));
        }

        // Calculate the solution from the previous decomposition and eigenvalues square root
        BDBtProductOperation(rMatrixSquareRoot, eigenvalues_matrix, eigenvectors_matrix);

        return is_converged;
    }

    /**
     * @brief Calculates the Factorial of a number k, Factorial = k!
     * @tparam Number The number of which the Factorial is computed
     */
    template<class TIntegerType>
    static inline TIntegerType Factorial(const TIntegerType Number)
    {
        if (Number == 0) {
            return 1;
        }
        TIntegerType k = Number;
        for (TIntegerType i = Number - 1; i > 0; --i){
            k *= i;
        }
        return k;
    }

    /**
     * @brief Calculates the exponential of a matrix
     * @brief see https://mathworld.wolfram.com/MatrixExponential.html
     * @tparam rMatrix: the matrix A of which exp is calculated
     * @tparam rExponentialMatrix: exp(A)
     */
    template<class TMatrixType>
    static inline void CalculateExponentialOfMatrix(
        const TMatrixType& rMatrix,
        TMatrixType& rExponentialMatrix,
        const double Tolerance = 1000.0*ZeroTolerance,
        const SizeType MaxTerms = 200
        )
    {
        SizeType series_term = 2;
        SizeType factorial = 1;
        const SizeType dimension = rMatrix.size1();

        noalias(rExponentialMatrix) = IdentityMatrix(dimension) + rMatrix;
        TMatrixType exponent_matrix = rMatrix;
        TMatrixType aux_matrix;

        while (series_term < MaxTerms) {
            noalias(aux_matrix) = prod(exponent_matrix, rMatrix);
            noalias(exponent_matrix) = aux_matrix;
            factorial = Factorial(series_term);
            noalias(rExponentialMatrix) += exponent_matrix / factorial;
            const double norm_series_term = std::abs(norm_frobenius(exponent_matrix) / factorial);
            if (norm_series_term < Tolerance)
                break;
            series_term++;
        }
    }


    static double DegreesToRadians(double AngleInDegrees)
    {
        return (AngleInDegrees * Globals::Pi) / 180.0;
    }

    ///@}

private:

    ///@name Unaccessible methods
    ///@{

    MathUtils(void);

    MathUtils(MathUtils& rSource);

    ///@}
}; /* Class MathUtils */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/
