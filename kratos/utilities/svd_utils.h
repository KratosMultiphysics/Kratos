//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class SVDUtils
 * @ingroup KratosCore
 * @brief Various mathematical utilities to compute SVD and the condition number of a matrix
 * @details Defines several utility functions:
 * - SingularValueDecomposition: Computes the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
 * - JacobiSingularValueDecomposition: Computes the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
 * - SVDConditionNumber: Computes the condition number of a given matrix using the SVD
 * @author Vicente Mataix Ferrandiz
 * @tparam TDataType Type of the data stored in the matrix.
 */
template<class TDataType>
class SVDUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the matrix type
    using MatrixType = Matrix;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of index type
    using IndexType = std::size_t;

    /// Definition of epsilon zero tolerance
    constexpr static TDataType ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function gives the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite inefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param rInputMatrix The matrix where perform the SVD
     * @param rUMatrix The unitary U matrix
     * @param rSMatrix The diagonal S matrix
     * @param rVMatrix The unitary V matrix
     * @param ThisParameters The configuration parameters
     * @return iter: The number of iterations
     * @tparam TInputMatrix The input matrix type
     * @tparam TUmatrix The U matrix type
     * @tparam TSMatrix The S matrix type
     * @tparam TVMatrix The V matrix type
     */
    template<class TInputMatrix, class TUmatrix, class TSMatrix, class TVMatrix>
    static inline std::size_t SingularValueDecomposition(
        const TInputMatrix& rInputMatrix,
        TUmatrix& rUMatrix,
        TSMatrix& rSMatrix,
        TVMatrix& rVMatrix,
        Parameters ThisParameters
        )
    {
        // Validating defaults
        Parameters default_parameters = Parameters(R"({
            "type_svd"             : "Jacobi",
            "tolerance"            : 0.0,
            "max_iter"             : 200
        })");
        default_parameters["tolerance"].SetDouble(ZeroTolerance);
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        const std::string& r_type_svd = ThisParameters["type_svd"].GetString();
        const double tolerance = ThisParameters["tolerance"].GetDouble();
        const double max_iter = ThisParameters["max_iter"].GetInt();
        return SingularValueDecomposition(rInputMatrix, rUMatrix, rSMatrix, rVMatrix, r_type_svd, tolerance, max_iter);
    }

    /**
     * @brief This function gives the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite inefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param rInputMatrix The matrix where perform the SVD
     * @param rUMatrix The unitary U matrix
     * @param rSMatrix The diagonal S matrix
     * @param rVMatrix The unitary V matrix
     * @param TypeSVD The type of SVD algorithm (Jacobi by default)
     * @param Tolerance The tolerance considered
     * @param MaxIter Maximum number of iterations
     * @return iter: The number of iterations
     * @tparam TInputMatrix The input matrix type
     * @tparam TUmatrix The U matrix type
     * @tparam TSMatrix The S matrix type
     * @tparam TVMatrix The V matrix type
     */
    template<class TInputMatrix, class TUmatrix, class TSMatrix, class TVMatrix>
    static inline std::size_t SingularValueDecomposition(
        const TInputMatrix& rInputMatrix,
        TUmatrix& rUMatrix,
        TSMatrix& rSMatrix,
        TVMatrix& rVMatrix,
        const std::string& TypeSVD = "Jacobi",
        const TDataType Tolerance = ZeroTolerance,
        const IndexType MaxIter = 200
        )
    {
        if (TypeSVD == "Jacobi") {
            return JacobiSingularValueDecomposition(rInputMatrix, rUMatrix, rSMatrix, rVMatrix, Tolerance, MaxIter);
        } else {
            KRATOS_ERROR << "SVD Type not implemented" << std::endl;
        }
    }

    /**
     * @brief This function gives the Jacobi SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite inefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param rInputMatrix The matrix where perform the SVD
     * @param rUMatrix The unitary U matrix
     * @param rSMatrix The diagonal S matrix
     * @param rVMatrix The unitary V matrix
     * @param Tolerance The tolerance considered
     * @param MaxIter Maximum number of iterations
     * @return iter: The number of iterations
     * @tparam TInputMatrix The input matrix type
     * @tparam TUmatrix The U matrix type
     * @tparam TSMatrix The S matrix type
     * @tparam TVMatrix The V matrix type
     */
    template<class TInputMatrix, class TUmatrix, class TSMatrix, class TVMatrix>
    static inline IndexType JacobiSingularValueDecomposition(
        const TInputMatrix& rInputMatrix,
        TUmatrix& rUMatrix,
        TSMatrix& rSMatrix,
        TVMatrix& rVMatrix,
        const TDataType Tolerance = ZeroTolerance,
        const IndexType MaxIter = 200
        )
    {
        // Define sizes
        const SizeType m = rInputMatrix.size1();
        const SizeType n = rInputMatrix.size2();
        KRATOS_ERROR_IF(m != n) << "Current Jacobi implementation only works with square matrices. Use \'LinearSolversApplication\' decompositions instead." << std::endl;

        // If dynamic matrix (we resize)
        if constexpr(std::is_same_v<TUmatrix, MatrixType>) {
            if (rUMatrix.size1() != m || rUMatrix.size2() != m) {
                rUMatrix.resize(m, m, false);
            }
        } else { // We will assume bounded matrices (we check the size)
            KRATOS_ERROR_IF(rUMatrix.size1() != m || rUMatrix.size2() != m) << "The size of the U matrix is not correct: " << rUMatrix.size1() << " vs " << m << " and " << rUMatrix.size2() << " vs " << m << std::endl;
        }

        // If dynamic matrix (we resize)
        if constexpr(std::is_same_v<TSMatrix, MatrixType>) {
            if (rSMatrix.size1() != m || rSMatrix.size2() != n) {
                rSMatrix.resize(m, n, false);
            }
        } else { // We will assume bounded matrices (we check the size)
            KRATOS_ERROR_IF(rSMatrix.size1() != m || rSMatrix.size2() != n) << "The size of the S matrix is not correct: " << rSMatrix.size1() << " vs " << m << " and " << rSMatrix.size2() << " vs " << n << std::endl;
        }

        // If dynamic matrix (we resize)
        if constexpr(std::is_same_v<TVMatrix, MatrixType>) {
            if (rVMatrix.size1() != n || rVMatrix.size2() != n) {
                rVMatrix.resize(n, n, false);
            }
        } else { // We will assume bounded matrices (we check the size)
            KRATOS_ERROR_IF(rVMatrix.size1() != n || rVMatrix.size2() != n) << "The size of the V matrix is not correct: " << rVMatrix.size1() << " vs " << n << " and " << rVMatrix.size2() << " vs " << n << std::endl;
        }

        // Assign values
        noalias(rSMatrix) = rInputMatrix;
        noalias(rUMatrix) = IdentityMatrix(m);
        noalias(rVMatrix) = IdentityMatrix(n);

        // If dynamic matrix
        if constexpr(std::is_same_v<TInputMatrix, MatrixType>) {
            // We create the auxiliary matrices (for aliased operations)
            MatrixType auxiliary_matrix_mn(m, n);
            MatrixType auxiliary_matrix_m(m, m);
            MatrixType auxiliary_matrix_n(n, n);

            // The Jacobi matrices
            MatrixType j1(m, m);
            MatrixType j2(n, n);

            return ComputeJacobiIterations(rInputMatrix, rUMatrix, rSMatrix, rVMatrix, j1, j2, auxiliary_matrix_mn, auxiliary_matrix_m, auxiliary_matrix_n, Tolerance, MaxIter);
        } else { // We will assume bounded matrices
            // We create the auxiliary matrices (for aliased operations)
            BoundedMatrix<TDataType, TInputMatrix::max_size1, TInputMatrix::max_size2> auxiliary_matrix_mn;
            BoundedMatrix<TDataType, TInputMatrix::max_size1, TInputMatrix::max_size1> auxiliary_matrix_m(m, m);
            BoundedMatrix<TDataType, TInputMatrix::max_size2, TInputMatrix::max_size2> auxiliary_matrix_n(n, n);

            // The Jacobi matrices
            BoundedMatrix<TDataType, TInputMatrix::max_size1, TInputMatrix::max_size1> j1(m, m);
            BoundedMatrix<TDataType, TInputMatrix::max_size2, TInputMatrix::max_size2> j2(n, n);

            return ComputeJacobiIterations(rInputMatrix, rUMatrix, rSMatrix, rVMatrix, j1, j2, auxiliary_matrix_mn, auxiliary_matrix_m, auxiliary_matrix_n, Tolerance, MaxIter);
        }
    }

    /**
     * @brief This method computes the condition number using the SVD
     * @details The condition number can be estimated as the ratio between the largest singular value and the smallest singular value
     * @param rInputMatrix The matrix to be evaluated
     * @param Tolerance The tolerance considered
     * @return condition_number: The ratio between the largest SV and the smallest SV
     */
    template<class TInputMatrix>
    static inline TDataType SVDConditionNumber(
        const TInputMatrix& rInputMatrix,
        const std::string TypeSVD = "Jacobi",
        const TDataType Tolerance = ZeroTolerance,
        const IndexType MaxIter = 200
        )
    {
        // We define the condition number
        TDataType condition_number = 0.0;

        // We get the size of the matrix S matrix, equal to the number of rows of the input matrix
        const SizeType size_s = rInputMatrix.size1();

        // If dynamic matrix
        if constexpr(std::is_same_v<TInputMatrix, MatrixType>) {
            MatrixType u_matrix, s_matrix, v_matrix;
            SingularValueDecomposition(rInputMatrix, u_matrix, s_matrix, v_matrix, TypeSVD, Tolerance, MaxIter);

            condition_number = s_matrix(0, 0)/s_matrix(size_s - 1, size_s - 1);
        } else { // We will assume bounded matrices
            BoundedMatrix<TDataType, TInputMatrix::max_size1, TInputMatrix::max_size2> s_matrix;
            BoundedMatrix<TDataType, TInputMatrix::max_size1, TInputMatrix::max_size1> u_matrix;
            BoundedMatrix<TDataType, TInputMatrix::max_size2, TInputMatrix::max_size2> v_matrix;
            SingularValueDecomposition(rInputMatrix, u_matrix, s_matrix, v_matrix, TypeSVD, Tolerance, MaxIter);

            condition_number = s_matrix(0, 0)/s_matrix(size_s - 1, size_s - 1);
        }


        return condition_number;
    }

    ///@}
private:
    ///@name Unaccessible methods
    ///@{

    SVDUtils(void);

    SVDUtils(SVDUtils& rSource);

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function gives the Jacobi SVD of a given 2x2 matrix, returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1
     * @param rInputMatrix The matrix where perform the SVD
     * @param rUMatrix The unitary U matrix
     * @param rSMatrix The diagonal S matrix
     * @param rVMatrix The unitary V matrix
     */
    static inline void SingularValueDecomposition2x2(
        const BoundedMatrix<TDataType, 2, 2>& rInputMatrix,
        BoundedMatrix<TDataType, 2, 2>& rUMatrix,
        BoundedMatrix<TDataType, 2, 2>& rSMatrix,
        BoundedMatrix<TDataType, 2, 2>& rVMatrix
        )
    {
        // If the trace (sum of diagonal elements) is approximately zero, handle as a special case
        if (std::abs(rInputMatrix(0, 0) + rInputMatrix(1, 1)) < ZeroTolerance) {
            noalias(rSMatrix) = rInputMatrix;      // S is just the original matrix
            noalias(rUMatrix) = IdentityMatrix(2); // U is identity
            noalias(rVMatrix) = rUMatrix;          // V is also identity
        } else {
            const TDataType t = (rInputMatrix(0, 1) - rInputMatrix(1, 0))/(rInputMatrix(0, 0) + rInputMatrix(1, 1));
            const TDataType c = 1.0/std::sqrt(1.0 + t*t);
            const TDataType s = t*c;
            BoundedMatrix<TDataType, 2, 2> r_matrix;
            r_matrix(0, 0) =  c;
            r_matrix(0, 1) = -s;
            r_matrix(1, 0) =  s;
            r_matrix(1, 1) =  c;

            BoundedMatrix<TDataType, 2, 2> m_matrix = prod(r_matrix, rInputMatrix);

            SingularValueDecomposition2x2Symmetric(m_matrix, rUMatrix, rSMatrix, rVMatrix);

            BoundedMatrix<TDataType, 2, 2> auxiliary_matrix_m;
            noalias(auxiliary_matrix_m) = prod(trans(r_matrix), rUMatrix);
            noalias(rUMatrix) = auxiliary_matrix_m;
        }
    }

    /**
     * @brief Computes the Jacobi Singular Value Decomposition (SVD) of a given 2x2 symmetric matrix.
     * @details This function decomposes the input matrix A into three matrices: U, S, and V such that:
     *      A = U * S * V^T
     * where:
     * - U and V are unitary (orthogonal for real matrices).
     * - S is a diagonal matrix containing singular values, which are always non-negative and sorted in descending order.
     * @tparam TDataType The numerical type used for matrix entries (e.g., double, float).
     * @param rInputMatrix The 2x2 symmetric matrix to be decomposed.
     * @param rUMatrix Output matrix representing the unitary U matrix.
     * @param rSMatrix Output matrix representing the diagonal S matrix with singular values.
     * @param rVMatrix Output matrix representing the unitary V matrix.
     */
    static inline void SingularValueDecomposition2x2Symmetric(
        const BoundedMatrix<TDataType, 2, 2>& rInputMatrix,
        BoundedMatrix<TDataType, 2, 2>& rUMatrix,
        BoundedMatrix<TDataType, 2, 2>& rSMatrix,
        BoundedMatrix<TDataType, 2, 2>& rVMatrix
        )
    {
        // If the off-diagonal element is approximately zero, the matrix is already diagonal
        if (std::abs(rInputMatrix(1, 0)) < ZeroTolerance) {
            noalias(rSMatrix) = rInputMatrix;      // S is just the original matrix
            noalias(rUMatrix) = IdentityMatrix(2); // U is identity
            noalias(rVMatrix) = rUMatrix;          // V is also identity
        } else {
            // Extract elements from the symmetric 2x2 matrix
            const TDataType w = rInputMatrix(0, 0); // Top-left element
            const TDataType y = rInputMatrix(1, 0); // Off-diagonal element
            const TDataType z = rInputMatrix(1, 1); // Bottom-right element

            // Compute the Jacobi rotation parameters
            const TDataType ro = (z - w) / (2.0 * y);
            const TDataType t = MathUtils<TDataType>::Sign(ro) / (std::abs(ro) + std::sqrt(1 + ro * ro));
            const TDataType c = 1.0 / std::sqrt(1.0 + t * t); // Cosine of rotation angle
            const TDataType s = t * c; // Sine of rotation angle

            // Construct the unitary rotation matrix U
            rUMatrix(0, 0) =  c;
            rUMatrix(0, 1) =  s;
            rUMatrix(1, 0) = -s;
            rUMatrix(1, 1) =  c;

            // Compute V as the transpose of U
            noalias(rVMatrix) = trans(rUMatrix);

            // Compute the diagonal matrix S: S = U^T * A * V
            noalias(rSMatrix) = prod(trans(rUMatrix), BoundedMatrix<TDataType, 2, 2>(prod(rInputMatrix, trans(rVMatrix))));
        }

        // Construct a diagonal sign matrix to ensure singular values are non-negative
        BoundedMatrix<TDataType, 2, 2> z_matrix;
        z_matrix(0, 0) = MathUtils<TDataType>::Sign(rSMatrix(0, 0)); // Sign of first singular value
        z_matrix(0, 1) = 0.0;
        z_matrix(1, 0) = 0.0;
        z_matrix(1, 1) = MathUtils<TDataType>::Sign(rSMatrix(1, 1)); // Sign of second singular value

        // Apply sign correction to U and S
        BoundedMatrix<TDataType, 2, 2> aux_2_2_matrix;
        noalias(aux_2_2_matrix) = prod(rUMatrix, z_matrix);
        noalias(rUMatrix) = aux_2_2_matrix;
        noalias(aux_2_2_matrix) = prod(z_matrix, rSMatrix);
        noalias(rSMatrix) = aux_2_2_matrix;

        // Ensure the singular values are sorted in descending order
        if (rSMatrix(0, 0) < rSMatrix(1, 1)) {
            // Define a permutation matrix to swap the singular values
            BoundedMatrix<TDataType, 2, 2> p_matrix;
            p_matrix(0, 0) = 0.0;
            p_matrix(0, 1) = 1.0;
            p_matrix(1, 0) = 1.0;
            p_matrix(1, 1) = 0.0;

            // Apply permutation to U, S, and V to reorder singular values
            noalias(aux_2_2_matrix) = prod(rUMatrix, p_matrix);
            noalias(rUMatrix) = aux_2_2_matrix;
            const BoundedMatrix<TDataType, 2, 2> aux_matrix = prod(rSMatrix, p_matrix);
            noalias(rSMatrix) = prod(p_matrix, aux_matrix);
            noalias(aux_2_2_matrix) = prod(p_matrix, rVMatrix);
            noalias(rVMatrix) = aux_2_2_matrix;
        }
    }

    /**
     * @brief This method computes the Jacobi rotation operation
     * @param J1 First Jacobi matrix
     * @param J2 Second Jacobi matrix
     * @param rInputMatrix The matrix to compute the Jacobi tolerance
     * @param Size1 The size of the matrix (number of rows)
     * @param Size2 The size of the matrix (number of columns)
     * @param Index1 The index to compute (row)
     * @param Index2 The index to compute (column)
     * @tparam TMatrixType1 The first Jacobi matrix type
     * @tparam TMatrixType2 The second Jacobi matrix type
     * @tparam TInputMatrix The input matrix type
     */
    template<class TMatrixType1, class TMatrixType2, class TInputMatrix>
    static inline void Jacobi(
        TMatrixType1& J1,
        TMatrixType2& J2,
        const TInputMatrix& rInputMatrix,
        const SizeType Size1,
        const SizeType Size2,
        const SizeType Index1,
        const SizeType Index2
        )
    {
        // We define the 2x2 matrix
        BoundedMatrix<TDataType, 2, 2> b_matrix;
        b_matrix(0, 0) = rInputMatrix(Index1, Index1);
        b_matrix(0, 1) = rInputMatrix(Index1, Index2);
        b_matrix(1, 0) = rInputMatrix(Index2, Index1);
        b_matrix(1, 1) = rInputMatrix(Index2, Index2);

        // We compute the singular value decomposition for a 2x2 matrix
        BoundedMatrix<TDataType, 2, 2> u_matrix, s_matrix, v_matrix;
        SingularValueDecomposition2x2(b_matrix, u_matrix, s_matrix, v_matrix);

        // We compute the Jacobi rotation
        J1 = IdentityMatrix(Size1);
        J1(Index1, Index1) = u_matrix(0, 0);
        J1(Index1, Index2) = u_matrix(1, 0);
        J1(Index2, Index1) = u_matrix(0, 1);
        J1(Index2, Index2) = u_matrix(1, 1);

        // We compute the Jacobi rotation
        J2 = IdentityMatrix(Size2);
        J2(Index1, Index1) = v_matrix(0, 0);
        J2(Index1, Index2) = v_matrix(1, 0);
        J2(Index2, Index1) = v_matrix(0, 1);
        J2(Index2, Index2) = v_matrix(1, 1);
    }

    /**
     * @brief This method computes the Jacobi rotation operation
     * @param J1 First Jacobi matrix
     * @param rInputMatrix The matrix to compute the Jacobi tolerance
     * @param Size1 The size of the matrix (number of rows)
     * @param Size2 The size of the matrix (number of columns)
     * @param Index1 The index to compute (row)
     * @param Index2 The index to compute (column)
     * @tparam TMatrixType1 The first Jacobi matrix type
     * @tparam TInputMatrix The input matrix type
     */
    template<class TMatrixType1, class TInputMatrix>
    static inline void Jacobi(
        TMatrixType1& J1,
        const TInputMatrix& rInputMatrix,
        const SizeType Size1,
        const SizeType Index1,
        const SizeType Index2
        )
    {
        // We define the 2x2 matrix
        BoundedMatrix<TDataType, 2, 2> b_matrix;
        b_matrix(0, 0) = rInputMatrix(Index1, Index1);
        b_matrix(0, 1) = 0.0;
        b_matrix(1, 0) = rInputMatrix(Index2, Index1);
        b_matrix(1, 1) = 0.0;

        // We compute the singular value decomposition for a 2x2 matrix
        BoundedMatrix<TDataType, 2, 2> u_matrix, s_matrix, v_matrix;
        SingularValueDecomposition2x2(b_matrix, u_matrix, s_matrix, v_matrix);

        // We compute the Jacobi rotation
        noalias(J1) = IdentityMatrix(Size1);
        J1(Index1, Index1) = u_matrix(0, 0);
        J1(Index1, Index2) = u_matrix(1, 0);
        J1(Index2, Index1) = u_matrix(0, 1);
        J1(Index2, Index2) = u_matrix(1, 1);
    }

    /**
     * @brief Computes the Singular Value Decomposition (SVD) of a given matrix using Jacobi rotations.
     * @details This function iteratively applies Jacobi rotations to transform the input matrix into a diagonal 
     * form, effectively computing its Singular Value Decomposition (SVD). It modifies the matrices U, S, and V such that:
     *      A = U * S * V^T
     * where:
     * - U is the left unitary (orthogonal for real matrices) matrix.
     * - S is a diagonal matrix containing the singular values.
     * - V is the right unitary matrix.
     *
     * The function continues iterating until the off-diagonal elements are below a specified tolerance 
     * or until a maximum number of iterations is reached.
     * @param rInputMatrix The input matrix for which SVD is to be computed.
     * @param rUMatrix The output matrix storing the left singular vectors (U).
     * @param rSMatrix The output matrix storing singular values in diagonal form (S).
     * @param rVMatrix The output matrix storing the right singular vectors (V).
     * @param rJ1 First Jacobi rotation matrix.
     * @param rJ2 Second Jacobi rotation matrix.
     * @param rAuxiliaryMatrixMN Auxiliary matrix for intermediate computations.
     * @param rAuxiliaryMatrixM Auxiliary matrix used for row transformations.
     * @param rAuxiliaryMatrixN Auxiliary matrix used for column transformations.
     * @param Tolerance The convergence tolerance for off-diagonal elements (default: ZeroTolerance).
     * @param MaxIter The maximum number of Jacobi iterations allowed (default: 200).
     * @tparam TInputMatrix Type of the input matrix.
     * @tparam TUmatrix Type of the left singular vector matrix (U).
     * @tparam TSMatrix Type of the diagonal matrix (S).
     * @tparam TVMatrix Type of the right singular vector matrix (V).
     * @tparam TJ1Matrix Type of the first Jacobi rotation matrix.
     * @tparam TJ2Matrix Type of the second Jacobi rotation matrix.
     * @tparam TAuxiliaryMatrixMN Type of the auxiliary matrix used for temporary calculations (MxN).
     * @tparam TAuxiliaryMatrixM Type of the auxiliary matrix for M-size operations.
     * @tparam TAuxiliaryMatrixN Type of the auxiliary matrix for N-size operations.
     * @return The number of iterations performed before convergence or reaching MaxIter.
     */
    template<class TInputMatrix, class TUmatrix, class TSMatrix, class TVMatrix,
            class TJ1Matrix, class TJ2Matrix, class TAuxiliaryMatrixMN,
            class TAuxiliaryMatrixM, class TAuxiliaryMatrixN>
    static inline IndexType ComputeJacobiIterations(
        const TInputMatrix& rInputMatrix,
        TUmatrix& rUMatrix,
        TSMatrix& rSMatrix,
        TVMatrix& rVMatrix,
        TJ1Matrix& rJ1,
        TJ2Matrix& rJ2,
        TAuxiliaryMatrixMN& rAuxiliaryMatrixMN,
        TAuxiliaryMatrixM& rAuxiliaryMatrixM,
        TAuxiliaryMatrixN& rAuxiliaryMatrixN,
        const TDataType Tolerance = ZeroTolerance,
        const IndexType MaxIter = 200
        )
    {
        // Get the dimensions of the input matrix
        const SizeType m = rInputMatrix.size1(); // Number of rows
        const SizeType n = rInputMatrix.size2(); // Number of columns

        // Compute the relative tolerance based on the norm of the input matrix
        const TDataType relative_tolerance = Tolerance * TwoNorm(rInputMatrix);

        // Iteration counter
        IndexType iter = 0;

        // Perform Jacobi SVD iterations until convergence (off-diagonal elements below threshold)
        while (JacobiNorm(rSMatrix) > relative_tolerance) {
            // Iterate over pairs of columns (i, j) for Jacobi transformations
            for (IndexType i = 0; i < n; i++) {
                for (IndexType j = i + 1; j < n; j++) {
                    // Compute Jacobi rotation matrices J1 and J2 for column pair (i, j)
                    Jacobi(rJ1, rJ2, rSMatrix, m, n, i, j);

                    // Apply Jacobi transformation to the singular values matrix S
                    noalias(rAuxiliaryMatrixMN) = prod(rSMatrix, rJ2);
                    noalias(rSMatrix) = prod(rJ1, rAuxiliaryMatrixMN);

                    // Update the left singular vector matrix U
                    noalias(rAuxiliaryMatrixM) = prod(rUMatrix, trans(rJ1));
                    noalias(rUMatrix) = rAuxiliaryMatrixM;

                    // Update the right singular vector matrix V
                    noalias(rAuxiliaryMatrixN) = prod(trans(rJ2), rVMatrix);
                    noalias(rVMatrix) = rAuxiliaryMatrixN;
                }

                // Iterate over rows to ensure full diagonalization (for rectangular matrices)
                for (IndexType j = n; j < m; j++) {
                    // Compute Jacobi rotation matrix J1 for row pair (i, j)
                    Jacobi(rJ1, rSMatrix, m, i, j);

                    // Apply transformation to the singular values matrix S
                    noalias(rAuxiliaryMatrixMN) = prod(rJ1, rSMatrix);
                    noalias(rSMatrix) = rAuxiliaryMatrixMN;

                    // Update the left singular vector matrix U
                    noalias(rAuxiliaryMatrixM) = prod(rUMatrix, trans(rJ1));
                    noalias(rUMatrix) = rAuxiliaryMatrixM;
                }
            }

            // Increment iteration counter
            ++iter;

            // Check if the maximum iteration limit is reached
            if (iter > MaxIter) {
                KRATOS_WARNING("JacobiSingularValueDecomposition")
                    << "Maximum number of iterations (" << MaxIter << ") reached without convergence."
                    << std::endl;
                break;
            }
        }

        return iter; // Return the total number of iterations performed
    }

    /**
     * @brief Computes the two-norm (Frobenious norm) of a given matrix.
     * @details This function calculates the two-norm of the input matrix @p rA by summing the squares of all its elements
     * and then taking the square root of the sum.
     * @param rA The input matrix whose two-norm is to be computed.
     * @return The two-norm of the input matrix.
     * @tparam TMatrixType Type of the input matrix.
     */
    template<class TMatrixType>
    static TDataType TwoNorm(const TMatrixType& rA)
    {
        const TDataType sum = IndexPartition<IndexType>(rA.size1()).for_each<SumReduction<TDataType>>([&](IndexType i) {
            TDataType aux_sum = TDataType();
            for (IndexType j = 0; j < rA.size2(); ++j) {
                aux_sum += std::pow(rA(i,j),2);
            }
            return aux_sum;
        });
        return std::sqrt(sum);
    }

    /**
     * @brief Computes the Jacobi norm of a given matrix.
     * @details This function calculates the Jacobi norm of the input matrix @p rA. The Jacobi norm is defined as the sum of the absolute values of the off-diagonal elements of the matrix.
     * @tparam TMatrixType Type of the input matrix.
     * @tparam TDataType Type of the data stored in the matrix.
     * @param rA The input matrix for which the Jacobi norm is to be computed.
     * @return The Jacobi norm of the input matrix.
     * @tparam TMatrixType Type of the input matrix.
     */
    template<class TMatrixType>
    static TDataType JacobiNorm(const TMatrixType& rA)
    {
        return IndexPartition<IndexType>(rA.size1()).for_each<SumReduction<TDataType>>([&](IndexType i) {
            TDataType aux_sum = TDataType();
            for (IndexType j = 0; j < rA.size2(); ++j) {
                if (i != j) {
                    aux_sum += std::abs(rA(i,j));
                }
            }
            return aux_sum;
        });
    }

    ///@}
}; /* Class SVDUtils */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/
