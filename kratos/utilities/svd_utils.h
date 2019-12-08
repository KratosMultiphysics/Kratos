//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SVD_UTILS )
#define  KRATOS_SVD_UTILS


/* System includes */


/* External includes */

/* Project includes */
#include "utilities/math_utils.h"
#include "spaces/ublas_space.h"


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
 * @details Defines several utility functions
 * @author Vicente Mataix Ferrandiz
 */
template<class TDataType>
class SVDUtils
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the matrix type
    typedef Matrix MatrixType;

    /// Definition of the vector type
    typedef Vector VectorType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of index type
    typedef std::size_t IndexType;

    /// Definition of local space
    typedef UblasSpace<TDataType, Matrix, Vector> LocalSpaceType;

    /// Definition of epsilon zero tolerance
    constexpr static TDataType ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /* Constructor */


    /** Destructor */

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function gives the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite innefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param InputMatrix The matrix where perform the SVD
     * @param UMatrix The unitary U matrix
     * @param SMatrix The diagonal S matrix
     * @param VMatrix The unitary V matrix
     * @param ThisParameters The configuration parameters
     * @return iter: The number of iterations
     */
    static inline std::size_t SingularValueDecomposition(
        const MatrixType& InputMatrix,
        MatrixType& UMatrix,
        MatrixType& SMatrix,
        MatrixType& VMatrix,
        Parameters ThisParameters)
    {
        // Validating defaults
        Parameters default_parameters = Parameters(R"(
        {
            "type_svd"             : "Jacobi",
            "tolerance"            : 0.0,
            "max_iter"             : 200
        })");
        default_parameters["tolerance"].SetDouble(std::numeric_limits<double>::epsilon());
        ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        const std::string& r_type_svd = ThisParameters["type_svd"].GetString();
        const double tolerance = ThisParameters["tolerance"].GetDouble();
        const double max_iter = ThisParameters["max_iter"].GetInt();
        return SingularValueDecomposition(InputMatrix, UMatrix, SMatrix, VMatrix, r_type_svd, tolerance, max_iter);
    }

    /**
     * @brief This function gives the SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite innefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param InputMatrix The matrix where perform the SVD
     * @param UMatrix The unitary U matrix
     * @param SMatrix The diagonal S matrix
     * @param VMatrix The unitary V matrix
     * @param TypeSVD The type of SVD algorithm (Jacobi by default)
     * @param Tolerance The tolerance considered
     * @param MaxIter Maximum number of iterations
     * @return iter: The number of iterations
     */
    static inline std::size_t SingularValueDecomposition(
        const MatrixType& InputMatrix,
        MatrixType& UMatrix,
        MatrixType& SMatrix,
        MatrixType& VMatrix,
        const std::string& TypeSVD = "Jacobi",
        const TDataType Tolerance = std::numeric_limits<double>::epsilon(),
        const IndexType MaxIter = 200)
    {
        if (TypeSVD == "Jacobi") {
            return JacobiSingularValueDecomposition(InputMatrix, UMatrix, SMatrix, VMatrix, Tolerance, MaxIter);
        } else {
            KRATOS_ERROR << "SVD Type not implemented" << std::endl;
        }
    }

    /**
     * @brief This function gives the Jacobi SVD of a given mxn matrix (m>=n), returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1 (which means that the biggest number is the first one and the smallest the last one)
     * @todo This version is quite innefficient, look for a real and mathematical implementation (not the algorithm found in Wikipedia!!)
     * @param InputMatrix The matrix where perform the SVD
     * @param UMatrix The unitary U matrix
     * @param SMatrix The diagonal S matrix
     * @param VMatrix The unitary V matrix
     * @param Tolerance The tolerance considered
     * @param MaxIter Maximum number of iterations
     * @return iter: The number of iterations
     */

    static inline std::size_t JacobiSingularValueDecomposition(
        const MatrixType& InputMatrix,
        MatrixType& UMatrix,
        MatrixType& SMatrix,
        MatrixType& VMatrix,
        const TDataType Tolerance = std::numeric_limits<double>::epsilon(),
        const IndexType MaxIter = 200)
    {
        const SizeType m = InputMatrix.size1();
        const SizeType n = InputMatrix.size2();

        if(SMatrix.size1() != m || SMatrix.size2() != n) {
            SMatrix.resize(m, n, false);
        }
        noalias(SMatrix) = InputMatrix;

        if(UMatrix.size1() != m || UMatrix.size2() != m) {
            UMatrix.resize(m, m, false);
        }
        noalias(UMatrix) = IdentityMatrix(m);

        if(VMatrix.size1() != n || VMatrix.size2() != n) {
            VMatrix.resize(n, n, false);
        }
        noalias(VMatrix) = IdentityMatrix(n);

        const TDataType relative_tolerance = Tolerance * LocalSpaceType::TwoNorm(InputMatrix);

        IndexType iter = 0;

        // We create the auxuliar matrices (for aliased operations)
        MatrixType auxiliar_matrix_mn(m, n);
        MatrixType auxiliar_matrix_m(m, m);
        MatrixType auxiliar_matrix_n(n, n);

        // We compute Jacobi
        while (LocalSpaceType::JacobiNorm(SMatrix) > relative_tolerance) {
            for (IndexType i = 0; i < n; i++) {
                for (IndexType j = i+1; j < n; j++) {
                    MatrixType j1(m, m);
                    MatrixType j2(n, n);

                    Jacobi(j1, j2, SMatrix, m, n, i, j);

                    const MatrixType aux_matrix = prod(SMatrix, j2);
                    noalias(SMatrix) = prod(j1, aux_matrix);
                    noalias(auxiliar_matrix_m) = prod(UMatrix, trans(j1));
                    noalias(UMatrix) = auxiliar_matrix_m;
                    noalias(auxiliar_matrix_n) = prod(trans(j2), VMatrix);
                    noalias(VMatrix) = auxiliar_matrix_n;
                }

                for (IndexType j = n; j < m; j++) {
                    MatrixType j1(m, m);

                    Jacobi(j1, SMatrix, m, i, j);

                    noalias(auxiliar_matrix_mn) = prod(j1, SMatrix);
                    noalias(SMatrix) = auxiliar_matrix_mn;
                    noalias(auxiliar_matrix_m) = prod(UMatrix, trans(j1));
                    noalias(UMatrix) = auxiliar_matrix_m;
                }
            }

            ++iter;
            if (iter > MaxIter) {
                KRATOS_WARNING("JacobiSingularValueDecomposition") << "Maximum number of iterations " << MaxIter << " reached." << std::endl;
                break;
            }
        }

        return iter;
    }

    /**
     * @brief This function gives the Jacobi SVD of a given 2x2 matrix, returns U,S; where A=U*S*V
     * @details U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1
     * @param InputMatrix The matrix where perform the SVD
     * @param UMatrix The unitary U matrix
     * @param SMatrix The diagonal S matrix
     * @param VMatrix The unitary V matrix
     */
    static inline void SingularValueDecomposition2x2(
        const MatrixType& InputMatrix,
        MatrixType& UMatrix,
        MatrixType& SMatrix,
        MatrixType& VMatrix
        )
    {
        const TDataType t = (InputMatrix(0, 1) - InputMatrix(1, 0))/(InputMatrix(0, 0) + InputMatrix(1, 1));
        const TDataType c = 1.0/std::sqrt(1.0 + t*t);
        const TDataType s = t*c;
        MatrixType r_matrix(2, 2);
        r_matrix(0, 0) =  c;
        r_matrix(0, 1) = -s;
        r_matrix(1, 0) =  s;
        r_matrix(1, 1) =  c;

        MatrixType m_matrix = prod(r_matrix, InputMatrix);

        SingularValueDecomposition2x2Symmetric(m_matrix, UMatrix, SMatrix, VMatrix);

        MatrixType auxiliar_matrix_m(UMatrix.size1(), UMatrix.size2());
        noalias(auxiliar_matrix_m) = prod(trans(r_matrix), UMatrix);
        noalias(UMatrix) = auxiliar_matrix_m;
    }

	/**
     * This function gives the Jacobi SVD of a given 2x2 matrix, returns U,S; where A=U*S*V
     * U and V are unitary, and S is a diagonal matrix.
     * Where s_i >= 0, and s_i >= s_i+1
     * @param InputMatrix The matrix where perform the SVD
     * @param UMatrix The unitary U matrix
     * @param SMatrix The diagonal S matrix
     * @param VMatrix The unitary V matrix
     */
    static inline void SingularValueDecomposition2x2Symmetric(
        const MatrixType& InputMatrix,
        MatrixType& UMatrix,
        MatrixType& SMatrix,
        MatrixType& VMatrix
        )
    {
        if(SMatrix.size1() != 2 || SMatrix.size2() != 2) {
            SMatrix.resize(2 ,2, false);
        }
        if(UMatrix.size1() != 2 || UMatrix.size2() != 2) {
            UMatrix.resize(2, 2, false);
        }
        if(VMatrix.size1() != 2 || VMatrix.size2() != 2) {
            VMatrix.resize(2, 2, false);
        }

        if (std::abs(InputMatrix(1, 0)) < ZeroTolerance) { // Already symmetric
            noalias(SMatrix) = InputMatrix;
            noalias(UMatrix) = IdentityMatrix(2);
            noalias(VMatrix) = UMatrix;
        } else {
            const TDataType w = InputMatrix(0, 0);
            const TDataType y = InputMatrix(1, 0);
            const TDataType z = InputMatrix(1, 1);
            const TDataType ro = (z - w)/(2.0 * y);
            const TDataType t = MathUtils<TDataType>::Sign(ro)/(std::abs(ro) + std::sqrt(1 + ro * ro));
            const TDataType c = 1.0/(std::sqrt(1.0 + t*t));
            const TDataType s = t*c;

            UMatrix(0, 0) =  c;
            UMatrix(0, 1) =  s;
            UMatrix(1, 0) = -s;
            UMatrix(1, 1) =  c;
            noalias(VMatrix) = trans(UMatrix);

            noalias(SMatrix) = prod(trans(UMatrix), MatrixType(prod(InputMatrix, trans(VMatrix))));
        }

        MatrixType z_matrix(2, 2);
        z_matrix(0, 0) = MathUtils<TDataType>::Sign(SMatrix(0, 0));
        z_matrix(0, 1) = 0.0;
        z_matrix(1, 0) = 0.0;
        z_matrix(1, 1) = MathUtils<TDataType>::Sign(SMatrix(1, 1));

        // Auxiliar matrix for alias operations
        MatrixType aux_2_2_matrix(2, 2);
        noalias(aux_2_2_matrix) = prod(UMatrix, z_matrix);
        noalias(UMatrix) = aux_2_2_matrix;
        noalias(aux_2_2_matrix) = prod(z_matrix, SMatrix);
        noalias(SMatrix) = aux_2_2_matrix;

        if (SMatrix(0, 0) < SMatrix(1, 1)) {
            MatrixType p_matrix(2, 2);
            p_matrix(0, 0) = 0.0;
            p_matrix(0, 1) = 1.0;
            p_matrix(1, 0) = 1.0;
            p_matrix(1, 1) = 0.0;

            noalias(aux_2_2_matrix) = prod(UMatrix, p_matrix);
            noalias(UMatrix) = aux_2_2_matrix;
            const MatrixType aux_matrix = prod(SMatrix, p_matrix);
            noalias(SMatrix) = prod(p_matrix, aux_matrix);
            noalias(aux_2_2_matrix) = prod(p_matrix, VMatrix);
            noalias(VMatrix) = aux_2_2_matrix;
        }
    }

    /**
     * @brief This method computes the Jacobi rotation operation
     * @param J1 First Jacobi matrix
     * @param J2 Second Jacobi matrix
     * @param InputMatrix The matrix to compute the Jacobi tolerance
     * @param Size1 The size of the matrix (number of rows)
     * @param Size2 The size of the matrix (number of columns)
     * @param Index1 The index to compute (row)
     * @param Index2 The index to compute (column)
     */
    static inline void Jacobi(
        MatrixType& J1,
        MatrixType& J2,
        const MatrixType& InputMatrix,
        const SizeType& Size1,
        const SizeType& Size2,
        const SizeType& Index1,
        const SizeType& Index2
        )
    {
        MatrixType b_matrix(2,2);
        b_matrix(0, 0) = InputMatrix(Index1, Index1);
        b_matrix(0, 1) = InputMatrix(Index1, Index2);
        b_matrix(1, 0) = InputMatrix(Index2, Index1);
        b_matrix(1, 1) = InputMatrix(Index2, Index2);

        MatrixType u_matrix, s_matrix, v_matrix;

        SingularValueDecomposition2x2(b_matrix, u_matrix, s_matrix, v_matrix);

        J1 = IdentityMatrix(Size1);
        J1(Index1, Index1) = u_matrix(0, 0);
        J1(Index1, Index2) = u_matrix(1, 0);
        J1(Index2, Index1) = u_matrix(0, 1);
        J1(Index2, Index2) = u_matrix(1, 1);

        J2 = IdentityMatrix(Size2);
        J2(Index1, Index1) = v_matrix(0, 0);
        J2(Index1, Index2) = v_matrix(1, 0);
        J2(Index2, Index1) = v_matrix(0, 1);
        J2(Index2, Index2) = v_matrix(1, 1);
    }

    /**
     * @brief This method computes the Jacobi rotation operation
     * @param J1 First Jacobi matrix
     * @param InputMatrix The matrix to compute the Jacobi tolerance
     * @param Size1 The size of the matrix (number of rows)
     * @param Size2 The size of the matrix (number of columns)
     * @param Index1 The index to compute (row)
     * @param Index2 The index to compute (column)
     */
    static inline void Jacobi(
        MatrixType& J1,
        const MatrixType& InputMatrix,
        const SizeType& Size1,
        const SizeType& Index1,
        const SizeType& Index2
        )
    {
        MatrixType b_matrix(2,2);
        b_matrix(0, 0) = InputMatrix(Index1, Index1);
        b_matrix(0, 1) = 0.0;
        b_matrix(1, 0) = InputMatrix(Index2, Index1);
        b_matrix(1, 1) = 0.0;

        MatrixType u_matrix, s_matrix, v_matrix;

        SingularValueDecomposition2x2(b_matrix, u_matrix, s_matrix, v_matrix);

        noalias(J1) = IdentityMatrix(Size1);
        J1(Index1, Index1) = u_matrix(0, 0);
        J1(Index1, Index2) = u_matrix(1, 0);
        J1(Index2, Index1) = u_matrix(0, 1);
        J1(Index2, Index2) = u_matrix(1, 1);
    }

    /**
     * @brief This method computes the condition number using the SVD
     * @details The condition number can be estimated as the ratio between the largest singular value and the smallest singular value
     * @param InputMatrix The matrix to be evaluated
     * @param Tolerance The tolerance considered
     * @return condition_number: The ratio between the largest SV and the smallest SV
     */
    static inline TDataType SVDConditionNumber(
        const MatrixType& InputMatrix,
        const std::string TypeSVD = "Jacobi",
        const TDataType Tolerance = std::numeric_limits<double>::epsilon(),
        const IndexType MaxIter = 200)
    {
        MatrixType u_matrix, s_matrix, v_matrix;
        SingularValueDecomposition(InputMatrix, u_matrix, s_matrix, v_matrix, TypeSVD, Tolerance, MaxIter);

        const SizeType size_s = s_matrix.size1();
        const TDataType condition_number = s_matrix(0, 0)/s_matrix(size_s - 1, size_s - 1);

        return condition_number;
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

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

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

    SVDUtils(void);

    SVDUtils(SVDUtils& rSource);

}; /* Class SVDUtils */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_SVD_UTILS  defined */
