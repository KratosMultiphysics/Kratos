//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

/* Trilinos includes */
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_FEVector.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_MatrixMatrix.h>

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "custom_utilities/trilinos_dof_updater.h"

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
 * @class TrilinosSpace
 * @ingroup TrilinosApplication
 * @brief The space adapted for Trilinos vectors and matrices
 * @author Riccardo Rossi
 * @tparam TMatrixType The matrix type considered
 * @tparam TVectorType the vector type considered
 */
template<class TMatrixType, class TVectorType>
class TrilinosSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosSpace
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosSpace);

    /// Definition of the data type
    using DataType = double;

    /// Definition of the matrix type
    using MatrixType = TMatrixType;

    /// Definition of the vector type
    using VectorType = TVectorType;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the pointer types
    using MatrixPointerType = typename Kratos::shared_ptr<TMatrixType>;
    using VectorPointerType = typename Kratos::shared_ptr<TVectorType>;

    /// Some other definitions
    using DofUpdaterType = TrilinosDofUpdater< TrilinosSpace<TMatrixType,TVectorType>>;
    using DofUpdaterPointerType = typename DofUpdater<TrilinosSpace<TMatrixType,TVectorType>>::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosSpace()
    {
    }

    /// Destructor.
    virtual ~TrilinosSpace()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an empty pointer to a matrix
     * @return The pointer to the matrix
     */
    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(nullptr);
    }

    /**
     * @brief This method creates an empty pointer to a vector
     * @return The pointer to the vector
     */
    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(nullptr);
    }

    /**
     * @brief This method creates an empty pointer to a matrix using epetra communicator
     * @param rComm The epetra communicator
     * @return The pointer to the matrix
     */
    static MatrixPointerType CreateEmptyMatrixPointer(Epetra_MpiComm& rComm)
    {
        const int global_elems = 0;
        Epetra_Map Map(global_elems, 0, rComm);
        return MatrixPointerType(new TMatrixType(::Copy, Map, 0));
    }

    /**
     * @brief This method creates an empty pointer to a vector using epetra communicator
     * @param rComm The epetra communicator
     * @return The pointer to the vector
     */
    static VectorPointerType CreateEmptyVectorPointer(Epetra_MpiComm& rComm)
    {
        const int global_elems = 0;
        Epetra_Map Map(global_elems, 0, rComm);
        return VectorPointerType(new TVectorType(Map));
    }

    /**
     * @brief Returns size of vector rV
     * @param rV The vector considered
     * @return The size of the vector
     */
    static IndexType Size(const VectorType& rV)
    {
        const int size = rV.GlobalLength();
        return size;
    }

    /**
     * @brief Returns number of rows of rM
     * @param rM The matrix considered
     * @return The number of rows of rM
     */
    static IndexType Size1(MatrixType const& rM)
    {
        const int size1 = rM.NumGlobalRows();
        return size1;
    }

    /**
     * @brief Returns number of columns of rM
     * @param rM The matrix considered
     * @return The number of columns of rM
     */
    static IndexType Size2(MatrixType const& rM)
    {
        const int size1 = rM.NumGlobalCols();
        return size1;
    }

    /**
     * @brief Returns the column of the matrix in the given position
     * @details rXi = rMij
     * @param j The position of the column
     * @param rM The matrix considered
     * @param rX The column considered
     * @todo Implement this method
     */
    static void GetColumn(
        const unsigned int j,
        const MatrixType& rM,
        VectorType& rX
        )
    {
        KRATOS_ERROR << "GetColumn method is not currently implemented" << std::endl;
    }

    /**
     * @brief Returns a copy of the matrix rX
     * @details rY = rX
     * @param rX The matrix considered
     * @param rY The copy of the matrix rX
     */
    static void Copy(
        const MatrixType& rX,
        MatrixType& rY
        )
    {
        rY = rX;
    }

    /**
     * @brief Returns a copy of the vector rX
     * @details rY = rX
     * @param rX The vector considered
     * @param rY The copy of the vector rX
     */
    static void Copy(
        const VectorType& rX,
        VectorType& rY
        )
    {
        rY = rX;
    }

    /**
     * @brief Returns the product of two vectors
     * @details rX * rY
     * @param rX The first vector considered
     * @param rY The second vector considered
     */
    static double Dot(
        const VectorType& rX,
        const VectorType& rY
        )
    {
        double value;
        const int sucess = rY.Dot(rX, &value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing dot product" <<  std::endl;
        return value;
    }

    /**
     * @brief Returns the maximum value of the vector rX
     * @param rX The vector considered
     * @return The maximum value of the vector rX
     */
    static double Max(const VectorType& rX)
    {
        double value;
        const int sucess = rX.MaxValue(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing maximum value" <<  std::endl;
        return value;
    }

    /**
     * @brief Returns the minimum value of the vector rX
     * @param rX The vector considered
     * @return The minimum value of the vector rX
     */
    static double Min(const VectorType& rX)
    {
        double value;
        const int sucess = rX.MinValue(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing minimum value" <<  std::endl;
        return value;
    }

    /**
     * @brief Returns the norm of the vector rX
     * @details ||rX||2
     * @param rX The vector considered
     * @return The norm of the vector rX
     */
    static double TwoNorm(const VectorType& rX)
    {
        double value;
        const int sucess = rX.Norm2(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing norm" <<  std::endl;
        return value;
    }

    /**
     * @brief Returns the Frobenius norm of the matrix rX
     * @details ||rA||2
     * @param rA The matrix considered
     * @return The Frobenius norm of the matrix rX
     */
    static double TwoNorm(const MatrixType& rA)
    {
        return rA.NormFrobenius();
    }

    /**
     * @brief Returns the multiplication of a matrix by a vector
     * @details y = A*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void Mult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        constexpr bool transpose_flag = false;
        const int ierr = rA.Multiply(transpose_flag, rX, rY);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra multiplication failure " << ierr << std::endl;
    }

    /**
     * @brief Returns the multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     */
    static void Mult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        constexpr bool transpose_flag = false;
        const int ierr = EpetraExt::MatrixMatrix::Multiply(rA, transpose_flag, rB, transpose_flag, rC, CallFillCompleteOnResult, KeepAllHardZeros);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra multiplication failure. This may result if A or B are not already Filled, or if errors occur in putting values into C, etc. " << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the transpose multiplication of a matrix by a vector
     * @details y = AT*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void TransposeMult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        constexpr bool transpose_flag = true;
        const int ierr = rA.Multiply(transpose_flag, rX, rY);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra multiplication failure " << ierr << std::endl;
    }

    /**
     * @brief Returns the transpose multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param TransposeFlag Flags to transpose the matrices
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     */
    static void TransposeMult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const std::pair<bool, bool> TransposeFlag = {false, false},
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        const int ierr = EpetraExt::MatrixMatrix::Multiply(rA, TransposeFlag.first, rB, TransposeFlag.second, rC, CallFillCompleteOnResult, KeepAllHardZeros);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra multiplication failure. This may result if A or B are not already Filled, or if errors occur in putting values into C, etc. " << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the product operation B'DB
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     * @tparam TEnforceInitialGraph If the initial graph is enforced, or a new one is generated
     */
    static void BtDBProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Define first auxiliary matrix
        std::vector<int> NumNz;
        MatrixType aux_1(::Copy, rA.RowMap(), NumNz.data());

        // First multiplication
        TransposeMult(rB, rD, aux_1, {true, false}, CallFillCompleteOnResult, KeepAllHardZeros);

        // If we enforce the initial connectivity
        if (EnforceInitialGraph) {
            // Second multiplication
            MatrixType aux_2(::Copy, rA.RowMap(), NumNz.data());
            Mult(aux_1, rB, aux_2, CallFillCompleteOnResult, KeepAllHardZeros);

            // Create an Epetra_Matrix
            MatrixType* aux_3 =  new MatrixType(::Copy, CombineMatricesGraphs(rA, aux_2));

            // Copy values
            CopyMatrixValues(*aux_3, aux_2);

            // Doing a swap
            std::swap(rA, *aux_3);

            // Delete the new matrix
            delete aux_3;
        } else { // A new matrix
            // Already existing matrix
            if (rA.NumGlobalNonzeros() > 0) {
                // Create an Epetra_Matrix
                MatrixType* aux_2 =  new MatrixType(::Copy, rB.RowMap(), NumNz.data());

                // Second multiplication
                Mult(aux_1, rB, *aux_2, CallFillCompleteOnResult, KeepAllHardZeros);

                // Doing a swap
                std::swap(rA, *aux_2);

                // Delete the new matrix
                delete aux_2;
            } else { // Empty matrix
                // Second multiplication
                Mult(aux_1, rB, rA, CallFillCompleteOnResult, KeepAllHardZeros);
            }
        }
    }

    /**
     * @brief Calculates the product operation BDB'
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     * @tparam TEnforceInitialGraph If the initial graph is enforced, or a new one is generated
     */
    static void BDBtProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Define first auxiliary matrix
        std::vector<int> NumNz;
        MatrixType aux_1(::Copy, rA.RowMap(), NumNz.data());

        // First multiplication
        Mult(rB, rD, aux_1, CallFillCompleteOnResult, KeepAllHardZeros);

        // If we enforce the initial connectivity
        if (EnforceInitialGraph) {
            // Second multiplication
            MatrixType aux_2(::Copy, rA.RowMap(), NumNz.data());
            TransposeMult(aux_1, rB, aux_2, {false, true}, CallFillCompleteOnResult, KeepAllHardZeros);

            // Create an Epetra_Matrix
            MatrixType* aux_3 =  new MatrixType(::Copy, CombineMatricesGraphs(rA, aux_2));

            // Copy values
            CopyMatrixValues(*aux_3, aux_2);

            // Doing a swap
            std::swap(rA, *aux_3);

            // Delete the new matrix
            delete aux_3;
        } else { // A new matrix
            // Already existing matrix
            if (rA.NumGlobalNonzeros() > 0) {
                // Create an Epetra_Matrix
                MatrixType* aux_2 =  new MatrixType(::Copy, rA.RowMap(), NumNz.data());

                // Second multiplication
                TransposeMult(aux_1, rB, *aux_2, {false, true}, CallFillCompleteOnResult, KeepAllHardZeros);

                // Doing a swap
                std::swap(rA, *aux_2);

                // Delete the new matrix
                delete aux_2;
            } else { // Empty matrix
                // Second multiplication
                TransposeMult(aux_1, rB, rA, {false, true}, CallFillCompleteOnResult, KeepAllHardZeros);
            }
        }
    }

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details y = A*x
     * Checks if a multiplication is needed and tries to do otherwise
     * @param rX The vector considered
     * @param A The scalar considered
     */
    static void InplaceMult(
        VectorType& rX,
        const double A
        )
    {
        if (A != 1.00) {
            const int ierr = rX.Scale(A);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra scaling failure " << ierr << std::endl;
        }
    }

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details x = A*y
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    static void Assign(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        if (A != 1.00) {
            const int ierr = rX.Scale(A, rY); //not sure
            KRATOS_ERROR_IF(ierr != 0) << "Epetra assign failure " << ierr << std::endl;
        } else {
            rX = rY;
        }
    }

    /**
     * @brief Returns the unaliased addition of a vector by a scalar times a vector
     * @details X += A*y;
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    static void UnaliasedAdd(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        const int ierr = rX.Update(A, rY, 1.0);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra unaliased add failure " << ierr << std::endl;
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rZ = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The second vector considered
     * @param rZ The resulting vector considered
     */
    static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        const VectorType& rY,
        VectorType& rZ
        )
    {
        const int ierr = rZ.Update(A, rX, B, rY, 0.0);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra scale and add failure " << ierr << std::endl;
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rY = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The resulting vector considered
     */
    static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        VectorType& rY
        )
    {
        const int ierr = rY.Update(A, rX, B);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra scale and add failure " << ierr << std::endl;
    }

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param value The value considered
     */
    static void SetValue(
        VectorType& rX,
        IndexType i,
        const double value
        )
    {
        Epetra_IntSerialDenseVector indices(1);
        Epetra_SerialDenseVector values(1);
        indices[0] = i;
        values[0] = value;

        int ierr = rX.ReplaceGlobalValues(indices, values);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;

        ierr = rX.GlobalAssemble(Insert,true); //Epetra_CombineMode mode=Add);
        KRATOS_ERROR_IF(ierr < 0) << "Epetra failure when attempting to insert value in function SetValue" << std::endl;
    }

    /**
     * @brief assigns a scalar to a vector
     * @details rX = A
     * @param rX The vector considered
     * @param A The scalar considered
     */
    static void Set(
        VectorType& rX,
        const DataType A
        )
    {
        const int ierr = rX.PutScalar(A);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra set failure " << ierr << std::endl;
    }

    /**
     * @brief Resizes a matrix
     * @param rA The matrix to be resized
     * @param m The new number of rows
     * @param n The new number of columns
     */
    static void Resize(
        MatrixType& rA,
        const SizeType m,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for Trilinos Sparse Matrix" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param rX The vector to be resized
     * @param n The new size
     */
    static void Resize(
        VectorType& rX,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param pA The pointer to the vector to be resized
     * @param n The new size
    */
    static void Resize(
        VectorPointerType& pX,
        const SizeType n
        )
    {
        //KRATOS_ERROR_IF(pX != NULL) << "Trying to resize a null pointer" << std::endl;
        int global_elems = n;
        Epetra_Map Map(global_elems, 0, pX->Comm());
        VectorPointerType pNewEmptyX = Kratos::make_shared<VectorType>(Map);
        pX.swap(pNewEmptyX);
    }

    /**
     * @brief Clears a matrix
     * @param pA The pointer to the matrix to be cleared
     */
    static void Clear(MatrixPointerType& pA)
    {
        if(pA != NULL) {
            int global_elems = 0;
            Epetra_Map Map(global_elems, 0, pA->Comm());
            MatrixPointerType pNewEmptyA = MatrixPointerType(new TMatrixType(::Copy, Map, 0));
            pA.swap(pNewEmptyA);
        }
    }

    /**
     * @brief Clears a vector
     * @param pX The pointer to the vector to be cleared
     */
    static void Clear(VectorPointerType& pX)
    {
        if(pX != NULL) {
            int global_elems = 0;
            Epetra_Map Map(global_elems, 0, pX->Comm());
            VectorPointerType pNewEmptyX = VectorPointerType(new VectorType(Map));
            pX.swap(pNewEmptyX);
        }
    }

    /**
     * @brief Sets a matrix to zero
     * @param rX The matrix to be set
     */
    inline static void SetToZero(MatrixType& rA)
    {
        const int ierr = rA.PutScalar(0.0);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra set to zero failure " << ierr << std::endl;
    }

    /**
     * @brief Sets a vector to zero
     * @param rX The vector to be set
     */
    inline static void SetToZero(VectorType& rX)
    {
        const int ierr = rX.PutScalar(0.0);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra set to zero failure " << ierr << std::endl;
    }

    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherMatrixType, class TEquationIdVectorType>

    /**
     * @brief Assembles the LHS of the system
     * @param rA The LHS matrix
     * @param rLHSContribution The contribution to the LHS
     * @param rEquationId The equation ids
     */
    inline static void AssembleLHS(
        MatrixType& rA,
        const Matrix& rLHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const unsigned int system_size = Size1(rA);

        // Count active indices
        unsigned int active_indices = 0;
        for (unsigned int i = 0; i < rEquationId.size(); i++)
            if (rEquationId[i] < system_size)
                ++active_indices;

        if (active_indices > 0) {
            // Size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseMatrix values(active_indices, active_indices);

            // Fill epetra vectors
            unsigned int loc_i = 0;
            for (unsigned int i = 0; i < rEquationId.size(); i++) {
                if (rEquationId[i] < system_size) {
                    indices[loc_i] = rEquationId[i];

                    unsigned int loc_j = 0;
                    for (unsigned int j = 0; j < rEquationId.size(); j++) {
                        if (rEquationId[j] < system_size) {
                            values(loc_i, loc_j) = rLHSContribution(i, j);
                            ++loc_j;
                        }
                    }
                    ++loc_i;
                }
            }

            const int ierr = rA.SumIntoGlobalValues(indices, values);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
        }
    }

    //***********************************************************************
    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherVectorType, class TEquationIdVectorType>

    /**
     * @brief Assembles the RHS of the system
     * @param rb The RHS vector
     * @param rRHSContribution The RHS contribution
     * @param rEquationId The equation ids
     */
    inline static void AssembleRHS(
        VectorType& rb,
        const Vector& rRHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const unsigned int system_size = Size(rb);

        // Count active indices
        unsigned int active_indices = 0;
        for (unsigned int i = 0; i < rEquationId.size(); i++)
            if (rEquationId[i] < system_size)
                ++active_indices;

        if (active_indices > 0) {
            // Size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseVector values(active_indices);

            // Fill epetra vectors
            unsigned int loc_i = 0;
            for (unsigned int i = 0; i < rEquationId.size(); i++) {
                if (rEquationId[i] < system_size) {
                    indices[loc_i] = rEquationId[i];
                    values[loc_i] = rRHSContribution[i];
                    ++loc_i;
                }
            }

            const int ierr = rb.SumIntoGlobalValues(indices, values);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
        }
    }

    /**
     * @brief This function returns if we are in a distributed system
     * @return True if we are in a distributed system, false otherwise (always true in this case)
     */
    inline static constexpr bool IsDistributed()
    {
        return true;
    }

    /**
     * @brief This function returns a value from a given vector according to a given index
     * @param rX The vector from which values are to be gathered
     * @param I The index of the value to be gathered
     * @return The value of the vector corresponding to the index I
     */
    inline static double GetValue(
        const VectorType& rX,
        const std::size_t I
        )
    {
        // index must be local to this proc
        KRATOS_ERROR_IF_NOT(rX.Map().MyGID(static_cast<int>(I))) << " non-local id: " << I << "." << std::endl;
        // Epetra_MultiVector::operator[] is used here to get the pointer to
        // the zeroth (only) local vector.
        return rX[0][rX.Map().LID(static_cast<int>(I))];
    }

    /**
     * @brief This function gathers the values of a given vector according to a given index array
     * @param rX The vector from which values are to be gathered
     * @param IndexArray The array containing the indices of the values to be gathered
     * @param pValues The array containing the gathered values
     */
    static void GatherValues(
        const VectorType& rX,
        const std::vector<int>& IndexArray,
        double* pValues
        )
    {
        KRATOS_TRY
        double tot_size = IndexArray.size();

        //defining a map as needed
        Epetra_Map dof_update_map(-1, tot_size, &(*(IndexArray.begin())), 0, rX.Comm());

        //defining the importer class
        Epetra_Import importer(dof_update_map, rX.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp(importer.TargetMap());

        //importing in the new temp vector the values
        int ierr = temp.Import(rX, importer, Insert);
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

        temp.ExtractCopy(&pValues);

        rX.Comm().Barrier();
        KRATOS_CATCH("")
    }

    /**
     * @brief Read a matrix from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param rComm The MPI communicator
     * @return The matrix read from the file
     */
    MatrixPointerType ReadMatrixMarket(
        const std::string FileName,
        Epetra_MpiComm& rComm
        )
    {
        KRATOS_TRY

        Epetra_CrsMatrix* pp = nullptr;

        int error_code = EpetraExt::MatrixMarketFileToCrsMatrix(FileName.c_str(), rComm, pp);

        KRATOS_ERROR_IF(error_code != 0) << "Eerror thrown while reading Matrix Market file "<<FileName<< " error code is : " << error_code;

        rComm.Barrier();

        const Epetra_CrsGraph& rGraph = pp->Graph();
        MatrixPointerType paux = Kratos::make_shared<Epetra_FECrsMatrix>( ::Copy, rGraph, false );

        IndexType NumMyRows = rGraph.RowMap().NumMyElements();

        int* MyGlobalElements = new int[NumMyRows];
        rGraph.RowMap().MyGlobalElements(MyGlobalElements);

        for(IndexType i = 0; i < NumMyRows; ++i) {
//             std::cout << pA->Comm().MyPID() << " : I=" << i << std::endl;
            IndexType GlobalRow = MyGlobalElements[i];

            int NumEntries;
            std::size_t Length = pp->NumGlobalEntries(GlobalRow);  // length of Values and Indices

            double* Values = new double[Length];     // extracted values for this row
            int* Indices = new int[Length];          // extracted global column indices for the corresponding values

            error_code = pp->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, Values, Indices);

            KRATOS_ERROR_IF(error_code != 0) << "Error thrown in ExtractGlobalRowCopy : " << error_code;

            error_code = paux->ReplaceGlobalValues(GlobalRow, Length, Values, Indices);

            KRATOS_ERROR_IF(error_code != 0) << "Error thrown in ReplaceGlobalValues : " << error_code;

            delete [] Values;
            delete [] Indices;
        }

        paux->GlobalAssemble();

        delete [] MyGlobalElements;
        delete pp;

        return paux;
        KRATOS_CATCH("");
    }

    /**
     * @brief Read a vector from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param rComm The MPI communicator
     * @param N The size of the vector
     */
    VectorPointerType ReadMatrixMarketVector(
        const std::string& rFileName,
        Epetra_MpiComm& rComm,
        const int N
        )
    {
        KRATOS_TRY

        Epetra_Map my_map(N, 0, rComm);
        Epetra_Vector* pv = nullptr;

        int error_code = EpetraExt::MatrixMarketFileToVector(rFileName.c_str(), my_map, pv);

        KRATOS_ERROR_IF(error_code != 0) << "error thrown while reading Matrix Market Vector file " << rFileName << " error code is: " << error_code;

        rComm.Barrier();

        IndexType num_my_rows = my_map.NumMyElements();
        std::vector<int> gids(num_my_rows);
        my_map.MyGlobalElements(gids.data());

        std::vector<double> values(num_my_rows);
        pv->ExtractCopy(values.data());

        VectorPointerType final_vector = Kratos::make_shared<VectorType>(my_map);
        int ierr = final_vector->ReplaceGlobalValues(gids.size(),gids.data(), values.data());
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found with code ierr = " << ierr << std::endl;

        final_vector->GlobalAssemble();

        delete pv;
        return final_vector;
        KRATOS_CATCH("");
    }

    /**
     * @brief Generates a graph combining the graphs of two matrices
     * @param rA The first matrix
     * @param rB The second matrix
     */
    static Epetra_CrsGraph CombineMatricesGraphs(
        const MatrixType& rA,
        const MatrixType& rB
        )
    {
        // Row maps must be the same
        KRATOS_ERROR_IF_NOT(rA.RowMap().SameAs(rB.RowMap())) << "Row maps are not compatible" << std::endl;

        // Getting the graphs
        const auto& r_graph_a = rA.Graph();
        const auto& r_graph_b = rB.Graph();

        // Assuming local indexes
        KRATOS_ERROR_IF_NOT(r_graph_a.IndicesAreLocal() && r_graph_b.IndicesAreLocal()) << "Graphs indexes must be local" << std::endl;

        // Some definitions
        int i, j, ierr;
        int num_entries; // Number of non-zero entries
        int* cols;       // Column indices of row non-zero values
        const bool same_col_map = rA.ColMap().SameAs(rB.ColMap());
        Epetra_CrsGraph graph = same_col_map ? Epetra_CrsGraph(::Copy, rA.RowMap(), rA.ColMap(), 1000) : Epetra_CrsGraph(::Copy, rA.RowMap(), 1000);

        // Same column map. Local indices, simpler and faster
        if (same_col_map) {
            for (i = 0; i < r_graph_a.NumMyRows(); i++) {
                std::unordered_set<int> combined_indexes;
                // First graph
                ierr = r_graph_a.ExtractMyRowView(i, num_entries, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting indices (I) with code ierr = " << ierr << std::endl;
                for (j = 0; j < num_entries; j++) {
                    combined_indexes.insert(cols[j]);
                }
                // Second graph
                ierr = r_graph_b.ExtractMyRowView(i, num_entries, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting indices (II) with code ierr = " << ierr << std::endl;
                for (j = 0; j < num_entries; j++) {
                    combined_indexes.insert(cols[j]);
                }
                // Vector equivalent
                std::vector<int> combined_indexes_vector(combined_indexes.begin(), combined_indexes.end());
                num_entries = combined_indexes_vector.size();
                // Adding to graph
                ierr = graph.InsertMyIndices(i, num_entries, combined_indexes_vector.data());
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure inserting indices with code ierr = " << ierr << std::endl;
            }
        } else { // Different column map, global indices
            for (i = 0; i < r_graph_a.NumMyRows(); i++) {
                const int global_row_index = r_graph_a.GRID(i);
                // First graph
                ierr = r_graph_a.ExtractMyRowView(i, num_entries, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting indices (I) with code ierr = " << ierr << std::endl;
                std::vector<int> combined_indexes_vector;
                combined_indexes_vector.reserve(num_entries);
                for (j = 0; j < num_entries; j++) {
                    combined_indexes_vector.push_back(r_graph_a.GCID(cols[j]));
                }
                // Adding to graph
                ierr = graph.InsertGlobalIndices(global_row_index, num_entries, combined_indexes_vector.data());
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure inserting indices with code ierr = " << ierr << std::endl;
            }
            for (i = 0; i < r_graph_b.NumMyRows(); i++) {
                const int global_row_index = r_graph_b.GRID(i);
                // Second graph
                ierr = r_graph_b.ExtractMyRowView(i, num_entries, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting indices (II) with code ierr = " << ierr << std::endl;
                std::vector<int> combined_indexes_vector;
                combined_indexes_vector.reserve(num_entries);
                for (j = 0; j < num_entries; j++) {
                    combined_indexes_vector.push_back(r_graph_b.GCID(cols[j]));
                }
                // Adding to graph
                ierr = graph.InsertGlobalIndices(global_row_index, num_entries, combined_indexes_vector.data());
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure inserting indices with code ierr = " << ierr << std::endl;
            }
        }

        // Finalizing graph construction
        ierr = graph.FillComplete();
        KRATOS_ERROR_IF(ierr < 0) << ": Epetra failure in Epetra_CrsGraph.FillComplete. Error code: " << ierr << std::endl;

        return graph;
    }

    /**
     * @brief Copy values from one matrix to another
     * @details It is assumed that the sparcity of both matrices is compatible
     * @param rA The matrix where assigning values
     * @param rB The matrix to be copied
     */
    static void CopyMatrixValues(
        MatrixType& rA,
        const MatrixType& rB
        )
    {   
        // Cleaning destination matrix
        SetToZero(rA);

        // The current process id
        const int rank = rA.Comm().MyPID();

        // Row maps must be the same
        const bool same_col_map = rA.ColMap().SameAs(rB.ColMap());

        // Getting the graphs
        const auto& r_graph_b = rB.Graph();

        // Copy values from rB to intermediate
        int i, ierr;
        int num_entries; // Number of non-zero entries (rB matrix)
        double* vals;    // Row non-zero values (rB matrix)
        int* cols;       // Column indices of row non-zero values (rB matrix)
        if (same_col_map) {
            for (i = 0; i < rB.NumMyRows(); i++) {
                ierr = rB.ExtractMyRowView(i, num_entries, vals, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting values in local row " << i << " in rank " << rank << " with code ierr = " << ierr << std::endl;
                ierr = rA.ReplaceMyValues(i, num_entries, vals, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found replacing values in local row " << i << " in rank " << rank << " with code ierr = " << ierr << std::endl;
            }
        } else {
            for (i = 0; i < rB.NumMyRows(); i++) {
                ierr = rB.ExtractMyRowView(i, num_entries, vals, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting values in local row " << i << " in rank " << rank << " with code ierr = " << ierr << std::endl;
                const int global_row_index = r_graph_b.GRID(i);
                for (int j = 0; j < num_entries; j++) {
                    cols[j] = r_graph_b.GCID(cols[j]);
                }
                ierr = rA.ReplaceGlobalValues(global_row_index, num_entries, vals, cols);
                KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found extracting values in global row " << global_row_index << " in rank " << rank << " with code ierr = " << ierr << std::endl;
            }
        }
    }

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal);

        for (int i = 0; i < rA.NumMyRows(); i++) {
            int numEntries; // Number of non-zero entries
            double* vals;   // Row non-zero values
            int* cols;      // Column indices of row non-zero values
            rA.ExtractMyRowView(i, numEntries, vals, cols);
            const int row_gid = rA.RowMap().GID(i);
            bool empty = true;
            int j;
            for (j = 0; j < numEntries; j++) {
                const int col_gid = rA.ColMap().GID(cols[j]);
                // Check diagonal value
                if (col_gid == row_gid) {
                    if(std::abs(vals[j]) > zero_tolerance) {
                        empty = false;
                    }
                    break;
                }
            }

            // If diagonal empty assign scale factor
            if (empty) {
                vals[j] = scale_factor;
                rb[0][i] = 0.0;
            }
        }

        // Global assembly
        rb.GlobalAssemble();
        rA.GlobalAssemble();

        return scale_factor;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(Size1(rA));
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
            default:
                return GetMaxDiagonal(rA);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    static double GetDiagonalNorm(const MatrixType& rA)
    {
        KRATOS_TRY

        Epetra_Vector diagonal(rA.RowMap());
        const int ierr = rA.ExtractDiagonalCopy(diagonal);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure extracting diagonal " << ierr << std::endl;

        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::TwoNorm(diagonal);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    static double GetMaxDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        Epetra_Vector diagonal(rA.RowMap());
        const int ierr = rA.ExtractDiagonalCopy(diagonal);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure extracting diagonal " << ierr << std::endl;
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Max(diagonal);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    static double GetMinDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        Epetra_Vector diagonal(rA.RowMap());
        const int ierr = rA.ExtractDiagonalCopy(diagonal);
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure extracting diagonal " << ierr << std::endl;
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Min(diagonal);

        KRATOS_CATCH("");
    }

   /**
    * @brief Check if the TrilinosSpace is distributed.
    * @details This static member function checks whether the TrilinosSpace is distributed or not.
    * @return True if the space is distributed, false otherwise.
    */
    static constexpr bool IsDistributedSpace()
    {
        return true;
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

    /**
     * @brief Turn back information as a string.
     * @return Info as a string.
     */
    virtual std::string Info() const
    {
        return "TrilinosSpace";
    }

    /**
     * @brief Print information about this object.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TrilinosSpace";
    }

    /**
     * @brief Print object's data.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    /**
     * @brief Writes a matrix to a file in MatrixMarket format
     * @param pFileName The name of the file to be written
     * @param rM The matrix to be written
     * @param Symmetric If the matrix is symmetric
     * @return True if the file was successfully written, false otherwise
     */
    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(
        const char* pFileName,
        const TOtherMatrixType& rM,
        const bool Symmetric
        )
    {
        // the argument "Symmetric" does not have an effect for Trilinos => needed for compatibility with other Spaces
        KRATOS_TRY;
        return EpetraExt::RowMatrixToMatrixMarketFile(pFileName, rM); // Returns 0 if no error, -1 if any problems with file system.
        KRATOS_CATCH("");
    }

    /**
     * @brief Writes a vector to a file in MatrixMarket format
     * @param pFileName The name of the file to be written
     * @param rV The vector to be written
     * @return True if the file was successfully written, false otherwise
     */
    template< class VectorType >
    static bool WriteMatrixMarketVector(
        const char* pFileName,
        const VectorType& rV
        )
    {
        KRATOS_TRY;
        return EpetraExt::MultiVectorToMatrixMarketFile(pFileName, rV);
        KRATOS_CATCH("");
    }

    /**
     * @brief Creates a new dof updater
     * @return The new dof updater
     */
    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
    }

    ///@}
private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosSpace & operator=(TrilinosSpace const& rOther);

    /// Copy constructor.
    TrilinosSpace(TrilinosSpace const& rOther);

    ///@}
}; // Class TrilinosSpace

///@}

} // namespace Kratos.
