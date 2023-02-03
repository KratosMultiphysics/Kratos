//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
 */
template<class TMatrixType, class TVectorType>
class TrilinosSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosSpace
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosSpace);

    typedef double DataType;

    typedef TMatrixType MatrixType;

    typedef TVectorType VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef typename Kratos::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename Kratos::shared_ptr< TVectorType > VectorPointerType;

    typedef TrilinosDofUpdater< TrilinosSpace<TMatrixType,TVectorType> > DofUpdaterType;
    typedef typename DofUpdater<TrilinosSpace<TMatrixType,TVectorType> >::UniquePointer DofUpdaterPointerType;

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
    static MatrixPointerType CreateEmptyMatrixPointer(Epetra_MpiComm& Comm)
    {
        int global_elems = 0;
        Epetra_Map Map(global_elems, 0, Comm);
        return MatrixPointerType(new TMatrixType(::Copy, Map, 0));
    }

    /**
     * @brief This method creates an empty pointer to a vector using epetra communicator
     * @param rComm The epetra communicator
     * @return The pointer to the vector
     */
    static VectorPointerType CreateEmptyVectorPointer(Epetra_MpiComm& Comm)
    {
        int global_elems = 0;
        Epetra_Map Map(global_elems, 0, Comm);
        return VectorPointerType(new TVectorType(Map));
    }

    /**
     * @brief Returns size of vector rV
     * @param rV The vector considered
     * @return The size of the vector
     */
    static IndexType Size(VectorType const& rV)
    {
        int size;
        size = rV.GlobalLength();
        return size;
    }

    /**
     * @brief Returns number of rows of rM
     * @param rM The matrix considered
     * @return The number of rows of rM
     */
    static IndexType Size1(MatrixType const& rM)
    {
        int size1;
        size1 = rM.NumGlobalRows();
        return size1;
    }

    /**
     * @brief Returns number of columns of rM
     * @param rM The matrix considered
     * @return The number of columns of rM
     */
    static IndexType Size2(MatrixType const& rM)
    {
        int size1;
        size1 = rM.NumGlobalCols();
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
    static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
    {
        KRATOS_ERROR << "GetColumn method is not currently implemented" << std::endl;
    }

    ///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
    /// rMij = rXi
    //      static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = row(rM, j);}

    /**
     * @brief Returns a copy of the matrix rX
     * @details rY = rX
     * @param rX The matrix considered
     * @param rY The copy of the matrix rX
     */
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY = rX;
    }

    /**
     * @brief Returns a copy of the vector rX
     * @details rY = rX
     * @param rX The vector considered
     * @param rY The copy of the vector rX
     */
    static void Copy(VectorType const& rX, VectorType& rY)
    {
        rY = rX;
    }

    /**
     * @brief Returns the product of two vectors
     * @details rX * rY
     * @param rX The first vector considered
     * @param rY The second vector considered
     */
    static double Dot(VectorType& rX, VectorType& rY)
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
    static double TwoNorm(VectorType const& rX)
    {
        double value;
        const int sucess = rX.Norm2(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing norm" <<  std::endl;
        return value;
    }

    /**
     * @brief Returns the multiplication of a matrix by a vector
     * @details y = A*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        //y = A*x
        bool transpose_flag = false;
        rA.Multiply(transpose_flag, rX, rY);
    }

    /**
     * @brief Returns the transpose multiplication of a matrix by a vector
     * @details y = AT*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        //y = A*x
        bool transpose_flag = true;
        rA.Multiply(transpose_flag, rX, rY);
    } // rY = rAT * rX

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details y = A*x
     * Checks if a multiplication is needed and tries to do otherwise
     * @param rX The vector considered
     * @param A The scalar considered
     */
    static void InplaceMult(VectorType& rX, const double A)
    {
        if (A != 1.00)
            rX.Scale(A);
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
    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        if (A != 1.00)
            rX.Scale(A, rY); //not sure
        else
            rX = rY;
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
    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        rX.Update(A, rY, 1.0);
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
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        rZ.Update(A, rX, B, rY, 0.0);
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar 
     * @details rY = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The resulting vector considered
     */
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        rY.Update(A, rX, B);
    }

    /// rA[i] * rX
    //       static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    // 	{
    // 	  return inner_prod(row(rA, i), rX);
    // 	}

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param value The value considered
     */
    static void SetValue(VectorType& rX, IndexType i, double value)
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
    static void Set(VectorType& rX, DataType A)
    {
        rX.PutScalar(A);
    }

    /**
     * @brief Resizes a matrix
     * @param rA The matrix to be resized
     * @param m The new number of rows
     * @param n The new number of columns
     */
    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        KRATOS_ERROR << "Resize is not defined for Trilinos Sparse Matrix" << std::endl;

    }

    /**
     * @brief Resizes a vector
     * @param rX The vector to be resized
     * @param n The new size
     */
    static void Resize(VectorType& rX, SizeType n)
    {
        KRATOS_ERROR << "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param pA The pointer to the vector to be resized
     * @param n The new size
     */
    static void Resize(VectorPointerType& pX, SizeType n)
    {
        //KRATOS_ERROR_IF(pX != NULL) << "trying to resize a null pointer" << std::endl;
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
        rA.PutScalar(0.0);
    }

    /**
     * @brief Sets a vector to zero
     * @param rX The vector to be set
     */
    inline static void SetToZero(VectorType& rX)
    {
        rX.PutScalar(0.0);
    }

    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherMatrixType, class TEquationIdVectorType>

    /**
     * @brief Assembles the LHS of the system
     * @param A The LHS matrix
     * @param LHS_Contribution The contribution to the LHS
     * @param EquationId The equation ids
     */
    inline static void AssembleLHS(
        MatrixType& A,
        Matrix& LHS_Contribution,
        std::vector<std::size_t>& EquationId
        )
    {
        unsigned int system_size = Size1(A);
        //unsigned int local_size = LHS_Contribution.size1();

        //count active indices
        unsigned int active_indices = 0;
        for (unsigned int i = 0; i < EquationId.size(); i++)
            if (EquationId[i] < system_size)
                active_indices += 1;

        if (active_indices > 0) {
            //size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseMatrix values(active_indices, active_indices);

            //fill epetra vectors
            int loc_i = 0;
            for (unsigned int i = 0; i < EquationId.size(); i++) {
                if (EquationId[i] < system_size) {
                    indices[loc_i] = EquationId[i];

                    int loc_j = 0;
                    for (unsigned int j = 0; j < EquationId.size(); j++) {
                        if (EquationId[j] < system_size) {
                            values(loc_i, loc_j) = LHS_Contribution(i, j);
                            loc_j += 1;
                        }
                    }
                    loc_i += 1;
                }
            }

            int ierr = A.SumIntoGlobalValues(indices, values);
            KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found" << std::endl;
        }
    }

    //***********************************************************************
    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherVectorType, class TEquationIdVectorType>

    /**
     * @brief Assembles the RHS of the system
     * @param b The RHS vector
     * @param RHS_Contribution The RHS contribution
     * @param EquationId The equation ids
     */
    inline static void AssembleRHS(
        VectorType& b,
        Vector& RHS_Contribution,
        std::vector<std::size_t>& EquationId
        )
    {
        unsigned int system_size = Size(b);
        //unsigned int local_size = RHS_Contribution.size();

        //count active indices
        int active_indices = 0;
        for (unsigned int i = 0; i < EquationId.size(); i++)
            if (EquationId[i] < system_size)
                active_indices += 1;

        if (active_indices > 0) {
            //size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseVector values(active_indices);

            //fill epetra vectors
            int loc_i = 0;
            for (unsigned int i = 0; i < EquationId.size(); i++) {
                if (EquationId[i] < system_size) {
                    indices[loc_i] = EquationId[i];
                    values[loc_i] = RHS_Contribution[i];
                    loc_i += 1;
                }
            }

            int ierr = b.SumIntoGlobalValues(indices, values);
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
    inline static double GetValue(const VectorType& x, std::size_t I)
    {
        // index must be local to this proc
        KRATOS_ERROR_IF_NOT(x.Map().MyGID(static_cast<int>(I))) << " non-local id: " << I << ".";
        // Epetra_MultiVector::operator[] is used here to get the pointer to
        // the zeroth (only) local vector.
        return x[0][x.Map().LID(static_cast<int>(I))];
    }

    /**
     * @brief This function gathers the values of a given vector according to a given index array
     * @param rX The vector from which values are to be gathered
     * @param IndexArray The array containing the indices of the values to be gathered
     * @param pValues The array containing the gathered values
     */
    static void GatherValues(const VectorType& rX, const std::vector<int>& IndexArray, double* pValues)
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
    MatrixPointerType ReadMatrixMarket(const std::string FileName,Epetra_MpiComm& Comm)
    {
        KRATOS_TRY

        Epetra_CrsMatrix* pp = nullptr;

        int error_code = EpetraExt::MatrixMarketFileToCrsMatrix(FileName.c_str(), Comm, pp);

        KRATOS_ERROR_IF(error_code != 0) << "Eerror thrown while reading Matrix Market file "<<FileName<< " error code is : " << error_code;

        Comm.Barrier();

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
    VectorPointerType ReadMatrixMarketVector(const std::string& rFileName, Epetra_MpiComm& rComm, int N)
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
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param rDataCommunicator The data communicator considered
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const DataCommunicator& rDataCommunicator,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, rDataCommunicator, ScalingDiagonal);

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
     * @param rDataCommunicator The data communicator considered
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator,
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
                return GetDiagonalNorm(rA, rDataCommunicator)/static_cast<double>(Size1(rA));
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA, rDataCommunicator);
            default:
                return GetMaxDiagonal(rA, rDataCommunicator);
        }

        KRATOS_CATCH("");
    }
    
    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    static double GetDiagonalNorm(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator
        )
    {
        KRATOS_TRY

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);

        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::TwoNorm(diagonal);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data communicator considered
     * @return The diagonal  max value
     */
    static double GetAveragevalueDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator
        )
    {
        KRATOS_TRY

        return 0.5 * (GetMaxDiagonal(rA, rDataCommunicator) + GetMinDiagonal(rA, rDataCommunicator));

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data communicator considered
     * @return The diagonal  max value
     */
    static double GetMaxDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator
        )
    {
        KRATOS_TRY

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Max(diagonal);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data communicator considered
     * @return The diagonal min value
     */
    static double GetMinDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator
        )
    {
        KRATOS_TRY

        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Min(diagonal);

        KRATOS_CATCH("");
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
    static bool WriteMatrixMarketMatrix(const char* pFileName, const TOtherMatrixType& rM, const bool Symmetric)
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
    static bool WriteMatrixMarketVector(const char* pFileName, const VectorType& rV)
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
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
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

