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
#include <EpetraExt_MatrixMatrix.h>

// NOTE: Epetra will be replaced with Tpetra in the future, there is an intermediate interface called Xpetra, which uses the interface of Toetra but allows to use Epetra

// Project includes
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "includes/parallel_environment.h"
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

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(nullptr);
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(nullptr);
    }

    static MatrixPointerType CreateEmptyMatrixPointer(Epetra_MpiComm& Comm)
    {
        int global_elems = 0;
        Epetra_Map Map(global_elems, 0, Comm);
        return MatrixPointerType(new TMatrixType(::Copy, Map, 0));
    }

    static VectorPointerType CreateEmptyVectorPointer(Epetra_MpiComm& Comm)
    {
        int global_elems = 0;
        Epetra_Map Map(global_elems, 0, Comm);
        return VectorPointerType(new TVectorType(Map));
    }

    /// return size of vector rV
    static IndexType Size(VectorType const& rV)
    {
        int size;
        size = rV.GlobalLength();
        return size;
    }

    /// return number of rows of rM
    static IndexType Size1(MatrixType const& rM)
    {
        int size1;
        size1 = rM.NumGlobalRows();
        return size1;
    }

    /// return number of columns of rM
    static IndexType Size2(MatrixType const& rM)
    {
        int size1;
        size1 = rM.NumGlobalCols();
        return size1;
    }

    /// rXi = rMij
    static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
    {
        KRATOS_ERROR << "GetColumn method is not currently implemented" << std::endl;
    }

    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY = rX;
    }

    /// rY = rX
    static void Copy(VectorType const& rX, VectorType& rY)
    {
        rY = rX;
    }

    /// rX * rY
    static double Dot(VectorType& rX, VectorType& rY)
    {
        double value;
        const int sucess = rY.Dot(rX, &value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing dot product" <<  std::endl;
        return value;
    }

    /// Maximum value
    static double Max(const VectorType& rX)
    {
        double value;
        const int sucess = rX.MaxValue(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing maximum value" <<  std::endl;
        return value;
    }

    /// Minimum value
    static double Min(const VectorType& rX)
    {
        double value;
        const int sucess = rX.MinValue(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing minimum value" <<  std::endl;
        return value;
    }

    /// ||rX||2
    static double TwoNorm(const VectorType& rX)
    {
        double value;
        const int sucess = rX.Norm2(&value); //it is prepared to handle vectors with multiple components
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing norm" <<  std::endl;
        return value;
    }

    /// ||rA||2
    static double TwoNorm(const MatrixType& rA)
    {
        double value;
        const int sucess = rA.Norm2(&value);
        KRATOS_ERROR_IF_NOT(sucess == 0) << "Error computing norm" <<  std::endl;
        return value;
    }

    static void Mult(
        const MatrixType& rA, 
        const VectorType& rX, 
        VectorType& rY
        )
    {
        //y = A*x
        const bool transpose_flag = false;
        rA.Multiply(transpose_flag, rX, rY);
    }

    static void Mult(
        const MatrixType& rA, 
        const MatrixType& rB, 
        MatrixType& rC
        )
    {
        //C = A*B
        const bool transpose_flag = false;
        EpetraExt::MatrixMatrix::Multiply(rA, transpose_flag, rB, transpose_flag, rC);
    }

    static void TransposeMult(
        const MatrixType& rA, 
        const VectorType& rX, 
        VectorType& rY
        )
    {
        //y = AT*x
        const bool transpose_flag = true;
        rA.Multiply(transpose_flag, rX, rY);
    }

    static void TransposeMult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const std::vector<bool> TransposeFlag = {false, false}
        )
    {
        //C = A*B
        KRATOS_ERROR_IF_NOT(TransposeFlag.size() > 1) << "Size of flags must be at least 2" << std::endl;
        EpetraExt::MatrixMatrix::Multiply(rA, TransposeFlag[0], rB, TransposeFlag[1], rC);
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    static void InplaceMult(VectorType& rX, const double A)
    {
        if (A != 1.00)
            rX.Scale(A);
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;

    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        if (A != 1.00)
            rX.Scale(A, rY); //not sure
        else
            rX = rY;
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;

    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        rX.Update(A, rY, 1.0);
    }

    //********************************************************************

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        rZ.Update(A, rX, B, rY, 0.0);
    }

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        rY.Update(A, rX, B);
    }

    /// rA[i] * rX
    //       static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    // 	{
    // 	  return inner_prod(row(rA, i), rX);
    // 	}
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

    /// rX = A

    static void Set(VectorType& rX, DataType A)
    {
        rX.PutScalar(A);
    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        KRATOS_ERROR << "Resize is not defined for Trilinos Sparse Matrix" << std::endl;
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        KRATOS_ERROR << "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer" << std::endl;
    }

    static void Resize(VectorPointerType& pX, SizeType n)
    {
        //KRATOS_ERROR_IF(pX != NULL) << "trying to resize a null pointer" << std::endl;
        int global_elems = n;
        Epetra_Map Map(global_elems, 0, pX->Comm());
        VectorPointerType pNewEmptyX = Kratos::make_shared<VectorType>(Map);
        pX.swap(pNewEmptyX);
    }

    static void Clear(MatrixPointerType& pA)
    {
        if(pA != NULL) {
            int global_elems = 0;
            Epetra_Map Map(global_elems, 0, pA->Comm());
            MatrixPointerType pNewEmptyA = MatrixPointerType(new TMatrixType(::Copy, Map, 0));
            pA.swap(pNewEmptyA);
        }
    }

    static void Clear(VectorPointerType& pX)
    {
        if(pX != NULL) {
            int global_elems = 0;
            Epetra_Map Map(global_elems, 0, pX->Comm());
            VectorPointerType pNewEmptyX = VectorPointerType(new VectorType(Map));
            pX.swap(pNewEmptyX);
        }
    }

    inline static void SetToZero(MatrixType& rA)
    {
        rA.PutScalar(0.0);
    }

    inline static void SetToZero(VectorType& rX)
    {
        rX.PutScalar(0.0);
    }

    /// TODO: creating the the calculating reaction version
    // 	template<class TOtherMatrixType, class TEquationIdVectorType>

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

    //***********************************************************************

    inline static constexpr bool IsDistributed()
    {
        return true;
    }

    //***********************************************************************

    inline static double GetValue(const VectorType& x, std::size_t I)
    {
        // index must be local to this proc
        KRATOS_ERROR_IF_NOT(x.Map().MyGID(static_cast<int>(I))) << " non-local id: " << I << ".";
        // Epetra_MultiVector::operator[] is used here to get the pointer to
        // the zeroth (only) local vector.
        return x[0][x.Map().LID(static_cast<int>(I))];
    }

    //***********************************************************************

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
     * @param ScalingDiagonal The type of caling diagonal considered
     * @param rDataCommunicator The data coomunicator considered
     * @return The scale norm
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal, rDataCommunicator);

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
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of caling diagonal considered
     * @param rDataCommunicator The data coomunicator considered
     * @return The scale norm
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
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
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    static double GetDiagonalNorm(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);

        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::TwoNorm(diagonal);
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data coomunicator considered
     * @return The diagonal  max value
     */
    static double GetAveragevalueDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
        return 0.5 * (GetMaxDiagonal(rA, rDataCommunicator) + GetMinDiagonal(rA, rDataCommunicator));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data coomunicator considered
     * @return The diagonal  max value
     */
    static double GetMaxDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Max(diagonal);
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @param rDataCommunicator The data coomunicator considered
     * @return The diagonal min value
     */
    static double GetMinDiagonal(
        const MatrixType& rA,
        const DataCommunicator& rDataCommunicator = ParallelEnvironment::GetDefaultDataCommunicator()
        )
    {
        // Generate Epetra communicator
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
        Epetra_MpiComm epetra_comm(raw_mpi_comm);

        const int global_elems = Size1(rA);
        Epetra_Map map(global_elems, 0, epetra_comm);
        Epetra_Vector diagonal(map);
        rA.ExtractDiagonalCopy(diagonal);
        return TrilinosSpace<Epetra_FECrsMatrix, Epetra_Vector>::Min(diagonal);
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

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        return "TrilinosSpace";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TrilinosSpace";
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char* pFileName, const TOtherMatrixType& rM, const bool Symmetric)
    {
        // the argument "Symmetric" does not have an effect for Trilinos => needed for compatibility with other Spaces
        KRATOS_TRY;
        return EpetraExt::RowMatrixToMatrixMarketFile(pFileName, rM); // Returns 0 if no error, -1 if any problems with file system.
        KRATOS_CATCH("");
    }

    template< class VectorType >
    static bool WriteMatrixMarketVector(const char* pFileName, const VectorType& rV)
    {
        KRATOS_TRY;
        return EpetraExt::MultiVectorToMatrixMarketFile(pFileName, rV);
        KRATOS_CATCH("");
    }


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

