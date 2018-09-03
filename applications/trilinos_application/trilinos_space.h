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
//

#if !defined(KRATOS_TRILINOS_SPACE_H_INCLUDED )
#define  KRATOS_TRILINOS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <sstream>

// External includes


//Trilinos includes
#include "mpi.h"
#include "Epetra_Import.h"

//#include "epetra_test_err.h"
#include "Epetra_MpiComm.h"
// #include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "EpetraExt_CrsMatrixIn.h"
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
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

/// Short class definition.

/** Detail class definition.
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
    //
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
        KRATOS_THROW_ERROR(std::logic_error, "GetColumn method is not currently implemented", "")
    }


    ///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
    /// rMij = rXi
    //      static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = row(rM, j);}

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
        double *temp = new double[1];
        rY.Dot(rX, temp); //it is prepared to handle vectors with multiple components
        value = temp[0];
        delete [] temp;
        return value;
    }

    /// ||rX||2

    static double TwoNorm(VectorType const& rX)
    {
        double value;
        double *temp = new double[1];
        rX.Norm2(temp); //it is prepared to handle vectors with multiple components
        value = temp[0];
        delete [] temp;
        return value;
    }

    static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        //y = A*x
        bool transpose_flag = false;
        rA.Multiply(transpose_flag, rX, rY);
    }

    static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        //y = A*x
        bool transpose_flag = true;
        rA.Multiply(transpose_flag, rX, rY);
    } // rY = rAT * rX


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
        KRATOS_THROW_ERROR(std::logic_error, "Resize is not defined for Trilinos Sparse Matrix", "")

    }

    static void Resize(VectorType& rX, SizeType n)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer", "")
    }

    static void Resize(VectorPointerType& pX, SizeType n)
    {
//         if(pX != NULL)
//             KRATOS_ERROR << "trying to resize a null pointer" ;
        int global_elems = n;
        Epetra_Map Map(global_elems, 0, pX->Comm());
        VectorPointerType pNewEmptyX = Kratos::make_shared<VectorType>(Map);
        pX.swap(pNewEmptyX);
    }

    // 	static void Clear(MatrixType& rA)
    // 	{
    // 		int global_elems = 0;
    // 		Epetra_Map Map(global_elems,0,rA.Comm() );
    // 		TMatrixType temp(::Copy,Map,0);
    //
    // // 		int global_elems = 0;
    // // 		Epetra_Map my_map(global_elems,0,rA.Comm() );
    // // 		Epetra_FECrsGraph Agraph(Copy, my_map, 0);
    // // 		Agraph.GlobalAssemble();
    // //    		TMatrixType temp(Copy,Agraph);
    //
    // // int NumGlobalElements = 10;
    // //  Epetra_Map Map(NumGlobalElements,0,rA.Comm() );
    // //
    // //   // create a diagonal FE crs matrix (one nonzero per row)
    // //   Epetra_FECrsMatrix temp(Copy,Map,1);
    //
    //   		rA = temp;
    // 	}

    static void Clear(MatrixPointerType& pA)
    {
        if(pA != NULL)
        {
            int global_elems = 0;
            Epetra_Map Map(global_elems, 0, pA->Comm());
            MatrixPointerType pNewEmptyA = MatrixPointerType(new TMatrixType(::Copy, Map, 0));
            pA.swap(pNewEmptyA);
        }
    }

    static void Clear(VectorPointerType& pX)
    {
        if(pX != NULL)
        {
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

        if (active_indices > 0)
        {
            //size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseMatrix values(active_indices, active_indices);

            //fill epetra vectors
            int loc_i = 0;
            for (unsigned int i = 0; i < EquationId.size(); i++)
            {
                if (EquationId[i] < system_size)
                {
                    indices[loc_i] = EquationId[i];

                    int loc_j = 0;
                    for (unsigned int j = 0; j < EquationId.size(); j++)
                    {
                        if (EquationId[j] < system_size)
                        {
                            values(loc_i, loc_j) = LHS_Contribution(i, j);
                            loc_j += 1;
                        }
                    }
                    loc_i += 1;
                }
            }

            int ierr = A.SumIntoGlobalValues(indices, values);
            if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");





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

        if (active_indices > 0)
        {
            //size Epetra vectors
            Epetra_IntSerialDenseVector indices(active_indices);
            Epetra_SerialDenseVector values(active_indices);

            //fill epetra vectors
            int loc_i = 0;
            for (unsigned int i = 0; i < EquationId.size(); i++)
            {
                if (EquationId[i] < system_size)
                {
                    indices[loc_i] = EquationId[i];
                    values[loc_i] = RHS_Contribution[i];
                    loc_i += 1;
                }
            }

            int ierr = b.SumIntoGlobalValues(indices, values);
            if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");

            //EPETRA_TEST_ERR( ierr );
        }
    }

    //***********************************************************************

    inline static bool IsDistributed()
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

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, double* pValues)
    {
        KRATOS_TRY
        double tot_size = IndexArray.size();

        //defining a map as needed
        Epetra_Map dof_update_map(-1, tot_size, &(*(IndexArray.begin())), 0, x.Comm());

        //defining the importer class
        Epetra_Import importer(dof_update_map, x.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp(importer.TargetMap());

        //importing in the new temp vector the values
        int ierr = temp.Import(x, importer, Insert);
        if(ierr != 0) KRATOS_THROW_ERROR(std::logic_error,"Epetra failure found","");


        temp.ExtractView(&pValues);

        x.Comm().Barrier();
        KRATOS_CATCH("")

    }

    MatrixPointerType ReadMatrixMarket(const std::string FileName,Epetra_MpiComm& Comm)
    {
        KRATOS_TRY

        Epetra_CrsMatrix* pp = nullptr;


        int error_code = EpetraExt::MatrixMarketFileToCrsMatrix(FileName.c_str(), Comm, pp);

        if(error_code != 0)
            KRATOS_ERROR << "error thrown while reading Matrix Market file "<<FileName<< " error code is : " << error_code;

        Comm.Barrier();

        const Epetra_CrsGraph& rGraph = pp->Graph();
        MatrixPointerType paux = Kratos::make_shared<Epetra_FECrsMatrix>( ::Copy, rGraph, false );

        IndexType NumMyRows = rGraph.RowMap().NumMyElements();

        int* MyGlobalElements = new int[NumMyRows];
        rGraph.RowMap().MyGlobalElements(MyGlobalElements);

        for(IndexType i = 0; i < NumMyRows; ++i)
        {
//             std::cout << pA->Comm().MyPID() << " : I=" << i << std::endl;
            IndexType GlobalRow = MyGlobalElements[i];

            int NumEntries;
            std::size_t Length = pp->NumGlobalEntries(GlobalRow);  // length of Values and Indices

            double* Values = new double[Length];     // extracted values for this row
            int* Indices = new int[Length];          // extracted global column indices for the corresponding values

            error_code = pp->ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, Values, Indices);

            if(error_code != 0)
                KRATOS_ERROR << "error thrown in ExtractGlobalRowCopy : " << error_code;

            error_code = paux->ReplaceGlobalValues(GlobalRow, Length, Values, Indices);

            if(error_code != 0)
                KRATOS_ERROR << "error thrown in ReplaceGlobalValues : " << error_code;

            delete [] Values;
            delete [] Indices;
        }

        paux->GlobalAssemble();

        delete [] MyGlobalElements;
        delete pp;

        return paux;
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

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    TrilinosSpace& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
// 				    const TrilinosSpace& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
///@}


} // namespace Kratos.

#endif // KRATOS_TRILINOS_SPACE_H_INCLUDED  defined
