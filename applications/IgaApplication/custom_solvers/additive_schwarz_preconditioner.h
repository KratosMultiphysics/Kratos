//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/preconditioner.h"

namespace Kratos
{

///@name  Preconditioners
///@{

/// AdditiveSchwarzPreconditioner class.
/**   */
template<class TSparseSpaceType, class TDenseSpaceType>
class AdditiveSchwarzPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AdditiveSchwarzPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION(AdditiveSchwarzPreconditioner);


    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdditiveSchwarzPreconditioner()
    {
        //@TODO
    }


    /// Copy constructor.
    AdditiveSchwarzPreconditioner(const AdditiveSchwarzPreconditioner& Other) {}


    /// Destructor.
    ~AdditiveSchwarzPreconditioner() override
    {
        //@TODO
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdditiveSchwarzPreconditioner& operator=(const AdditiveSchwarzPreconditioner& Other)
    {

    }

    ///@}
    ///@name Operations
    ///@{

    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {
        VectorType z = rX;
        TSparseSpaceType::Mult(rA,z, rY);
        ApplyLeft(rY);
    }

    void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {
        VectorType z = rX;
        ApplyTransposeLeft(z);
        TSparseSpaceType::TransposeMult(rA,z, rY);
    }

    /** multiply first rX by L^-1 and store result in temp
        then multiply temp by U^-1 and store result in rX
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        // const int size = TSparseSpaceType::Size(rX);
        // VectorType temp(size);
        // double sum;
        // int i, indexj;
        // for (i=0; i<size; i++)
        // {
        //     sum=rX[i];
        //     for (indexj=iL[i]; indexj<iL[i+1]; indexj++)
        //     {
        //         sum=sum-L[indexj]*temp[jL[indexj]];
        //     }
        //     temp[i]=sum;
        // }
        // for (i=size-1; i>=0; i--)
        // {
        //     sum=temp[i];
        //     for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++)
        //     {
        //         sum=sum-U[indexj]*rX[jU[indexj]];
        //     }
        //     rX[i]=sum/U[iU[i]];
        // }
        // return rX;
    }

    /** Multiply first rX by U^-T and store result in temp
        then multiply temp by L^-T and store result in rX
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyTransposeLeft(VectorType& rX) override
    {
        // const int size = TSparseSpaceType::Size(rX);
        // VectorType temp(size);
        // int i, indexj;
        // double tempi, rxi;
        // for (i=0; i<size; i++) temp[i]=rX[i];
        // for (i=0; i<size; i++)
        // {
        //     temp[i]=temp[i]/U[iU[i]];
        //     tempi=temp[i];
        //     for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++)
        //     {
        //         temp[jU[indexj]]=temp[jU[indexj]]-tempi*U[indexj];
        //     }
        // }
        // for (i=0; i<size; i++) rX[i]=temp[i];
        // for (i=size-1; i>=0; i--)
        // {
        //     rxi=rX[i];
        //     for (indexj=iL[i]; indexj<iL[i+1]; indexj++)
        //     {
        //         rX[jL[indexj]]=rX[jL[indexj]]-rxi*L[indexj];
        //     }
        // }
        // return rX;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Return information about this object.
    std::string Info() const override
    {
        return "AdditiveSchwarzPreconditioner";
    }


    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "AdditiveSchwarzPreconditioner";
    }

    void PrintData(std::ostream& OStream) const override
    {
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    unsigned int mILUSize;
    int *iL, *jL, *iU, *jU;
    double *L, *U;

    ///@}

}; // Class AdditiveSchwarzPreconditioner

///@}

///@}

///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,
                                  AdditiveSchwarzPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}


/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const AdditiveSchwarzPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);


    return OStream;
}
///@}


}  // namespace Kratos.


#endif // KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED  defined

