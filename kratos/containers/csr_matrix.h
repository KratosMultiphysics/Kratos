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

#if !defined(KRATOS_CSR_MATRIX_H_INCLUDED )
#define  KRATOS_CSR_MATRIX_H_INCLUDED


// System includes
#include <iostream>
#include "containers/sparse_contiguous_row_graph.h"

// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

/// This function implements "serial" CSR matrix, including capabilities for FEM assembly
template< class TDataType=double, class TIndexType=std::size_t>
class CsrMatrix
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;

    /// Pointer definition of CsrMatrix
    KRATOS_CLASS_POINTER_DEFINITION(CsrMatrix);

    ///@}
    ///@name Life Cycle
    ///@{
    CsrMatrix()
    {
    }

    /// Default constructor.
    template<class TGraphType>
    CsrMatrix(const TGraphType& rSparseGraph)
    {
        rSparseGraph.ExportCSRArrays(mRowIndices,mColIndices);
        mNrows = rSparseGraph.Size();

        //TODO: parallelize
        IndexType max_col = 0;
        for(IndexType i=0; i<mColIndices.size(); ++i)
            max_col = std::max(max_col,mColIndices[i]);
        mNcols = max_col;

        //initialize mValuesVector to zero
        mValuesVector.resize(mColIndices.size(),false);
        SetValue(0.0);
    }

    /// Destructor.
    virtual ~CsrMatrix(){}

    /// Assignment operator. TODO: decide if we do want to allow it
    CsrMatrix& operator=(CsrMatrix const& rOther)=delete;
    // {
    //     this->AddEntries(rOther.GetGraph());
    //     return *this;
    // }

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        mRowIndices.clear();
        mColIndices.clear();
        mValuesVector.clear();
    }

    void SetValue(const TDataType value)
    {
        #pragma omp parallel for
        for(IndexType i=0; i<mValuesVector.size(); ++i)
            mValuesVector[i] = value;
    }

    IndexType size1()
    {
        return mRowIndices.size()-1;
    }

    IndexType size2()
    {
        return mNcols;
    }

    inline IndexType nnz() const{
        return index2_data().size();
    }

    inline DenseVector<IndexType>& index1_data(){
        return mRowIndices;
    }
    inline DenseVector<IndexType>& index2_data(){
        return mColIndices;
    }
    inline DenseVector<TDataType>& value_data(){
        return mValuesVector;
    }

    inline const DenseVector<IndexType>& index1_data() const{
        return mRowIndices;
    }
    inline const DenseVector<IndexType>& index2_data() const{
        return mColIndices;
    }
    inline const DenseVector<TDataType>& value_data() const{
        return mValuesVector;
    }

    TDataType& operator()(IndexType I, IndexType J){
        const IndexType row_begin = index1_data()[I];
        const IndexType row_end = index1_data()[I+1];
        IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
        KRATOS_DEBUG_ERROR_IF(k<0) << "local indices I,J : " << I << " " << J << " not found in matrix" << std::endl;
        return value_data()[k];
    }

    bool Has(IndexType I, IndexType J) const {
        const IndexType row_begin = index1_data()[I];
        const IndexType row_end = index1_data()[I+1];
        IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
        if(k >= 0)
            return true;
        else
            return false;
    }

    ///@}
    ///@name Operations
    ///@{
    //TODO
    //+=
    //-=
    //+
    //-
    //*

    template<class TOutputVector, class TInputVector>
    static void MultAndAdd(
        TOutputVector& rOutputVector,
        const CsrMatrix& rA,
        const TInputVector& rInputVector
    )
    {
        #pragma omp parallel for
        for(int i=0; i<rA.index1_data().size()-1; ++i)
        {
            for(int k=rA.index1_data()[i]; k<rA.index1_data()[i+1]; ++k)
            {
                auto j = rA.index2_data()[k];
                rOutputVector[i] += rA.value_data()[k]*rInputVector[j];
            }
        }
    }

    void reserve(IndexType nrows, IndexType nnz){
        index1_data().resize(nrows+1,false);
        if(nrows > 0)
            index1_data()[0] = 0;
        index2_data().resize(nnz,false);
        value_data().resize(nnz,false);
        mNrows = nrows;
    }


    void BeginAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case

    void FinalizeAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case





    template<class TMatrixType, class TIndexVectorType >
    void Assemble(
        const TMatrixType& rMatrixInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size1() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size2() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;

        unsigned int local_size = rMatrixInput.size1();

        for (unsigned int i_local = 0; i_local < local_size; ++i_local)
        {
            const IndexType I = EquationId[i_local];
            const IndexType row_begin = index1_data()[I];
            const IndexType row_end = index1_data()[I+1];

            //find first entry (note that we know it exists since local_size > 0)
            IndexType J = EquationId[0];
            IndexType k = BinarySearch(index2_data(), row_begin, row_end, EquationId[0]);
            IndexType lastJ = J;

            #pragma omp atomic
            value_data()[k] += rMatrixInput(i_local,0);

            //now find other entries. note that we assume that it is probably that next entries immediately follow in the ordering
            for(unsigned int j_local=1; j_local<local_size; ++j_local){
                J = EquationId[j_local];

                if(k+1<index2_data().size() && index2_data()[k+1] == J){
                    k = k+1;
                }
                else if(J > lastJ){ //note that the case k+2 >= index2_data().size() should be impossible
                    k = BinarySearch(index2_data(), k+2, row_end, J);
                }
                else if(J < lastJ){
                    k = BinarySearch(index2_data(), row_begin, k-1, J);
                }
                //the last missing case is J == lastJ, which should never happen in FEM. If that happens we can reuse k

                #pragma omp atomic
                value_data()[k] += rMatrixInput(i_local,j_local);

                lastJ = J;
            }
        }
    }

    //TODO
    // LeftScaling
    // RightScaling
    // SymmetricScaling

    //TODO
    //void ApplyDirichlet

    //TODO
    //NormFrobenius

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
        std::stringstream buffer;
        buffer << "CsrMatrix" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "CsrMatrix";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
    // A non-recursive binary search function. It returns 
    // location of x in given array arr[l..r] is present, 
    // otherwise -1
    template< class TVectorType > 
    inline IndexType BinarySearch(const TVectorType& arr, 
                        IndexType l, IndexType r, IndexType x) const
    { 
        while (l <= r) { 
            int m = l + (r - l) / 2; 
    
            // Check if x is present at mid 
            if (arr[m] == x) 
                return m; 
    
            // If x greater, ignore left half 
            if (arr[m] < x) 
                l = m + 1; 
    
            // If x is smaller, ignore right half 
            else
                r = m - 1; 
        } 
    
        // if we reach here, then element was not present 
        return -1; 
    } 

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
    vector<IndexType> mRowIndices;
    vector<IndexType> mColIndices;
    vector<TDataType> mValuesVector;
    IndexType mNrows=0;
    IndexType mNcols=0;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("rows",mRowIndices);
        rSerializer.save("cols",mColIndices);
        rSerializer.save("values",mValuesVector);
        rSerializer.save("Nrow",mNrows);
        rSerializer.save("Ncol",mNcols);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("rows",mRowIndices);
        rSerializer.load("cols",mColIndices);
        rSerializer.load("values",mValuesVector);
        rSerializer.load("Nrow",mNrows);
        rSerializer.load("Ncol",mNcols);
    }


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



    ///@}

}; // Class CsrMatrix

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDataType>
inline std::istream& operator >> (std::istream& rIStream,
                CsrMatrix<TDataType>& rThis){
                    return rIStream;
                }

/// output stream function
template< class TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                const CsrMatrix<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CSR_MATRIX_H_INCLUDED  defined


