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
#include <limits>
#include <span/span.hpp>

#include "containers/sparse_contiguous_row_graph.h"
#include "containers/system_vector.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"
#include "includes/key_hash.h"
#include "includes/parallel_environment.h"

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

/// This class implements "serial" CSR matrix, including capabilities for FEM assembly
template< class TDataType=double, class TIndexType=std::size_t>
class CsrMatrix
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef std::unordered_map<std::pair<IndexType, IndexType>,
                          double,
                          PairHasher<IndexType, IndexType>,
                          PairComparor<IndexType, IndexType>
                          > MatrixMapType;

    /// Pointer definition of CsrMatrix
    KRATOS_CLASS_POINTER_DEFINITION(CsrMatrix);

    ///@}
    ///@name Life Cycle
    ///@{
    CsrMatrix() //needs to be public, since one could use the low level API to construc the CSR matrix
    {
        mpComm = &ParallelEnvironment::GetDataCommunicator("Serial");
    }

    CsrMatrix(const DataCommunicator& rComm) //needs to be public, since one could use the low level API to construc the CSR matrix
    {
        if(rComm.IsDistributed())
            KRATOS_ERROR << "Attempting to construct a serial CsrMatrix with a distributed communicator" << std::endl;

        mpComm = &rComm;
    }
    

    /// constructor.
    template<class TGraphType>
    CsrMatrix(const TGraphType& rSparseGraph)
    {
        mpComm = rSparseGraph.pGetComm();
        TIndexType row_data_size=0;
        TIndexType col_data_size=0;
        rSparseGraph.ExportCSRArrays(mpRowIndicesData,row_data_size,mpColIndicesData, col_data_size);
        mRowIndices = Kratos::span<TIndexType>(mpRowIndicesData, row_data_size); //no copying of data happening here
        mColIndices = Kratos::span<TIndexType>(mpColIndicesData, col_data_size);

        mNrows = size1();

        ComputeColSize();

        //initialize mValuesVector to zero
        ResizeValueData(mColIndices.size());
        SetValue(0.0);
    }

    explicit CsrMatrix(const CsrMatrix<TDataType,TIndexType>& rOtherMatrix)
    {
        mpComm = rOtherMatrix.mpComm;
        ResizeIndex1Data(rOtherMatrix.mRowIndices.size());
        ResizeIndex2Data(rOtherMatrix.mColIndices.size());
        ResizeValueData(rOtherMatrix.mValuesVector.size());

        mNrows = rOtherMatrix.mNrows;
        mNcols = rOtherMatrix.mNcols;

        IndexPartition<IndexType>(mRowIndices.size()).for_each( [&](IndexType i){
            mRowIndices[i] = rOtherMatrix.mRowIndices[i];
        });

        IndexPartition<IndexType>(mColIndices.size()).for_each( [&](IndexType i){
            mColIndices[i] = rOtherMatrix.mColIndices[i];
        });

        IndexPartition<IndexType>(mValuesVector.size()).for_each( [&](IndexType i){
            mValuesVector[i] = rOtherMatrix.mValuesVector[i];
        });
    }

    //move constructor
    CsrMatrix(CsrMatrix<TDataType,TIndexType>&& rOtherMatrix)
    {
        mIsOwnerOfData=rOtherMatrix.mIsOwnerOfData;
        rOtherMatrix.mIsOwnerOfData=false;

        //swap the pointers to take owership of data
        mpRowIndicesData = rOtherMatrix.mpRowIndicesData;
        mpColIndicesData = rOtherMatrix.mpColIndicesData;
        mpValuesVectorData = rOtherMatrix.mpValuesVectorData;

        //here we assign the span
        mRowIndices = rOtherMatrix.mRowIndices;
        mColIndices = rOtherMatrix.mColIndices;
        mValuesVector = rOtherMatrix.mValuesVector;
        
        mNrows = rOtherMatrix.mNrows;
        mNcols = rOtherMatrix.mNcols;

    }

    /// Destructor.
    virtual ~CsrMatrix(){
        AssignIndex1Data(nullptr,0);
        AssignIndex2Data(nullptr,0);
        AssignValueData(nullptr,0);        
    }

    /// Assignment operator. 
    CsrMatrix& operator=(CsrMatrix const& rOtherMatrix) = delete; //i really think this should not be allowed, too risky
 
    //move assignement operator
    CsrMatrix& operator=(CsrMatrix&& rOtherMatrix)
    {
        mpComm = rOtherMatrix.mpComm;
        mIsOwnerOfData=rOtherMatrix.mIsOwnerOfData;
        rOtherMatrix.mIsOwnerOfData=false;

        //swap the pointers to take owership of data
        mpRowIndicesData = rOtherMatrix.mpRowIndicesData;
        mpColIndicesData = rOtherMatrix.mpColIndicesData;
        mpValuesVectorData = rOtherMatrix.mpValuesVectorData;

        //here we assign the span
        mRowIndices = rOtherMatrix.mRowIndices;
        mColIndices = rOtherMatrix.mColIndices;
        mValuesVector = rOtherMatrix.mValuesVector;
        
        mNrows = rOtherMatrix.mNrows;
        mNcols = rOtherMatrix.mNcols;
        return *this;
    }
    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        AssignIndex1Data(nullptr,0);
        AssignIndex2Data(nullptr,0);
        AssignValueData(nullptr,0);
        mNrows=0;
        mNcols=0;
    }

    const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    void SetValue(const TDataType value)
    {
        IndexPartition<IndexType>(mValuesVector.size()).for_each([&](IndexType i){
            mValuesVector[i] = value;
        });
    }

    IndexType size1() const
    {
        return mRowIndices.size()-1;
    }

    IndexType size2() const
    {
        return mNcols;
    }

    inline IndexType nnz() const{
        return index2_data().size();
    }

    bool IsOwnerOfData() const{
        return mIsOwnerOfData;
    }

    void SetIsOwnerOfData(bool IsOwner){
        mIsOwnerOfData = IsOwner;
    }

    inline Kratos::span<IndexType>& index1_data(){
        return mRowIndices;
    }
    inline Kratos::span<IndexType>& index2_data(){
        return mColIndices;
    }
    inline Kratos::span<TDataType>& value_data(){
        return mValuesVector;
    }

    inline const Kratos::span<IndexType>& index1_data() const{
        return mRowIndices;
    }
    inline const Kratos::span<IndexType>& index2_data() const{
        return mColIndices;
    }
    inline const Kratos::span<TDataType>& value_data() const{
        return mValuesVector;
    }

    void SetColSize(IndexType Ncols){
        mNcols = Ncols;
    }

    void SetRowSize(IndexType Nrows){
        mNrows = Nrows;
    }

    void ComputeColSize()
    {
        //compute the max id 
        IndexType max_col = IndexPartition<IndexType>(mColIndices.size()).template for_each< MaxReduction<double> >([&](IndexType i){
                return mColIndices[i];
            });

        mNcols = max_col+1; //note that we must add 1 to the greatest column id
    }

    void AssignIndex1Data(TIndexType* pExternalData, TIndexType DataSize){
        if(IsOwnerOfData() && mpRowIndicesData != nullptr)
            delete [] mpRowIndicesData;
        mpRowIndicesData = pExternalData;
        if(DataSize!=0)
            mRowIndices = Kratos::span<TIndexType>(mpRowIndicesData, DataSize);
        else
            mRowIndices = Kratos::span<TIndexType>();
    }

    void AssignIndex2Data(TIndexType* pExternalData, TIndexType DataSize){
        if(IsOwnerOfData() && mpColIndicesData != nullptr)
            delete [] mpColIndicesData;
        mpColIndicesData = pExternalData;
        if(DataSize!=0)
            mColIndices = Kratos::span<TIndexType>(mpColIndicesData, DataSize);
        else
            mColIndices = Kratos::span<TIndexType>();
    }

    void AssignValueData(TDataType* pExternalData, TIndexType DataSize){
        if(IsOwnerOfData() && mpValuesVectorData != nullptr)
            delete [] mpValuesVectorData;
        mpValuesVectorData = pExternalData;
        if(DataSize!=0)
            mValuesVector = Kratos::span<TDataType>(mpValuesVectorData, DataSize);
        else
            mValuesVector = Kratos::span<TDataType>();
    }

    void ResizeIndex1Data(TIndexType DataSize){
        KRATOS_ERROR_IF_NOT(IsOwnerOfData()) << "ResizeIndex1Data is only allowed if the data are locally owned" << std::endl;
        if(mpRowIndicesData != nullptr)
            delete [] mpRowIndicesData;
        mpRowIndicesData = new TIndexType[DataSize];
        mRowIndices = Kratos::span<TIndexType>(mpRowIndicesData, DataSize);
    }

    void ResizeIndex2Data(TIndexType DataSize){
        KRATOS_ERROR_IF_NOT(IsOwnerOfData()) << "ResizeIndex2Data is only allowed if the data are locally owned" << std::endl;
        if(mpColIndicesData != nullptr)
            delete [] mpColIndicesData;
        mpColIndicesData = new TIndexType[DataSize];
        mColIndices = Kratos::span<TIndexType>(mpColIndicesData, DataSize);
    }

    void ResizeValueData(TIndexType DataSize){
        KRATOS_ERROR_IF_NOT(IsOwnerOfData()) << "ResizeValueData is only allowed if the data are locally owned" << std::endl;
        if(mpValuesVectorData != nullptr)
            delete [] mpValuesVectorData;
        mpValuesVectorData = new TDataType[DataSize];
        mValuesVector = Kratos::span<TDataType>(mpValuesVectorData, DataSize);
    }



    void CheckColSize()
    {
        IndexType max_col = 0;
        for(IndexType i=0; i<mColIndices.size(); ++i)
            max_col = std::max(max_col,mColIndices[i]);
        if(max_col > mNcols)
            KRATOS_ERROR << " max column index : " << max_col << " exceeds mNcols :" << mNcols << std::endl;
    }

    TDataType& operator()(IndexType I, IndexType J){
        const IndexType row_begin = index1_data()[I];
        const IndexType row_end = index1_data()[I+1];
        IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
        KRATOS_DEBUG_ERROR_IF(k==std::numeric_limits<IndexType>::max()) << "local indices I,J : " << I << " " << J << " not found in matrix" << std::endl;
        return value_data()[k];
    }

    const TDataType& operator()(IndexType I, IndexType J) const{
        const IndexType row_begin = index1_data()[I];
        const IndexType row_end = index1_data()[I+1];
        IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
        KRATOS_DEBUG_ERROR_IF(k==std::numeric_limits<IndexType>::max()) << "local indices I,J : " << I << " " << J << " not found in matrix" << std::endl;
        return value_data()[k];
    }

    bool Has(IndexType I, IndexType J) const {
        const IndexType row_begin = index1_data()[I];
        const IndexType row_end = index1_data()[I+1];
        IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
        return k != std::numeric_limits<IndexType>::max();
    }

    // y += A*x  -- where A is *this
    template<class TInputVectorType, class TOutputVectorType>
    void SpMV(const TInputVectorType& x, TOutputVectorType& y) const
    {
        KRATOS_ERROR_IF(size1() != y.size() ) << "SpMV: mismatch between row sizes : " << size1()  << " and destination vector size " << y.size() << std::endl;
        KRATOS_ERROR_IF(size2() != x.size() ) << "SpmV: mismatch between col sizes : " << size2()  << " and input vector size " << x.size() << std::endl;
        if(nnz() != 0)
        {
            IndexPartition<IndexType>(y.size()).for_each( [&](IndexType i){
                IndexType row_begin = index1_data()[i];
                IndexType row_end   = index1_data()[i+1];
                for(IndexType k = row_begin; k < row_end; ++k){
                    IndexType col = index2_data()[k];
                    y(i) += value_data()[k] * x(col);
                }  
            });
        }
    }

    //y = alpha*y + beta*A*x
    template<class TInputVectorType, class TOutputVectorType>
    void SpMV(const TDataType alpha,
              const TInputVectorType& x, 
              const TDataType beta,
              TOutputVectorType& y) const
    {
        KRATOS_ERROR_IF(size1() != y.size() ) << "SpMV: mismatch between matrix sizes : " << size1() << " " <<size2() << " and destination vector size " << y.size() << std::endl;
        KRATOS_ERROR_IF(size2() != x.size() ) << "SpmV: mismatch between matrix sizes : " << size1() << " " <<size2() << " and input vector size " << x.size() << std::endl;
        IndexPartition<IndexType>(y.size()).for_each( [&](IndexType i){
            IndexType row_begin = index1_data()[i];
            IndexType row_end   = index1_data()[i+1];
            TDataType aux = TDataType();
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType col = index2_data()[k];
                aux += value_data()[k] * x(col);
            }  
            y(i) = beta*y(i) + alpha*aux;
        });
    }

    // y += A^t*x  -- where A is *this    
    template<class TInputVectorType, class TOutputVectorType>
    void TransposeSpMV(const TInputVectorType& x, TOutputVectorType& y) const
    {
        KRATOS_ERROR_IF(size2() != y.size() ) << "TransposeSpMV: mismatch between transpose matrix sizes : " << size2() << " " <<size1() << " and destination vector size " << y.size() << std::endl;
        KRATOS_ERROR_IF(size1() != x.size() ) << "TransposeSpMV: mismatch between transpose matrix sizes : " << size2() << " " <<size1() << " and input vector size " << x.size() << std::endl;
        IndexPartition<IndexType>(size1()).for_each( [&](IndexType i){
            IndexType row_begin = index1_data()[i];
            IndexType row_end   = index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType j = index2_data()[k];
                AtomicAdd(y(j), value_data()[k] * x(i) );
            }  
        });
    }

    //y = alpha*y + beta*A^t*x
    template<class TInputVectorType, class TOutputVectorType>
    void TransposeSpMV(const TDataType alpha,
              const TInputVectorType& x, 
              const TDataType beta,
              TOutputVectorType& y) const
    {
        KRATOS_ERROR_IF(size2() != y.size() ) << "TransposeSpMV: mismatch between transpose matrix sizes : " << size2() << " " <<size1() << " and destination vector size " << y.size() << std::endl;
        KRATOS_ERROR_IF(size1() != x.size() ) << "TransposeSpMV: mismatch between transpose matrix sizes : " << size2() << " " <<size1() << " and input vector size " << x.size() << std::endl;
        y *= beta;
        IndexPartition<IndexType>(size1()).for_each( [&](IndexType i)
        {
            IndexType row_begin = index1_data()[i];
            IndexType row_end   = index1_data()[i+1];
            TDataType aux = alpha*x(i);
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType j = index2_data()[k];
                AtomicAdd(y(j), value_data()[k] * x(i) );
            }  
        });
    }

    TDataType NormFrobenius() const
    {
        auto sum2 = IndexPartition<TIndexType>(this->value_data().size()).template for_each< SumReduction<TDataType> >( [this](TIndexType i){
                return std::pow(this->value_data()[i],2);
            });
        return std::sqrt(sum2);
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

    void reserve(IndexType NRows, IndexType nnz){
        ResizeIndex1Data(NRows+1);

        if(NRows > 0)
            index1_data()[0] = 0;

        ResizeIndex2Data(nnz);
        ResizeValueData(nnz);
        mNrows = NRows;
    }

    MatrixMapType ToMap() const
    {
        MatrixMapType value_map;
        for(unsigned int i=0; i<size1(); ++i)
        {
            IndexType row_begin = index1_data()[i];
            IndexType row_end   = index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType j = index2_data()[k];
                TDataType v = value_data()[k];
                value_map[{i,j}] = v;
            }  
        }
        return value_map;
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

            AtomicAdd(value_data()[k], rMatrixInput(i_local,0));

            //now find other entries. note that we assume that it is probably that next entries immediately follow in the ordering
            for(unsigned int j_local=1; j_local<local_size; ++j_local){
                J = EquationId[j_local];

                if(k+1<row_end && index2_data()[k+1] == J){
                    k = k+1;
                }
                else if(J > lastJ){ //note that the case k+2 >= index2_data().size() should be impossible
                    k = BinarySearch(index2_data(), k+2, row_end, J);
                }
                else if(J < lastJ){
                    k = BinarySearch(index2_data(), row_begin, k-1, J);
                }
                //the last missing case is J == lastJ, which should never happen in FEM. If that happens we can reuse k

                AtomicAdd(value_data()[k] , rMatrixInput(i_local,j_local));

                lastJ = J;
            }
        }
    }

    template<class TMatrixType, class TIndexVectorType >
    void Assemble(
        const TMatrixType& rMatrixInput,
        const TIndexVectorType& RowEquationId,
        const TIndexVectorType& ColEquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size1() != RowEquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size2() != ColEquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;

        unsigned int local_size = rMatrixInput.size1();
        unsigned int col_size = rMatrixInput.size2();

        for (unsigned int i_local = 0; i_local < local_size; ++i_local)
        {
            const IndexType I = RowEquationId[i_local];
            const IndexType row_begin = index1_data()[I];
            const IndexType row_end = index1_data()[I+1];

            //find first entry (note that we know it exists since local_size > 0)
            IndexType J = ColEquationId[0];
            IndexType k = BinarySearch(index2_data(), row_begin, row_end, J);
            IndexType lastJ = J;

            AtomicAdd(value_data()[k], rMatrixInput(i_local,0));

            //now find other entries. note that we assume that it is probably that next entries immediately follow in the ordering
            for(unsigned int j_local=1; j_local<col_size; ++j_local){
                J = ColEquationId[j_local];

                if(k+1<row_end && index2_data()[k+1] == J){
                    k = k+1;
                }
                else if(J > lastJ){ //note that the case k+2 >= index2_data().size() should be impossible
                    k = BinarySearch(index2_data(), k+2, row_end, J);
                }
                else if(J < lastJ){
                    k = BinarySearch(index2_data(), row_begin, k-1, J);
                }
                //the last case is J == lastJ, which should never happen in FEM. If that happens we can reuse k
                AtomicAdd(value_data()[k] , rMatrixInput(i_local,j_local));

                lastJ = J;
            }
        }
    }

    void AssembleEntry(const TDataType Value, const IndexType GlobalI, const IndexType GlobalJ)
    {
        IndexType k = BinarySearch(index2_data(), index1_data()[GlobalI], index1_data()[GlobalI+1], GlobalJ);
        AtomicAdd(value_data()[k], Value);
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
    virtual void PrintData(std::ostream& rOStream) const {
        rOStream << "size1 : " << size1() <<std::endl;
        rOStream << "size2 : " << size2() <<std::endl;
        rOStream << "nnz : " << nnz() <<std::endl; 
        rOStream << "index1_data : " << std::endl;
        for(auto item : index1_data())
            rOStream << item << ",";
        rOStream << std::endl;
        rOStream << "index2_data : " << std::endl;
        for(auto item : index2_data())
            rOStream << item << ",";
        rOStream << std::endl;        
        rOStream << "value_data  : " << std::endl;
        for(auto item : value_data())
            rOStream << item << ",";
        rOStream << std::endl;
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
    // A non-recursive binary search function. It returns 
    // location of x in given array arr[l..r] is present, 
    // otherwise std::dec << std::numeric_limits<IndexType>::max()
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
        return std::numeric_limits<IndexType>::max(); 
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
    const DataCommunicator* mpComm;
    bool mIsOwnerOfData = true;
    IndexType* mpRowIndicesData = nullptr;
    IndexType* mpColIndicesData = nullptr;
    TDataType* mpValuesVectorData = nullptr;
    Kratos::span<IndexType> mRowIndices;
    Kratos::span<IndexType> mColIndices;
    Kratos::span<TDataType> mValuesVector;
    IndexType mNrows=0;
    IndexType mNcols=0;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("IsOwnerOfData",mIsOwnerOfData);

        rSerializer.save("nrows",mRowIndices.size());
        for(IndexType i=0; i<mRowIndices.size(); ++i)
            rSerializer.save("i",mRowIndices[i]);

        rSerializer.save("cols_size",mColIndices.size());
        for(IndexType i=0; i<mColIndices.size(); ++i)
            rSerializer.save("i",mColIndices[i]);

        rSerializer.save("val_size",mValuesVector.size());
        for(IndexType i=0; i<mValuesVector.size(); ++i)
            rSerializer.save("d",mValuesVector[i]);
        
        rSerializer.save("Nrow",mNrows);
        rSerializer.save("Ncol",mNcols);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("IsOwnerOfData",mIsOwnerOfData);
        if(mIsOwnerOfData == false)
        {
            mIsOwnerOfData=true;
            KRATOS_WARNING("csr_matrix becomes owner of a copy of data after serialization");
        }

        IndexType rows_size;
        rSerializer.load("nrows",rows_size);
        mpRowIndicesData = new TIndexType[rows_size];
        mRowIndices = Kratos::span<TIndexType>(mpRowIndicesData, rows_size);
        for(IndexType i=0; i<rows_size; ++i)
            rSerializer.load("i",mRowIndices[i]);

        IndexType cols_size;
        rSerializer.load("cols_size",cols_size);
        mpColIndicesData = new TIndexType[cols_size];
        mColIndices = Kratos::span<TIndexType>(mpRowIndicesData, cols_size);
        for(IndexType i=0; i<mColIndices.size(); ++i)
            rSerializer.load("i",mColIndices[i]);

        IndexType vals_size;
        rSerializer.load("val_size",vals_size);
        mpValuesVectorData = new TDataType[vals_size];
        mValuesVector = Kratos::span<TDataType>(mpValuesVectorData, vals_size);
        for(IndexType i=0; i<mValuesVector.size(); ++i)
            rSerializer.load("d",mValuesVector[i]);
        
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
