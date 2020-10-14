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

#if !defined(KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED


// System includes
#include <iostream>
#include "containers/csr_matrix.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_vector_importer.h"
#include "includes/key_hash.h"
#include "utilities/atomic_utilities.h"


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
class DistributedCsrMatrix
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;
    typedef typename CsrMatrix<double,IndexType>::MatrixMapType MatrixMapType;

    /// Pointer definition of DistributedCsrMatrix
    KRATOS_CLASS_POINTER_DEFINITION(DistributedCsrMatrix);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistributedCsrMatrix(const DistributedSparseGraph& rSparseGraph)
        :
            mrComm(rSparseGraph.GetComm())
    {
        mpRowNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >( rSparseGraph.GetRowNumbering());

        //compute the columns size
        IndexType max_col_index = rSparseGraph.ComputeMaxGlobalColumnIndex();
        IndexType tot_col_size = max_col_index+1;
        mpColNumbering = Kratos::make_unique<DistributedNumbering<IndexType>>(mrComm, tot_col_size, mrComm.Size()); 

        mOffDiagonalLocalIds.clear(); //this is the map that allows to transform from global_ids to local Ids for entries in the non_diag block

        //count entries in diagonal and off_diag blocks
        const auto& local_graph = rSparseGraph.GetLocalGraph();
        const IndexType nlocal_rows = local_graph.Size();

        IndexType diag_nnz = 0, offdiag_nnz=0;
        for(const auto& entries : local_graph){
            for(const auto& global_j : entries){
                if(GetColNumbering().IsLocal(global_j)){ 
                    diag_nnz += 1;
                }
                else{
                    offdiag_nnz += 1;
                    mOffDiagonalLocalIds[global_j] = 0;
                }
            }
        }

        //begin by computing local Ids for non diagonal block
        IndexType counter = 0;
        for(auto& item : mOffDiagonalLocalIds)
            item.second = counter++;

        mOffDiagonalGlobalIds.resize(mOffDiagonalLocalIds.size(), false);
        counter = 0;
        for(auto& item : mOffDiagonalLocalIds)
        {
            IndexType global_i = item.first;
            mOffDiagonalGlobalIds[counter++] = global_i;
        }

        

        //*******************************
        //construct diagonal block
        mDiagBlock.reserve(nlocal_rows, diag_nnz);
        counter = 0;
        for(const auto& entries : local_graph){
            unsigned int k = 0;
            IndexType row_begin = mDiagBlock.index1_data()[counter];
            for(auto global_j : entries){
                
                if(GetColNumbering().IsLocal(global_j)){
                    IndexType local_j = GetColNumbering().LocalId(global_j);
                    mDiagBlock.index2_data()[row_begin+k] = local_j;
                    mDiagBlock.value_data()[row_begin+k] = 0.0;
                    k++;
                }
            }
            mDiagBlock.index1_data()[counter+1] = row_begin + k;
            counter++;
        }
        mDiagBlock.SetColSize(GetColNumbering().LocalSize());
        
#ifdef KRATOS_DEBUG
        mDiagBlock.CheckColSize();
#endif

        //ensure columns are ordered
        for(IndexType i = 0; i<nlocal_rows; ++i){
            IndexType row_begin = mDiagBlock.index1_data()[i];
            IndexType row_end = mDiagBlock.index1_data()[i+1];
            std::sort(mDiagBlock.index2_data().begin() + row_begin ,mDiagBlock.index2_data().begin()+row_end); 
        }


        //*******************************
        //construct offdiagonal block
        
        //store off-diagonal block
        mOffDiagBlock.reserve(nlocal_rows, offdiag_nnz); //TODO: merge with the previous loop
        counter = 0;
        for(const auto& entries : local_graph){
            unsigned int k = 0;
            IndexType row_begin = mOffDiagBlock.index1_data()[counter];
            for(auto global_j : entries){
                if( ! GetColNumbering().IsLocal(global_j)){
                    IndexType local_j = GetOffDiagLocalId(global_j);
                    mOffDiagBlock.index2_data()[row_begin+k] = local_j;
                    mOffDiagBlock.value_data()[row_begin+k] = 0.0;
                    k++;
                }
            }
            mOffDiagBlock.index1_data()[counter+1] = row_begin + k;
            counter++;
        }
        mOffDiagBlock.SetColSize(mOffDiagonalLocalIds.size());

#ifdef KRATOS_DEBUG
        mOffDiagBlock.CheckColSize();
#endif

        //ensure columns are ordered
        for(IndexType i = 0; i<nlocal_rows; ++i){
            IndexType row_begin = mOffDiagBlock.index1_data()[i];
            IndexType row_end = mOffDiagBlock.index1_data()[i+1];
            std::sort(mOffDiagBlock.index2_data().begin() + row_begin ,mOffDiagBlock.index2_data().begin()+row_end); 
        }

        PrepareNonLocalCommunications(rSparseGraph);

        //mount importer for SpMV calculations

        auto pimporter = Kratos::make_unique<DistributedVectorImporter<double,IndexType>>(GetComm(),mOffDiagonalGlobalIds, GetColNumbering()); 
        mpVectorImporter.swap(pimporter);
    }

    /// Destructor.
    virtual ~DistributedCsrMatrix(){}

    /// Assignment operator. TODO: decide if we do want to allow it
    DistributedCsrMatrix& operator=(DistributedCsrMatrix const& rOther)=delete;
    // {
    //     this->AddEntries(rOther.GetGraph());
    //     return *this;
    // }

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {

    }

    inline const DistributedNumbering<IndexType>& GetRowNumbering() const
    {
        return *mpRowNumbering;
    }
    
    inline const DistributedNumbering<IndexType>& GetColNumbering() const
    {
        return *mpColNumbering;
    }

    void SetValue(const TDataType value)
    {
        mDiagBlock.SetValue(value);
        mOffDiagBlock.SetValue(value);
    }

    IndexType local_size1() const
    {
        return mDiagBlock.size1();
        //TODO decide if we should give back the local or globale sizes
    }

    IndexType size2() const
    {
        return -1; 
        //TODO decide if we should give back the local or globale sizes
    }

    inline IndexType local_nnz() const{
        return mDiagBlock.nnz();
    }

    const DataCommunicator& GetComm(){
        return mrComm;
    }

    inline CsrMatrix<double,IndexType>& GetDiagBlock(){
        return mDiagBlock;
    }
    inline CsrMatrix<double,IndexType>& GetOffDiagBlock(){
        return mOffDiagBlock;
    }

    inline const CsrMatrix<double,IndexType>& GetDiagBlock() const{
        return mDiagBlock;
    }
    inline const CsrMatrix<double,IndexType>& GetOffDiagBlock() const{
        return mOffDiagBlock;
    }



    

    IndexType GetOffDiagLocalId(IndexType GlobalJ) const{
        auto it = mOffDiagonalLocalIds.find(GlobalJ);
        KRATOS_DEBUG_ERROR_IF( it == mOffDiagonalLocalIds.end() ) << "GlobalJ is not in the nonlocal list" << std::endl;
        return it->second;
    }

    TDataType& GetLocalDataByGlobalId(IndexType GlobalI, IndexType GlobalJ){
        KRATOS_DEBUG_ERROR_IF(  ! GetRowNumbering().IsLocal(GlobalI) ) << "non local row access for GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << std::endl;

        IndexType LocalI = GetRowNumbering().LocalId(GlobalI);
        if(GetColNumbering().IsLocal(GlobalJ)){ 
            return mDiagBlock( LocalI, GetColNumbering().LocalId(GlobalJ) );
        }
        else{
            return mOffDiagBlock( LocalI, GetOffDiagLocalId(GlobalJ) );
        }
    }

    TDataType& GetNonLocalDataByGlobalId(IndexType GlobalI, IndexType GlobalJ){
        KRATOS_DEBUG_ERROR_IF(  GetRowNumbering().IsLocal(GlobalI) ) << " local row access for GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << " expected to be nonlocal" << std::endl;
        auto it = mNonLocalData.find(std::make_pair(GlobalI,GlobalJ));
        KRATOS_DEBUG_ERROR_IF(it == mNonLocalData.end()) << " entry GlobalI,GlobalJ = " << GlobalI << " " << GlobalJ << " not found in NonLocalData" << std::endl;
        return it->second;
    }


    // TDataType& operator()(IndexType I, IndexType J){
    // }

    ///@}
    ///@name Operations
    ///@{
    //TODO
    //+=
    //-=
    //+
    //-
    //*

   void SpMV(DistributedSystemVector<TDataType,TIndexType>& y,
              const DistributedSystemVector<TDataType,TIndexType>& x) const
    {
        //get off diagonal terms (requires communication)
        auto off_diag_x = mpVectorImporter->ImportData(x);
        mOffDiagBlock.SpMV(y.GetLocalData(),off_diag_x);
        mDiagBlock.SpMV(y.GetLocalData(),x.GetLocalData());
    }




    void BeginAssemble(){
        //set to zero non local data 
        for(auto& item : mNonLocalData)
            item.second = 0.0;

    } 

    void FinalizeAssemble(){
        //communicate data to finalize the assembly
        auto& rComm = GetComm();

        //sendrecv data
        for(auto color : mfem_assemble_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& direct_senddata_access = mPointersToSendValues[color];
                const auto& direct_recvdata_access = mPointersToRecvValues[color];

                auto& send_data = msend_buffers[color];
                
                for(IndexType i=0; i<send_data.size(); ++i)
                {
                    send_data[i] = *(direct_senddata_access[i]);

                    //slower but more understandable version TOTO remove
                    // IndexType I = mIndicesToSendValues[color][i].first;
                    // IndexType J = mIndicesToSendValues[color][i].second;
                    // send_data[i] = GetNonLocalDataByGlobalId(I,J);
                    // if(I == 33 && J == 33) std::cout << "send entry 33,33 :" << send_data[i] << std::endl;


                }
                //NOTE: this can be made nonblocking 
                auto& recv_data = mrecv_buffers[color];
                rComm.SendRecv(send_data, color, 0, recv_data, color, 0); //TODO, we know all the sizes, we shall use that!

                for(IndexType i=0; i<recv_data.size(); ++i)
                {
                    *(direct_recvdata_access[i]) += recv_data[i]; //here we assemble the nonlocal contribution to the local data

                    //slower but more understandable version TODO remove
                    // IndexType I = mIndicesToRecvValues[color][i].first;
                    // IndexType J = mIndicesToRecvValues[color][i].second;
                    // if(I == 33 && J == 33) std::cout << "recv entry 33,33 :" << recv_data[i] << std::endl;
                    // GetLocalDataByGlobalId(I,J) += recv_data[i];

                    
                }
            }
        }


    } 





    template<class TMatrixType, class TIndexVectorType >
    void Assemble(
        const TMatrixType& rMatrixInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size1() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;
        KRATOS_DEBUG_ERROR_IF(rMatrixInput.size2() != EquationId.size()) << "sizes of matrix and equation id do not match in Assemble" << std::endl;

        for(unsigned int i=0; i<EquationId.size(); ++i){
            const IndexType global_i = EquationId[i];
            if(GetRowNumbering().IsLocal(global_i))
            {
               for(unsigned int j = 0; j<EquationId.size(); ++j)
                {
                    const IndexType global_j = EquationId[j];
                    TDataType& value = GetLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }
            }
            else
            {
                for(unsigned int j = 0; j<EquationId.size(); ++j)
                {
                    const IndexType global_j = EquationId[j];
                    TDataType& value = GetNonLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }                

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

        for(unsigned int i=0; i<RowEquationId.size(); ++i){
            const IndexType global_i = RowEquationId[i];

            if(GetRowNumbering().IsLocal(global_i))
            {
               for(unsigned int j = 0; j<ColEquationId.size(); ++j)
                {
                    const IndexType global_j = ColEquationId[j];
                    TDataType& value = GetLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }
            }
            else
            {
                for(unsigned int j = 0; j<ColEquationId.size(); ++j)
                {
                    const IndexType global_j = ColEquationId[j];
                    TDataType& value = GetNonLocalDataByGlobalId(global_i,global_j);
                    AtomicAdd(value, rMatrixInput(i,j));
                }                

            }
        }
    }

    MatrixMapType ToMap() const
    {
        MatrixMapType value_map;
        for(unsigned int i=0; i<local_size1(); ++i)
        {
            IndexType row_begin = mDiagBlock.index1_data()[i];
            IndexType row_end   = mDiagBlock.index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType j = mDiagBlock.index2_data()[k];
                TDataType v = mDiagBlock.value_data()[k];
                value_map[{GetRowNumbering().GlobalId(i),GetColNumbering().GlobalId(j)}] = v;
            }  
        }

        for(unsigned int i=0; i<local_size1(); ++i)
        {
            IndexType row_begin = mOffDiagBlock.index1_data()[i];
            IndexType row_end   = mOffDiagBlock.index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType j = mOffDiagBlock.index2_data()[k];
                TDataType v = mOffDiagBlock.value_data()[k];
                value_map[{GetRowNumbering().GlobalId(i),mOffDiagonalGlobalIds[j]}] = v;
            }  
        }

        return value_map;
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
        buffer << "DistributedCsrMatrix" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedCsrMatrix";}

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
                        IndexType l, IndexType r, IndexType x) 
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
    void PrepareNonLocalCommunications(const DistributedSparseGraph& rSparseGraph)
    {        
        auto& rComm = GetComm();

        // mPointersToReceivedValues.resize(world_size);
        // mPointersToReceivedValues.resize(world_size);

        const auto& nonlocal_graphs = rSparseGraph.GetNonLocalGraphs();
        std::vector<int> send_list;
        for(unsigned int id = 0; id<nonlocal_graphs.size(); ++id)
                if( !nonlocal_graphs[id].IsEmpty())
                    send_list.push_back(id);

        mfem_assemble_colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, rComm);

        //sendrecv data
        for(auto color : mfem_assemble_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& send_graph = nonlocal_graphs[color];
                auto& direct_senddata_access = mPointersToSendValues[color];
                std::vector<IndexType> send_ij;

                for(auto row_it=send_graph.begin(); row_it!=send_graph.end(); ++row_it)
                {
                    const auto remote_local_I = row_it.GetRowIndex();
                    const auto remote_global_I = GetRowNumbering().RemoteGlobalId(remote_local_I, color);

                    for(auto J : *row_it){
                        TDataType& value = mNonLocalData[std::make_pair(remote_global_I,J)]; //here we create the I,J entry in the nonlocal data (entry was there in the graph!)
                        direct_senddata_access.push_back(&value); //storing a direct pointer to the value contained in the data structure
                    //    indices_senddata_access.push_back(std::make_pair(remote_global_I,J)); //TODO: remove, for debug

                        send_ij.push_back(remote_global_I);
                        send_ij.push_back(J);
                    }
                }
                
                //NOTE: this can be made nonblocking 
                const auto recv_ij = rComm.SendRecv(send_ij, color, color);

                auto& direct_recvdata_access = mPointersToRecvValues[color];

                for(IndexType k=0; k<recv_ij.size(); k+=2)
                {
                    IndexType I = recv_ij[k];
                    IndexType J = recv_ij[k+1];
                    auto& value = GetLocalDataByGlobalId(I,J);
                    //      indices_recvdata_access.push_back(std::make_pair(I,J)); //TODO: remove, for debug
                    direct_recvdata_access.push_back(&value);
                }

                //resizing buffers to be later used for sending and receiving
                msend_buffers[color].resize(direct_senddata_access.size());
                mrecv_buffers[color].resize(direct_recvdata_access.size());

            }
        }
    }




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
    const DataCommunicator& mrComm;

    typename DistributedNumbering<IndexType>::UniquePointer mpRowNumbering;
    typename DistributedNumbering<IndexType>::UniquePointer mpColNumbering;

    CsrMatrix<double,IndexType> mDiagBlock;
    CsrMatrix<double,IndexType> mOffDiagBlock;
    MatrixMapType mNonLocalData; //data which is assembled locally and needs to be communicated to the owner

    //this map tells for an index J which does not belong to the local diagonal block which is the corresponding localJ
    std::map<IndexType, IndexType> mOffDiagonalLocalIds; //usage: mOffDiagonalLocalIds[global_id] contains the local_id associated to that global_id (for a off diagonal block entry)
    DenseVector<IndexType> mOffDiagonalGlobalIds; //usage: mOffDiagonalGlobalIds[local_id] contains the global_id associated

    std::vector<int> mfem_assemble_colors; //coloring of communication
    std::unordered_map< unsigned int, std::vector<TDataType*> > mPointersToRecvValues; //this contains direct pointers into the data contained in mNonLocalData, prepared so to speed up communications
    std::unordered_map< unsigned int, std::vector<TDataType*> > mPointersToSendValues; //this contains direct pointers into mDiagBlock and mOffDiagBlock, prepared so to speed up communication
    
    // std::unordered_map< unsigned int, std::vector<std::pair<IndexType,IndexType>>> mIndicesToSendValues; //TODO remove, for debug
    // std::unordered_map< unsigned int, std::vector<std::pair<IndexType,IndexType>>> mIndicesToRecvValues; //TODO remove, for debug

    std::unordered_map< unsigned int, std::vector<TDataType> > msend_buffers;
    std::unordered_map< unsigned int, std::vector<TDataType> > mrecv_buffers;

    std::unique_ptr<DistributedVectorImporter<double,IndexType>> mpVectorImporter;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        //TODO
    }

    void load(Serializer& rSerializer)
    {
        //TODO
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
    DistributedCsrMatrix(){};



    ///@}

}; // Class DistributedCsrMatrix

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDataType>
inline std::istream& operator >> (std::istream& rIStream,
                DistributedCsrMatrix<TDataType>& rThis){
                    return rIStream;
                }

/// output stream function
template< class TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedCsrMatrix<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_CSR_MATRIX_H_INCLUDED  defined


