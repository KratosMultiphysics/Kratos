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

#pragma once

// System includes
#include <iostream>
#include <mutex>
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "includes/parallel_environment.h"
#include "containers/distributed_numbering.h"
#include "utilities/communication_coloring_utilities.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "utilities/parallel_utilities.h"

// External includes
#include <unordered_map>
#include <unordered_set>

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

/// Short class definition.
//class to construct and store a matrix graph. Can be used to construct efficiently a CSR matrix (or other sparse matrix types)

/** This class is designed to store a matrix graph, aimed at the fast construction of other
 * sparse matrix formats (particularly CSR)
 * IMPORTANT NOTE: it is BY DESIGN NOT threadsafe! (a graph should be computed in each thread and then merged)
*/
template< class TIndexType=std::size_t >
class DistributedSparseGraph final
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;
    typedef SparseContiguousRowGraph<IndexType> LocalGraphType; //using a map since we need it ordered
    typedef SparseGraph<IndexType> NonLocalGraphType; //using a map since we need it ordered

    /// Pointer definition of DistributedSparseGraph
    KRATOS_CLASS_POINTER_DEFINITION(DistributedSparseGraph);

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor.
    DistributedSparseGraph(const IndexType LocalSize,
                           DataCommunicator& rComm)
    :
      mpComm(&rComm),
      mLocalGraph(LocalSize)
    {
        mNonLocalGraphs.resize(mpComm->Size(),false);
        mNonLocalLocks = decltype(mNonLocalLocks)(mpComm->Size());

        mpRowNumbering = Kratos::make_unique<DistributedNumbering<IndexType>>(*mpComm,LocalSize);
    }


    /// Destructor.
    ~DistributedSparseGraph() {}

    inline const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    inline const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    inline const DistributedNumbering<IndexType>& GetRowNumbering() const
    {
        return *mpRowNumbering;
    }

    inline IndexType Size() const
    {
        return mpRowNumbering->Size();
    }

    inline IndexType LocalSize() const
    {
        return mpRowNumbering->LocalSize();
    }

    bool Has(const IndexType GlobalI, const IndexType GlobalJ) const
    {
        return mLocalGraph.Has(GetRowNumbering().LocalId(GlobalI),GlobalJ);
    }

    //this function detects the maximum and minimum globalJ found in the local graph
    void ComputeLocalMinMaxColumnIndex(IndexType& rMinJ, IndexType& rMaxJ) const
    {
        rMaxJ = 0;
        rMinJ = 0;
        for(IndexType local_i = 0; local_i<mLocalGraph.Size(); ++local_i)
        {
            for(auto J : mLocalGraph[local_i] ) //J here is the global index
            {
                rMaxJ = std::max(rMaxJ, J);
                rMinJ = std::min(rMinJ, J);
            }
        }
    }

    IndexType ComputeMaxGlobalColumnIndex() const
    {
        IndexType MinJ, MaxJ;
        ComputeLocalMinMaxColumnIndex(MinJ,MaxJ);
        return GetComm().MaxAll(MaxJ);
    }


    ///@}
    ///@name Operators
    ///@{
    const typename LocalGraphType::GraphType::value_type& operator[](const IndexType& LocalPosition) const
    {
        return mLocalGraph[LocalPosition];
    }

    void Clear()
    {
        mLocalGraph.Clear();
        mNonLocalGraphs.clear();
    }

    void AddEntry(const IndexType RowIndex, const IndexType ColIndex)
    {
        if(GetRowNumbering().IsLocal(RowIndex))
        {
            mLocalGraph.AddEntry(GetRowNumbering().LocalId(RowIndex), ColIndex);
        }
        else
        {
            IndexType owner = GetRowNumbering().OwnerRank(RowIndex);
            const std::lock_guard<LockObject> scope_lock(mNonLocalLocks[owner]);
            mNonLocalGraphs[owner].AddEntry(GetRowNumbering().RemoteLocalId(RowIndex,owner), ColIndex);
        }
    }

    template<class TContainerType>
    void AddEntries(const IndexType RowIndex, const TContainerType& rColIndices)
    {
        if(GetRowNumbering().IsLocal(RowIndex))
        {
            mLocalGraph.AddEntries(GetRowNumbering().LocalId(RowIndex), rColIndices);
        }
        else
        {
            IndexType owner = GetRowNumbering().OwnerRank(RowIndex);
            mNonLocalLocks[owner].lock();
            mNonLocalGraphs[owner].AddEntries(GetRowNumbering().RemoteLocalId(RowIndex,owner), rColIndices);
            mNonLocalLocks[owner].unlock();
        }
    }

    template<class TIteratorType>
    void AddEntries(const IndexType RowIndex,
                    const TIteratorType& rColBegin,
                    const TIteratorType& rColEnd
                   )
    {
        if(GetRowNumbering().IsLocal(RowIndex))
        {
            mLocalGraph.AddEntries(GetRowNumbering().LocalId(RowIndex), rColBegin, rColEnd);
        }
        else
        {
            IndexType owner = GetRowNumbering().OwnerRank(RowIndex);
            mNonLocalLocks[owner].lock();
            mNonLocalGraphs[owner].AddEntries(GetRowNumbering().RemoteLocalId(RowIndex,owner), rColBegin, rColEnd);
            mNonLocalLocks[owner].unlock();
        }
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rIndices)
    {
        for(auto I : rIndices)
        {
            if(GetRowNumbering().IsLocal(I))
            {
                mLocalGraph.AddEntries(GetRowNumbering().LocalId(I), rIndices);
            }
            else
            {
                IndexType owner = GetRowNumbering().OwnerRank(I);
                mNonLocalLocks[owner].lock();
                mNonLocalGraphs[owner].AddEntries(GetRowNumbering().RemoteLocalId(I,owner), rIndices);;
                mNonLocalLocks[owner].unlock();
            }
        }
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rRowIndices, const TContainerType& rColIndices)
    {
        for(auto I : rRowIndices)
        {
            AddEntries(I, rColIndices);
        }
    }

    // void AddEntries(DistributedSparseGraph& rOtherGraph)
    // {
    //     //TODO
    // }

    void Finalize()
    {
        std::vector<MpiIndexType> send_list;

        for(unsigned int id = 0; id<mNonLocalGraphs.size(); ++id)
            if( !mNonLocalGraphs[id].IsEmpty())
                send_list.push_back(id);

        auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, *mpComm);

        //sendrecv data
        for(auto color : colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                //TODO: this can be made nonblocking

                //using serialization
                //const auto recv_graph = mpComm->SendRecv(mNonLocalGraphs[color], color, color);
                // for(auto row_it=recv_graph.begin(); row_it!=recv_graph.end(); ++row_it)
                // {
                //     auto I = row_it.GetRowIndex();
                //     mLocalGraph.AddEntries(I,*row_it);
                // }

                //using native calls
                auto send_single_vector_repr = mNonLocalGraphs[color].ExportSingleVectorRepresentation();
                const auto recv_single_vector_repr = mpComm->SendRecv(send_single_vector_repr, color, color);
                mLocalGraph.AddFromSingleVectorRepresentation(recv_single_vector_repr);

            }
        }

    }

    const LocalGraphType& GetLocalGraph() const
    {
        return mLocalGraph;
    }

    const NonLocalGraphType& GetNonLocalGraph(IndexType Rank) const
    {
        return mNonLocalGraphs[Rank];
    }

    const DenseVector<NonLocalGraphType>& GetNonLocalGraphs() const
    {
        return mNonLocalGraphs;
    }


    ///@}
    ///@name Operations
    ///@{


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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DistributedSparseGraph" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DistributedSparseGraph";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

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
    typename DistributedNumbering<IndexType>::UniquePointer mpRowNumbering = nullptr;
    const DataCommunicator* mpComm;

    LocalGraphType mLocalGraph;
    DenseVector<NonLocalGraphType> mNonLocalGraphs;
    std::vector<LockObject> mNonLocalLocks;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
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

    /// Assignment operator.
    DistributedSparseGraph& operator=(DistributedSparseGraph const& rOther) = delete;

    /// Copy constructor.
    DistributedSparseGraph(DistributedSparseGraph const& rOther) = delete;

    ///@}

}; // Class DistributedSparseGraph

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                                  DistributedSparseGraph<TIndexType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DistributedSparseGraph<TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
