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

#if !defined(KRATOS_DISTRIBUTED_SPARSE_GRAPH_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_SPARSE_GRAPH_H_INCLUDED


// System includes
#include <iostream>
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "includes/parallel_environment.h"
#include "utilities/communication_coloring_utilities.h"

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

class DistributedSparseGraph //: public SparseGraph
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType; //note that this could be different from the one in the basetype
    typedef SparseGraph LocalGraphType; //using a map since we need it ordered

    /// Pointer definition of DistributedSparseGraph
    KRATOS_CLASS_POINTER_DEFINITION(DistributedSparseGraph);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistributedSparseGraph(const std::vector<IndexType> limits,
                           DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator())
    : mLocalBounds(limits),
      mrComm(rComm),
      mLocalGraph() //limits[1]-limits[0]),
    {

        auto all_limits = mrComm.AllGather(mLocalBounds);

        //construct and sort a list containing the limits
        std::vector< std::array<IndexType,3> > tmp;
        IndexType cpu_id = 0;
        for(IndexType i=0; i<all_limits.size(); i+=2)
        {
            tmp.push_back({all_limits[i], all_limits[i+1], cpu_id++});
        }
        std::sort(tmp.begin(),
                  tmp.end(),
                  [](const std::array<IndexType,3>& a,
                     const std::array<IndexType,3>& b)
                     {return (a[0]<b[0]);}
                );

        //do consistency check
        for(unsigned int i=0; i<tmp.size()-1; ++i)
        {
            //check that no gaps are present
            if(tmp[i][1]!=tmp[i+1][0])
                KRATOS_ERROR << "upper bound of cpu" << i << " : " << tmp[i][1] <<
                            "is not consistent with the lower bound of cpu " <<
                            i+1 << " : " << tmp[i+1][0] << std::endl;

            //check that the first cpu has the lowest ids and progressively higher
            if(tmp[i][2] > tmp[i+1][2])
                KRATOS_ERROR << "cpu bounds are not ordered correctly" << std::endl;
        }
        mCpuBounds.push_back(tmp[0][0]);
        for(unsigned int i=0; i<tmp.size(); ++i)
        {
            //save to the final list
            mCpuBounds.push_back(tmp[i][1]);
        }

    }

    /// Destructor.
    virtual ~DistributedSparseGraph(){}

    bool IsLocal(const IndexType I)
    {
        return (I>=mLocalBounds[0] && I<mLocalBounds[1]);
    }

    IndexType LocalId(const IndexType rGlobalId) const
    {
        return rGlobalId-mLocalBounds[0];
    }

    IndexType GlobalId(const IndexType rLocalId) const
    {
        return rLocalId+mLocalBounds[0];
    }

    IndexType RemoteLocalId(const IndexType rGlobalId, const IndexType rOwnerRank) const
    {
        return rGlobalId-mCpuBounds[rOwnerRank];
    }

    IndexType OwnerRank(const IndexType RowIndex)
    {
        //NOTE: here we assume that CPUs with lower rank get the rows with lower rank
        //this leads to efficient checks but may well be too restrictive
        //TODO: decide if we want to relax this limitation

        //position of element just larger than RowIndex in mCpuBounds
        auto it = std::upper_bound(mCpuBounds.begin(), mCpuBounds.end(), RowIndex);

        KRATOS_DEBUG_ERROR_IF(it == mCpuBounds.end()) <<
            "row RowIndex " << RowIndex <<
            " is not owned by any processor " << std::endl;

        IndexType owner_rank = (it-mCpuBounds.begin()-1);

        KRATOS_DEBUG_ERROR_IF(owner_rank < 0) <<
            "row RowIndex " << RowIndex <<
            " is not owned by any processor " << std::endl;

        return owner_rank;

    }



    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        mLocalGraph.Clear();
        mNonLocalGraphs.clear();
    }

    void AddEntry(const IndexType RowIndex, const IndexType ColIndex)
    {
        if(IsLocal(RowIndex)){
            mLocalGraph.AddEntry(LocalId(RowIndex), ColIndex);
        }
        else{
            IndexType owner = OwnerRank(RowIndex);
            mNonLocalGraphs[owner].AddEntry(RemoteLocalId(RowIndex,owner), ColIndex);
        }
    }

    template<class TContainerType>
    void AddEntries(const IndexType RowIndex, const TContainerType& rColIndices)
    {
        if(IsLocal(RowIndex)){
            mLocalGraph.AddEntries(LocalId(RowIndex), rColIndices);
        }
        else{
            IndexType owner = OwnerRank(RowIndex);
            mNonLocalGraphs[owner].AddEntries(RemoteLocalId(RowIndex,owner), rColIndices);
        }
    }

    template<class TIteratorType>
    void AddEntries(const IndexType RowIndex,
                    const TIteratorType& rColBegin,
                    const TIteratorType& rColEnd
                    )
    {
        if(IsLocal(RowIndex)){
            mLocalGraph.AddEntries(LocalId(RowIndex), rColBegin, rColEnd);
        }
        else{
            IndexType owner = OwnerRank(RowIndex);
            mNonLocalGraphs[owner].AddEntries(RemoteLocalId(RowIndex,owner), rColBegin, rColEnd);
        }
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rIndices)
    {
        for(auto I : rIndices)
        {
            if(IsLocal(I)){
                mLocalGraph.AddEntries(LocalId(I), rIndices);
            }
            else{
                //TODO
                mNonLocalGraphs[OwnerRank(I)].AddEntries(I, rIndices);;
            }
        }
    }

    void AddEntries(DistributedSparseGraph& rOtherGraph)
    {
        //TODO
    }

    void AddEntries(const SparseGraph& rOtherGraph)
    {
        for(auto it = rOtherGraph.begin(); it!=rOtherGraph.end(); ++it)
        {
            AddEntries(it.GetRowIndex(), *it);
        }
    }

    void Finalize()
    {
        std::vector<int> send_list;
        for(auto& item : mNonLocalGraphs)
            send_list.push_back(item.first);

        auto colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, mrComm);

        //sendrecv data
        for(auto color : colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                //TODO: this can be made nonblocking
                auto recv_graph = mrComm.SendRecv(mNonLocalGraphs[color], color, color);

                // SparseGraph recv_graph;
                mLocalGraph.AddEntries(recv_graph);
            }
        }
    }

    const LocalGraphType& GetLocalGraph() const{
        return mLocalGraph;
    }

    const LocalGraphType& GetNonLocalGraph(IndexType Rank) const{
        return mNonLocalGraphs.at(Rank);
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DistributedSparseGraph" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedSparseGraph";}

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
    std::vector<IndexType> mLocalBounds;
    DataCommunicator& mrComm;

    std::vector<IndexType> mCpuBounds;
    LocalGraphType mLocalGraph;
    std::unordered_map<IndexType,LocalGraphType> mNonLocalGraphs;

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
inline std::istream& operator >> (std::istream& rIStream,
                DistributedSparseGraph& rThis){
                    return rIStream;
                }

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedSparseGraph& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_SPARSE_GRAPH_H_INCLUDED  defined


