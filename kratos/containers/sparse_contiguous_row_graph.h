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

#if !defined(KRATOS_SPARSE_CONTIGUOUS_ROW_GRAPH_H_INCLUDED )
#define  KRATOS_SPARSE_CONTIGUOUS_ROW_GRAPH_H_INCLUDED


// System includes
#include <iostream>
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "includes/lock_object.h"
#include "utilities/parallel_utilities.h"

// External includes
#include <unordered_map>
#include <unordered_set>

// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup KratosCore
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
 * IMPORTANT NOTE: AddEntries IS threasafe
*/

template<typename TIndexType=std::size_t>
class SparseContiguousRowGraph
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef DenseVector<std::unordered_set<IndexType> > GraphType; 
    typedef typename GraphType::const_iterator const_row_iterator;

    /// Pointer definition of SparseContiguousRowGraph
    KRATOS_CLASS_POINTER_DEFINITION(SparseContiguousRowGraph);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. - needs to be public for communicator, but it will fail if used in any other mode
    SparseContiguousRowGraph()
    {
    }

    SparseContiguousRowGraph(IndexType GraphSize)
    {
        mGraph.resize(GraphSize,false);
        mLocks.resize(GraphSize);

        // @RiccardoRossi why is this needed? isn't this done with resizing?
        //doing first touching
        IndexPartition<IndexType>(GraphSize).for_each([&](IndexType i){
            // mLocks[i] = LockObject();
            mGraph[i] = std::unordered_set<IndexType>();
        });
    }

    /// Destructor.
    virtual ~SparseContiguousRowGraph(){}

    /// Assignment operator. TODO: decide if we do want to allow it
    SparseContiguousRowGraph& operator=(SparseContiguousRowGraph const& rOther)=delete;
    // {
    //     this->AddEntries(rOther.GetGraph());
    //     return *this;
    // }

    /// Copy constructor. 
    SparseContiguousRowGraph(const SparseContiguousRowGraph& rOther)
    {
        this->AddEntries(rOther);
    }

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        mGraph.clear();
        mLocks.clear();
    }

    inline IndexType Size() const{
        return mGraph.size();
    }

    bool Has(const IndexType I, const IndexType J) const
    {
        const auto& row = mGraph[I];
        return (row.find(J) != row.end());
    }

    const typename GraphType::value_type& operator[](const IndexType& Key) const
    {
		return mGraph[Key];
    }

    inline IndexType LocalIndex(const IndexType GlobalIndex) const{
        return GlobalIndex;
    }

    inline IndexType GlobalIndex(const IndexType LocalIndex) const{
         return LocalIndex;
    }

    void AddEntry(const IndexType RowIndex, const IndexType ColIndex)
    {
        mLocks[RowIndex].SetLock();
        mGraph[RowIndex].insert(ColIndex);
        mLocks[RowIndex].UnSetLock();
    }

    template<class TContainerType>
    void AddEntries(const IndexType RowIndex, const TContainerType& rColIndices)
    {
        mLocks[RowIndex].SetLock();
        mGraph[RowIndex].insert(rColIndices.begin(), rColIndices.end());
        mLocks[RowIndex].UnSetLock();
    }

    template<class TIteratorType>
    void AddEntries(const IndexType RowIndex,
                    const TIteratorType& rColBegin,
                    const TIteratorType& rColEnd
                    )
    {
        mLocks[RowIndex].SetLock();
        mGraph[RowIndex].insert(rColBegin, rColEnd);
        mLocks[RowIndex].UnSetLock();
    }

    //adds a square FEM matrix, identified by rIndices
    template<class TContainerType>
    void AddEntries(const TContainerType& rIndices)
    {
        for(auto I : rIndices){
            KRATOS_DEBUG_ERROR_IF(I > this->Size()) << "Index : " << I
                << " exceeds the graph size : " << Size() << std::endl;
            mLocks[I].SetLock();
            mGraph[I].insert(rIndices.begin(), rIndices.end());
            mLocks[I].UnSetLock();
        }
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rRowIndices, const TContainerType& rColIndices)
    {
        for(auto I : rRowIndices){
            AddEntries(I, rColIndices);
        }
    }


    void AddEntries(const SparseContiguousRowGraph& rOtherGraph)
    {
        for(IndexType i=0; i<rOtherGraph.Size(); ++i)
        {
            AddEntries(i, rOtherGraph.GetGraph()[i]);
        }
    }

    void Finalize()
    {

    }

    const GraphType& GetGraph() const{
        return mGraph;
    }

    template<class TVectorType=DenseVector<IndexType>>
    IndexType ExportCSRArrays(
        TVectorType& rRowIndices,
        TVectorType& rColIndices
    ) const
    {
        //need to detect the number of rows this way since there may be gaps
        IndexType nrows=Size();

        if(rRowIndices.size() != nrows+1)
        {
            rRowIndices.resize(nrows+1, false);
        }
        //set it to zero in parallel to allow first touching
        IndexPartition<IndexType>(nrows+1).for_each([&](IndexType i){
            rRowIndices[i] = 0;
        });

        //count the entries 
        IndexPartition<IndexType>(nrows).for_each([&](IndexType i){
            rRowIndices[i+1] = mGraph[i].size();
        });
        
        //sum entries
        for(IndexType i = 1; i<static_cast<IndexType>(rRowIndices.size()); ++i){
            rRowIndices[i] += rRowIndices[i-1];
        }

        IndexType nnz = rRowIndices[nrows];
        if(rColIndices.size() != nnz){
            rColIndices.resize(nnz, false);
        }
        //set it to zero in parallel to allow first touching
        IndexPartition<IndexType>(rColIndices.size()).for_each([&](IndexType i){
            rColIndices[i] = 0;
        });

        //count the entries 
        IndexPartition<IndexType>(nrows).for_each([&](IndexType i){

            IndexType start = rRowIndices[i];

            IndexType counter = 0;
            for(auto index : mGraph[i]){
                rColIndices[start+counter] = index;
                counter++;
            }
        });

        //reorder columns
        IndexPartition<IndexType>(rRowIndices.size()-1).for_each([&](IndexType i){
            std::sort(rColIndices.begin()+rRowIndices[i], rColIndices.begin()+rRowIndices[i+1]);
        });

        return nrows;
    }

    //this function returns the Graph as a single vector
    //in the form of 
    //  RowIndex NumberOfEntriesInTheRow .... list of all Indices in the row 
    //every row is pushed back one after the other 
    std::vector<IndexType> ExportSingleVectorRepresentation()
    {
        std::vector< IndexType > IJ;
        IJ.push_back(GetGraph().size()); //number of rows
        for(unsigned int I=0; I<GetGraph().size(); ++I)
        {
            IJ.push_back(I); //id of the row
            IJ.push_back(mGraph[I].size()); //number of Js in the row
            for(auto J : mGraph[I])
                IJ.push_back(J); //J
        }
        return IJ;
    }

    void AddFromSingleVectorRepresentation(const std::vector<IndexType>& rSingleVectorRepresentation)
    {
        auto graph_size = rSingleVectorRepresentation[0];
        KRATOS_ERROR_IF(graph_size > GetGraph().size() ) << "mismatching size - attempting to add a graph with more rows than the ones allowed in graph" << std::endl;
        IndexType counter = 1;
        while(counter < rSingleVectorRepresentation.size())
        {
            auto I = rSingleVectorRepresentation[counter++];
            auto nrow = rSingleVectorRepresentation[counter++];
            auto begin = &rSingleVectorRepresentation[counter];
            AddEntries(I, begin, begin+nrow);
            counter += nrow;
        }
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{
    class const_iterator_adaptor : public std::iterator<
        std::forward_iterator_tag,
        typename GraphType::value_type
        >
	{
		const_row_iterator map_iterator;
        const_row_iterator mbegin;
	public:
		const_iterator_adaptor(const_row_iterator it) :map_iterator(it),mbegin(it) {}
		const_iterator_adaptor(const const_iterator_adaptor& it)
            : map_iterator(it.map_iterator),mbegin(it.mbegin) {}
		const_iterator_adaptor& operator++() { map_iterator++; return *this; }
		const_iterator_adaptor operator++(int) { const_iterator_adaptor tmp(*this); operator++(); return tmp; }
		bool operator==(const const_iterator_adaptor& rhs) const
            { return map_iterator == rhs.map_iterator; }
		bool operator!=(const const_iterator_adaptor& rhs) const
            { return map_iterator != rhs.map_iterator; }
        //TODO: is it correct that the two following operators are the same?
		const typename GraphType::value_type& operator*() const { return *map_iterator; }
		const typename GraphType::value_type& operator->() const { return *map_iterator; }
		const_row_iterator& base() { return map_iterator; }
		const_row_iterator const& base() const { return map_iterator; }
        IndexType GetRowIndex() const{
            return map_iterator-mbegin;
        }
	};

    const_iterator_adaptor begin() const
    {
        return const_iterator_adaptor( mGraph.begin() );
    }
    const_iterator_adaptor end() const
    {
        return const_iterator_adaptor( mGraph.end() );
    }


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
        buffer << "SparseContiguousRowGraph" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SparseContiguousRowGraph";}

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
    GraphType mGraph;
    std::vector<LockObject> mLocks;

    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        const IndexType N = this->Size();
        rSerializer.save("GraphSize",N);
        for(IndexType I=0; I<N; ++I)
        {
            IndexType row_size = mGraph[I].size();
            rSerializer.save("row_size",row_size);
            for(auto J : mGraph[I]){
                rSerializer.save("J",J);
            }
        }
    }

    void load(Serializer& rSerializer)
    {
        IndexType size;
        rSerializer.load("GraphSize",size);

        mLocks.resize(size);
        mGraph.resize(size);

        for(IndexType I=0; I<size; ++I)
        {
            IndexType row_size;
            rSerializer.load("row_size",row_size);
            for(IndexType k=0; k<row_size; ++k){
                IndexType J;
                rSerializer.load("J",J);
                AddEntry(I,J);
            }
        }
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

}; // Class SparseContiguousRowGraph

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TIndexType=std::size_t>
inline std::istream& operator >> (std::istream& rIStream,
                SparseContiguousRowGraph<TIndexType>& rThis){
                    return rIStream;
                }

/// output stream function
template<class TIndexType=std::size_t>
inline std::ostream& operator << (std::ostream& rOStream,
                const SparseContiguousRowGraph<TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPARSE_CONTIGUOUS_ROW_GRAPH_H_INCLUDED  defined


