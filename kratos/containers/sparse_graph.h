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

#if !defined(KRATOS_SPARSE_GRAPH_H_INCLUDED )
#define  KRATOS_SPARSE_GRAPH_H_INCLUDED


// System includes
#include <iostream>
#include <unordered_map>
#include <unordered_set>

// External includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "includes/parallel_environment.h"

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

template<class TIndexType=std::size_t>
class SparseGraph final
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef std::map<IndexType, std::unordered_set<IndexType> > GraphType; //using a map since we need it ordered
    typedef typename GraphType::const_iterator const_row_iterator;

    /// Pointer definition of SparseGraph
    KRATOS_CLASS_POINTER_DEFINITION(SparseGraph);

    ///@}
    ///@name Life Cycle
    ///@{

    SparseGraph(IndexType N)
    {
        mpComm = &ParallelEnvironment::GetDataCommunicator("Serial"); //that's the only option
    }

    /// Default constructor.
    SparseGraph()
    {
        mpComm = &ParallelEnvironment::GetDataCommunicator("Serial"); //thats the only option
    }

    SparseGraph(DataCommunicator& rComm)
    {
        if(rComm.IsDistributed())
            KRATOS_ERROR << "Attempting to construct a serial CsrMatrix with a distributed communicator" << std::endl;

        mpComm = &rComm;
    }

    /// Destructor.
    ~SparseGraph(){}

    /// Copy constructor. 
    SparseGraph(const SparseGraph& rOther)
    :
        mpComm(rOther.mpComm),
        mGraph(rOther.mGraph)
    {
    }

    SparseGraph(const std::vector<IndexType>& rSingleVectorRepresentation)
    {
        AddFromSingleVectorRepresentation(rSingleVectorRepresentation);
    }

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm() const
    {
        return *mpComm;
    }

    const DataCommunicator* pGetComm() const
    {
        return mpComm;
    }

    IndexType Size() const
    {
        if(!this->GetGraph().empty())
            return (--(this->GetGraph().end()))->first + 1; //last key in the map + 1
        return 0;
    }

    bool IsEmpty() const
    {
        return mGraph.empty();
    }

    bool Has(const IndexType I, const IndexType J) const
    {
        const auto& row_it = mGraph.find(I);
        if(row_it != mGraph.end() ) {
            if((row_it->second).find(J) != (row_it->second).end())
                return true;
        }
        return false;
    }

    const typename GraphType::mapped_type& operator[](const IndexType& Key) const
    {
		return (mGraph.find(Key))->second;
    }

    void Clear()
    {
        mGraph.clear();
    }

    void AddEntry(const IndexType RowIndex, const IndexType ColIndex)
    {
        mGraph[RowIndex].insert(ColIndex);
    }

    template<class TContainerType>
    void AddEntries(const IndexType RowIndex, const TContainerType& rColIndices)
    {
        mGraph[RowIndex].insert(rColIndices.begin(), rColIndices.end());
    }

    template<class TIteratorType>
    void AddEntries(const IndexType RowIndex,
                    const TIteratorType& rColBegin,
                    const TIteratorType& rColEnd
                    )
    {
        mGraph[RowIndex].insert(rColBegin, rColEnd);
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rIndices)
    {
        for(auto I : rIndices)
            mGraph[I].insert(rIndices.begin(), rIndices.end());
    }

    template<class TContainerType>
    void AddEntries(const TContainerType& rRowIndices, const TContainerType& rColIndices)
    {
        for(auto I : rRowIndices)
            AddEntries(I, rColIndices);
    }


    void AddEntries(SparseGraph& rOtherGraph)
    {
        for(const auto& item: rOtherGraph.GetGraph())
        {
            AddEntries(item.first, item.second);
        }
    }

    void Finalize()
    {
    }

    const GraphType& GetGraph() const{
        return mGraph;
    }


    //note that this function is slightly less efficient than the span version below
    template<class TVectorType=DenseVector<IndexType>>
    IndexType ExportCSRArrays(
        TVectorType& rRowIndices,
        TVectorType& rColIndices
    ) const
    {
        IndexType* pRowIndicesData=nullptr;
        IndexType RowIndicesDataSize=0;
        IndexType* pColIndicesData=nullptr;
        IndexType ColIndicesDataSize=0;
        ExportCSRArrays(pRowIndicesData,RowIndicesDataSize,pColIndicesData,ColIndicesDataSize);
        if(rRowIndices.size() != RowIndicesDataSize)
            rRowIndices.resize(RowIndicesDataSize);
        IndexPartition<IndexType>(RowIndicesDataSize).for_each( 
            [&](IndexType i){rRowIndices[i] = pRowIndicesData[i];}
        );
        
        delete [] pRowIndicesData;
        if(rColIndices.size() != ColIndicesDataSize)
            rColIndices.resize(ColIndicesDataSize);
        IndexPartition<IndexType>(ColIndicesDataSize).for_each( 
            [&](IndexType i){rColIndices[i] = pColIndicesData[i];}
        );
        delete [] pColIndicesData;

        return rRowIndices.size();
    }

    IndexType ExportCSRArrays(
        Kratos::span<IndexType>& rRowIndices,
        Kratos::span<IndexType>& rColIndices
    ) const = delete;

    //NOTE this function will transfer ownership of pRowIndicesData and pColIndicesData to the caller
    //hence the caller will be in charge of deleting that array
    IndexType ExportCSRArrays(
        IndexType*& pRowIndicesData,
        IndexType& rRowDataSize,
        IndexType*& pColIndicesData,
        IndexType& rColDataSize
    ) const
    {
        //need to detect the number of rows this way since there may be gaps
        IndexType nrows=this->Size();

        pRowIndicesData = new IndexType[nrows+1];
        rRowDataSize = nrows+1;
        Kratos::span<IndexType> row_indices(pRowIndicesData, nrows+1);
        
        //set it to zero in parallel to allow first touching
        IndexPartition<IndexType>(row_indices.size()).for_each([&](IndexType i){
                    row_indices[i] = 0;
                });            

        //count the entries TODO: do the loop in parallel if possible
        for(const auto& item : this->GetGraph())
        {
            row_indices[item.first+1] = item.second.size();
        }

        //sum entries
        for(int i = 1; i<static_cast<int>(row_indices.size()); ++i){
            row_indices[i] += row_indices[i-1];
        }


        IndexType nnz = row_indices[nrows];
        rColDataSize = nnz;
        pColIndicesData = new IndexType[nnz];
        Kratos::span<IndexType> col_indices(pColIndicesData, nnz);
        
        //set it to zero in parallel to allow first touching
        IndexPartition<IndexType>(col_indices.size()).for_each([&](IndexType i){
                    col_indices[i] = 0;
                });            

        //count the entries TODO: do the loop in parallel if possible
        for(const auto& item : this->GetGraph()){
            IndexType start = row_indices[item.first];

            IndexType counter = 0;
            for(auto index : item.second){
                col_indices[start+counter] = index;
                counter++;
            }
        }

        //reorder columns
        IndexPartition<IndexType>(row_indices.size()-1).for_each([&](IndexType i){
            std::sort(col_indices.begin()+row_indices[i], col_indices.begin()+row_indices[i+1]);
        });
        return nrows;
    }

    //this function returns the Graph as a single vector
    //in the form of 
    //  RowIndex NumberOfEntriesInTheRow .... list of all Indices in the row 
    //every row is pushed back one after the other 
    std::vector<IndexType> ExportSingleVectorRepresentation() const
    {
        std::vector< IndexType > single_vector_representation;

        IndexType nrows = this->Size(); //note that the number of non zero rows may be different from this one
        single_vector_representation.push_back(nrows);

        for(const auto& item : this->GetGraph()){
            IndexType I = item.first;
            single_vector_representation.push_back(I); //we store the index of the rows
            single_vector_representation.push_back(item.second.size()); //the number of items in the row 
            for(auto J : item.second)
                single_vector_representation.push_back(J); //the columns
        }
        return single_vector_representation;
    }

    void AddFromSingleVectorRepresentation(const std::vector<IndexType>& rSingleVectorRepresentation)
    {
        //IndexType graph_size = rSingleVectorRepresentation[0]; 
        //we actually do not need the graph_size since it will be reconstructed, 
        //however it is important that the graph_size is stored to make the single_vector
        //representation also compatible with the sparse_contiguous_row_graph
        
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
    class const_iterator_adaptor
	{
		const_row_iterator map_iterator;
	public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = typename GraphType::value_type;
        using pointer           = typename GraphType::value_type*;
        using reference         = typename GraphType::value_type&;

		const_iterator_adaptor(const_row_iterator it) :map_iterator(it) {}
		const_iterator_adaptor(const const_iterator_adaptor& it)
            : map_iterator(it.map_iterator) {}
		const_iterator_adaptor& operator++() { map_iterator++; return *this; }
		const_iterator_adaptor operator++(int) { const_iterator_adaptor tmp(*this); operator++(); return tmp; }
		bool operator==(const const_iterator_adaptor& rhs) const
            { return map_iterator == rhs.map_iterator; }
		bool operator!=(const const_iterator_adaptor& rhs) const
            { return map_iterator != rhs.map_iterator; }
        //TODO: is it correct that the two following operators are the same?
		const typename GraphType::mapped_type& operator*() const { return (map_iterator->second); }
		const typename GraphType::mapped_type& operator->() const { return map_iterator->second; }
		const_row_iterator& base() { return map_iterator; }
		const_row_iterator const& base() const { return map_iterator; }
        IndexType GetRowIndex() const{
            return map_iterator->first;
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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SparseGraph" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {rOStream << "SparseGraph";}

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
    DataCommunicator* mpComm;
    GraphType mGraph;


    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {   
        std::vector< IndexType > IJ = ExportSingleVectorRepresentation();
        rSerializer.save("IJ",IJ);
    }

    void load(Serializer& rSerializer)
    {
        std::vector< IndexType > IJ;
        rSerializer.load("IJ",IJ);
        AddFromSingleVectorRepresentation(IJ);
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
    SparseGraph& operator=(SparseGraph const& rOther) = delete;


    ///@}

}; // Class SparseGraph

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TIndexType=std::size_t>
inline std::istream& operator >> (std::istream& rIStream,
                SparseGraph<TIndexType>& rThis){
                    return rIStream;
                }

/// output stream function
template<class TIndexType=std::size_t>
inline std::ostream& operator << (std::ostream& rOStream,
                const SparseGraph<TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPARSE_GRAPH_H_INCLUDED  defined


