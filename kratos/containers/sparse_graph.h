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
#include "includes/ublas_interface.h"
#include "includes/serializer.h"

// External includes
#include <unordered_map>
#include <unordered_set>

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

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
class SparseGraph
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
    : mGraphSize(N)
    {
        mSizeIsAvailable = true;
    }

    /// Default constructor.
    SparseGraph()
    : mGraphSize(0)
    {
        mSizeIsAvailable = false;
    }

    /// Destructor.
    virtual ~SparseGraph(){}

    /// Copy constructor. TODO: we need it otherwise the sendrecv does not work...
    ///but i don't know why :-(
    SparseGraph(const SparseGraph& rOther)
    {
        this->mGraph = rOther.mGraph;
        mGraphSize = rOther.mGraphSize;
        mSizeIsAvailable = rOther.mSizeIsAvailable;
    }

    ///@}
    ///@name Operators
    ///@{

    IndexType Size() const{
        KRATOS_ERROR_IF(mSizeIsAvailable == false) << "SparseGraph Size method can only be called after Finalize, or if provided in the constructor" << std::endl;
        return mGraphSize; //note that this is only valid after Finalize has been called
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
        mSizeIsAvailable = false;
        mGraphSize = 0;
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

    void AddEntries(SparseGraph& rOtherGraph)
    {
        for(const auto& item: rOtherGraph.GetGraph())
        {
            AddEntries(item.first, item.second);
        }
    }

    void Finalize()
    {
        if(!mSizeIsAvailable){
            mGraphSize=0;
            for(const auto& item : this->GetGraph())
            {
                mGraphSize = std::max(mGraphSize,item.first+1);
            }
            mSizeIsAvailable = true;
        }
    }

    const GraphType& GetGraph() const{
        return mGraph;
    }

    IndexType ExportCSRArrays(
        vector<IndexType>& rRowIndices,
        vector<IndexType>& rColIndices
    ) const
    {
        //need to detect the number of rows this way since there may be gaps
        IndexType nrows=this->Size();

        if(rRowIndices.size() != nrows+1)
        {
            rRowIndices.resize(nrows+1, false);
        }
        //set it to zero in parallel to allow first touching
        IndexPartition<IndexType>(rRowIndices.size()).for_each([&](IndexType i){
                    rRowIndices[i] = 0;
                });            

        //count the entries TODO: do the loop in parallel if possible
        for(const auto& item : this->GetGraph())
        {
            rRowIndices[item.first+1] = item.second.size();
        }
        //sum entries
        for(int i = 1; i<static_cast<int>(rRowIndices.size()); ++i){
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

        //count the entries TODO: do the loop in parallel if possible
        for(const auto& item : this->GetGraph()){
            IndexType start = rRowIndices[item.first];

            IndexType counter = 0;
            for(auto index : item.second){
                rColIndices[start+counter] = index;
                counter++;
            }
        }

        //reorder columns
        IndexPartition<IndexType>(rRowIndices.size()-1).for_each([&](IndexType i){
            std::sort(rColIndices.begin()+rRowIndices[i], rColIndices.begin()+rRowIndices[i+1]);
        });
        return nrows;
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
	public:
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SparseGraph" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SparseGraph";}

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
    IndexType mGraphSize = 0;
    GraphType mGraph;
    bool mSizeIsAvailable = false;


    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        std::vector< IndexType > IJ;
        for(const auto& item : this->GetGraph()){
            IndexType I = item.first;
            IJ.push_back(I);
            IJ.push_back(item.second.size());
            for(auto J : item.second)
                IJ.push_back(J);
        }
        rSerializer.save("IJ",IJ);
        rSerializer.save("size",mGraphSize);
        rSerializer.save("mSizeIsAvailable",mSizeIsAvailable);
    }

    void load(Serializer& rSerializer)
    {
        std::vector< IndexType > IJ;
        rSerializer.load("IJ",IJ);

        if(IJ.size() != 0)
        {
            IndexType counter = 0;
            while(counter < IJ.size())
            {
                auto I = IJ[counter++];
                auto nrow = IJ[counter++];
                auto begin = &IJ[counter];s
                AddEntries(I, begin, begin+nrow);
                counter += nrow;
            }
        }

        rSerializer.load("size",mGraphSize);
        rSerializer.load("mSizeIsAvailable",mSizeIsAvailable);

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


