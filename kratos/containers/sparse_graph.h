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

class SparseGraph
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::map<IndexType, std::unordered_set<IndexType> > GraphType; //using a map since we need it ordered

    /// Pointer definition of SparseGraph
    KRATOS_CLASS_POINTER_DEFINITION(SparseGraph);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SparseGraph(){}

    /// Destructor.
    virtual ~SparseGraph(){}

    ///@}
    ///@name Operators
    ///@{
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

    IndexType ExportCSRArrays(
        vector<IndexType>& rRowIndices,
        vector<IndexType>& rColIndices
    )
    {
        //need to detect the number of rows this way since there may be gaps
        IndexType nrows=0;
        for(const auto& item : this->GetGraph())
        {
            nrows = std::max(nrows,item.first+1);
        }

        if(rRowIndices.size() != nrows+1)
        {
            rRowIndices.resize(nrows+1, false);
        }
        //set it to zero in parallel to allow first touching
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(nrows+1); ++i)
            rRowIndices[i] = 0;

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
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rColIndices.size());++i)
            rColIndices[i] = 0;

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
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rRowIndices.size()-1);++i){
            std::sort(rColIndices.begin()+rRowIndices[i], rColIndices.begin()+rRowIndices[i+1]);
        }
        return nrows;
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
    GraphType mGraph;

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
                AddEntries(I, &IJ[counter],&IJ[counter+nrow]);
                counter += nrow;
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

    /// Assignment operator.
    SparseGraph& operator=(SparseGraph const& rOther) = delete;

    /// Copy constructor.
    SparseGraph(SparseGraph const& rOther) = delete;

    ///@}

}; // Class SparseGraph

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                SparseGraph& rThis){
                    return rIStream;
                }

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const SparseGraph& rThis)
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


