//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <vector>
#include <limits>

// External includes

// Project includes
#include "utilities/pointer_communicator.h"
#include "spatial_containers/spatial_search_result.h"
#include "includes/indexed_object.h"
#include "containers/array_1d.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@}
///@name  Enum's
///@{

/**
 * @enum SpatialSearchCommunication
 * @brief Enum class for defining types of spatial search communication.
 * @details This enum class is used to specify the type of communication used to communicate search results between partitions.
 */
enum class SpatialSearchCommunication
{
    SYNCHRONOUS,   ///< Synchronous where all partitions know everything.
    ASYNCHRONOUS   ///< Asynchronous communication. Not implemented yet.
};

///@}
///@name Kratos Classes
///@{

class DataCommunicator;  // forward declaration

/**
 * @class SpatialSearchResultContainer
 * @brief Spatial search result container.
 * @details This class is used to store the results of a spatial search for one point.
 * @tparam TObjectType The type of the object.
 * @tparam TSpatialSearchCommunication The type of spatial search communication considered.
 * @ingroup KratosCore
 * @author Vicente Mataix Ferrandiz
 */
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication = SpatialSearchCommunication::SYNCHRONOUS>
class KRATOS_API(KRATOS_CORE) SpatialSearchResultContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResultContainer
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResultContainer);

    /// The index type
    using IndexType = std::size_t;

    /// Defining signed index type
    using SignedIndexType = std::ptrdiff_t;

    /// Spatial search result type
    using SpatialSearchResultType = SpatialSearchResult<TObjectType>;

    /// The spatial search result reference type
    using SpatialSearchResultReferenceType = SpatialSearchResultType&;

    /// The spatial search result const reference type
    using SpatialSearchResultConstReferenceType = const SpatialSearchResultType&;

    /// The spatial search result pointer type
    using SpatialSearchResultPointerType = typename SpatialSearchResultType::Pointer;

    /// Local vector of SpatialSearchResult
    using LocalResultsVector = std::vector<SpatialSearchResultType>;

    /// The global pointer communicator
    using GlobalPointerCommunicatorType = GlobalPointerCommunicator<SpatialSearchResultType>;

    /// The global pointer communicator pointer
    using GlobalPointerCommunicatorPointerType = typename GlobalPointerCommunicatorType::Pointer;

    /// The global pointer result type
    using GlobalPointerResultType = GlobalPointer<SpatialSearchResult<TObjectType>>;

    /// The global vector of SpatialSearchResult
    using GlobalResultsVector = GlobalPointersVector<SpatialSearchResultType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     */
    SpatialSearchResultContainer();

    /// Destructor.
    virtual ~SpatialSearchResultContainer() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Operator <<
     * @param os The output stream
     * @return The output stream
     */
    std::ostream& operator<<(std::ostream& os)
    {
        this->PrintData(os);
        return os;
    }

    /**
     * @brief Operator []
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultReferenceType operator[](const std::size_t Index)
    {
        return mLocalResults[Index];
    }

    /**
     * @brief Operator [] const version
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultConstReferenceType operator[](const std::size_t Index) const
    {
        return mLocalResults[Index];
    }

    /**
     * @brief Operator ()
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultReferenceType operator()(const std::size_t Index)
    {
        // Check if the communicator has been created
        KRATOS_ERROR_IF_NOT(mIsSynchronized) << "The data has not been synchronized" << std::endl;
        return mGlobalResults[Index];
    }

    /**
     * @brief Operator () const version
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultConstReferenceType operator()(const std::size_t Index) const
    {
        // Check if data is synchronized
        KRATOS_ERROR_IF_NOT(mIsSynchronized) << "The data has not been synchronized" << std::endl;
        return mGlobalResults[Index];
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns an iterator pointing to the beginning of the container.
     * @return Iterator pointing to the beginning of the container.
     */
    typename GlobalResultsVector::iterator begin()
    {
        return mGlobalResults.begin();
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the container.
     * @return Constant iterator pointing to the beginning of the container.
     */
    typename GlobalResultsVector::const_iterator begin() const
    {
        return mGlobalResults.begin();
    }

    /**
     * @brief Returns an iterator pointing to the end of the container.
     * @return Iterator pointing to the end of the container.
     */
    typename GlobalResultsVector::iterator end()
    {
        return mGlobalResults.end();
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the container.
     * @return Constant iterator pointing to the end of the container.
     */
    typename GlobalResultsVector::const_iterator end() const
    {
        return mGlobalResults.end();
    }

    /**
     * @brief Returns a reverse iterator pointing to the last element of the container.
     * @return Reverse iterator pointing to the last element of the container.
     */
    typename GlobalResultsVector::reverse_iterator rbegin()
    {
        return mGlobalResults.rbegin();
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the last element of the container.
     * @return Constant reverse iterator pointing to the last element of the container.
     */
    typename GlobalResultsVector::const_reverse_iterator rbegin() const
    {
        return mGlobalResults.rbegin();
    }

    /**
     * @brief Returns a reverse iterator pointing to the theoretical element preceding the first element of the container.
     * @return Reverse iterator pointing to the theoretical element preceding the first element.
     */
    typename GlobalResultsVector::reverse_iterator rend()
    {
        return mGlobalResults.rend();
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the theoretical element preceding the first element of the container.
     * @return Constant reverse iterator pointing to the theoretical element preceding the first element.
     */
    typename GlobalResultsVector::const_reverse_iterator rend() const
    {
        return mGlobalResults.rend();
    }

    /**
     * @brief Returns a pointer iterator pointing to the beginning of the container.
     * @return Pointer iterator pointing to the beginning of the container.
     */
    typename GlobalResultsVector::ptr_iterator ptr_begin()
    {
        return mGlobalResults.ptr_begin();
    }

    /**
     * @brief Returns a constant pointer iterator pointing to the beginning of the container.
     * @return Constant pointer iterator pointing to the beginning of the container.
     */
    typename GlobalResultsVector::ptr_const_iterator ptr_begin() const
    {
        return mGlobalResults.ptr_begin();
    }

    /**
     * @brief Returns a pointer iterator pointing to the end of the container.
     * @return Pointer iterator pointing to the end of the container.
     */
    typename GlobalResultsVector::ptr_iterator ptr_end()
    {
        return mGlobalResults.ptr_end();
    }

    /**
     * @brief Returns a constant pointer iterator pointing to the end of the container.
     * @return Constant pointer iterator pointing to the end of the container.
     */
    typename GlobalResultsVector::ptr_const_iterator ptr_end() const
    {
        return mGlobalResults.ptr_end();
    }

    /**
     * @brief Returns a reverse pointer iterator pointing to the last element of the container.
     * @return Reverse pointer iterator pointing to the last element of the container.
     */
    typename GlobalResultsVector::ptr_reverse_iterator ptr_rbegin()
    {
        return mGlobalResults.ptr_rbegin();
    }

    /**
     * @brief Returns a constant reverse pointer iterator pointing to the last element of the container.
     * @return Constant reverse pointer iterator pointing to the last element of the container.
     */
    typename GlobalResultsVector::ptr_const_reverse_iterator ptr_rbegin() const
    {
        return mGlobalResults.ptr_rbegin();
    }

    /**
     * @brief Returns a reverse pointer iterator pointing to the theoretical element preceding the first element of the container.
     * @return Reverse pointer iterator pointing to the theoretical element preceding the first element.
     */
    typename GlobalResultsVector::ptr_reverse_iterator ptr_rend()
    {
        return mGlobalResults.ptr_rend();
    }

    /**
     * @brief Returns a constant reverse pointer iterator pointing to the theoretical element preceding the first element of the container.
     * @return Constant reverse pointer iterator pointing to the theoretical element preceding the first element.
     */
    typename GlobalResultsVector::ptr_const_reverse_iterator ptr_rend() const
    {
        return mGlobalResults.ptr_rend();
    }

    /**
     * @brief Returns a reference to the first element of the container.
     * @return Reference to the first element of the container.
     */
    typename GlobalResultsVector::reference front() /* nothrow */
    {
        return mGlobalResults.front();
    }

    /**
     * @brief Returns a constant reference to the first element of the container.
     * @return Constant reference to the first element of the container.
     */
    typename GlobalResultsVector::const_reference front() const /* nothrow */
    {
        return mGlobalResults.front();
    }

    /**
     * @brief Returns a reference to the last element of the container.
     * @return Reference to the last element of the container.
     */
    typename GlobalResultsVector::reference back() /* nothrow */
    {
        return mGlobalResults.back();
    }

    /**
     * @brief Returns a constant reference to the last element of the container.
     * @return Constant reference to the last element of the container.
     */
    typename GlobalResultsVector::const_reference back() const /* nothrow */
    {
        return mGlobalResults.back();
    }

    /**
     * @brief Returns if at least one result is found
     * @return If at least one result is found
     */
    bool IsObjectFound() const
    {
        return static_cast<bool>(mGlobalResults.size());
    }

    /**
     * @brief Returns the local pointers size
     * @return The local pointers size
     */
    std::size_t NumberOfLocalResults() const
    {
        return mLocalResults.size();
    }

    /**
     * @brief Returns the global pointers size
     * @return The global pointers size
     */
    std::size_t NumberOfGlobalResults() const
    {
        return mGlobalResults.size();
    }

    /**
     * @brief Reserves the container
     * @details Only local
     * @param Size The size of the container
     */
    void Reserve(const std::size_t Size)
    {
        // Only local
        mLocalResults.reserve(Size);
    }

    /**
     * @brief Add a result to the container
     * @param rResult The result to be added
     */
    void AddResult(SpatialSearchResultReferenceType rResult);

    /**
     * @brief Pushes back a result to the container
     * @param rResult The result to be added
     */
    void push_back(SpatialSearchResultReferenceType rResult)
    {
        AddResult(rResult);
    }

    /**
     * @brief Add a result to the container
     * @param pResult The result to be added
     */
    void AddResult(
        TObjectType* pResult,
        const int Rank = 0
        );

    /**
     * @brief Add a result to the container
     * @param pResult The result to be added
     * @param Distance The distance to be added
     */
    void AddResult(
        TObjectType* pResult,
        const double Distance,
        const int Rank = 0
        );

    /**
     * @brief Pushes back a result to the container
     * @param pResult The result to be added
     */
    void push_back(TObjectType* pResult)
    {
        AddResult(pResult);
    }

    /**
     * @brief Clear the containers
     * @details This method clears the containers
     */
    void Clear();

    /**
     * @brief Removes elements from the given ranks.
     * @details This function takes a list of ranks and removes the elements at those ranks from the list.
     * @param rRanks A constant reference to a std::vector<int> containing the ranks where no local solution is expected.
     * @param rAllRanks A constant reference to a std::vector<int> containing all ranks.
     * @param rDataCommunicator The data communicator.
     */
    void RemoveResultsFromRanksList(
        const std::vector<int>& rRanks,
        const std::vector<int>& rAllRanks,
        const DataCommunicator& rDataCommunicator
        );

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Accessor for mLocalResults.
     * @details This method returns a reference to the LocalResultsVector mLocalResults.
     * @return A reference to the LocalResultsVector mLocalResults.
     */
    LocalResultsVector& GetLocalResults()
    {
        return mLocalResults;
    }

    /**
     * @brief Accessor for mGlobalResults.
     * @details This method returns a reference to the GlobalResultsVector mGlobalResults.
     * @return A reference to the GlobalResultsVector mGlobalResults.
     */
    GlobalResultsVector& GetGlobalResults()
    {
        return mGlobalResults;
    }

    /**
    * @brief Retrieve the global index value.
    * @return The global index value.
    */
    IndexType GetGlobalIndex() const
    {
        return mGlobalIndex;
    }

    /**
    * @brief Update the global index value.
    * @param Index The global index value to assign.
    */
    void SetGlobalIndex(const IndexType GlobalIndex)
    {
        // Assign the global index
        mGlobalIndex = GlobalIndex;
    }

    /**
    * @brief Get the local index value.
    * @return The local index value.
    */
    SignedIndexType GetLocalIndex() const
    {
        return mLocalIndex;
    }

    /**
    * @brief Set the local index value.
    * @param Index The local index value to set.
    */
    void SetLocalIndex(const SignedIndexType LocalIndex)
    {
        // Assign index
        mLocalIndex = LocalIndex;
    }

    /**
     * @brief Sets if the data is synchronized
     * @param IsSynchronized true if the data is synchronized, false otherwise
     */
    void SetIsSynchronized(const bool IsSynchronized)
    {
        mIsSynchronized = IsSynchronized;
    }

    /**
     * @brief Returns if the data is synchronized
     * @return true if the data is synchronized, false otherwise
     */
    bool GetIsSynchronized() const
    {
        return mIsSynchronized;
    }

    /**
     * @brief Check if the data is synchronized
     * @return true if the data is synchronized, false otherwise
     */
    bool IsSynchronized() const
    {
        return mIsSynchronized;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the point of the search is local.
     * @return true if the ranks match, false otherwise.
     */
    bool IsLocalPoint() const
    {
        return mLocalIndex >= 0;
    }

    /**
     * @brief Check if the search rank is the same as the rank of the data communicator.
     * @param rDataCommunicator The data communicator.
     * @param rResultsRank The ranks of the results.
     * @param GlobalIndex The global index of the search. Default is 0.
     * @return true if the ranks match, false otherwise.
     */
    bool IsLocalSearch(
        const DataCommunicator& rDataCommunicator,
        const std::vector<int>& rResultsRank,
        const IndexType GlobalIndex = 0
        )
    {
        return rResultsRank[GlobalIndex] == rDataCommunicator.Rank();
    }

    /**
     * @brief Check if the search rank is the same as the rank of the data communicator.
     * @param rDataCommunicator The data communicator.
     * @param rResultRank The rank of the result.
     * @return true if the ranks match, false otherwise.
     */
    bool IsLocalSearch(
        const DataCommunicator& rDataCommunicator,
        const int ResultRank
        )
    {
        return ResultRank == rDataCommunicator.Rank();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    ///@}
private:
    ///@name Member Variables
    ///@{

    LocalResultsVector mLocalResults;                                           /// Local results
    GlobalResultsVector mGlobalResults;                                         /// Global results

    SignedIndexType mLocalIndex = -1;                                           /// Some index considered for identification (local)
    IndexType mGlobalIndex = 0;                                                 /// Some index considered for identification (global)

    bool mIsSynchronized = false;                                               /// If the container is synchronized

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    /**
     * @brief Save the object in a serializer
     * @param rSerializer The serializer
     */
    void save(Serializer& rSerializer) const;

    /**
     * @brief Load the object in a serializer
     * @param rSerializer The serializer
     */
    void load(Serializer& rSerializer);

    ///@}
}; // Class SpatialSearchResultContainer

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
inline std::istream& operator>>(std::istream& rIStream,
                                SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.