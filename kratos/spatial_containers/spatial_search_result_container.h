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

///@name Kratos Classes
///@{

class DataCommunicator;  // forward declaration

/**
 * @class SpatialSearchResultContainer
 * @brief Spatial search result container
 * @details This class is used to store the results of a spatial search
 * @tparam TObjectType The type of the object
 * @ingroup KratosCore
 * @author Vicente Mataix Ferrandiz
 */
template <class TObjectType>
class KRATOS_API(KRATOS_CORE) SpatialSearchResultContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResultContainer
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResultContainer);

    /// The index type
    using IndexType = std::size_t;

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
        KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created. Therefore is not synchronized" << std::endl;
        return mGlobalResults[Index];
    }

    /**
     * @brief Operator () const version
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultConstReferenceType operator()(const std::size_t Index) const
    {
        // Check if the communicator has been created
        KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created. Therefore is not synchronized" << std::endl;
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
     * @brief Generate the global pointer communicator
     * @param rDataCommunicator The data communicator
     */
    void GenerateGlobalPointerCommunicator(const DataCommunicator& rDataCommunicator);

    /**
     * @brief Synchronize the container between partitions
     * @details This method synchronizes the container between partitions
     * @param rDataCommunicator The data communicator
     */
    void SynchronizeAll(const DataCommunicator& rDataCommunicator);

    /**
     * @brief Applies a user-provided function to the global pointers and return a proxy to the results.
     * @tparam TFunctorType Functor type.
     * @param UserFunctor The user-provided function.
     * @return A proxy to the results.
     */
    template<class TFunctorType>
    ResultsProxy<
    SpatialSearchResultType,
    TFunctorType // TODO: Unfortunately this is deprecated in c++17, so we will have to change this call in the future
    > Apply(TFunctorType&& UserFunctor)
    {
        // Check if the communicator has been created
        KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

        // Apply the user-provided function
        return mpGlobalPointerCommunicator->Apply(std::forward<TFunctorType>(UserFunctor));
    }

    /**
     * @brief Retrieves the global distances
     * @return A vector containing all the distances
     */
    std::vector<double> GetDistances();

    /**
     * @brief Retrieves if is local the entity
     * @return A vector containing all the booleans showing is local the entity
     */
    std::vector<bool> GetResultIsLocal();

    /**
     * @brief Retrieves the rank of the entity
     * @return A vector containing all the ranks of the entity
     */
    std::vector<int> GetResultRank();

    /**
     * @brief Retrieves if is active the entity
     * @return A vector containing all the booleans showing is active the entity
     */
    std::vector<bool> GetResultIsActive();

    /**
     * @brief Retrieves if inside the geometry
     * @param rPoint The point coordinates
     * @param Tolerance The tolerance considered
     * @return A vector containing all the booleans showing is inside the geometry
     */
    std::vector<bool> GetResultIsInside(
        const array_1d<double, 3>& rPoint,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        );

    /**
     * @brief Considers the global pointer communicator to get the shape functions of the resulting object
     * @param rPoint The point coordinates
     * @return A vector containing all the shape functions
     */
    std::vector<Vector> GetResultShapeFunctions(const array_1d<double, 3>& rPoint);

    /**
     * @brief Considers the global pointer communicator to get the indices of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<IndexType> GetResultIndices();

    /**
     * @brief Considers the global pointer communicator to get the indices of the nodes of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<std::vector<IndexType>> GetResultNodeIndices();

    /**
     * @brief Considers the global pointer communicator to get the partition indices of the nodes of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<std::vector<int>> GetResultPartitionIndices();

    /**
     * @brief Considers the global pointer communicator to get the coordinates of the resulting object
     * @return A vector containing all the coordinates
     */
    std::vector<std::vector<array_1d<double, 3>>> GetResultCoordinates();

    /**
     * @brief Removes elements at specified indexes from a list.
     * @details This function takes a list of indexes and removes the elements at those indexes from the list.
     * @param rIndexes A constant reference to a std::vector<IndexType> containing the indexes of elements to be removed.
     */
    void RemoveResultsFromIndexesList(const std::vector<IndexType>& rIndexes);

    /**
    * @brief Generates a vector of indexes where the elements in the input vector are greater than zero.
    * @details This function takes an input vector and iterates through its elements. It adds the indexes of elements
    * that are greater than zero to a new vector and returns that vector.
    * @param rInputVector The input vector to process.
    * @return A vector of indexes where elements in the input vector are greater than zero.
    * @note Static method that can be used without an instance of this class. Will be called in the vector class.
    */
    static std::vector<int> GenerateGreaterThanZeroIndexes(const std::vector<int>& rInputVector);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Accessor for mLocalResults.
     * This method returns a reference to the LocalResultsVector mLocalResults.
     * @return A reference to the LocalResultsVector mLocalResults.
     */
    LocalResultsVector& GetLocalResults()
    {
        return mLocalResults;
    }

    /**
     * @brief Accessor for mGlobalResults.
     * This method returns a reference to the GlobalResultsVector mGlobalResults.
     * @return A reference to the GlobalResultsVector mGlobalResults.
     */
    GlobalResultsVector& GetGlobalResults()
    {
        return mGlobalResults;
    }

    /**
     * @brief Accessor for mpGlobalPointerCommunicator.
     * This method returns the GlobalPointerCommunicatorPointer mpGlobalPointerCommunicator.
     * @return The GlobalPointerCommunicatorPointer mpGlobalPointerCommunicator.
     */
    GlobalPointerCommunicatorPointerType GetGlobalPointerCommunicator()
    {
        return mpGlobalPointerCommunicator;
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

    GlobalPointerCommunicatorPointerType mpGlobalPointerCommunicator = nullptr; /// Global pointer to the communicator

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

/// input stream function
template <class TObjectType>
inline std::istream& operator>>(std::istream& rIStream,
                                SpatialSearchResultContainer<TObjectType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TObjectType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SpatialSearchResultContainer<TObjectType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.