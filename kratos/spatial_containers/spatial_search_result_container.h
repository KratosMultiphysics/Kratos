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
#include <unordered_map>

// External includes

// Project includes
#include "utilities/pointer_communicator.h"
#include "spatial_containers/spatial_search_result.h"
#include "containers/pointer_vector_set.h"
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

    /// Local vector of SpatialSearchResult
    using LocalResultsVector = PointerVectorSet<SpatialSearchResultType,
                               IndexedObject,
                               std::less<typename IndexedObject::result_type>,
                               std::equal_to<typename IndexedObject::result_type>,
                               typename SpatialSearchResultType::Pointer,
                               std::vector<typename SpatialSearchResultType::Pointer>
                               >;

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
    SpatialSearchResultType& operator[](const std::size_t Index)
    {
        return *(mLocalResults.begin() + Index);
    }

    /**
     * @brief Operator [] const version
     * @param Index The index
     * @return The result container
     */
    const SpatialSearchResultType& operator[](const std::size_t Index) const 
    {
        return *(mLocalResults.begin() + Index);
    }

    /**
     * @brief Operator ()
     * @param Index The index
     * @return The result container
     */
    SpatialSearchResultType& operator()(const std::size_t Index)
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
    const SpatialSearchResultType& operator()(const std::size_t Index) const 
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
    void AddResult(SpatialSearchResultType& rResult);

    /**
     * @brief Pushes back a result to the container
     * @param rResult The result to be added
     */
    void push_back(SpatialSearchResultType& rResult)
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
     * @brief Considers the global pointer communicator to get the shape functions of the resulting object
     * @param rPoint The point coordinates
     * @return A vector containing all the shape functions
     */
    std::vector<Vector> GetResultShapeFunctions(const array_1d<double, 3>& rPoint);

    /**
     * @brief Considers the global pointer communicator to get the indices of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<std::size_t> GetResultIndices();

    /**
     * @brief Considers the global pointer communicator to get the indices of the nodes of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<std::vector<std::size_t>> GetResultNodeIndices();

    /**
     * @brief Considers the global pointer communicator to get the partition indices of the nodes of the resulting object
     * @return A vector containing all the indices
     */
    std::vector<std::vector<std::size_t>> GetResultPartitionIndices();

    /**
     * @brief Considers the global pointer communicator to get the coordinates of the resulting object
     * @return A vector containing all the coordinates
     */
    std::vector<std::vector<array_1d<double, 3>>> GetResultCoordinates();

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

/**
 * @class SpatialSearchResultContainerMap
 * @brief Spatial search result container map
 * @details This class is used to store the results of a spatial search, in a map to be identify results for a given coordinates
 * @tparam TObjectType The type of the object
 * @ingroup KratosCore
 * @author Vicente Mataix Ferrandiz
 */
template <class TObjectType>
class KRATOS_API(KRATOS_CORE) SpatialSearchResultContainerMap
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResultContainerMap
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResultContainerMap);

    /// The hash type
    using HashType = std::size_t;

    /// The container type
    using ContainerType = std::unordered_map<HashType, SpatialSearchResultContainer<TObjectType>>;

    // Define the iterator class
    class iterator {
    public:
        /// The category of the iterator, indicating forward iteration.
        using iterator_category = std::forward_iterator_tag;

        /// The type of the value pointed to by the iterator.
        using value_type = std::pair<const HashType, SpatialSearchResultContainer<TObjectType>>;

        /// The difference type between two iterators.
        using difference_type = std::ptrdiff_t;

        /// A pointer to the value type.
        using pointer = value_type*;

        /// A reference to the value type.
        using reference = value_type&;

        /**
         * @brief Constructs an iterator pointing to the specified position.
         * @param iter The iterator to initialize from.
         */
        iterator(typename ContainerType::iterator iter) : iter_(iter) {}

        /**
         * @brief Prefix increment operator.
         * @return Reference to the incremented iterator.
         */
        iterator& operator++() {
            ++iter_;
            return *this;
        }

        /**
         * @brief Postfix increment operator.
         * @return An iterator before increment.
         */
        iterator operator++(int) {
            iterator temp = *this;
            ++(*this);
            return temp;
        }

        /**
         * @brief Equality comparison operator.
         * @param other The iterator to compare with.
         * @return True if the iterators are equal, false otherwise.
         */
        bool operator==(const iterator& other) const {
            return iter_ == other.iter_;
        }

        /**
         * @brief Inequality comparison operator.
         * @param other The iterator to compare with.
         * @return True if the iterators are not equal, false otherwise.
         */
        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }

        /**
         * @brief Dereference operator.
         * @return Reference to the value pointed to by the iterator.
         */
        reference operator*() const {
            return *iter_;
        }

        /**
         * @brief Member access operator.
         * @return Pointer to the value pointed to by the iterator.
         */
        pointer operator->() const {
            return &(*iter_);
        }

    private:
        typename ContainerType::iterator iter_;
    };

    // Define the const_iterator class
    class const_iterator {
    public:
        /// The category of the iterator, indicating forward iteration.
        using iterator_category = std::forward_iterator_tag;

        /// The type of the value pointed to by the iterator.
        using value_type = std::pair<const HashType, SpatialSearchResultContainer<TObjectType>>;

        /// The difference type between two iterators.
        using difference_type = std::ptrdiff_t;

        /// A pointer to the value type.
        using pointer = const value_type*;

        /// A reference to the value type.
        using reference = const value_type&;

        /**
         * @brief Constructs a constant iterator pointing to the specified position.
         * @param iter The constant iterator to initialize from.
         */
        const_iterator(typename ContainerType::const_iterator iter) : iter_(iter) {}

        /**
         * @brief Prefix increment operator.
         * @return Reference to the incremented constant iterator.
         */
        const_iterator& operator++() {
            ++iter_;
            return *this;
        }

        /**
         * @brief Postfix increment operator.
         * @return A constant iterator before increment.
         */
        const_iterator operator++(int) {
            const_iterator temp = *this;
            ++(*this);
            return temp;
        }

        /**
         * @brief Equality comparison operator.
         * @param other The constant iterator to compare with.
         * @return True if the constant iterators are equal, false otherwise.
         */
        bool operator==(const const_iterator& other) const {
            return iter_ == other.iter_;
        }

        /**
         * @brief Inequality comparison operator.
         * @param other The constant iterator to compare with.
         * @return True if the constant iterators are not equal, false otherwise.
         */
        bool operator!=(const const_iterator& other) const {
            return !(*this == other);
        }

        /**
         * @brief Dereference operator.
         * @return Reference to the value pointed to by the constant iterator.
         */
        reference operator*() const {
            return *iter_;
        }

        /**
         * @brief Member access operator.
         * @return Pointer to the value pointed to by the constant iterator.
         */
        pointer operator->() const {
            return &(*iter_);
        }

    private:
        typename ContainerType::const_iterator iter_;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SpatialSearchResultContainerMap() = default;

    /// Destructor.
    virtual ~SpatialSearchResultContainerMap() = default;

    ///@}
    ///@name Operators
    ///@{

   /**
     * @brief Operator []
     * @param Index The index to be initialized
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator[](const IndexType Index)
    {
        const auto it = mPointResults.find(Index);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for index: " << Index << std::endl;
        return it->second;
    }

    /**
     * @brief Operator []
     * @param Index The index to be initialized
     * @return The result container
     */
    const SpatialSearchResultContainer<TObjectType>& operator[](const IndexType Index) const
    {
        const auto it = mPointResults.find(Index);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for index: " << Index << std::endl;
        return it->second;
    }

    /**
     * @brief Operator ()
     * @param Index The index to be initialized
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator()(const IndexType Index)
    {
        const auto it = mPointResults.find(Index);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for index: " << Index << std::endl;
        return it->second;
    }

    /**
     * @brief Operator ()
     * @param Index The index to be initialized
     * @return The result container
     */
    const SpatialSearchResultContainer<TObjectType>& operator()(const IndexType Index) const
    {
        const auto it = mPointResults.find(Index);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for index: " << Index << std::endl;
        return it->second;
    }

    /**
     * @brief Operator []
     * @param rCoordinates The coordinates
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator[](const array_1d<double, 3>& rCoordinates)
    {
        const HashType hash = Hash(rCoordinates);
        const auto it = mPointResults.find(hash);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for point: " << rCoordinates[0] << ", " << rCoordinates[1] << ", " << rCoordinates[2] << std::endl;
        return it->second;
    }

    /**
     * @brief Operator []
     * @param rCoordinates The coordinates
     * @return The result container
     */
    const SpatialSearchResultContainer<TObjectType>& operator[](const array_1d<double, 3>& rCoordinates) const
    {
        const HashType hash = Hash(rCoordinates);
        const auto it = mPointResults.find(hash);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for point: " << rCoordinates[0] << ", " << rCoordinates[1] << ", " << rCoordinates[2] << std::endl;
        return it->second;
    }

    /**
     * @brief Operator ()
     * @param rCoordinates The coordinates
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator()(const array_1d<double, 3>& rCoordinates)
    {
        const HashType hash = Hash(rCoordinates);
        const auto it = mPointResults.find(hash);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for point: " << rCoordinates[0] << ", " << rCoordinates[1] << ", " << rCoordinates[2] << std::endl;
        return it->second;
    }

    /**
     * @brief Operator ()
     * @param rCoordinates The coordinates
     * @return The result container
     */
    const SpatialSearchResultContainer<TObjectType>& operator()(const array_1d<double, 3>& rCoordinates) const
    {
        const HashType hash = Hash(rCoordinates);
        const auto it = mPointResults.find(hash);
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist for point: " << rCoordinates[0] << ", " << rCoordinates[1] << ", " << rCoordinates[2] << std::endl;
        return it->second;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns an iterator pointing to the beginning of the container.
     * @return An iterator pointing to the beginning of the container.
     */
    iterator begin() 
    {
        return iterator(mPointResults.begin());
    }

    /**
     * @brief Returns an iterator pointing to the end of the container.
     * @return An iterator pointing to the end of the container.
     */
    iterator end() 
    {
        return iterator(mPointResults.end());
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the container.
     * @return A constant iterator pointing to the beginning of the container.
     */
    const_iterator begin() const 
    {
        return const_iterator(mPointResults.begin());
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the container.
     * @return A constant iterator pointing to the end of the container.
     */
    const_iterator end() const 
    {
        return const_iterator(mPointResults.end());
    }

    /**
     * @brief Returns the number of points results
     * @return The number of points results
     */
    std::size_t NumberOfSearchResults() const
    {
        return mPointResults.size();
    }

    /**
     * @brief Initialize the container
     * @param Index The index to be initialized
     */
    SpatialSearchResultContainer<TObjectType>& InitializeResult(const IndexType Index);

    /**
     * @brief Initialize the container
     * @param rCoordinates The coordinates
     */
    SpatialSearchResultContainer<TObjectType>& InitializeResult(const array_1d<double, 3>& rCoordinates);

    /**
     * @brief Check if coordinates are initialized
     * @param Index The index to be initialized
     * @return True if hash is initialized, false otherwise
    */
    bool HasResult(const IndexType Index) const;

    /**
     * @brief Check if coordinates are initialized
     * @param rCoordinates The coordinates
     * @return True if coordinates are initialized, false otherwise
    */
    bool HasResult(const array_1d<double, 3>& rCoordinates) const;

    /**
     * @brief Clear the containers
     * @details This method clears the containers
     */
    void Clear();

    /**
     * @brief Synchronize all container between partitions
     * @details This method synchronizes all the container between partitions
     * @param rDataCommunicator The data communicator
     */
    void SynchronizeAll(const DataCommunicator& rDataCommunicator);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the container
     * @return The container
     */
    ContainerType& GetContainer() 
    {
        return mPointResults;
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
    
    ContainerType mPointResults;  /// The results of each point

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Hash the coordinates
     * @param rCoordinates The coordinates
     * @return The hash
     */
    HashType Hash(const array_1d<double, 3>& rCoordinates) const;

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
}; // Class SpatialSearchResultContainerMap

/// input stream function
template <class TObjectType>
inline std::istream& operator>>(std::istream& rIStream,
                                SpatialSearchResultContainerMap<TObjectType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TObjectType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SpatialSearchResultContainerMap<TObjectType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.