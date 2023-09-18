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
#include "spatial_containers/spatial_search_result_container.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

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