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
 * @class SpatialSearchResultContainerVector
 * @brief Spatial search result container map
 * @details This class is used to store the results of a spatial search, in a map to be identify results for a given coordinates
 * @tparam TObjectType The type of the object
 * @ingroup KratosCore
 * @author Vicente Mataix Ferrandiz
 */
template <class TObjectType>
class KRATOS_API(KRATOS_CORE) SpatialSearchResultContainerVector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResultContainerVector
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResultContainerVector);

    /// The spatial search result container type
    using SpatialSearchResultContainerType = SpatialSearchResultContainer<TObjectType>;

    /// The spatial search result container reference type
    using SpatialSearchResultContainerReferenceType = SpatialSearchResultContainerType&;

    /// The spatial search result container pointer type
    using SpatialSearchResultContainerPointerType = SpatialSearchResultContainerType*;

    /// The container type
    using ContainerType = std::vector<SpatialSearchResultContainerPointerType>;

    // Define the iterator class
    class iterator {
    public:
        /// The category of the iterator, indicating forward iteration.
        using iterator_category = std::forward_iterator_tag;

        /// The type of the value pointed to by the iterator.
        using value_type = SpatialSearchResultContainerType;

        /// The difference type between two iterators.
        using difference_type = std::ptrdiff_t;

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
        SpatialSearchResultContainerReferenceType operator*() const {
            return **iter_;
        }

        /**
        * @brief Member access operator.
        * @return Pointer to the value pointed to by the iterator.
        */
        SpatialSearchResultContainerPointerType operator->() const {
            return *iter_;
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
        using value_type = SpatialSearchResultContainerType;

        /// The difference type between two iterators.
        using difference_type = std::ptrdiff_t;

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
        * @return Reference to the value pointed to by the iterator.
        */
        const SpatialSearchResultContainerReferenceType operator*() const {
            return **iter_;
        }

        /**
        * @brief Member access operator.
        * @return Pointer to the value pointed to by the iterator.
        */
        const SpatialSearchResultContainerType* operator->() const {
            return *iter_;
        }

    private:
        typename ContainerType::const_iterator iter_;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SpatialSearchResultContainerVector() = default;

    /// Destructor.
    virtual ~SpatialSearchResultContainerVector() {
        // Make sure to delete the pointers stored in the container
        for (auto pResult : mPointResults) {
            delete pResult;
        }
    }

    ///@}
    ///@name Operators
    ///@{

   /**
     * @brief Operator []
     * @param Index The index to be initialized
     * @return The result container
     */
    SpatialSearchResultContainerReferenceType operator[](const IndexType Index)
    {
        KRATOS_ERROR_IF_NOT(this->HasResult(Index)) << "The result container does not exist for index: " << Index << std::endl;
        return *mPointResults[Index];
    }

    /**
     * @brief Operator []
     * @param Index The index to be initialized
     * @return The result container
     */
    const SpatialSearchResultContainerType& operator[](const IndexType Index) const
    {
        KRATOS_ERROR_IF_NOT(this->HasResult(Index)) << "The result container does not exist for index: " << Index << std::endl;
        return *mPointResults[Index];
    }

    /**
     * @brief Operator ()
     * @param Index The index to be initialized
     * @return The result container
     */
    SpatialSearchResultContainerReferenceType operator()(const IndexType Index)
    {
        KRATOS_ERROR_IF_NOT(this->HasResult(Index)) << "The result container does not exist for index: " << Index << std::endl;
        return *mPointResults[Index];
    }

    /**
     * @brief Operator ()
     * @param Index The index to be initialized
     * @return The result container
     */
    const SpatialSearchResultContainerType& operator()(const IndexType Index) const
    {
        KRATOS_ERROR_IF_NOT(this->HasResult(Index)) << "The result container does not exist for index: " << Index << std::endl;
        return *mPointResults[Index];
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
    std::size_t NumberOfSearchResults() const;

    /**
     * @brief Initialize the container
     * @param Index The index to be initialized
     */
    SpatialSearchResultContainerReferenceType InitializeResult(const IndexType Index);

    /**
     * @brief Check if coordinates are initialized
     * @param Index The index to be initialized
     * @return True if hash is initialized, false otherwise
    */
    bool HasResult(const IndexType Index) const;

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
}; // Class SpatialSearchResultContainerVector

/// input stream function
template <class TObjectType>
inline std::istream& operator>>(std::istream& rIStream,
                                SpatialSearchResultContainerVector<TObjectType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TObjectType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SpatialSearchResultContainerVector<TObjectType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.