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

///@}
///@name Kratos Classes
///@{

/**
 * @class SpatialSearchResultContainerVector
 * @brief Spatial search result container vector.
 * @details This class is used to store the results of a spatial search of several points, in a vector which assumes the order.
 * @tparam TObjectType The type of the object.
 * @tparam TSpatialSearchCommunication The type of spatial search communication considered.
 * @ingroup KratosCore
 * @author Vicente Mataix Ferrandiz
 */
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication = SpatialSearchCommunication::SYNCHRONOUS>
class KRATOS_API(KRATOS_CORE) SpatialSearchResultContainerVector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResultContainerVector
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResultContainerVector);

    /// The index map type
    using IndexMapType = std::unordered_map<IndexType, IndexType>;

    /// The spatial search result container type
    using SpatialSearchResultContainerType = SpatialSearchResultContainer<TObjectType, TSpatialSearchCommunication>;

    /// Spatial search result type
    using SpatialSearchResultType = typename SpatialSearchResultContainerType::SpatialSearchResultType;

    /// The global pointer result type
    using GlobalPointerResultType = typename SpatialSearchResultContainerType::GlobalPointerResultType;

    /// The global results vector
    using GlobalResultsVector = typename SpatialSearchResultContainerType::GlobalResultsVector;

    /// The global pointer communicator type
    using GlobalPointerCommunicatorType = typename SpatialSearchResultContainerType::GlobalPointerCommunicatorType;

    /// The global pointer communicator pointer type
    using GlobalPointerCommunicatorPointerType = typename SpatialSearchResultContainerType::GlobalPointerCommunicatorPointerType;

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
        const SpatialSearchResultContainerType& operator*() const {
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

    /// Definition of the max value
    static constexpr double MaxValue = std::numeric_limits<double>::max();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SpatialSearchResultContainerVector() = default;

    /// Destructor.
    virtual ~SpatialSearchResultContainerVector();

    ///@}
    ///@name Operators
    ///@{

   /**
     * @brief Operator []
     * @param Index The index to be retrieved.
     * @return The result container.
     */
    SpatialSearchResultContainerReferenceType operator[](const IndexType Index)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(HasResult(Index)) << "Index " << Index << " not available. Size: " << mPointResults.size() << std::endl;
        return *mPointResults[Index];
    }

    /**
     * @brief Operator []
     * @param Index The index to be retrieved.
     * @return The result container.
     */
    const SpatialSearchResultContainerType& operator[](const IndexType Index) const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(HasResult(Index)) << "Index " << Index << " not available. Size: " << mPointResults.size() << std::endl;
        return *mPointResults[Index];
    }

    /**
     * @brief Operator ()
     * @param Index The index to be retrieved.
     * @return The result container pointer.
     */
    SpatialSearchResultContainerPointerType operator()(const IndexType Index)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(HasResult(Index)) << "Index " << Index << " not available. Size: " << mPointResults.size() << std::endl;
        return mPointResults[Index];
    }

    /**
     * @brief Operator ()
     * @param Index The index to be retrieved.
     * @return The result container pointer.
     */
    const SpatialSearchResultContainerType* operator()(const IndexType Index) const
    {
        KRATOS_DEBUG_ERROR_IF_NOT(HasResult(Index)) << "Index " << Index << " not available. Size: " << mPointResults.size() << std::endl;
        return mPointResults[Index];
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
     * @brief Returns the number of points results.
     * @return The number of points results.
     */
    std::size_t NumberOfSearchResults() const;

    /**
     * @brief Initialize the container.
     * @return The result container.
     */
    SpatialSearchResultContainerReferenceType InitializeResult();

    /**
     * @brief Initialize the container.
     * @param NumberOfResults The number of results to be initialized.
     */
    void InitializeResults(const std::size_t NumberOfResults);

    /**
     * @brief Check if coordinates are initialized.
     * @param Index The index to be initialized.
     * @return True if hash is initialized, false otherwise.
    */
    bool HasResult(const IndexType Index) const;

    /**
     * @brief Clear the containers.
     * @details This method clears the containers.
     */
    void Clear();

    /**
     * @brief Generate the global pointer communicator.
     * @param rDataCommunicator The data communicator considered.
     */
    void GenerateGlobalPointerCommunicator(const DataCommunicator& rDataCommunicator);

    /**
     * @brief Generate the global pointer communicator.
     * @param rAllGlobalResults The vector containing all the global results.
     * @param rDataCommunicator The data communicator considered.
     */
    void GenerateGlobalPointerCommunicator(
        const GlobalResultsVector& rAllGlobalResults,
        const DataCommunicator& rDataCommunicator
        );

    /**
     * @brief Synchronize all container between partitions
     * @details This method synchronizes all the container between partitions
     * @param rDataCommunicator The data communicator considered
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
     * @brief Retrieves the global distances.
     * @param rResults The vector containing all the distances.
     */
    void GetDistances(std::vector<std::vector<double>>& rResults);

    /**
     * @brief Retrieves if is local the entity.
     * @param rResults The vector containing all the booleans showing is local the entity
     * @param rDataCommunicator The data communicator.
     */
    void GetResultIsLocal(
        std::vector<std::vector<bool>>& rResults,
        const DataCommunicator& rDataCommunicator
        );

    /**
     * @brief Retrieves the rank of the entity.
     * @param rResults The vector containing all the ranks of the entity.
     */
    void GetResultRank(std::vector<std::vector<int>>& rResults);

    /**
     * @brief Retrieves if is active the entity.
     * @param rResults The vector containing all the booleans showing is active the entity
     * @param rDataCommunicator The data communicator.
     */
    void GetResultIsActive(
        std::vector<std::vector<bool>>& rResults,
        const DataCommunicator& rDataCommunicator
        );

    /**
     * @brief Retrieves if inside the geometry.
     * @param rResults The vector containing all the booleans showing is inside the geometry.
     * @param rPoints The points container.
     * @param rDataCommunicator The data communicator.
     * @param Tolerance The tolerance considered.
     */
    template<typename TPointContainerType>
    void GetResultIsInside(
        std::vector<std::vector<bool>>& rResults,
        const TPointContainerType& rPoints,
        const DataCommunicator& rDataCommunicator,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        )
    {
        // Get the iterator to the points
        auto it_point_begin = rPoints.begin();

        // Define the coordinates vector
        const std::size_t number_of_global_solutions = mPointResults.size();
        if (rResults.size() != number_of_global_solutions) {
            rResults.resize(number_of_global_solutions);
        }

        // Get the is inside
        Point point;
        for(std::size_t i = 0; i<number_of_global_solutions; ++i) {
            auto& r_partial_results = mPointResults[i];
            auto& r_global_results = r_partial_results->GetGlobalResults();

            // If local point
            const bool is_local_point = r_partial_results->IsLocalPoint();

            if (is_local_point) {
                auto it_point = it_point_begin + r_partial_results->GetLocalIndex();
                noalias(point) = it_point->Coordinates();
            } else {
                for (unsigned int i = 0; i < 3; i++) {
                    point[i] = -MaxValue;
                }
            }
            for (unsigned int i = 0; i < 3; i++) {
                point[i] = rDataCommunicator.MaxAll(point[i]);
            }

            // Call Apply to get the proxy
            auto proxy = this->Apply([&Tolerance, &point](GlobalPointerResultType& rGP) -> bool {
                auto p_object = rGP->Get();
                if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
                    auto& r_geometry = p_object->GetGeometry();
                    Point::CoordinatesArrayType aux_coords;
                    return r_geometry.IsInside(point, aux_coords, Tolerance);
                } else if constexpr (std::is_same<TObjectType, Node>::value) {
                    KRATOS_ERROR << "Nodes do not provide is inside. Not possible to compute is inside for point: " << point[0]<< "\t" << point[1] << "\t" << point[2] << " with tolerance " << Tolerance << std::endl;
                    return false;
                } else {
                    KRATOS_ERROR << "Not implemented yet. Not possible to compute is inside for point: " << point[0]<< "\t" << point[1] << "\t" << point[2] << " with tolerance " << Tolerance << std::endl;
                    return false;
                }
            });

            // Get the is inside
            const std::size_t number_of_gp = r_global_results.size();
            auto& r_is_inside = rResults[i];
            r_is_inside.resize(number_of_gp);
            for(std::size_t j = 0; j < number_of_gp; ++j) {
                auto& r_gp = r_global_results(j);
                r_is_inside[j] = rDataCommunicator.MaxAll(proxy.Get(r_gp));
            }
        }
    }

    /**
     * @brief Considers the global pointer communicator to get the shape functions of the resulting object.
     * @param rResults A vector to store the shape functions of the resulting object.
     * @param rPoints The points container.
     */
    template<typename TPointContainerType>
    void GetResultShapeFunctions(
        std::vector<std::vector<Vector>>& rResults,
        const TPointContainerType& rPoints,
        const DataCommunicator& rDataCommunicator
        )
    {
        // Get the iterator to the points
        auto it_point_begin = rPoints.begin();

        // Define the coordinates vector
        const std::size_t number_of_global_solutions = mPointResults.size();
        if (rResults.size() != number_of_global_solutions) {
            rResults.resize(number_of_global_solutions);
        }

        // Get the shape functions
        Point point;
        for(std::size_t i = 0; i<number_of_global_solutions; ++i) {
            auto& r_partial_results = mPointResults[i];
            auto& r_global_results = r_partial_results->GetGlobalResults();

            // If local point
            const bool is_local_point = r_partial_results->IsLocalPoint();

            if (is_local_point) {
                auto it_point = it_point_begin + r_partial_results->GetLocalIndex();
                noalias(point) = it_point->Coordinates();
            } else {
                for (unsigned int i = 0; i < 3; i++) {
                    point[i] = -MaxValue;
                }
            }
            for (unsigned int i = 0; i < 3; i++) {
                point[i] = rDataCommunicator.MaxAll(point[i]);
            }

            // Call Apply to get the proxy
            auto proxy = this->Apply([&point](GlobalPointerResultType& rGP) -> Vector {
                auto p_object = rGP->Get();
                if constexpr (std::is_same<TObjectType, GeometricalObject>::value) {
                    auto& r_geometry = p_object->GetGeometry();
                    Vector N(r_geometry.size());
                    array_1d<double, 3> local_coordinates;
                    r_geometry.PointLocalCoordinates(local_coordinates, point);
                    r_geometry.ShapeFunctionsValues(N, local_coordinates);
                    return N;
                } else if constexpr (std::is_same<TObjectType, Node>::value) {
                    KRATOS_ERROR << "Nodes do not provide shape functions. Not possible to compute shape functions for point: " << point[0]<< "\t" << point[1] << "\t" << point[2] << std::endl;
                    Vector N;
                    return N;
                } else {
                    KRATOS_ERROR << "Not implemented yet. Not possible to compute shape functions for point: " << point[0]<< "\t" << point[1] << "\t" << point[2] << std::endl;
                    Vector N;
                    return N;
                }
            });

            // Get the shape functions
            const std::size_t number_of_gp = r_global_results.size();
            auto& r_shape_functions = rResults[i];
            r_shape_functions.resize(number_of_gp);
            for(std::size_t j = 0; j < number_of_gp; ++j) {
                auto& r_gp = r_global_results(j);
                r_shape_functions[j] = proxy.Get(r_gp);
            }
        }
    }

    /**
     * @brief Considers the global pointer communicator to get the indices of the resulting object.
     * @param rResults A vector to store the indices of the resulting object.
     */
    void GetResultIndices(std::vector<std::vector<IndexType>>& rResults);

    /**
     * @brief Considers the global pointer communicator to get the indices of the nodes of the resulting object.
     * @param rResults A vector to store the indices of the nodes of the resulting object.
     */
    void GetResultNodeIndices(std::vector<std::vector<std::vector<IndexType>>>& rResults);

    /**
     * @brief Considers the global pointer communicator to get the partition indices of the nodes of the resulting object.
     * @param rResults A vector to store the partition indices of the resulting object.
     */
    void GetResultPartitionIndices(std::vector<std::vector<std::vector<int>>>& rResults);

    /**
     * @brief Considers the global pointer communicator to get the coordinates of the resulting object.
     * @param rResults A vector to store the coordinates of the resulting object.
     */
    void GetResultCoordinates(std::vector<std::vector<std::vector<array_1d<double, 3>>>>& rResults);

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

    /**
     * @brief Accessor for mpGlobalPointerCommunicator.
     * @details This method returns the GlobalPointerCommunicatorPointer mpGlobalPointerCommunicator.
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

    ContainerType mPointResults;                                                /// The results of each point

    GlobalPointerCommunicatorPointerType mpGlobalPointerCommunicator = nullptr; /// Global pointer to the communicator

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Copies values from a vector of GlobalPointerResultType to global results vectors.
     * @details This function copies values from a vector of GlobalPointerResultType to global results vectors based on the provided active results and result global sizes. It iterates through the input vector of GlobalPointerResultType and assigns values to the corresponding global results vectors.
     * @param rGlobalResults The vector of GlobalPointerResultType to copy values from.
     * @param rActiveResults The vector of active results indices.
     * @param rResultGlobalSize The vector of result global sizes.
     */
    void CopyingValuesToGlobalResultsVector(
        const std::vector<GlobalPointerResultType>& rGlobalResults,
        const std::vector<std::size_t>& rActiveResults,
        const std::vector<std::size_t>& rResultGlobalSize
        );

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

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
inline std::istream& operator>>(std::istream& rIStream,
                                SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TObjectType, SpatialSearchCommunication TSpatialSearchCommunication>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SpatialSearchResultContainerVector<TObjectType, TSpatialSearchCommunication>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.