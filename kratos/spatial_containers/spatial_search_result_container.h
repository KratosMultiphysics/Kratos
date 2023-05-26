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

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Add a result to the container
     * @param rResult The result
     */
    void AddResult(SpatialSearchResult<TObjectType>& rResult);

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
    TObjectType,
    TFunctorType // TODO: Unfortunately this is deprecated in c++17, so we will have to change this call in the future
    > Apply(TFunctorType&& UserFunctor)
    {
        // Check if the communicator has been created
        KRATOS_ERROR_IF(mpGlobalPointerCommunicator == nullptr) << "The communicator has not been created." << std::endl;

        // Apply the user-provided function
        return mpGlobalPointerCommunicator->Apply(std::forward<TFunctorType>(UserFunctor));
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
    
    PointerVectorSet<TObjectType,
                     IndexedObject,
                     std::less<typename IndexedObject::result_type>,
                     std::equal_to<typename IndexedObject::result_type>,
                     typename TObjectType::Pointer,
                     std::vector< typename TObjectType::Pointer >
                     > mLocalPointers;                                                              /// Local pointers of the container
    GlobalPointersVector<TObjectType> mGlobalPointers;                                              /// Global pointers of the container
    std::unordered_map<IndexType, double> mLocalDistances;                                          /// The local distances 
    typename GlobalPointerCommunicator<TObjectType>::Pointer mpGlobalPointerCommunicator = nullptr; /// Global pointer to the communicator 

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
     * @param rCoordinates The coordinates
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator[](const array_1d<double, 3>& rCoordinates)
    {
        return mPointResults[Hash(rCoordinates)];
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
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist." << std::endl;
        return it->second;
    }

    /**
     * @brief Operator ()
     * @param rCoordinates The coordinates
     * @return The result container
     */
    SpatialSearchResultContainer<TObjectType>& operator()(const array_1d<double, 3>& rCoordinates)
    {
        return mPointResults[Hash(rCoordinates)];
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
        KRATOS_ERROR_IF(it == mPointResults.end()) << "The result container does not exist." << std::endl;
        return it->second;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the container
     * @param rCoordinates The coordinates
     */
    SpatialSearchResultContainer<TObjectType>& InitializeResult(const array_1d<double, 3>& rCoordinates);

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
    
    std::unordered_map<HashType, SpatialSearchResultContainer<TObjectType>> mPointResults;  /// The results of each point

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
}; // Class SpatialSearchResultContainer

///@}

///@} addtogroup block

}  // namespace Kratos.