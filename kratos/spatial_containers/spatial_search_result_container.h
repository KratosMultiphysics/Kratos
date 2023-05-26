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

// External includes

// Project includes
#include "utilities/pointer_communicator.h"
#include "spatial_containers/spatial_search_result.h"
#include "containers/pointer_vector_set.h"
#include "includes/indexed_object.h"

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

    /// Global pointer definition of TObjectType
    using TPointerType = GlobalPointer<TObjectType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rDataCommunicator The data communicator
     */
    SpatialSearchResultContainer(const DataCommunicator& rDataCommunicator);

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
     */
    void SynchronizeAll();

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
    
    const DataCommunicator& mrDataCommunicator;                                                     /// The data communicator
    PointerVectorSet<TObjectType,
                     IndexedObject,
                     std::less<typename IndexedObject::result_type>,
                     std::equal_to<typename IndexedObject::result_type>,
                     typename TObjectType::Pointer,
                     std::vector< typename TObjectType::Pointer >
                     > mLocalPointers;                                                              /// Local pointers of the container
    GlobalPointersVector<TObjectType> mGlobalPointers;                                              /// Global pointers of the container
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

///@}

///@} addtogroup block

}  // namespace Kratos.