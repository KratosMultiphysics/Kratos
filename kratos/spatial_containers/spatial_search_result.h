//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/indexed_object.h"
#include "input_output/logger.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/** 
 * @class SpatialSearchResult
 * @brief This class is the result of the spatial searches.
 * @details It provides:
 *  - A global pointer to the object found.
 *  - Distance to the object if IsDistanceCalculated() is true
 *  - IsObjectFound if for example search nearest fails or not
 * @ingroup KratosCore
 * @author Pooyan Dadvand
 */
template <typename TObjectType>
class SpatialSearchResult
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialSearchResult
    KRATOS_CLASS_POINTER_DEFINITION(SpatialSearchResult);

    /// Global pointer definition of TObjectType
    using TPointerType = GlobalPointer<TObjectType>;

    /// Define index type
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialSearchResult()
        : mpObject(nullptr),
          mDistance(0.0),
          mIsObjectFound(false),
          mIsDistanceCalculated(false)
    {
    }

    /// Constructor with the resulted object
    SpatialSearchResult(
        TObjectType* pObject,
        const int Rank = 0
        ) : mpObject(pObject, Rank),
            mDistance(0.0),
            mIsObjectFound(false),
            mIsDistanceCalculated(false)
    {
        if (mpObject.get() != nullptr)
            mIsObjectFound = true;
    }

    /// Copy constructor.
    SpatialSearchResult(SpatialSearchResult const& /* Other */) = default;

    /// Move constructor.
    SpatialSearchResult(SpatialSearchResult&& /* Other */) = default;

    /// Destructor.
    virtual ~SpatialSearchResult(){}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SpatialSearchResult& operator=(SpatialSearchResult const& /*Other*/) = default;

    ///@}
    ///@name Operations
    ///@{

    /// Reset the result
    void Reset()
    {
        mpObject = nullptr;
        mDistance = 0.0;
        mIsObjectFound = false;
        mIsDistanceCalculated = false;
    }

    ///@}
    ///@name Access
    ///@{

    /**
    * @brief Get the ID of the object.
    * @details If the object is derived from IndexedObject, retrieves its ID. Otherwise, returns 0.
    * @return IndexType - The ID of the object or 0 if not derived from IndexedObject.
    */
    IndexType Id() const
    {
        // Get Id if the object is derived from IndexedObject
        if constexpr (std::is_base_of_v<IndexedObject, TObjectType>) {
            return mpObject->Id();
        } else { // Otherwise 0
            KRATOS_ERROR << "Object does not provide an ID" << std::endl;
            return 0;
        }
    }

    /**
    * @brief Get the ID of the object.
    * @details This function essentially performs the same operation as the Id() function. If the object is derived from IndexedObject, retrieves its ID. Otherwise, returns 0.
    * @note Consider refactoring to avoid redundancy with Id() function.
    * @return IndexType - The ID of the object or 0 if not derived from IndexedObject.
    */
    IndexType GetId() const
    {
        // Get Id if the object is derived from IndexedObject
        if constexpr (std::is_base_of_v<IndexedObject, TObjectType>) {
            return mpObject->Id();
        } else { // Otherwise 0
            KRATOS_ERROR << "Object does not provide an ID" << std::endl;
            return 0;
        }
    }

    /// Returns the global pointer to the object
    TPointerType Get() {
        return mpObject;
    }

    /// Returns a const global pointer to the object
    TPointerType const Get() const {
        return mpObject;
    }

    /// Set the object to be pointed
    void Set(TObjectType* pObject) {
        mpObject = pObject;
        mIsObjectFound = true;
    }

    /// Getting the result distance
    double GetDistance() const {
        return mDistance;
    }

    /// Setting the result distance
    void SetDistance(const double TheDistance) {
        mDistance = TheDistance;
        mIsDistanceCalculated = true;
    }

    /// Getting if the object is found
    bool GetIsObjectFound() const {
        return mIsObjectFound;
    }

    /// Getting if the ibject is found
    void SetIsObjectFound(const bool IsObjectFound) {
        mIsObjectFound = IsObjectFound;
    }

    /// Getting if the distance is calculated
    bool GetIsDistanceCalculated() const {
        return mIsDistanceCalculated;
    }

    /// Setting if the distance is calculated
    void SetIsDistanceCalculated(const bool IsDistanceCalculated) {
        mIsDistanceCalculated = IsDistanceCalculated;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Returns true if the object is set
    bool IsObjectFound() const {
        return mIsObjectFound;
    }

    /// Returns true if the distance is set
    bool IsDistanceCalculated() const {
        return mIsDistanceCalculated;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "SpatialSearchResult";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SpatialSearchResult";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const 
    {

    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    TPointerType mpObject;      /// The object found
    double mDistance;           /// The distance to the object
    bool mIsObjectFound;        /// If the object is found
    bool mIsDistanceCalculated; /// If the distance is calculated

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    /**
     * @brief Saves the object data to a serializer.
     * @param rSerializer The serializer to save the data to.
     */
    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Object", mpObject);
        rSerializer.save("Distance", mDistance);
        rSerializer.save("Is Object Found", mIsObjectFound);
        rSerializer.save("Is Distance Calculated", mIsDistanceCalculated);
    }

    /**
     * @brief Loads the data of an object from a Serializer.
     * @param rSerializer The Serializer to load the data from.
     */
    void load(Serializer& rSerializer)
    {
        rSerializer.load("Object", mpObject);
        rSerializer.load("Distance", mDistance);
        rSerializer.load("Is Object Found", mIsObjectFound);
        rSerializer.load("Is Distance Calculated", mIsDistanceCalculated);
    }

    ///@}

}; // Class SpatialSearchResult

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <typename TObjectType>
inline std::istream& operator >> (std::istream& rIStream,
                SpatialSearchResult<TObjectType>& rThis){
                    return rIStream;
                }

/// output stream function
template <typename TObjectType>
inline std::ostream& operator << (std::ostream& rOStream,
                const SpatialSearchResult<TObjectType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.