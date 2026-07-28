//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo Crescenzio
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup MPMApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class BruteForceMaterialPointLocator
 * @ingroup MPMApplication
 * @brief Utility class to find a material point entity based on a location
 * @details Based on the location of a point, the corresponding material point
 * entity (element or condition) is found and its id is returned
 */
class KRATOS_API(MPM_APPLICATION) BruteForceMaterialPointLocator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BruteForceMaterialPointLocator
    KRATOS_CLASS_POINTER_DEFINITION(BruteForceMaterialPointLocator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit BruteForceMaterialPointLocator(
        ModelPart& rModelPart
        ) : mrModelPart(rModelPart)
    {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function finds a material point element based on a location
     * @param rThePoint the location to search
     * @param tolerance tolerance for finding closest material point element
     * @return Id of the found element. -1 if no element was found
     */
    int FindElement(
        const Point& rThePoint,
        const double tolerance = 1e-6
        ) const;

    /**
     * @brief This function finds a material point condition based on a location
     * @param rThePoint the location to search
     * @param tolerance tolerance for finding closest material point condition
     * @return Id of the found condition. -1 if no condition was found
     */
    int FindCondition(
        const Point& rThePoint,
        const double tolerance = 1e-6
        ) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function finds an object based on a location
     * @param rObjects the objects to search
     * @param rObjectType type of the object: "Element"/"Condition"
     * @param rThePoint the location to search
     * @param rObjectId Id of the found object
     * @param tolerance tolerance for finding closest material point condition
     */
    template<typename TObjectType>
    void FindObject(
        const TObjectType& rObjects,
        const std::string& rObjectType,
        const Point& rThePoint,
        int& rObjectId,
        const double tolerance
        ) const;

    ///@}

    class KRATOS_API(MPM_APPLICATION) MPMMinDistanceReduction
    {
    public:
        using value_type = std::tuple<double,int>;
        using return_type = int;

        value_type mValue = {std::numeric_limits<double>::max(), -1};

        /// access to reduced value
        return_type GetValue() const
        {
            return get<1>(mValue);
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type value){
            mValue = get<0>(value) < get<0>(mValue) ? value : mValue;
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const MPMMinDistanceReduction& rOther)
        {
            KRATOS_CRITICAL_SECTION
            LocalReduce(rOther.mValue);
        }
    };

}; // Class BruteForceMaterialPointLocator

///@}

///@} addtogroup block

}  // namespace Kratos.
