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
    explicit BruteForceMaterialPointLocator(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

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
        const double tolerance = 1e-6) const;

    /**
     * @brief This function finds a material point condition based on a location
     * @param rThePoint the location to search
     * @param tolerance tolerance for finding closest material point condition
     * @return Id of the found condition. -1 if no condition was found
     */
    int FindCondition(
        const Point& rThePoint,
        const double tolerance = 1e-6) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "BruteForceMaterialPointLocator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "BruteForceMaterialPointLocator";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
     * @param rObjectType type of the object => "Element"/"Condition"
     * @param rThePoint the location to search
     * @param rObjectId Id of the found condition. -1 if no object was found
     * @param tolerance tolerance local-coordinates for IsInside
     */
    template<typename TObjectType>
    void FindObject(
        const TObjectType& rObjects,
        const std::string& rObjectType,
        const Point& rThePoint,
        int& rObjectId,
        const double tolerance) const;

    ///@}

}; // Class BruteForceMaterialPointLocator

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const BruteForceMaterialPointLocator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
