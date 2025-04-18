//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos Roig
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/master_slave_constraint.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @ingroup KratosCore
 * @class CompareElementsAndConditionsUtility
 * @brief Utility class to compare elements and conditions in Kratos.
 * @details This utility provides static methods to retrieve the registered names of geometries,
 * elements, and conditions. It supports both direct references and pointer-based access.
 * @author Carlos Roig and Vicente Mataix Ferrandiz
 */
class CompareElementsAndConditionsUtility
{
public:
    ///@name Operations
    ///@{

    /**
     * @brief Retrieves the registered name of a geometry.
     * @details This method retrieves the registered name associated with the given geometry object.
     * @param rGeometry The geometry object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Geometry<Node>& rGeometry,
        std::string& rName
        );

    /**
     * @brief Retrieves the registered name of a geometry pointer.
     * @details This method retrieves the registered name associated with the given geometry pointer.
     * @param pGeometry A pointer to the geometry object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void GetRegisteredName(
        const Geometry<Node>* pGeometry,
        std::string& rName)
    {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pGeometry, rName);
    }

    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Element& rElement,
        std::string& rName
        );

    /**
     * @brief Retrieves the registered name of an element pointer.
     * @details This method retrieves the registered name associated with the given element pointer.
     * @param pElement A pointer to the element object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void GetRegisteredName(
        const Element* pElement,
        std::string& rName
        )
    {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pElement, rName);
    }

    /**
     * @brief Retrieves the registered name of a condition.
     * @details This method retrieves the registered name associated with the given condition object.
     * @param rCondition The condition object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const Condition& rCondition,
        std::string& rName
        );

    /**
     * @brief Retrieves the registered name of a condition pointer.
     * @details This method retrieves the registered name associated with the given condition pointer.
     * @param pCondition A pointer to the condition object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void GetRegisteredName(
        const Condition* pCondition,
        std::string& rName
        )
    {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pCondition, rName);
    }

    /**
     * @brief Retrieves the registered name of a master-slave constraint.
     * @details This method retrieves the registered name associated with the given master-slave constraint object.
     * @param rConstraint The master-slave constraint object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
     static void KRATOS_API(KRATOS_CORE) GetRegisteredName(
        const MasterSlaveConstraint& rConstraint,
        std::string& rName
        );

    /**
     * @brief Retrieves the registered name of a master-slave constraint pointer.
     * @details This method retrieves the registered name associated with the given master-slave constraint pointer.
     * @param pConstraint A pointer to the master-slave constraint object whose registered name is to be retrieved.
     * @param rName A reference to a string where the registered name will be stored.
     */
    static void GetRegisteredName(
        const MasterSlaveConstraint* pConstraint,
        std::string& rName
        )
    {
        CompareElementsAndConditionsUtility::GetRegisteredName(*pConstraint, rName);
    }

    ///@}
private:
    ///@name Private Operations
    ///@{

    ///@}
};
///@} // Kratos Classes
} // namespace Kratos.