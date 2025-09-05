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
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PointObject
 * @ingroup KratosCore
 * @brief Custom Point container to be used by the search
 * @details It stores the pointer of a certain object
 * @author Vicente Mataix Ferrandiz
 */
template<class TObject>
class KRATOS_API(KRATOS_CORE) PointObject
    : public Point
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    using BaseType = Point;

    /// Definition of the object type
    using ObjectType = TObject;

    /// Counted pointer of PointObject
    KRATOS_CLASS_POINTER_DEFINITION( PointObject );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointObject();

    /**
     * @brief Constructor that takes the coordinates of the point.
     * @param Coords The coordinates of the point.
     */
    PointObject(const array_1d<double, 3>& Coords);

    /**
     * @brief Constructor with object
     * @param pObject The pointer to the object
     */
    PointObject(typename TObject::Pointer pObject);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint();

    ///@}
    ///@name Acess
    ///@{

    /**
     * @brief Returns the geometry associated to the point
     * @return mrGeometry The reference to the geometry associated to the point
     */
    typename TObject::Pointer pGetObject() const;

    /**
     * @brief Sets the object associated to the point
     * @param pObject The pointer to the object
     */
    void pSetObject(typename TObject::Pointer pObject);

    /**
     * @brief This method checks everything is right
     */
    void Check() const;

    ///@}
private:
    ///@}
    ///@name Member Variables
    ///@{

    typename TObject::Pointer mpObject = nullptr;

    ///@}

}; // Class PointObject

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@}addtogroup block

}  // namespace Kratos.