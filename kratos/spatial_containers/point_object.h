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
class PointObject
    : public Point
{
public:

    ///@name Type Definitions
    ///@{

    /// Base class definition
    using BaseType = Point;

    /// Counted pointer of PointObject
    KRATOS_CLASS_POINTER_DEFINITION( PointObject );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointObject():
        BaseType()
    {
    }

    /**
     * @brief Constructor that takes the coordinates of the point.
     * @param Coords The coordinates of the point.
     */
    PointObject(const array_1d<double, 3>& Coords)
        :BaseType(Coords)
    {}


    /**
     * @brief Constructor with object
     * @param pObject The pointer to the object
     */
    PointObject(typename TObject::Pointer pObject):
        mpObject(pObject)
    {
        UpdatePoint();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint()
    {
        if constexpr (std::is_same<TObject, Node>::value) {
            noalias(this->Coordinates()) = mpObject->Coordinates();
        } else if constexpr (std::is_same<TObject, GeometricalObject>::value || std::is_same<TObject, Condition>::value || std::is_same<TObject, Element>::value) {
            noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
        } else {
            static_assert((std::is_same<TObject, Node>::value || std::is_same<TObject, GeometricalObject>::value || std::is_same<TObject, Condition>::value || std::is_same<TObject, Element>::value), "PointObject is implemented for Node, GeometricalObject, Condition and Element");
        }
    }

    ///@}
    ///@name Acess
    ///@{

    /**
     * @brief Returns the geometry associated to the point
     * @return mrGeometry The reference to the geometry associated to the point
     */
    typename TObject::Pointer pGetObject()
    {
        return mpObject;
    }

    /**
     * @brief Sets the object associated to the point
     * @param pObject The pointer to the object
     */
    void pSetObject(typename TObject::Pointer pObject)
    {
        mpObject = pObject;
        UpdatePoint();
    }

    /**
     * @brief This method checks everything is right
     */
    void Check()
    {
        KRATOS_TRY;

        auto aux_coord = std::make_shared<array_1d<double, 3>>(this->Coordinates());
        KRATOS_ERROR_IF(!aux_coord) << "Coordinates no initialized in the PointObject class" << std::endl;
        KRATOS_ERROR_IF(mpObject.get() == nullptr) << "TEntity no initialized in the PointObject class" << std::endl;

        KRATOS_CATCH("Error checking the PointObject");
    }

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