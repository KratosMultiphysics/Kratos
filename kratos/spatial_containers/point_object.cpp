//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |  (   | |  (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "spatial_containers/point_object.h"

namespace Kratos
{

template<class TObject>
PointObject<TObject>::PointObject():
    BaseType()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
PointObject<TObject>::PointObject(const array_1d<double, 3>& Coords)
    :BaseType(Coords)
{

}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
PointObject<TObject>::PointObject(typename TObject::Pointer pObject):
    mpObject(pObject)
{
    UpdatePoint();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
void PointObject<TObject>::UpdatePoint()
{
    if constexpr (std::is_same<TObject, Node>::value) {
        noalias(this->Coordinates()) = mpObject->Coordinates();
    } else if constexpr (std::is_same<TObject, GeometricalObject>::value || std::is_same<TObject, Condition>::value || std::is_same<TObject, Element>::value) {
        noalias(this->Coordinates()) = mpObject->GetGeometry().Center().Coordinates();
    } else if constexpr (std::is_same<TObject, Geometry<Node>>::value || std::is_same<TObject, Condition>::value || std::is_same<TObject, Element>::value) {
        noalias(this->Coordinates()) = mpObject->Center().Coordinates();
    } else {
        static_assert((std::is_same<TObject, Node>::value || std::is_same<TObject, GeometricalObject>::value || std::is_same<TObject, Condition>::value || std::is_same<TObject, Element>::value), "PointObject is implemented for Node, GeometricalObject, Condition and Element");
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
typename TObject::Pointer PointObject<TObject>::pGetObject() const
{
    return mpObject;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
void PointObject<TObject>::pSetObject(typename TObject::Pointer pObject)
{
    mpObject = pObject;
    UpdatePoint();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TObject>
void PointObject<TObject>::Check() const
{
    KRATOS_TRY;

    auto aux_coord = std::make_shared<array_1d<double, 3>>(this->Coordinates());
    KRATOS_ERROR_IF(!aux_coord) << "Coordinates no initialized in the PointObject class" << std::endl;
    KRATOS_ERROR_IF(mpObject.get() == nullptr) << "TEntity no initialized in the PointObject class" << std::endl;

    KRATOS_CATCH("Error checking the PointObject");
}

/***********************************************************************************/
/***********************************************************************************/

template class PointObject<Node>;
template class PointObject<GeometricalObject>;
template class PointObject<Condition>;
template class PointObject<Element>;
template class PointObject<Geometry<Node>>;

} // namespace Kratos.