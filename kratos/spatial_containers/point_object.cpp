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

template class PointObject<Node>;
template class PointObject<GeometricalObject>;
template class PointObject<Condition>;
template class PointObject<Element>;

} // namespace Kratos.