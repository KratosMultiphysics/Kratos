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

#if !defined(KRATOS_REGISTER_KRATOS_COMPONENTS_FOR_GEOMETRY_H_INCLUDED )
#define  KRATOS_REGISTER_KRATOS_COMPONENTS_FOR_GEOMETRY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<Geometry<Node<3>>>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, Geometry<Node<3>> const& ThisComponent);

}  // namespace Kratos.

#endif // KRATOS_REGISTER_KRATOS_COMPONENTS_FOR_GEOMETRY_H_INCLUDED  defined
