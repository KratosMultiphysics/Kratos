// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

// System includes

// External includes

// Project includes
#include "shell_7p_element.hpp"

namespace Kratos
{

Shell7pElement::Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

Shell7pElement::Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer Shell7pElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_intrusive<Shell7pElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer Shell7pElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_intrusive<Shell7pElement>(NewId, pGeom, pProperties);
}

} // namespace Kratos