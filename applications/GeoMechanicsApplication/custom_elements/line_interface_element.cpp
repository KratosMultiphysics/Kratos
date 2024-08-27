// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "line_interface_element.h"

namespace Kratos
{

LineInterfaceElement::LineInterfaceElement() = default;

Element::Pointer LineInterfaceElement::Create(IndexType               NewId,
                                              const NodesArrayType&   rNodes,
                                              PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

Element::Pointer LineInterfaceElement::Create(IndexType               NewId,
                                              GeometryType::Pointer   pGeometry,
                                              PropertiesType::Pointer pProperties) const
{
    return {make_intrusive<LineInterfaceElement>(NewId, pGeometry, pProperties)};
}

} // namespace Kratos
