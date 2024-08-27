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
#include "custom_utilities/dof_utilities.h"

namespace Kratos
{

LineInterfaceElement::LineInterfaceElement() = default;

void LineInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    DofsVectorType dofs;
    GetDofList(dofs, rCurrentProcessInfo);
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(dofs);
}

void LineInterfaceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    // At this point we only look at the U dofs, so we leave the water pressure nodes empty.
    rElementalDofList = Geo::DofUtilities::ExtractUPwDofsFromNodes(
        GetGeometry(), Geometry<Node>(), GetGeometry().WorkingSpaceDimension());
}

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
