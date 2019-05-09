// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/rocco_ring_element.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
RoccoRingElement::RoccoRingElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

RoccoRingElement::RoccoRingElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
RoccoRingElement::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_intrusive<RoccoRingElement>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
RoccoRingElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<RoccoRingElement>(NewId, pGeom,
                                               pProperties);
}

RoccoRingElement::~RoccoRingElement() {}

void RoccoRingElement::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != msLocalSize) rResult.resize(msLocalSize);

    for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * 3;
    rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index + 2] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

void RoccoRingElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != msLocalSize) rElementalDofList.resize(msLocalSize);


    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * 3;
        rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] =
            this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] =
            this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }
}

} // namespace Kratos.
