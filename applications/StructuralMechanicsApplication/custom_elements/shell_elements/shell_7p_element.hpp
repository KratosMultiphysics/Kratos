// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) Shell7pElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell7pElement);

    // A private default constructor necessary for serialization
    Shell7pElement() = default;

    Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry);

    Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~Shell7pElement() override = default;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;
};

} // namespace Kratos