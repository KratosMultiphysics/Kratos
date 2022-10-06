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
#include "custom_elements/brake_element_3D2N.hpp"
#include "includes/define.h"
#include "cable_net_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
BrakeElement3D2N::BrakeElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

BrakeElement3D2N::BrakeElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
BrakeElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                         PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<BrakeElement3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
BrakeElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<BrakeElement3D2N>(NewId, pGeom,
            pProperties);
}

BrakeElement3D2N::~BrakeElement3D2N() {}


double BrakeElement3D2N::ReturnReferenceLength() const
{
    return GetProperties()[EXPERIMENTAL_LENGTH];   // defined to fit experimental curves
}

double BrakeElement3D2N::ReturnCurrentLength() const
{
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    Vector current_displacement_global = ZeroVector(msLocalSize);
    GetValuesVector(current_displacement_global, 0);

    Vector current_displacement_local =  prod(trans(transformation_matrix), current_displacement_global);

    return ReturnReferenceLength() + current_displacement_local[3] - current_displacement_local[0];
}


void BrakeElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
}
void BrakeElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
}
} // namespace Kratos.
