//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "circular_convection_element_data.h"

namespace Kratos
{
namespace StabilizationValidationElementData
{
CircularConvectionElementData::CircularConvectionElementData(
    const GeometryType& rGeometry,
    const Properties& rProperties,
    const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo)
{
    KRATOS_TRY

    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_coordinates = rGeometry[i].Coordinates();
        for (IndexType j = 0; j < 2; ++j) {
            mNodalCoordinates(j, i) = r_coordinates[j];
        }
    }

    KRATOS_CATCH("");
}

const Variable<double>& CircularConvectionElementData::GetScalarVariable()
{
    return VELOCITY_POTENTIAL;
}

void CircularConvectionElementData::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CIRCULAR_CONVECTION_ROTATION_CLOCKWISE))
        << "CIRCULAR_CONVECTION_ROTATION_CLOCKWISE is not found int process info.\n";

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CIRCULAR_CONVECTION_ROTATION_CENTER))
        << "CIRCULAR_CONVECTION_ROTATION_CENTER is not found int process info.\n";

    const auto& r_geometry = rElement.GetGeometry();

    for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_POTENTIAL, r_node);
    }

    KRATOS_CATCH("");
}

void CircularConvectionElementData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mRotationFactor = (rCurrentProcessInfo[CIRCULAR_CONVECTION_ROTATION_CLOCKWISE]) ? 1.0 : -1.0;
    mRotationCenter = rCurrentProcessInfo[CIRCULAR_CONVECTION_ROTATION_CENTER];
    mEffectiveKinematicViscosity = 0.0;
    mReactionTerm = 0.0;
    mSourceTerm = 0.0;
}

void CircularConvectionElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    const ArrayD& gauss_coordinates = prod(mNodalCoordinates, rShapeFunctions);
    mEffectiveVelocity[0] = (-mRotationCenter[1] + gauss_coordinates[1]) * mRotationFactor;
    mEffectiveVelocity[1] = ( mRotationCenter[0] - gauss_coordinates[0]) * mRotationFactor;

    KRATOS_CATCH("");
}

} // namespace StabilizationValidationElementData

} // namespace Kratos