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

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "circular_convection_element_data_derivatives.h"

namespace Kratos
{

namespace StabilizationValidationElementData
{

/***************************************************************/
/************************ Element Data *************************/
/***************************************************************/

const Variable<double>& CircularConvectionElementDataDerivatives::Data::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

void CircularConvectionElementDataDerivatives::Data::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CIRCULAR_CONVECTION_ROTATION_CLOCKWISE))
        << "CIRCULAR_CONVECTION_ROTATION_CLOCKWISE is not found int process info.\n";

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CIRCULAR_CONVECTION_ROTATION_CENTER))
        << "CIRCULAR_CONVECTION_ROTATION_CENTER is not found int process info.\n";

    for (const auto& r_node : rElement.GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_1_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

/***************************************************************/
/*********************** Phi Derivative  ***********************/
/***************************************************************/

CircularConvectionElementDataDerivatives::PhiDerivative::PhiDerivative(
    const Data& rData,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
    : BaseType(NodeIndex, DirectionIndex, rData.GetGeometry(), W, rN, rdNdX, WDerivative, DetJDerivative, rdNdXDerivative),
      mrData(rData)
{

}

const Variable<double>& CircularConvectionElementDataDerivatives::PhiDerivative::GetDerivativeVariable() const
{
    return VELOCITY_POTENTIAL;
}

array_1d<double, 2> CircularConvectionElementDataDerivatives::PhiDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(2);
}

double CircularConvectionElementDataDerivatives::PhiDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    return 0.0;
}

double CircularConvectionElementDataDerivatives::PhiDerivative::CalculateReactionTermDerivative() const
{
    return 0.0;
}

double CircularConvectionElementDataDerivatives::PhiDerivative::CalculateSourceTermDerivative() const
{
    return 0.0;
}

/***************************************************************/
/*********************** Shape Derivative **********************/
/***************************************************************/

const Variable<double>& CircularConvectionElementDataDerivatives::ShapeDerivative::GetDerivativeVariable() const
{
    switch (this->mDirectionIndex) {
    case 0:
        return SHAPE_SENSITIVITY_X;
        break;
    case 1:
        return SHAPE_SENSITIVITY_Y;
        break;
    case 2:
        return SHAPE_SENSITIVITY_Z;
        break;
    default:
        return Variable<double>::StaticObject();
    };
}

array_1d<double, 2> CircularConvectionElementDataDerivatives::ShapeDerivative::CalculateEffectiveVelocityDerivative() const
{
    const double N_value = this->mrN[this->mNodeIndex] * mrData.mRotationFactor;
    return ArrayD({N_value * (this->mDirectionIndex == 1),
                   -N_value * (this->mDirectionIndex == 0)});
}

double CircularConvectionElementDataDerivatives::ShapeDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    return 0.0;
}

double CircularConvectionElementDataDerivatives::ShapeDerivative::CalculateReactionTermDerivative() const
{
    return 0.0;
}

double CircularConvectionElementDataDerivatives::ShapeDerivative::CalculateSourceTermDerivative() const
{
    return 0.0;
}

} // namespace StabilizationValidationElementData

} // namespace Kratos
