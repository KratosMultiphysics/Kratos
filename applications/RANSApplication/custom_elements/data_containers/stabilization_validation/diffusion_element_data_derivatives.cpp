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
#include "includes/cfd_variables.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "diffusion_element_data_derivatives.h"

namespace Kratos
{

namespace StabilizationValidationElementData
{

/***************************************************************/
/************************ Element Data *************************/
/***************************************************************/

const Variable<double>& DiffusionElementDataDerivatives::Data::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

void DiffusionElementDataDerivatives::Data::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_properties = rElement.GetProperties();

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in properties with id "
        << r_properties.Id() << " in element with id " << rElement.Id() << ".\n";

    KRATOS_ERROR_IF(r_properties[DYNAMIC_VISCOSITY] <= 0.0)
        << "DYNAMIC_VISCOSITY is needs to be greater than zero in properties "
           "with id "
        << r_properties.Id() << " in element with id " << rElement.Id()
        << " [ DYNAMIC_VISCOSITY = " << r_properties[DYNAMIC_VISCOSITY] << " ].\n";

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

DiffusionElementDataDerivatives::PhiDerivative::PhiDerivative(
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

const Variable<double>& DiffusionElementDataDerivatives::PhiDerivative::GetDerivativeVariable() const
{
    return VELOCITY_POTENTIAL;
}

array_1d<double, 2> DiffusionElementDataDerivatives::PhiDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(2);
}

double DiffusionElementDataDerivatives::PhiDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    return 0.0;
}

double DiffusionElementDataDerivatives::PhiDerivative::CalculateReactionTermDerivative() const
{
    return 0.0;
}

double DiffusionElementDataDerivatives::PhiDerivative::CalculateSourceTermDerivative() const
{
    return 0.0;
}

/***************************************************************/
/*********************** Shape Derivative **********************/
/***************************************************************/

const Variable<double>& DiffusionElementDataDerivatives::ShapeDerivative::GetDerivativeVariable() const
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

array_1d<double, 2> DiffusionElementDataDerivatives::ShapeDerivative::CalculateEffectiveVelocityDerivative() const
{
    return ZeroVector(2);
}

double DiffusionElementDataDerivatives::ShapeDerivative::CalculateEffectiveKinematicViscosityDerivative() const
{
    return 0.0;
}

double DiffusionElementDataDerivatives::ShapeDerivative::CalculateReactionTermDerivative() const
{
    return 0.0;
}

double DiffusionElementDataDerivatives::ShapeDerivative::CalculateSourceTermDerivative() const
{
    return 0.0;
}

} // namespace StabilizationValidationElementData

} // namespace Kratos
