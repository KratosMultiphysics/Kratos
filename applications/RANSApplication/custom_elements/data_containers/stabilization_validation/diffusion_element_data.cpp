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
#include "includes/cfd_variables.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "diffusion_element_data.h"

namespace Kratos
{
namespace StabilizationValidationElementData
{
DiffusionElementData::DiffusionElementData(
    const GeometryType& rGeometry,
    const Properties& rProperties,
    const ProcessInfo& rProcessInfo)
    : BaseType(rGeometry, rProperties, rProcessInfo)
{
}

const Variable<double>& DiffusionElementData::GetScalarVariable()
{
    return VELOCITY_POTENTIAL;
}

void DiffusionElementData::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();

    const auto& r_properties = rElement.GetProperties();

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in properties with id "
        << r_properties.Id() << " in element with id " << rElement.Id() << ".\n";

    KRATOS_ERROR_IF(r_properties[DYNAMIC_VISCOSITY] <= 0.0)
        << "DYNAMIC_VISCOSITY is needs to be greater than zero in properties "
           "with id "
        << r_properties.Id() << " in element with id " << rElement.Id()
        << " [ DYNAMIC_VISCOSITY = " << r_properties[DYNAMIC_VISCOSITY] << " ].\n";

    for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_POTENTIAL, r_node);
    }

    KRATOS_CATCH("");
}

void DiffusionElementData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mEffectiveKinematicViscosity = this->GetProperties()[DYNAMIC_VISCOSITY];
    mReactionTerm = 0.0;
    mSourceTerm = 0.0;
    mEffectiveVelocity.clear();
}

void DiffusionElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
}

} // namespace StabilizationValidationElementData

} // namespace Kratos