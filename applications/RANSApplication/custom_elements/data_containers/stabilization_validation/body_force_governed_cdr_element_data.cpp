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
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "body_force_governed_cdr_element_data.h"

namespace Kratos
{
namespace StabilizationValidationElementData
{
BodyForceGovernedCDRElementData::BodyForceGovernedCDRElementData(
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

const Variable<double>& BodyForceGovernedCDRElementData::GetScalarVariable()
{
    return VELOCITY_POTENTIAL;
}

void BodyForceGovernedCDRElementData::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();

    KRATOS_ERROR_IF_NOT(r_properties.Has(VELOCITY))
        << "VELOCITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(REACTION_COEFFICIENT))
        << "REACTION_COEFFICIENT is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_POTENTIAL, r_node);
    }

    KRATOS_CATCH("");
}

void BodyForceGovernedCDRElementData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_properties = this->GetProperties();
    const auto& r_velocity = r_properties[VELOCITY];
    mEffectiveVelocity[0] = r_velocity[0];
    mEffectiveVelocity[1] = r_velocity[1];

    mEffectiveKinematicViscosity = r_properties[DYNAMIC_VISCOSITY];
    mReactionTerm = r_properties[REACTION_COEFFICIENT];
}

void BodyForceGovernedCDRElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    const ArrayD& gauss_coordinates = prod(mNodalCoordinates, rShapeFunctions);

    // the source term is derrived to obtain following solution
    // phi = 100 * (x**2) * ((1-x)**2) * (y**2) * ((1-y)**2)

    mSourceTerm =
        100 * std::pow(gauss_coordinates[0], 2) * std::pow(gauss_coordinates[1], 2) *
            mReactionTerm * std::pow(gauss_coordinates[0] - 1, 2) *
            std::pow(gauss_coordinates[1] - 1, 2) +
        200 * std::pow(gauss_coordinates[0], 2) * gauss_coordinates[1] *
            mEffectiveVelocity[1] * std::pow(gauss_coordinates[0] - 1, 2) *
            (gauss_coordinates[1] - 1) * (2 * gauss_coordinates[1] - 1) +
        200 * gauss_coordinates[0] * std::pow(gauss_coordinates[1], 2) *
            mEffectiveVelocity[0] * (gauss_coordinates[0] - 1) *
            (2 * gauss_coordinates[0] - 1) * std::pow(gauss_coordinates[1] - 1, 2) -
        200 * mEffectiveKinematicViscosity *
            (std::pow(gauss_coordinates[0], 2) * std::pow(gauss_coordinates[1], 2) *
                 std::pow(gauss_coordinates[0] - 1, 2) +
             std::pow(gauss_coordinates[0], 2) * std::pow(gauss_coordinates[1], 2) *
                 std::pow(gauss_coordinates[1] - 1, 2) +
             4 * std::pow(gauss_coordinates[0], 2) * gauss_coordinates[1] *
                 std::pow(gauss_coordinates[0] - 1, 2) * (gauss_coordinates[1] - 1) +
             std::pow(gauss_coordinates[0], 2) * std::pow(gauss_coordinates[0] - 1, 2) *
                 std::pow(gauss_coordinates[1] - 1, 2) +
             4 * gauss_coordinates[0] * std::pow(gauss_coordinates[1], 2) *
                 (gauss_coordinates[0] - 1) * std::pow(gauss_coordinates[1] - 1, 2) +
             std::pow(gauss_coordinates[1], 2) * std::pow(gauss_coordinates[0] - 1, 2) *
                 std::pow(gauss_coordinates[1] - 1, 2));

    KRATOS_CATCH("");
}

} // namespace StabilizationValidationElementData

} // namespace Kratos