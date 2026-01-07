// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/steady_state_Pw_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                     pProperties, this->GetStressStatePolicy().Clone(),
                                                     this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwElement(NewId, pGeom, pProperties,
                                                     this->GetStressStatePolicy().Clone(),
                                                     this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();
    const auto& r_geometry   = this->GetGeometry();

    CheckUtilities::CheckDomainSize(r_geometry.DomainSize(), this->Id());

    CheckUtilities::CheckHasNodalSolutionStepData(
        r_geometry, {std::cref(WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)});
    CheckUtilities::CheckHasDofs(r_geometry, {std::cref(WATER_PRESSURE)});

    const CheckProperties check_properties(r_properties, "material properties",
                                           CheckProperties::Bounds::AllExclusive);
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllInclusive).Check(DENSITY_WATER);
    constexpr auto max_value_porosity = 1.0;
    check_properties.Check(POROSITY, max_value_porosity);
    check_properties.Check(DYNAMIC_VISCOSITY);

    if constexpr (TDim == 2) CheckUtilities::CheckForNonZeroZCoordinateIn2D(r_geometry);

    check_properties.CheckPermeabilityProperties(TDim);

    // Verify that the constitutive law has the correct dimension

    return RetentionLaw::Check(mRetentionLawVector, r_properties, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                         VectorType&        rRightHandSideVector,
                                                         const ProcessInfo& rCurrentProcessInfo,
                                                         bool CalculateStiffnessMatrixFlag,
                                                         bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const GeometryType&                             r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    const auto fluid_pressures = GeoTransportEquationUtilities::CalculateFluidPressures(
        Variables.NContainer, Variables.PressureVector);
    const auto relative_permeability_values = this->CalculateRelativePermeabilityValues(fluid_pressures);
    const auto bishop_coefficients = this->CalculateBishopCoefficients(fluid_pressures);
    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, Variables.detJContainer);
    const auto degrees_of_saturation = this->CalculateDegreesOfSaturation(fluid_pressures);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++) {
        // Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        Variables.RelativePermeability = relative_permeability_values[GPoint];
        Variables.DegreeOfSaturation   = degrees_of_saturation[GPoint];
        Variables.BishopCoefficient    = bishop_coefficients[GPoint];

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                               ElementVariables& rVariables)
{
    noalias(rLeftHandSideMatrix) += GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
        rVariables.GradNpT, rVariables.DynamicViscosityInverse, rVariables.PermeabilityMatrix,
        rVariables.RelativePermeability, rVariables.IntegrationCoefficient);
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               ElementVariables& rVariables,
                                                               unsigned int      GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template class SteadyStatePwElement<2, 3>;
template class SteadyStatePwElement<2, 4>;
template class SteadyStatePwElement<3, 4>;
template class SteadyStatePwElement<3, 8>;

template class SteadyStatePwElement<2, 6>;
template class SteadyStatePwElement<2, 8>;
template class SteadyStatePwElement<2, 9>;
template class SteadyStatePwElement<2, 10>;
template class SteadyStatePwElement<2, 15>;
template class SteadyStatePwElement<3, 10>;
template class SteadyStatePwElement<3, 20>;
template class SteadyStatePwElement<3, 27>;

} // Namespace Kratos
