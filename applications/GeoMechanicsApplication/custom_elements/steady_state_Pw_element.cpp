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

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geometry   = this->GetGeometry();

    CheckUtilities::CheckDomainSize(r_geometry.DomainSize(), this->Id());

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (r_geometry[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << r_geometry[i].Id() << std::endl;

        if (r_geometry[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "missing VOLUME_ACCELERATION variable on node " << r_geometry[i].Id() << std::endl;

        if (r_geometry[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE "
                            "on node "
                         << r_geometry[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if (!r_properties.Has(DENSITY_WATER) || r_properties[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(POROSITY) || r_properties[POROSITY] < 0.0 || r_properties[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (r_geometry[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << r_geometry[i].Id()
                             << std::endl;
        }
    }

    // Verify specific properties
    if (!r_properties.Has(DYNAMIC_VISCOSITY) || r_properties[DYNAMIC_VISCOSITY] < 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(PERMEABILITY_XX) || r_properties[PERMEABILITY_XX] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(PERMEABILITY_YY) || r_properties[PERMEABILITY_YY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_YY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(PERMEABILITY_XY) || r_properties[PERMEABILITY_XY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if constexpr (TDim > 2) {
        if (!r_properties.Has(PERMEABILITY_ZZ) || r_properties[PERMEABILITY_ZZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (!r_properties.Has(PERMEABILITY_YZ) || r_properties[PERMEABILITY_YZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_YZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (!r_properties.Has(PERMEABILITY_ZX) || r_properties[PERMEABILITY_ZX] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZX does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;
    }

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
