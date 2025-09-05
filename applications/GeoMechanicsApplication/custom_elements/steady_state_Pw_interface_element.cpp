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
#include "custom_elements/steady_state_Pw_interface_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                        NodesArrayType const& ThisNodes,
                                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwInterfaceElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwInterfaceElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                        GeometryType::Pointer pGeom,
                                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwInterfaceElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwInterfaceElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geometry   = this->GetGeometry();

    if (this->Id() < 1)
        KRATOS_ERROR << "Element found with Id 0 or negative, element: " << this->Id() << std::endl;

    CheckUtilities::CheckHasNodalSolutionStepData(
        r_geometry, {std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)});
    CheckUtilities::CheckHasDofs(r_geometry, {std::cref(WATER_PRESSURE)});

    // Verify specific properties
    if (!r_properties.Has(MINIMUM_JOINT_WIDTH) || r_properties[MINIMUM_JOINT_WIDTH] <= 0.0)
        KRATOS_ERROR << "MINIMUM_JOINT_WIDTH has Key zero, is not defined or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(TRANSVERSAL_PERMEABILITY) || r_properties[TRANSVERSAL_PERMEABILITY] < 0.0)
        KRATOS_ERROR << "TRANSVERSAL_PERMEABILITY has Key zero, is not defined "
                        "or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(DYNAMIC_VISCOSITY) || r_properties[DYNAMIC_VISCOSITY] <= 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY has Key zero, is not defined or has "
                        "an invalid value at element"
                     << this->Id() << std::endl;

    // Verify properties
    if (!r_properties.Has(DENSITY_WATER) || r_properties[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (!r_properties.Has(POROSITY) || r_properties[POROSITY] < 0.0 || r_properties[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    // Verify the constitutive law
    if (!r_properties.Has(CONSTITUTIVE_LAW))
        KRATOS_ERROR << "CONSTITUTIVE_LAW has Key zero or is not defined at "
                        "element "
                     << this->Id() << std::endl;
    if (!r_properties[CONSTITUTIVE_LAW])
        KRATOS_ERROR << "A constitutive law needs to be specified for the "
                        "element "
                     << this->Id() << std::endl;

    // Check constitutive law
    ierr = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, this->GetGeometry(), rCurrentProcessInfo);

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                  VectorType& rRightHandSideVector,
                                                                  const ProcessInfo& CurrentProcessInfo,
                                                                  bool CalculateStiffnessMatrixFlag,
                                                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           r_properties = this->GetProperties();
    const GeometryType&                             r_geometry   = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geometry.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        r_geometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    r_geometry.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    r_geometry.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, r_geometry, r_properties, CurrentProcessInfo);

    // VG: TODO
    // Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = r_properties[MINIMUM_JOINT_WIDTH];

    // Auxiliary variables
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, detJContainer);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer, GPoint);

        this->template CalculateShapeFunctionsGradients<Matrix>(
            Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
            DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

        // Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

        InterfaceElementUtilities::FillPermeabilityMatrix(
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, r_properties[TRANSVERSAL_PERMEABILITY]);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

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
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                        InterfaceElementVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                        InterfaceElementVariables& rVariables,
                                                                        unsigned int GPoint)
{
    KRATOS_TRY

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template class SteadyStatePwInterfaceElement<2, 4>;
template class SteadyStatePwInterfaceElement<3, 6>;
template class SteadyStatePwInterfaceElement<3, 8>;

} // Namespace Kratos
