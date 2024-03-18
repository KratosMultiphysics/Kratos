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

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwInterfaceElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                        NodesArrayType const& ThisNodes,
                                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwInterfaceElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwInterfaceElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                        GeometryType::Pointer pGeom,
                                                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwInterfaceElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwInterfaceElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType&   Geom = this->GetGeometry();

    if (this->Id() < 1)
        KRATOS_ERROR << "Element found with Id 0 or negative, element: " << this->Id() << std::endl;

    // Verify dof variables
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (Geom[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if (Geom[i].SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable DT_WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if (Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "missing variable VOLUME_ACCELERATION on node " << Geom[i].Id() << std::endl;

        if (Geom[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;
    }

    // Verify specific properties
    if (Prop.Has(MINIMUM_JOINT_WIDTH) == false || Prop[MINIMUM_JOINT_WIDTH] <= 0.0)
        KRATOS_ERROR << "MINIMUM_JOINT_WIDTH has Key zero, is not defined or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(TRANSVERSAL_PERMEABILITY) == false || Prop[TRANSVERSAL_PERMEABILITY] < 0.0)
        KRATOS_ERROR << "TRANSVERSAL_PERMEABILITY has Key zero, is not defined "
                        "or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] <= 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY has Key zero, is not defined or has "
                        "an invalid value at element"
                     << this->Id() << std::endl;

    // Verify properties
    if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    // Verify the constitutive law
    if (Prop.Has(CONSTITUTIVE_LAW) == false)
        KRATOS_ERROR << "CONSTITUTIVE_LAW has Key zero or is not defined at "
                        "element "
                     << this->Id() << std::endl;

    if (Prop[CONSTITUTIVE_LAW] != NULL) {
        // Check constitutive law
        ierr = Prop[CONSTITUTIVE_LAW]->Check(Prop, this->GetGeometry(), rCurrentProcessInfo);
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the "
                        "element "
                     << this->Id() << std::endl;

    return ierr;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                                  VectorType& rRightHandSideVector,
                                                                  const ProcessInfo& CurrentProcessInfo,
                                                                  const bool CalculateStiffnessMatrixFlag,
                                                                  const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        Geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Element variables
    InterfaceElementVariables Variables;
    this->InitializeElementVariables(Variables, Geom, Prop, CurrentProcessInfo);

    // VG: TODO
    // Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = Prop[MINIMUM_JOINT_WIDTH];

    // Auxiliary variables
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), CurrentProcessInfo);

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
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, Prop[TRANSVERSAL_PERMEABILITY]);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, detJContainer[GPoint]);

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                        InterfaceElementVariables& rVariables)
{
    KRATOS_TRY;

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwInterfaceElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                        InterfaceElementVariables& rVariables,
                                                                        unsigned int GPoint)
{
    KRATOS_TRY;

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
template class SteadyStatePwInterfaceElement<2, 4>;
template class SteadyStatePwInterfaceElement<3, 6>;
template class SteadyStatePwInterfaceElement<3, 8>;

} // Namespace Kratos
