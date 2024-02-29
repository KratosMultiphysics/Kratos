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

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               NodesArrayType const& ThisNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwElement(NewId, this->GetGeometry().Create(ThisNodes),
                                                     pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                               GeometryType::Pointer pGeom,
                                                               PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(
        new SteadyStatePwElement(NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType&   Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        if (Geom[i].SolutionStepsDataHas(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing variable WATER_PRESSURE on node " << Geom[i].Id() << std::endl;

        if (Geom[i].SolutionStepsDataHas(VOLUME_ACCELERATION) == false)
            KRATOS_ERROR << "missing VOLUME_ACCELERATION variable on node " << Geom[i].Id() << std::endl;

        if (Geom[i].HasDofFor(WATER_PRESSURE) == false)
            KRATOS_ERROR << "missing the dof for the variable WATER_PRESSURE "
                            "on node "
                         << Geom[i].Id() << std::endl;
    }

    // Verify ProcessInfo variables

    // Verify properties
    if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
        KRATOS_ERROR << "DENSITY_WATER does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
        KRATOS_ERROR << "POROSITY does not exist in the material properties or "
                        "has an invalid value at element"
                     << this->Id() << std::endl;

    if (TDim == 2) {
        // If this is a 2D problem, nodes must be in XY plane
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (Geom[i].Z() != 0.0)
                KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    // Verify specific properties
    if (Prop.Has(DYNAMIC_VISCOSITY) == false || Prop[DYNAMIC_VISCOSITY] < 0.0)
        KRATOS_ERROR << "DYNAMIC_VISCOSITY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_XX) == false || Prop[PERMEABILITY_XX] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XX does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_YY) == false || Prop[PERMEABILITY_YY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_YY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if (Prop.Has(PERMEABILITY_XY) == false || Prop[PERMEABILITY_XY] < 0.0)
        KRATOS_ERROR << "PERMEABILITY_XY does not exist in the material "
                        "properties or has an invalid value at element"
                     << this->Id() << std::endl;

    if constexpr (TDim > 2) {
        if (Prop.Has(PERMEABILITY_ZZ) == false || Prop[PERMEABILITY_ZZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (Prop.Has(PERMEABILITY_YZ) == false || Prop[PERMEABILITY_YZ] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_YZ does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;

        if (Prop.Has(PERMEABILITY_ZX) == false || Prop[PERMEABILITY_ZX] < 0.0)
            KRATOS_ERROR << "PERMEABILITY_ZX does not exist in the material "
                            "properties or has an invalid value at element"
                         << this->Id() << std::endl;
    }

    // Verify that the constitutive law has the correct dimension

    // Check constitutive law
    if (mRetentionLawVector.size() > 0) {
        return mRetentionLawVector[0]->Check(Prop, rCurrentProcessInfo);
    }

    return 0;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAll(MatrixType&        rLeftHandSideMatrix,
                                                         VectorType&        rRightHandSideVector,
                                                         const ProcessInfo& rCurrentProcessInfo,
                                                         const bool CalculateStiffnessMatrixFlag,
                                                         const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Previous definitions
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables, rCurrentProcessInfo);

    // create general parameters of retention law
    RetentionLaw::Parameters RetentionParameters(this->GetProperties(), rCurrentProcessInfo);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++) {
        // Compute GradNpT, B and StrainVector
        this->CalculateKinematics(Variables, GPoint);

        // Compute Nu and BodyAcceleration
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.Nu, Variables.NContainer, GPoint);
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, Variables.NContainer, Variables.VolumeAcceleration, GPoint);

        CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);

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
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                               ElementVariables& rVariables)
{
    KRATOS_TRY;

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               ElementVariables& rVariables,
                                                               unsigned int      GPoint)
{
    KRATOS_TRY;

    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
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
