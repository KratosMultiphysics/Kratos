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
#include "geo_mechanics_application_variables.h"
#include "custom_elements/geo_curved_beam_element.hpp"
#include "custom_utilities/element_utilities.hpp"

#include <math.h>

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoCurvedBeamElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
{
    KRATOS_INFO("0-GeoCurvedBeamElement::Create") << this->Id() << std::endl;

    return Element::Pointer( new GeoCurvedBeamElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer GeoCurvedBeamElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties ) const
{
    KRATOS_INFO("0-GeoCurvedBeamElement::Create") << this->Id() << std::endl;

    return Element::Pointer( new GeoCurvedBeamElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_ERROR << "GeoCurvedBeamElement::CalculateLeftHandSide not implemented" << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
int GeoCurvedBeamElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    KRATOS_INFO("0-GeoCurvedBeamElement::Check") << this->Id() << std::endl;

    // Base class checks for positive area and Id > 0
    int ierr = GeoStructuralBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& rProp = this->GetProperties();

    if ( I33.Key() == 0 ||
         rProp.Has( I33 ) == false ||
         rProp[I33] < 0.0 )
        KRATOS_ERROR << "I33 has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if ( CROSS_AREA.Key() == 0 ||
         rProp.Has( CROSS_AREA ) == false ||
         rProp[CROSS_AREA] < 0.0 )
        KRATOS_ERROR << "CROSS_AREA has Key zero, is not defined or has an invalid value at element: "
                     << this->Id()
                     << std::endl;

    if (TDim > 2) {
        if ( TORSIONAL_INERTIA.Key() == 0 ||
             rProp.Has( TORSIONAL_INERTIA ) == false ||
             rProp[TORSIONAL_INERTIA] < 0.0 )
             KRATOS_ERROR << "TORSIONAL_INERTIA has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;

        if ( I22.Key() == 0 ||
             rProp.Has( I22 ) == false ||
             rProp[I22] < 0.0 )
             KRATOS_ERROR << "I22 has Key zero, is not defined or has an invalid value at element: "
                          << this->Id()
                          << std::endl;
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::Check") << std::endl;

    return 0;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    unsigned int index = 0;
    rRotationalInertia[index++] = rProp[TORSIONAL_INERTIA];
    rRotationalInertia[index++] = rProp[I22];
    rRotationalInertia[index++] = rProp[I33];

    KRATOS_INFO("1-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    SetRotationalInertiaVector(const PropertiesType& rProp, Vector& rRotationalInertia) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    if ( rRotationalInertia.size() != N_DOF_NODE_ROT )
        rRotationalInertia.resize( N_DOF_NODE_ROT, false );

    rRotationalInertia[0] = rProp[I33];

    KRATOS_INFO("1-GeoCurvedBeamElement::SetRotationalInertiaVector") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateMassMatrix( MatrixType& rMassMatrix,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateMassMatrix") << std::endl;

    //Resizing mass matrix
    if ( rMassMatrix.size1() != N_DOF_ELEMENT )
        rMassMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rMassMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Defining shape functions and the determinant of the jacobian at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJContainer(IntegrationPoints.size());
    rGeom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    //Defining necessary variables
    const double Density = rProp[DENSITY];

    Vector RotationalInertia;
    this->SetRotationalInertiaVector(rProp, RotationalInertia);

     unsigned int index = 0;
    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); GPoint++ ) {
        //calculating weighting coefficient for integration
        double IntegrationCoefficient = IntegrationPoints[GPoint].Weight() * detJContainer[GPoint];

        //Adding contribution to Mass matrix
        // loop over nodes
        for (unsigned int node = 0; node < TNumNodes; ++node) {
            // displacement degrees of freedom
            for (unsigned int dof = 0; dof < N_DOF_NODE_DISP; ++dof) {
                const unsigned int i = index++;
                rMassMatrix(i, i) += Density * NContainer(GPoint, node) * IntegrationCoefficient;
            }

            // rotational degrees of freedom
            for (unsigned int dof = 0; dof < N_DOF_NODE_ROT; ++dof) {
                const unsigned int i = index++;
                rMassMatrix(i, i) += RotationalInertia[dof] * NContainer(GPoint, node) * IntegrationCoefficient;
            }
        }
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateMassMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateStiffnessMatrix( MatrixType& rStiffnessMatrix,
                              const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateStiffnessMatrix") << std::endl;

    //Resizing mass matrix
    if ( rStiffnessMatrix.size1() != N_DOF_ELEMENT )
        rStiffnessMatrix.resize( N_DOF_ELEMENT, N_DOF_ELEMENT, false );
    noalias( rStiffnessMatrix ) = ZeroMatrix( N_DOF_ELEMENT, N_DOF_ELEMENT );

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPoints.size());
    Vector detJContainer(IntegrationPoints.size());
    //rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);
    this->ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,
                                     ConstitutiveParameters,
                                     rGeom,
                                     rProp,
                                     CurrentProcessInfo);

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Nu, GradNu, B and StrainVector
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);
        noalias(Variables.GradNu) = DN_DXContainer[GPointAlong];

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNu );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;

            this->CalculateDeterminantJacobian( GPointCross,
                                                Variables,
                                                DetJacobianMatrix );

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());

            //Compute stiffness matrix
            this->CalculateAndAddLHS(rStiffnessMatrix, Variables);
        }
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateStiffnessMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAll") << this->Id() << std::endl;

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType&
        IntegrationPointsAlong = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPointsAlong.size());

    Vector detJContainer(IntegrationPointsAlong.size());

    // rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,
    //                                                detJContainer,
    //                                                mThisIntegrationMethod);
    this->ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,
                                                   detJContainer,
                                                   mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      rGeom,
                                      rProp,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPointsAlong.size(); ++GPointAlong) {
        KRATOS_INFO("GPointAlong") << GPointAlong << std::endl;

        //Compute Nu, GradNu, B and StrainVector
        noalias(Variables.GradNu) = DN_DXContainer[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);
        GeoElementUtilities::
            CalculateNuMatrix<TDim, TNumNodes>(Variables.NuTot,
                                               NContainer,
                                               GPointAlong);

        KRATOS_INFO("Variables.GradNu") << Variables.GradNu << std::endl;

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNu );
        KRATOS_INFO("Variables.TransformationMatrix") << Variables.TransformationMatrix << std::endl;

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < GetCrossNumberIntegrationPoints(); ++GPointCross) {
            KRATOS_INFO("GPointCross") << GPointCross << std::endl;

            int GPoint= GPointAlong * GetCrossNumberIntegrationPoints() + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;

            this->CalculateDeterminantJacobian( GPointCross,
                                                Variables,
                                                DetJacobianMatrix );
            KRATOS_INFO("DetJacobianMatrix") << DetJacobianMatrix << std::endl;

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);
            KRATOS_INFO("InvertDetJacobianMatrix") << InvertDetJacobianMatrix << std::endl;

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            KRATOS_INFO("Variables") << Variables.StrainVector << std::endl;

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPointsAlong[GPointAlong].Weight());

            KRATOS_INFO("08-GeoCurvedBeamElement::CalculateAll") << std::endl;

            //Contributions to the left hand side
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAll") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    InitializeElementVariables( ElementVariables& rVariables,
                                ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                const GeometryType& rGeom,
                                const PropertiesType& rProp,
                                const ProcessInfo& CurrentProcessInfo ) const
{
    KRATOS_TRY

    KRATOS_INFO("0-GeoCurvedBeamElement::InitializeElementVariables") << std::endl;

    GeoStructuralBaseElement<TDim,TNumNodes>::
        InitializeElementVariables( rVariables,
                                    rConstitutiveParameters,
                                    rGeom,
                                    rProp,
                                    CurrentProcessInfo );

    rVariables.HalfThickness = 0.5 * std::sqrt(12.0 * (rProp[I33] / rProp[CROSS_AREA]));

    KRATOS_INFO("1-GeoCurvedBeamElement::InitializeElementVariables") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateRHS( VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateRHS") << std::endl;

    //Previous definitions
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( mThisIntegrationMethod );

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(IntegrationPoints.size());
    Vector detJContainer(IntegrationPoints.size());
    // rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);
    this->ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, mThisIntegrationMethod);

    //Constitutive Law parameters
    ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom,rProp,CurrentProcessInfo);
    ConstitutiveParameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveParameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables( Variables,
                                      ConstitutiveParameters,
                                      rGeom,
                                      rProp,
                                      CurrentProcessInfo );

    //Loop over integration points
    for (unsigned int GPointAlong = 0; GPointAlong < IntegrationPoints.size(); ++GPointAlong) {
        //Compute Nu, GradNu, B and StrainVector
        noalias(Variables.GradNu) = DN_DXContainer[GPointAlong];
        noalias(Variables.Nu)     = row(NContainer, GPointAlong);
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Variables.NuTot, NContainer, GPointAlong);

        this->CalculateTransformationMatrix( Variables.TransformationMatrix,
                                             Variables.GradNu );

        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim,TNumNodes>( Variables.BodyAcceleration,
                                                               NContainer,
                                                               Variables.VolumeAcceleration,
                                                               GPointAlong );

        for (unsigned int GPointCross = 0; GPointCross < N_POINT_CROSS; ++GPointCross) {
            int GPoint= GPointAlong * N_POINT_CROSS + GPointCross;

            BoundedMatrix<double,TDim, TDim> DetJacobianMatrix;
            this->CalculateDeterminantJacobian(GPointCross,
                                               Variables,
                                               DetJacobianMatrix);

            double detJacobian;
            BoundedMatrix<double,TDim, TDim> InvertDetJacobianMatrix;
            MathUtils<double>::InvertMatrix(DetJacobianMatrix,
                                            InvertDetJacobianMatrix,
                                            detJacobian);

            this->CalculateBMatrix(Variables.B, GPointCross, InvertDetJacobianMatrix, Variables);

            this->CalculateStrainVector(Variables);

            //Compute constitutive tensor
            ConstitutiveParameters.SetStressVector(mStressVector[GPoint]);
            mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient = this->CalculateIntegrationCoefficient(GPointCross,
                                                                                     detJacobian,
                                                                                     IntegrationPoints[GPointAlong].Weight());

            //Contributions to the right hand side
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        }
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateRHS") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, 
                       ElementVariables& rVariables ) const
{
    KRATOS_TRY

    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddLHS") << std::endl;

    noalias(rVariables.UVoigtMatrix) = prod(trans(rVariables.B),rVariables.ConstitutiveMatrix);
    rLeftHandSideMatrix =  prod(rVariables.UVoigtMatrix,rVariables.B)
                         * rVariables.IntegrationCoefficient;

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddLHS") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddRHS(VectorType& rRightHandSideVector,
                       ElementVariables& rVariables) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddRHS") << std::endl;

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddBodyForce(rRightHandSideVector, rVariables);

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddRHS") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddBodyForce(VectorType& rRightHandSideVector,
                             ElementVariables& rVariables) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddBodyForce") << std::endl;

    const PropertiesType& rProp = this->GetProperties();
    const double &density = rProp[DENSITY];

    KRATOS_INFO("density") << density << std::endl;
    KRATOS_INFO("rVariables.NuTot") << rVariables.NuTot << std::endl;
    KRATOS_INFO("rVariables.BodyAcceleration") << rVariables.BodyAcceleration << std::endl;
    KRATOS_INFO("rVariables.IntegrationCoefficient") << rVariables.IntegrationCoefficient << std::endl;
    KRATOS_INFO("rRightHandSideVector") << rRightHandSideVector << std::endl;

    //Distribute body force block vector into elemental vector
    rRightHandSideVector +=  density
                           * prod(trans(rVariables.NuTot), rVariables.BodyAcceleration)
                           * rVariables.IntegrationCoefficient;

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddBodyForce") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector,
                                  ElementVariables& rVariables) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAndAddStiffnessForce") << std::endl;

    //Distribute stiffness block vector into elemental vector
    rRightHandSideVector -=  prod(trans(rVariables.B),rVariables.StressVector)
                           * rVariables.IntegrationCoefficient;

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAndAddStiffnessForce") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAngleAtNode(unsigned int GPoint,
                         const BoundedMatrix<double,TNumNodes, TNumNodes> & DN_DXContainer) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAngleAtNode") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += DN_DXContainer(GPoint, node) * rGeom[node].X0();
        dy += DN_DXContainer(GPoint, node) * rGeom[node].Y0();
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAngleAtNode") << std::endl;

    return atan2(dx, -dy);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateAngleAtGaussPoint(const Matrix &GradNu) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateAngleAtGaussPoint") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();

    double dx = 0;
    double dy = 0;

    // loop over nodes
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        dx += GradNu(node, 0) * rGeom[node].X0();
        dy += GradNu(node, 0) * rGeom[node].Y0();
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateAngleAtGaussPoint, dy, dx") << dy << ", " << dx << std::endl;

    return atan2(dy, dx);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<2,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNu ) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateTransformationMatrix") << std::endl;

    const double phi = CalculateAngleAtGaussPoint(GradNu);
    const double cosPhi = cos(phi);
    const double sinPhi = sin(phi);

    TransformationMatrix(0,0) = cosPhi * cosPhi;
    TransformationMatrix(0,1) = sinPhi * sinPhi;
    TransformationMatrix(0,2) = -2.0 * sinPhi * cosPhi;

    TransformationMatrix(1,0) = TransformationMatrix(0,1);
    TransformationMatrix(1,1) = TransformationMatrix(0,0);
    TransformationMatrix(1,2) = -TransformationMatrix(0,2);

    TransformationMatrix(2,0) = sinPhi * cosPhi;
    TransformationMatrix(2,1) = -TransformationMatrix(2,0);
    TransformationMatrix(2,2) = TransformationMatrix(0,0) - TransformationMatrix(0,1);

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateTransformationMatrix") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< >
void GeoCurvedBeamElement<3,3>::
    CalculateTransformationMatrix( Matrix &TransformationMatrix,
                                   const Matrix &GradNu ) const
{
    KRATOS_TRY

    KRATOS_ERROR << "Undefined dimension in CalculateTransformationMatrix" << std::endl;;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateNodalCrossDirection(Matrix& NodalCrossDirection) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateNodalCrossDirection") << std::endl;

    BoundedMatrix<double, TNumNodes, TNumNodes> DN_DXContainer;
    GeoElementUtilities::CalculateNewtonCotesShapeFunctionsGradients(DN_DXContainer);
    KRATOS_INFO("NewtonCotes_DN_DXContainer") << DN_DXContainer << std::endl;

    for (unsigned int node = 0; node < TNumNodes; ++node) {
        double phi = CalculateAngleAtNode(node, DN_DXContainer);
        KRATOS_INFO("phi") << phi << std::endl;

        NodalCrossDirection(0, node) = cos(phi);
        NodalCrossDirection(1, node) = sin(phi);
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateNodalCrossDirection") << std::endl;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateDeterminantJacobian(unsigned int GPointCross,
                                  const ElementVariables &rVariables,
                                  BoundedMatrix<double,TDim, TDim> &DeterminantJacobian) const
{
    KRATOS_TRY
    // details can be found in
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"
    // Eqs.(20) - (22)
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateDeterminantJacobian") << std::endl;

    const GeometryType& rGeom = this->GetGeometry();
    const double &t = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0/sqrt(3), 1.0/sqrt(3)};

    noalias(DeterminantJacobian) = ZeroMatrix(TDim, TDim);

    KRATOS_INFO("rVariables.GradNu")         << rVariables.GradNu << std::endl;
    KRATOS_INFO("rVariables.NodalCrossDirection") << rVariables.NodalCrossDirection << std::endl;

    for (unsigned int node=0; node < TNumNodes; ++node) {
        DeterminantJacobian(0, 0) += rVariables.GradNu(node, 0) * rGeom[node].X0() + t * CrossEta[GPointCross] * rVariables.NodalCrossDirection(0, node);
        DeterminantJacobian(1, 0) += rVariables.GradNu(node, 0) * rGeom[node].Y0() + t * CrossEta[GPointCross] * rVariables.NodalCrossDirection(1, node);

        DeterminantJacobian(0, 1) += rVariables.Nu(node) * t * rVariables.NodalCrossDirection(0, node);
        DeterminantJacobian(1, 1) += rVariables.Nu(node) * t * rVariables.NodalCrossDirection(1, node);
    }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateDeterminantJacobian") << std::endl;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateBMatrix( Matrix &BTransformed,
                      unsigned int GPointCross,
                      const BoundedMatrix<double,TDim, TDim> &InvJ,
                      ElementVariables &rVariables ) const
{
    KRATOS_TRY
    // Details of derivation of linear part of B-Matrix can be found in:
    // "1. Geometrically non-linear formulation for the axisymmetric shell elements"

    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateBMatrix") << std::endl;

    const double &t = rVariables.HalfThickness;
    const std::vector<double> CrossEta{-1.0/sqrt(3), 1.0/sqrt(3)};

    Matrix B;
    B.resize(VoigtSize, N_DOF_ELEMENT, false);
    noalias(B) = ZeroMatrix(VoigtSize, N_DOF_ELEMENT);

    for (unsigned int node=0; node < TNumNodes; ++node) {
        const unsigned int index = N_DOF_NODE * node;
        const unsigned int index_x = index + INDEX_2D_BEAM_X;
        const unsigned int index_y = index + INDEX_2D_BEAM_Y;

        // Parts of B-matrix due to u_x and u_y (u, v)
        B(INDEX_2D_BEAM_XX, index_x) = rVariables.GradNu(node, 0) * InvJ(0, 0);
        B(INDEX_2D_BEAM_YY, index_y) = rVariables.GradNu(node, 0) * InvJ(1, 0);
        B(INDEX_2D_BEAM_XY, index_x) = B(INDEX_2D_BEAM_YY, index_y);
        B(INDEX_2D_BEAM_XY, index_y) = B(INDEX_2D_BEAM_XX, index_x);

        // Parts of B-matrix due to rotation (theta)
        const unsigned int index_t = index + INDEX_2D_BEAM_T;

        // eps_xx due to theta
        const double term_xx = InvJ(0,1)*rVariables.Nu(node) + CrossEta[GPointCross] * B(INDEX_2D_BEAM_XX, index_x);
        B(INDEX_2D_BEAM_XX, index_t) = -t * rVariables.NodalCrossDirection(1, node) * term_xx;

        // eps_yy due to theta
        const double term_yy = InvJ(1,1)*rVariables.Nu(node) + CrossEta[GPointCross] * B(INDEX_2D_BEAM_YY, index_y);
        B(INDEX_2D_BEAM_YY, index_t) =  t * rVariables.NodalCrossDirection(0, node) * term_yy;

        // eps_xy due to theta
        B(INDEX_2D_BEAM_XY, index_t) = - t * rVariables.NodalCrossDirection(1, node) * term_yy
                                       + t * rVariables.NodalCrossDirection(0, node) * term_xx;
    }

    // KRATOS_INFO("B") << B << std::endl;
    // KRATOS_INFO("T") << rVariables.TransformationMatrix << std::endl;

    if (BTransformed.size1() != rVariables.TransformationMatrix.size1() ||
        BTransformed.size2() != B.size2()     )
        BTransformed.resize( rVariables.TransformationMatrix.size1(), B.size2(), false );

    // KRATOS_INFO("0-BTransformed") << BTransformed << std::endl;

    noalias(BTransformed) = prod(rVariables.TransformationMatrix, B);

    KRATOS_INFO("1-BTransformed") << BTransformed << std::endl;

    // noalias(BTransformed) = ZeroMatrix(BTransformed.size1(), BTransformed.size2());
    // for (unsigned int idim = 0; idim < TDim; ++idim) {
    //     for (unsigned int idof=0; idof < N_DOF_ELEMENT; ++idof) {
    //         for (unsigned int jvoigt=0; jvoigt < VoigtSize; ++jvoigt) {
    //             BTransformed(idim, idof) += rVariables.TransformationMatrix(idim, jvoigt) * B(jvoigt,idof);
    //         }
    //     }
    // }

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateBMatrix") << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
     CalculateStrainVector(ElementVariables &rVariables) const
{
    KRATOS_TRY
    KRATOS_INFO("0-GeoCurvedBeamElement::CalculateStrainVector") << std::endl;

    noalias(rVariables.StrainVector) = prod(rVariables.B, rVariables.DofValuesVector);

    KRATOS_INFO("1-GeoCurvedBeamElement::CalculateStrainVector") << std::endl;
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
double GeoCurvedBeamElement<TDim,TNumNodes>::
    CalculateIntegrationCoefficient(unsigned int GPointCross,
                                    double detJ,
                                    double weight) const
{
    KRATOS_INFO("weight") << weight << std::endl;
    KRATOS_INFO("CrossWeight[GPointCross]") << CrossWeight[GPointCross] << std::endl;
    KRATOS_INFO("detJ") << detJ << std::endl;

    return weight * CrossWeight[GPointCross] * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
SizeType GeoCurvedBeamElement<TDim,TNumNodes>::
    GetAlongNumberIntegrationPoints() const
{
    return this->GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
SizeType GeoCurvedBeamElement<TDim,TNumNodes>::
    GetCrossNumberIntegrationPoints() const
{
    return N_POINT_CROSS;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void GeoCurvedBeamElement<TDim,TNumNodes>::
    ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult,
                                             Vector& determinants_of_jacobian,
                                             const GeometryData::IntegrationMethod& ThisMethod ) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int integration_points_number = rGeom.IntegrationPointsNumber( ThisMethod );

    // KRATOS_INFO("integration_points_number") << integration_points_number << std::endl;

    if ( integration_points_number == 0 )
        KRATOS_ERROR << "This integration method is not supported " << *this << std::endl;

    if ( rResult.size() != integration_points_number )
        rResult.resize(  rGeom.IntegrationPointsNumber( ThisMethod ), false  );
    if ( determinants_of_jacobian.size() != integration_points_number )
        determinants_of_jacobian.resize(  rGeom.IntegrationPointsNumber( ThisMethod ), false  );
    // KRATOS_INFO("determinants_of_jacobian") << determinants_of_jacobian << std::endl;

    //calculating the local gradients
    const ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( ThisMethod );
    // KRATOS_INFO("DN_De") << DN_De << std::endl;

    //loop over all integration points
    Matrix J(rGeom.WorkingSpaceDimension(),rGeom.LocalSpaceDimension());
    // KRATOS_INFO("J") << J << std::endl;

    Matrix Jinv(rGeom.WorkingSpaceDimension(),rGeom.LocalSpaceDimension());
    // KRATOS_INFO("Jinv") << Jinv << std::endl;

    // KRATOS_INFO("rGeom.WorkingSpaceDimension()") << rGeom.WorkingSpaceDimension() << std::endl;
    // KRATOS_INFO("rGeom.size()")                  << rGeom.size() << std::endl;
    // KRATOS_INFO("rGeom.LocalSpaceDimension()")   << rGeom.LocalSpaceDimension() << std::endl;

    double DetJ;
    for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
    {
        if( rResult[pnt].size1() != rGeom.size() ||
            rResult[pnt].size2() != rGeom.WorkingSpaceDimension()     )
            rResult[pnt].resize( rGeom.size(), rGeom.WorkingSpaceDimension(), false );

        // KRATOS_INFO("0-rResult[pnt]") << pnt << ",  "<< rResult[pnt] << std::endl;

        rGeom.Jacobian(J,pnt, ThisMethod);
        // KRATOS_INFO("1-J") << J << std::endl;

        MathUtils<double>::GeneralizedInvertMatrix( J, Jinv, DetJ );
        // KRATOS_INFO("1-Jinv") << Jinv << std::endl;

        noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
        // KRATOS_INFO("1-rResult[pnt]") << rResult[pnt] << std::endl;

        determinants_of_jacobian[pnt] = DetJ;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class GeoCurvedBeamElement<2,3>;
template class GeoCurvedBeamElement<3,3>;

} // Namespace Kratos
